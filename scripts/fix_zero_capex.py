# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  Open Energy Transition gGmbH
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import os
import sys
sys.path.append(os.path.abspath(os.path.join(__file__ ,"../../")))
import numpy as np
import pandas as pd
import logging
import warnings
warnings.filterwarnings("ignore")
from plots._helpers import mock_snakemake, load_network, load_unsolved_network, \
                    save_network, save_unsolved_network, change_path_to_pypsa_eur, \
                    change_path_to_base


# nominal attribute by component
NOMINAL_ATTRS = {
    "Generator": "p_nom",
    "Link": "p_nom",
    "Store": "e_nom",
    "StorageUnit": "p_nom",
}


def fix_zero_capex_extendables(n, tolerance=1e-3, exclude=None):
    # extendables without capital cost have degenerate optimal capacity, fix them to infinite capacity
    pattern = "|".join(exclude) if exclude else None
    fixed = []
    for c, attr in NOMINAL_ATTRS.items():
        df = n.static(c)
        if df.empty:
            continue
        mask = df[f"{attr}_extendable"] & (df.capital_cost.abs() <= tolerance)
        if pattern:
            mask &= ~df.index.str.contains(pattern)
        # skip both directions of bidirectional links (their p_nom are synced by constraint)
        if c == "Link" and "reversed" in df.columns:
            reversed_i = df.index[df.reversed.fillna(False).astype(bool)]
            mask &= ~df.index.isin(reversed_i.union(reversed_i.str.replace("-reversed", "")))
        idx = df.index[mask]
        if idx.empty:
            continue
        fixed.append(pd.DataFrame({"component": c, "name": idx, "nom": df.loc[idx, attr].values}))
        df.loc[idx, attr] = np.inf
        df.loc[idx, f"{attr}_extendable"] = False
        logging.info(f"Fixed {len(idx)} zero capital cost extendable {c}s to {attr}=inf "
                     f"(carriers: {', '.join(sorted(df.loc[idx, 'carrier'].unique()))})")

    if not fixed:
        return pd.DataFrame(columns=["component", "name", "nom"])
    return pd.concat(fixed, ignore_index=True)


def restore_zero_capex_extendables(n, fixed):
    # make fixed components extendable again with optimal capacity equal to peak dispatch
    for c, group in fixed.groupby("component"):
        attr = NOMINAL_ATTRS[c]
        idx = pd.Index(group.name)
        dynamic = n.dynamic(c)
        if c == "Link":
            peak = dynamic["p0"].reindex(columns=idx).abs().max()
        elif c == "Store":
            peak = dynamic["e"].reindex(columns=idx).abs().max()
        elif c == "StorageUnit":
            peak = np.maximum(dynamic["p_dispatch"].reindex(columns=idx).abs().max(),
                              dynamic["p_store"].reindex(columns=idx).abs().max())
        else:
            peak = dynamic["p"].reindex(columns=idx).abs().max()
        df = n.static(c)
        df.loc[idx, attr] = group.nom.values
        df.loc[idx, f"{attr}_extendable"] = True
        df.loc[idx, f"{attr}_opt"] = peak.fillna(0.0).values
        logging.info(f"Restored {len(idx)} {c}s as extendable with {attr}_opt set to peak dispatch")
    return n


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "fix_zero_capex",
            clusters="48",
            planning_horizon="2030",
            scenario="flexible",
            step="fix",
        )
    logging.basicConfig(level=logging.INFO)

    # network parameters
    config = snakemake.config
    opts = config.get("fix_zero_capex_extendables", {})
    planning_horizon = snakemake.wildcards.planning_horizon
    scenario = snakemake.wildcards.scenario
    step = snakemake.wildcards.step

    # list of fixed components is kept in base directory
    fixed_path = os.path.abspath(
        f"scripts/logs/zero_capex_fixed_{snakemake.wildcards.clusters}_{planning_horizon}_{scenario}.csv"
    )

    # move to pypsa-eur directory
    change_path_to_pypsa_eur()

    if step == "fix":
        n = load_unsolved_network(scenario, planning_horizon)
        if n is None:
            change_path_to_base()
            raise FileNotFoundError(f"Missing composed network for {scenario} {planning_horizon}")
        fixed = pd.DataFrame(columns=["component", "name", "nom"])
        if opts.get("enable", False):
            fixed = fix_zero_capex_extendables(n, opts.get("tolerance", 1e-3), opts.get("exclude", []))
            save_unsolved_network(n, scenario, planning_horizon)
        fixed.to_csv(fixed_path, index=False)
    elif step == "restore":
        fixed = pd.read_csv(fixed_path) if os.path.exists(fixed_path) else pd.DataFrame()
        if not fixed.empty:
            n = load_network(scenario, planning_horizon)
            if n is None:
                change_path_to_base()
                raise FileNotFoundError(f"Missing solved network for {scenario} {planning_horizon}")
            n = restore_zero_capex_extendables(n, fixed)
            save_network(n, scenario, planning_horizon)
        # avoid restoring twice
        if os.path.exists(fixed_path):
            os.remove(fixed_path)
    else:
        change_path_to_base()
        raise ValueError(f"Unknown step {step}, use 'fix' or 'restore'")

    # move to base directory
    change_path_to_base()

    # write logs
    with open(snakemake.output.logs, 'w') as f:
        f.write(f"""Scenario: {scenario}
                \nPlanning horizon: {planning_horizon}
                \nStep: {step}
                \nComponents: {len(fixed)}
                \nSuccess: True""")
