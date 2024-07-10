# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  Open Energy Transition gGmbH
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import os
import sys
sys.path.append("../submodules/pypsa-eur")
import pypsa
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
from _helpers import mock_snakemake, update_config_from_wildcards, load_network, \
                     change_path_to_pypsa_eur, change_path_to_base, \
                     LINE_LIMITS, CO2L_LIMITS, BAU_HORIZON, replace_multiindex_values


def define_res_share_dataframe():
    # Define index
    idx = ['RES generation (TWh)', 'Total load (TWh)', "RES share (%)"]

    # Define column levels
    col_level_0 = ["2030"]*4 + ["2040"]*4 + ["2050"]*4
    col_level_1 = ["Optimal Renovation and Cost-Optimal Heating", "Optimal Renovation and Electric Heating", 
                   "Limited Renovation and Cost-Optimal Heating", "No Renovation and Electric Heating"] + \
                   ["Optimal Renovation and Cost-Optimal Heating", "Optimal Renovation and Electric Heating", 
                   "Limited Renovation and Cost-Optimal Heating", "No Renovation and Electric Heating"]*2

    # Create a MultiColumns
    multi_cols = pd.MultiIndex.from_arrays([col_level_0, col_level_1], names=['Year', 'Scenario'])

    # initialize DataFrame for storing heat pumps
    df_res_share = pd.DataFrame(index=idx, columns=multi_cols)

    return df_res_share


def compute_load(n):
    # get electricity load
    ac_load_buses = n.loads.query("carrier in 'electricity'").index
    elec_load = n.loads_t.p_set.multiply(n.snapshot_weightings.objective, axis=0)[ac_load_buses].sum().sum()

    # get net heat load (initial heat load - heat saved by retrofitting)
    heat_load = n.loads_t.p_set.multiply(n.snapshot_weightings.objective, axis=0).filter(like="heat").sum().sum()
    retrofitting = n.generators_t.p.multiply(n.snapshot_weightings.objective, axis=0).filter(like="retrofitting").sum().sum()
    net_heat_load = heat_load - retrofitting

    # get land transport EV load
    ev_load_buses = n.loads.query("carrier in 'land transport EV'").index
    ev_load = n.loads_t.p_set.multiply(n.snapshot_weightings.objective, axis=0)[ev_load_buses].sum().sum()

    # total electricity + heat load
    total_load = elec_load + net_heat_load + ev_load

    return total_load


def compute_RES_generation(n):
    # carriers for RES generation
    transmission_carrier = ['offwind-ac', 'onwind', 'ror', 'solar', 'offwind-dc']
    transmission_discharger_carrier = ['battery discharger', 'H2 Fuel Cell']
    transmission_charger_carrier = ['battery charger', 'H2 Electrolysis']
    distribution_carrier = ['solar rooftop']
    distribution_discharger_carrier = ['home battery discharger']
    distribution_charger_carrier = ['home battery charger']
    heat_carrier = ['residential rural solar thermal', 'services rural solar thermal',
                    'residential urban decentral solar thermal', 'services urban decentral solar thermal',
                    'urban central solar thermal']
    
    # electricity RES generation in transmission level
    trans_gens = n.generators.query("carrier in @transmission_carrier").index
    trans_generation = n.generators_t.p.multiply(n.snapshot_weightings.objective, axis=0)[trans_gens].sum().sum()
    
    # electricity discharge of transmission level stores
    trans_dischargers = n.links.query("carrier in @transmission_discharger_carrier").index
    trans_discharge = -n.links_t.p1.multiply(n.snapshot_weightings.objective, axis=0)[trans_dischargers].sum().sum()

    # electricity charge of transmission level store
    trans_chargers = n.links.query("carrier in @transmission_charger_carrier").index
    trans_charge = n.links_t.p0.multiply(n.snapshot_weightings.objective, axis=0)[trans_chargers].sum().sum()

    # electricity RES generation in distribution level
    dist_gens = n.generators.query("carrier in @distribution_carrier").index
    dist_generation = n.generators_t.p.multiply(n.snapshot_weightings.objective, axis=0)[dist_gens].sum().sum()

    # electricity discharge of distribution level store
    dist_dischargers = n.links.query("carrier in @distribution_discharger_carrier").index
    dist_discharge = -n.links_t.p1.multiply(n.snapshot_weightings.objective, axis=0)[dist_dischargers].sum().sum()

    # electricity charge of distribution level store
    dist_chargers = n.links.query("carrier in @distribution_charger_carrier").index
    dist_charge = n.links_t.p0.multiply(n.snapshot_weightings.objective, axis=0)[dist_chargers].sum().sum()

    # heat RES generation
    heat_gens = n.generators.query("carrier in @heat_carrier").index
    heat_generation = n.generators_t.p[heat_gens].sum().sum()

    # total res consumption
    total_res_gen = (trans_generation + trans_discharge - trans_charge) + (dist_generation + dist_discharge - dist_charge) + heat_generation

    return total_res_gen


def compute_res_ratio(total_load, total_res_consumed):
    res_ratio = 100 * total_res_consumed / total_load
    return res_ratio


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "get_res_share", 
            clusters="48",
        )
    # update config based on wildcards
    config = update_config_from_wildcards(snakemake.config, snakemake.wildcards)


    # move to submodules/pypsa-eur
    change_path_to_pypsa_eur()
    # network parameters
    co2l_limits = CO2L_LIMITS
    line_limits = LINE_LIMITS
    clusters = config["plotting"]["clusters"]
    planning_horizons = config["plotting"]["planning_horizon"]
    planning_horizons = [str(x) for x in planning_horizons if not str(x) == BAU_HORIZON]
    opts = config["plotting"]["sector_opts"]

    # define scenario namings
    scenarios = {"flexible": "Optimal Renovation and Cost-Optimal Heating", 
                 "retro_tes": "Optimal Renovation and Electric Heating", 
                 "flexible-moderate": "Limited Renovation and Cost-Optimal Heating", 
                 "rigid": "No Renovation and Electric Heating"}

    # define heat pumps dataframe
    df_res_share = define_res_share_dataframe()

    # heat pumps estimation
    for planning_horizon in planning_horizons:
        lineex = line_limits[planning_horizon]
        sector_opts = f"Co2L{co2l_limits[planning_horizon]}-{opts}"

        for scenario, nice_name in scenarios.items():
            # load networks
            n = load_network(lineex, clusters, sector_opts, planning_horizon, scenario)

            if n is None:
                # Skip further computation for this scenario if network is not loaded
                print(f"Network is not found for scenario '{scenario}', planning year '{planning_horizon}'. Skipping...")
                continue

            # compute total load
            total_load = compute_load(n)

            # compute res consumption
            total_res_consumed = compute_RES_generation(n)

            # compute RES share
            res_share = compute_res_ratio(total_load=total_load, total_res_consumed=total_res_consumed)

            # store into table
            df_res_share.loc[:, (planning_horizon, nice_name)] = [round(total_res_consumed/1e6, 2), round(total_load/1e6, 2), round(res_share, 2)]


    # add BAU scenario
    BAU_horizon = BAU_HORIZON
    scenario = "BAU"
    lineex = line_limits[BAU_horizon]
    sector_opts = f"Co2L{co2l_limits[BAU_horizon]}-{opts}"

    # load BAU network
    n = load_network(lineex, clusters, sector_opts, BAU_horizon, scenario)

    # move to base directory
    change_path_to_base()

    if n is None:
        # Skip further computation for this scenario if network is not loaded
        print(f"Network is not found for scenario '{scenario}', planning year '{BAU_horizon}'. Skipping...")
    else:
        # compute total load
        total_load = compute_load(n)

        # compute res consumption
        total_res_consumed = compute_RES_generation(n)

        # compute RES share
        res_share = compute_res_ratio(total_load=total_load, total_res_consumed=total_res_consumed)

        # store into table
        df_res_share.loc[:, (2020, "BAU")] = [round(total_res_consumed/1e6, 2), round(total_load/1e6, 2), round(res_share, 2)]


    # move to base directory
    change_path_to_base()

    # save the heat pumps data in Excel format
    df_res_share.columns = replace_multiindex_values(df_res_share.columns, 
                                                     ("2040", "Limited Renovation and Cost-Optimal Heating"),
                                                     ("2040", "Limited Renovation and Electric Heating"))
    df_res_share.columns = replace_multiindex_values(df_res_share.columns, 
                                                     ("2050", "Limited Renovation and Cost-Optimal Heating"),
                                                     ("2050", "Limited Renovation and Electric Heating"))
    df_res_share.to_csv(snakemake.output.table)
