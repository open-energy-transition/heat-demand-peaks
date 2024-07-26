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

def DSR(n, style, sector):

    if sector == "heating":
        techs = ["residential rural heat", "residential urban decentral heat", "urban central heat"]
        query = "carrier in @techs"
    elif sector == "transport":
        query = "carrier == 'Li ion'"

    if style == "upward":
        DSR = n.stores_t.p.multiply(n.snapshot_weightings.objective,axis=0)[n.stores.query(query).index].clip(lower=0).sum().sum()
    elif style == "downward":
        DSR = n.stores_t.p.multiply(n.snapshot_weightings.objective,axis=0)[n.stores.query(query).index].clip(upper=0).sum().sum()
    return DSR

if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "get_heat_pump", 
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
    df_DSR_heat = pd.DataFrame(
            index = [scenario + dsr for scenario in scenarios.keys() for dsr in [" upward", " downward"]],
            columns = planning_horizons,
        )
    df_DSR_transport = pd.DataFrame(
            index = [scenario + dsr for scenario in scenarios.keys() for dsr in [" upward", " downward"]],
            columns = planning_horizons,
        )

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

            # store into table
            df_DSR_heat.loc[scenario + " upward", planning_horizon] = DSR(n, "upward", "heating")/1e6
            df_DSR_heat.loc[scenario + " downward", planning_horizon] = -1*DSR(n, "downward", "heating")/1e6
    
            df_DSR_transport.loc[scenario + " upward", planning_horizon] = DSR(n, "upward", "transport")/1e6
            df_DSR_transport.loc[scenario + " downward", planning_horizon] = DSR(n, "downward", "transport")/1e6


    df_DSR_heat.index = df_DSR_heat.index + " heating"
    df_DSR_transport.index = df_DSR_transport.index + " transport"
    df_DSR = pd.concat([df_DSR_heat, df_DSR_transport])
    print(df_DSR)

    # move to base directory
    change_path_to_base()
        
    df_DSR.to_csv(snakemake.output.table)
