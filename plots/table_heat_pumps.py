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


def define_heat_pump_dataframe():
    # Define index levels
    index_level_0 = ['Optimal Heat Pump Capacity [GW_el]', 
                     'Optimal Heat Pump Capacity [GW_el]',
                     'Optimal Heat Pump Capacity [GW_el]', 
                     'Optimal Heat Pump Capacity [GW_el]',
                     'Number of heat pumps [millions]',
                     'Number of heat pumps [millions]', 
                     'Number of heat pumps [millions]',  
                     'Number of heat pumps [millions]']
    index_level_1 = ['Residential decentral', 
                     'Services', 
                     'Central',
                     'Total',
                     'Residential decentral', 
                     'Services', 
                     'Central',
                     'Total']
    # Create a MultiIndex
    multi_index = pd.MultiIndex.from_arrays([index_level_0, index_level_1])

    # Define column levels
    col_level_0 = ["2030"]*5 + ["2040"]*4 + ["2050"]*4
    col_level_1 = ["Widespread Renovation", "Widespread Renovation and Electrification", 
                   "Limited Renovation", "Business as Usual and Electrification", 
                   "EU action plan (Announced Pledges Scenario) [2]"] + \
                   ["Widespread Renovation", "Widespread Renovation and Electrification", 
                   "Limited Renovation", "Business as Usual and Electrification"]*2

    # Create a MultiColumns
    multi_cols = pd.MultiIndex.from_arrays([col_level_0, col_level_1], names=['Year', 'Scenario'])

    # initialize DataFrame for storing heat pumps
    df_heat_pumps = pd.DataFrame(index=multi_index, columns=multi_cols)

    return df_heat_pumps


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

    # heat pump techs
    heat_pumps = {"residential": "Residential decentral", 
                  "services": "Services",
                  "urban central": "Central"}
    # heat pump capacities for different building types
    hp_capacity = {"residential": config["heat_pumps"]["residential_decentral"],
                   "services": config["heat_pumps"]["services"],
                   "urban central": config["heat_pumps"]["central"]}

    # define scenario namings
    scenarios = {"flexible": "Widespread Renovation",
                 "retro_tes": "Widespread Renovation and Electrification",
                 "flexible-moderate": "Limited Renovation",
                 "rigid": "Business as Usual and Electrification"}

    # define heat pumps dataframe
    df_heat_pumps = define_heat_pump_dataframe()

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

            # compute heat pump capacities
            for h, h_name in heat_pumps.items():
                p_nom_opt = n.links.filter(like=h, axis=0).filter(like="heat pump", axis=0).p_nom_opt.sum()
                cop = n.links_t.efficiency.filter(like=h).filter(like="heat pump").mean().mean()
                
                df_heat_pumps.loc[("Optimal Heat Pump Capacity [GW_el]", h_name), (planning_horizon, nice_name)] = (p_nom_opt / 1e3).round(2)
                df_heat_pumps.loc[("Number of heat pumps [millions]", h_name), (planning_horizon, nice_name)] = (p_nom_opt * cop / (hp_capacity[h] * 1e3)).round(2)

            # estimate heat pump 
            df_heat_pumps.loc[("Number of heat pumps [millions]", "Total"), (planning_horizon, nice_name)] = df_heat_pumps.loc[("Number of heat pumps [millions]", heat_pumps.values()), (planning_horizon, nice_name)].sum()
            df_heat_pumps.loc[("Optimal Heat Pump Capacity [GW_el]", "Total"), (planning_horizon, nice_name)] = df_heat_pumps.loc[("Optimal Heat Pump Capacity [GW_el]", heat_pumps.values()), (planning_horizon, nice_name)].sum()

    # set heat plan amount based on EU action plan
    df_heat_pumps.loc[('Number of heat pumps [millions]', 'Total'), ("2030", "EU action plan (Announced Pledges Scenario) [2]")] = 45
        
    # move to base directory
    change_path_to_base()

    # save the heat pumps data in Excel format
    df_heat_pumps.columns = replace_multiindex_values(df_heat_pumps.columns, 
                                                      ("2040", "Limited Renovation"),
                                                      ("2040", "Limited Renovation and Electrification"))
    df_heat_pumps.columns = replace_multiindex_values(df_heat_pumps.columns, 
                                                      ("2050", "Limited Renovation"),
                                                      ("2050", "Limited Renovation and Electrification"))
    df_heat_pumps.to_csv(snakemake.output.table)
