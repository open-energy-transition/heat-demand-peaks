# SPDX-FileCopyrightText:  Open Energy Transition gGmbH
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import os
import sys
sys.path.append("../submodules/pypsa-eur")
import pypsa
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import geopandas as gpd
import cartopy.crs as ccrs
import logging
import colors as c
import warnings
warnings.filterwarnings("ignore")
from _helpers import mock_snakemake, update_config_from_wildcards, load_network, \
                     change_path_to_pypsa_eur, change_path_to_base, \
                     LINE_LIMITS, CO2L_LIMITS, BAU_HORIZON, replace_multiindex_values

logger = logging.getLogger(__name__)

RESULTS_DIR = "plots/results"


def get_curtailment(n, nice_name):
    curtailments = n.statistics()[["Curtailment"]]
    techs = {
             "Solar Rooftop PV": ["solar rooftop"],
             "Solar Utility PV": ["Solar"],
             "Onshore Wind": ["Onshore Wind"],
             "Offshore Wind": ["Offshore Wind (AC)", "Offshore Wind (DC)"]
             }
    
    curtailment_dict = {}

    for name, tech in techs.items():
        curtailment_dict[name] = curtailments[curtailments.index.get_level_values(1).isin(tech)].sum().item()

    curtailment_dict["Total"] = sum(curtailment_dict.values())

    return curtailment_dict


def plot_curtailment(df_curtailment):
    # color codes for legend
    color_codes = {"OROH":"purple", 
                   "OREH":"limegreen", 
                   "LROH":"royalblue", 
                   "NREH":"#f4b609",
                   "BAU": "grey"}
    
    # MWh to TWh
    df_curtailment = df_curtailment / 1e6
    # get BAU year
    BAU_year = df_curtailment.filter(like="BAU", axis=1).columns.get_level_values(0)

    fig, ax = plt.subplots(figsize=(7, 3))
    for nice_name, color_code in color_codes.items():
        # set name for Limited retrofitting for 2040 and 2050
        if planning_horizon in ["2040", "2050"] and nice_name == 'LROH':
            label_name = "LROH/LREH"
        else:
            label_name = nice_name

        if not nice_name == "BAU":
            df_curtailment.loc["Total", (slice(None), nice_name)].plot(ax=ax, color=color_code, 
                                                                       linewidth=2, marker='o', label=label_name, zorder=5)
        elif nice_name == "BAU" and not BAU_year.empty:
            ax.axhline(y=df_curtailment.loc["Total", (BAU_year, nice_name)].values, 
                       color=color_code, linestyle='--', label=nice_name, zorder=1)

    unique_years = sorted(set(df_curtailment.columns.get_level_values(0)) - set(BAU_year))
    ax.set_xticks(range(len(unique_years)))  # Set the tick locations
    ax.set_xticklabels(unique_years)  # Set the tick labels
    ax.set_ylabel("Curtailment [TWh]")
    ax.set_xlabel(None)
    ax.legend(loc="upper left", facecolor="white", fontsize='xx-small')
    plt.savefig(snakemake.output.figure, dpi=600, bbox_inches = 'tight')
    

def define_table_df(scenarios):
    # Define column levels
    col_level_0 = ["2030"]*4 + ["2040"]*4 + ["2050"]*4
    col_level_1 = list(scenarios.values()) * 3
    # Create a MultiColumns
    multi_cols = pd.MultiIndex.from_arrays([col_level_0, col_level_1], names=['Year', 'Scenario'])
    df = pd.DataFrame(columns=multi_cols, index=["Solar Rooftop PV", "Solar Utility PV", 
                                                 "Onshore Wind", "Offshore Wind", "Total"])
    return df


def fill_table_df(df, planning_horizon, scenarios, values):
    for scenario in scenarios.values():
        for tech_name, _ in values.iterrows():
            df.loc[tech_name, (planning_horizon, scenario)] = values.loc[tech_name, scenario]
    return df


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "plot_curtailment", 
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
    opts = config["plotting"]["sector_opts"]
    planning_horizons = config["plotting"]["planning_horizon"]
    planning_horizons = [str(x) for x in planning_horizons if not str(x) == BAU_HORIZON]

    # define scenario namings
    scenarios = {"flexible": "OROH", 
                "retro_tes": "OREH", 
                "flexible-moderate": "LROH", 
                "rigid": "NREH"}


    # initialize df for storing curtailment information
    curtailment_df = define_table_df(scenarios)

    for planning_horizon in planning_horizons:
        lineex = line_limits[planning_horizon]
        sector_opts = f"Co2L{co2l_limits[planning_horizon]}-{opts}"
    
        # load networks
        for scenario, nice_name in scenarios.items():
            n = load_network(lineex, clusters, sector_opts, planning_horizon, scenario)

            if n is None:
                # Skip further computation for this scenario if network is not loaded
                print(f"Network is not found for scenario '{scenario}', planning year '{planning_horizon}'. Skipping...")
                continue
            
            curtailment_dict = get_curtailment(n, nice_name)
            curtailment_df.loc[:, (planning_horizon, nice_name)] = pd.Series(curtailment_dict)
    
    # Add BAU scenario
    BAU_horizon = BAU_HORIZON
    scenario = "BAU"
    lineex = line_limits[BAU_horizon]
    sector_opts = f"Co2L{co2l_limits[BAU_horizon]}-{opts}"

    n = load_network(lineex, clusters, sector_opts, BAU_horizon, scenario)

    if n is None:
        # Skip further computation for this scenario if network is not loaded
        print(f"Network is not found for scenario '{scenario}', planning year '{BAU_horizon}'. Skipping...")
    else:
        curtailment_dict = get_curtailment(n, scenario)
        curtailment_df.loc[:, (BAU_horizon, scenario)] = pd.Series(curtailment_dict)

    # move to base directory
    change_path_to_base()

    # store to csv and png
    if not curtailment_df.empty:
        # make plot
        plot_curtailment(curtailment_df)
        # save to csv
        curtailment_df.columns = replace_multiindex_values(curtailment_df.columns, 
                                                           ("2040", "LROH"),
                                                           ("2040", "LREH"))
        curtailment_df.columns = replace_multiindex_values(curtailment_df.columns, 
                                                           ("2050", "LROH"),
                                                           ("2050", "LREH"))
        curtailment_df.to_csv(snakemake.output.table)
