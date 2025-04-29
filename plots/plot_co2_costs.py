# SPDX-FileCopyrightText:  Open Energy Transition gGmbH
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import os
import sys
sys.path.append("../submodules/pypsa-eur")
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import logging
import warnings
warnings.filterwarnings("ignore")
from _helpers import mock_snakemake, update_config_from_wildcards, load_network, \
                     change_path_to_pypsa_eur, change_path_to_base, \
                     LINE_LIMITS, CO2L_LIMITS, BAU_HORIZON, replace_multiindex_values

logger = logging.getLogger(__name__)

RESULTS_DIR = "plots/results"


def get_co2_costs(n, nice_name):
    # get co2 costs (shadow price)
    co2_cost = -n.global_constraints.loc["CO2Limit", "mu"]
    df = pd.DataFrame([co2_cost], columns=[nice_name], index=["CO2 costs [EUR/tCO2_eq]"])
    return df


def plot_co2_costs(co2_df, clusters, planning_horizon, plot_width=7):
    # color codes for legend
    color_codes = {"WIDE":"purple",
                   "WIDE\n+ELEC":"limegreen",
                   "LIMIT":"royalblue",
                   "BAU\n+ELEC":"#f4b609",
                   "Baseline\n2023": "grey"}

    # plot co2 costs
    fig, ax = plt.subplots(figsize=(plot_width,9))
    co2_df.iloc[0].plot.bar(ax=ax, legend=False, color=[color_codes[x] for x in co2_df.columns])
    # configure plot
    plt.xticks(rotation=0, fontsize=14)
    ax.set_ylabel("CO$_2$ costs [EUR/tCO$_{2-eq}$]", fontsize=14)
    ax.set_xlabel("")
    ax.set_yticks(np.arange(0, 500, 50))
    ax.set_ylim([0,450])

    x_ticks = list(co2_df.columns)
    if planning_horizon in ["2040", "2050"] and "LIMIT" in x_ticks:
        # replace name for Limited Renovation scenario for 2030 to be LROH
        x_ticks[x_ticks.index("LIMIT")] = "LIMIT\n+ELEC"

    ax.set_xticklabels(x_ticks)

    # Turn off both horizontal and vertical grid lines
    ax.grid(False, which='both')

    if planning_horizon == BAU_HORIZON:
        ax.set_title("2023", fontsize=15)
    else:
        ax.set_title(planning_horizon, fontsize=15)

    # Calculate x-coordinates for the groups
    unique_x_coords = sorted(list(set([bar.get_x() + bar.get_width() / 2 for bar in ax.patches])))
    
    # Annotate each group
    for i, x in enumerate(unique_x_coords):
        plt.annotate(
            f'{co2_df.iloc[0, i]:.2f}', 
            xy=(x, co2_df[co2_df>0].iloc[:,i].sum()), 
            xytext=(x, co2_df[co2_df>0].iloc[:,i].sum()+2), 
            fontsize=15, 
            ha='center'
        )

    ax.set_facecolor('white')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    ax.grid(axis='y', linestyle='--', linewidth=0.5, color='gray', zorder=0)
    plt.savefig(f"{RESULTS_DIR}/plot_co2_costs_{clusters}_{planning_horizon}.png", dpi=600, bbox_inches='tight')


def define_table_df(scenarios):
    # Define column levels
    col_level_0 = ["2030"]*4 + ["2040"]*4 + ["2050"]*4
    col_level_1 = list(scenarios.values()) * 3
    # Create a MultiColumns
    multi_cols = pd.MultiIndex.from_arrays([col_level_0, col_level_1], names=['Year', 'Scenario'])
    df = pd.DataFrame(columns=multi_cols)
    return df


def fill_table_df(df, planning_horizon, scenarios, values):
    for scenario in scenarios.values():
        for tech_name, _ in values.iterrows():
            if scenario in values.columns:
                df.loc[tech_name, (planning_horizon, scenario)] = values.loc[tech_name, scenario]
    return df


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "plot_co2_cost", 
            clusters="48",
        )
    # update config based on wildcards
    config = update_config_from_wildcards(snakemake.config, snakemake.wildcards)


    # network parameters
    co2l_limits = CO2L_LIMITS
    line_limits = LINE_LIMITS
    clusters = config["plotting"]["clusters"]
    opts = config["plotting"]["sector_opts"]
    planning_horizons = config["plotting"]["planning_horizon"]
    planning_horizons = [str(x) for x in planning_horizons if not str(x) == BAU_HORIZON]

    # define scenario namings
    scenarios = {"flexible": "WIDE",
                "retro_tes": "WIDE\n+ELEC",
                "flexible-moderate": "LIMIT",
                "rigid": "BAU\n+ELEC"}


    # initialize df for storing co2 costs information
    table_co2_df = define_table_df(scenarios)

    for planning_horizon in planning_horizons:
        lineex = line_limits[planning_horizon]
        sector_opts = f"Co2L{co2l_limits[planning_horizon]}-{opts}"

        # move to submodules/pypsa-eur
        change_path_to_pypsa_eur()
    
        # load networks
        co2_df = pd.DataFrame()
        for scenario, nice_name in scenarios.items():
            n = load_network(lineex, clusters, sector_opts, planning_horizon, scenario)

            if n is None:
                # Skip further computation for this scenario if network is not loaded
                print(f"Network is not found for scenario '{scenario}', planning year '{planning_horizon}'. Skipping...")
                continue
            
            co2_costs = get_co2_costs(n, nice_name)
            co2_df = co2_df.join(co2_costs, how="outer").fillna(0)

        # move to base directory
        change_path_to_base()

        # plot co2 balance
        if not co2_df.empty:
            plot_co2_costs(co2_df, clusters, planning_horizon)
            table_co2_df = fill_table_df(table_co2_df, planning_horizon, scenarios, co2_df)
    
    # Add BAU scenario
    BAU_horizon = BAU_HORIZON
    scenario = "BAU"
    lineex = line_limits[BAU_horizon]
    sector_opts = f"Co2L{co2l_limits[BAU_horizon]}-{opts}"

    # move to submodules/pypsa-eur
    change_path_to_pypsa_eur()

    n = load_network(lineex, clusters, sector_opts, BAU_horizon, scenario)

    # move to base directory
    change_path_to_base()

    if n is None:
        # Skip further computation for this scenario if network is not loaded
        print(f"Network is not found for scenario '{scenario}', planning year '{BAU_horizon}'. Skipping...")
    else:
        # get co2 balance for BAU and group technologies
        co2_BAU = get_co2_costs(n, "Baseline\n2023")
        if not table_co2_df.empty and not co2_BAU.empty:
            plot_co2_costs(co2_BAU, clusters, BAU_horizon, plot_width=1.5)
            table_co2_df = fill_table_df(table_co2_df, "2023", {"BAU":"Baseline\n2023"}, co2_BAU)
            table_co2_df = table_co2_df.fillna(0)

    # move to base directory
    change_path_to_base()

    # store to csv
    if not table_co2_df.empty:
        # save to csv
        table_co2_df.columns = replace_multiindex_values(table_co2_df.columns, 
                                                         ("2040", "LIMIT"),
                                                         ("2040","LIMIT\n+ELEC"))
        table_co2_df.columns = replace_multiindex_values(table_co2_df.columns, 
                                                         ("2050", "LIMIT"),
                                                         ("2050","LIMIT\n+ELEC"))
        table_co2_df.to_csv(snakemake.output.table)
