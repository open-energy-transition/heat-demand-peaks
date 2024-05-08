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
                     LINE_LIMITS, CO2L_LIMITS, BAU_HORIZON

logger = logging.getLogger(__name__)

RESULTS_DIR = "plots/results"


def get_co2(n, nice_name):
    co2_carrier = 'co2'
    co2_bus = n.stores.query("carrier in @co2_carrier").index
    # get co2 level in MtCO2_eq
    co2_level = n.stores_t.e[co2_bus].iloc[-1]
    return co2_level.values


def plot_co2(df_co2):
    # color codes for legend
    color_codes = {"No Renovation and Green Heating":"#f4b609",
                   "BAU": "grey"}
    
    # tCO2_eq to MtCO2_eq
    df_co2 = df_co2 / 1e6
    # get BAU year
    BAU_year = df_co2.filter(like="BAU", axis=1).columns.get_level_values(0)

    fig, ax = plt.subplots(figsize=(7, 3))
    for nice_name, color_code in color_codes.items():
        if not nice_name == "BAU":
            df_co2.loc["Total", (slice(None), nice_name)].plot(ax=ax, color=color_code, 
                                                                       linewidth=2, marker='o', label="Scenarios", zorder=5)
        elif nice_name == "BAU" and not BAU_year.empty:
            ax.axhline(y=df_co2.loc["Total", (BAU_year, nice_name)].values, 
                       color=color_code, linestyle='--', label=nice_name, zorder=1)

    unique_years = sorted(set(df_co2.columns.get_level_values(0)) - set(BAU_year))
    ax.set_xticks(range(len(unique_years)))  # Set the tick locations
    ax.set_xticklabels(unique_years)  # Set the tick labels
    ax.set_ylabel(r"CO$_2$ emissions [MtCO$_{2-eq}$]")
    ax.set_xlabel(None)
    ax.legend(facecolor="white", fontsize='x-small')
    plt.savefig(snakemake.output.figure, dpi=600, bbox_inches = 'tight')
    

def define_table_df(scenarios):
    # Define column levels
    col_level_0 = ["2030"]*4 + ["2040"]*4 + ["2050"]*4
    col_level_1 = list(scenarios.values()) * 3
    # Create a MultiColumns
    multi_cols = pd.MultiIndex.from_arrays([col_level_0, col_level_1], names=['Year', 'Scenario'])
    df = pd.DataFrame(columns=multi_cols, index=["Total"])
    return df


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "plot_co2_level", 
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
    time_resolution = config["plotting"]["time_resolution"]

    # define scenario namings
    scenarios = {"flexible": "Optimal Renovation and Heating", 
                "retro_tes": "Optimal Renovation and Green Heating", 
                "flexible-moderate": "Limited Renovation and Optimal Heating", 
                "rigid": "No Renovation and Green Heating"}


    # initialize df for storing curtailment information
    co2_df = define_table_df(scenarios)

    for planning_horizon in planning_horizons:
        lineex = line_limits[planning_horizon]
        sector_opts = f"Co2L{co2l_limits[planning_horizon]}-{time_resolution}-{opts}"
    
        # load networks
        for scenario, nice_name in scenarios.items():
            n = load_network(lineex, clusters, sector_opts, planning_horizon, scenario)

            if n is None:
                # Skip further computation for this scenario if network is not loaded
                print(f"Network is not found for scenario '{scenario}', planning year '{planning_horizon}', and time resolution of '{time_resolution}'. Skipping...")
                continue
            
            co2_level = get_co2(n, nice_name)
            co2_df.loc["Total", (planning_horizon, nice_name)] = co2_level
    
    # Add BAU scenario
    BAU_horizon = BAU_HORIZON
    scenario = "BAU"
    lineex = line_limits[BAU_horizon]
    sector_opts = f"Co2L{co2l_limits[BAU_horizon]}-{time_resolution}-{opts}"

    n = load_network(lineex, clusters, sector_opts, BAU_horizon, scenario)

    if n is None:
        # Skip further computation for this scenario if network is not loaded
        print(f"Network is not found for scenario '{scenario}', planning year '{BAU_horizon}', and time resolution of '{time_resolution}'. Skipping...")
    else:
        co2_level = get_co2(n, scenario)
        co2_df.loc["Total", (BAU_horizon, scenario)] = co2_level

    # move to base directory
    change_path_to_base()

    # store to csv and png
    if not co2_df.empty:
        # save to csv
        co2_df.to_csv(snakemake.output.table)
        co2_df.index.name = "CO2 emissions [tCO2_eq]"
        # make plot
        plot_co2(co2_df)
