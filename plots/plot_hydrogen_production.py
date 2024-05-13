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


def get_hydrogen_production(n, param_dict):
    tech = 'H2 Electrolysis'
    hydrogen_links = n.links.query("carrier in @tech").index

    hydrogen_dict = {}

    for name, param in param_dict.items():
        sign = -1 if param == "p1" else 1
        hydrogen_dict[name] = n.links_t[param][hydrogen_links].multiply(n.snapshot_weightings.stores, axis=0).sum().sum() * sign / 1e6

    return hydrogen_dict


def plot_hydrogen_production(df_hydrogen, figure_dict):
    # color codes for legend
    color_codes = {"Optimal Renovation and Heating":"purple", 
                   "Optimal Renovation and Green Heating":"limegreen", 
                   "Limited Renovation and Optimal Heating":"royalblue", 
                   "No Renovation and Green Heating":"#f4b609",
                   "BAU": "grey"}
    
    # get BAU year
    BAU_year = df_hydrogen.filter(like="BAU", axis=1).columns.get_level_values(0)

    for idx, output in figure_dict.items():
        _, ax = plt.subplots(figsize=(7, 3))
        for nice_name, color_code in color_codes.items():
            if not nice_name == "BAU":
                df_hydrogen.loc[idx, (slice(None), nice_name)].plot(ax=ax, color=color_code, 
                                                                        linewidth=2, marker='o', label=nice_name, zorder=5)
            elif nice_name == "BAU" and not BAU_year.empty:
                ax.axhline(y=df_hydrogen.loc[idx, (BAU_year, nice_name)].values, 
                        color=color_code, linestyle='--', label=nice_name, zorder=1)

        unique_years = sorted(set(df_hydrogen.columns.get_level_values(0)) - set(BAU_year))
        ax.set_xticks(range(len(unique_years)))  # Set the tick locations
        ax.set_xticklabels(unique_years)  # Set the tick labels
        ax.set_ylabel(idx)
        ax.set_xlabel(None)
        ax.legend(loc="upper left", facecolor="white", fontsize='xx-small')
        plt.savefig(output, dpi=600, bbox_inches = 'tight')
        

def define_table_df(scenarios, idx):
    # Define column levels
    col_level_0 = ["2030"]*4 + ["2040"]*4 + ["2050"]*4
    col_level_1 = list(scenarios.values()) * 3
    # Create a MultiColumns
    multi_cols = pd.MultiIndex.from_arrays([col_level_0, col_level_1], names=['Year', 'Scenario'])
    df = pd.DataFrame(columns=multi_cols, index=idx)
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
    time_resolution = config["plotting"]["time_resolution"]

    # define scenario namings
    scenarios = {"flexible": "Optimal Renovation and Heating", 
                "retro_tes": "Optimal Renovation and Green Heating", 
                "flexible-moderate": "Limited Renovation and Optimal Heating", 
                "rigid": "No Renovation and Green Heating"}

    # define dictionary to map to links buses
    param_dict = {"Electricity used [TWh]": "p0",
                  "Hydrogen produced [TWh]": "p1"}
    
    # define dictionary that contain snakemake outputs
    figure_dict = {"Electricity used [TWh]": snakemake.output.figure_elec,
                   "Hydrogen produced [TWh]": snakemake.output.figure_prod}

    # initialize df for storing hydrogen production information
    hydrogen_df = define_table_df(scenarios, idx=param_dict.keys())

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
            
            hydrogen_dict = get_hydrogen_production(n, param_dict=param_dict)
            hydrogen_df.loc[:, (planning_horizon, nice_name)] = pd.Series(hydrogen_dict)
    
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
        hydrogen_dict = get_hydrogen_production(n, param_dict=param_dict)
        hydrogen_df.loc[:, (BAU_horizon, scenario)] = pd.Series(hydrogen_dict)

    # move to base directory
    change_path_to_base()

    # store to csv and png
    if not hydrogen_df.empty:
        # save to csv
        hydrogen_df.to_csv(snakemake.output.table)
        # make plot
        plot_hydrogen_production(hydrogen_df, figure_dict=figure_dict)
