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
                     change_path_to_pypsa_eur, change_path_to_base

logger = logging.getLogger(__name__)

RESULTS_DIR = "plots/results"


def get_curtailment(n, nice_name):
    curtailments = n.statistics()[["Curtailment"]]
    techs = {"Solar": ["solar rooftop", "Solar"],
             "Wind": ["Offshore Wind (AC)", "Offshore Wind (DC)", "Onshore Wind"]
             }
    
    curtailment_dict = {}

    for name, tech in techs.items():
        curtailment_dict[name] = curtailments[curtailments.index.get_level_values(1).isin(tech)].sum().item()

    curtailment_dict["Total"] = sum(curtailment_dict.values())

    return curtailment_dict


def plot_curtailment(df_curtailment):
    pass
    

def define_table_df(scenarios):
    # Define column levels
    col_level_0 = ["2030"]*4 + ["2040"]*4 + ["2050"]*4
    col_level_1 = list(scenarios.values()) * 3
    # Create a MultiColumns
    multi_cols = pd.MultiIndex.from_arrays([col_level_0, col_level_1], names=['Year', 'Scenario'])
    df = pd.DataFrame(columns=multi_cols, index=["Solar", "Wind", "Total"])
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
    co2l_limits = {"2030":"0.45", "2040":"0.1", "2050":"0.0"}
    line_limits = {"2030":"v1.15", "2040":"v1.3", "2050":"v1.5"}
    clusters = config["plotting"]["clusters"]
    planning_horizons = ["2030", "2040", "2050"]
    time_resolution = config["plotting"]["time_resolution"]

    # define scenario namings
    scenarios = {"flexible": "Optimal Renovation and Heating", 
                "retro_tes": "Optimal Renovation and Green Heating", 
                "flexible-moderate": "Limited Renovation and Optimal Heating", 
                "rigid": "No Renovation and Green Heating"}


    # initialize df for storing curtailment information
    curtailment_df = define_table_df(scenarios)

    for planning_horizon in planning_horizons:
        lineex = line_limits[planning_horizon]
        sector_opts = f"Co2L{co2l_limits[planning_horizon]}-{time_resolution}-T-H-B-I"
    
        # load networks
        for scenario, nice_name in scenarios.items():
            n = load_network(lineex, clusters, sector_opts, planning_horizon, scenario)

            if n is None:
                # Skip further computation for this scenario if network is not loaded
                print(f"Network is not found for scenario '{scenario}', planning year '{planning_horizon}', and time resolution of '{time_resolution}'. Skipping...")
                continue
            
            curtailment_dict = get_curtailment(n, nice_name)
            print(curtailment_dict)
            
    #         # calculate capital costs for scenario
    #         cap_costs = compute_costs(n, nice_name, "Capital")
    #         cap_cost_df = cap_cost_df.join(cap_costs, how="outer").fillna(0)

    #         # calculate operational costs for scenario
    #         op_costs = compute_costs(n, nice_name, "Operational")
    #         op_cost_df = op_cost_df.join(op_costs, how="outer").fillna(0)

    #         # calculate capacities for scenario
    #         capacities = compute_capacities(n, nice_name)
    #         capacities_df = capacities_df.join(capacities, how="outer").fillna(0)

    #         # get p_nom_opt for scenario
    #         p_nom_opt = get_p_nom_opt(n, nice_name)
    #         p_nom_opt_df = p_nom_opt_df.join(p_nom_opt, how="outer").fillna(0)

    #     # update capital costs based on previous horizons
    #     cap_costs_dict[planning_horizon] = cap_cost_df
    #     p_nom_opt_dict[planning_horizon] = p_nom_opt_df

    #     # if capital costs data is present (not empty)
    #     if not any([x.empty for x in cap_costs_dict.values()]):
    #         updated_caps_df = update_capital_cost(cap_costs_dict, p_nom_opt_dict, planning_horizon)

    #         # add capital and operational costs
    #         cost_df = sum_costs(updated_caps_df, op_cost_df)
    #         # reorder scenarios
    #         reorder_columns = [s for s in scenarios.values() if s in cost_df.columns]
    #         cost_df = cost_df[reorder_columns]

    #     # move to base directory
    #     change_path_to_base()

    #     # plot costs
    #     if not cost_df.empty:
    #         processed_cost_df = plot_costs(cost_df, clusters, planning_horizon)
    #         table_cost_df = fill_table_df(table_cost_df, planning_horizon, scenarios, processed_cost_df)

    #     # plot capacities
    #     if not capacities_df.empty:
    #         processed_capacities_df = plot_capacities(capacities_df, clusters, planning_horizon)
    #         table_cap_df = fill_table_df(table_cap_df, planning_horizon, scenarios, processed_capacities_df)

        
    # # save all costs to csv
    # if not table_cost_df.empty:
    #     table_cost_df.index.name = "System cost [EUR billion per year]"
    #     table_cost_df.to_csv(snakemake.output.costs)

    # # save all capacities to csv
    # if not table_cap_df.empty:
    #     table_cap_df.index.name = "Installed capacity [GW]"
    #     table_cap_df.to_csv(snakemake.output.capacities) 
