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
                     change_path_to_pypsa_eur, change_path_to_base


def define_heat_pump_dataframe():
    # Define index levels
    index_level_0 = ['Optimal Heat Pump Capacity [MW]', 
                     'Optimal Heat Pump Capacity [MW]',
                     'Optimal Heat Pump Capacity [MW]', 
                     'Approximate number of heat pumps [millions]', 
                     'Approximate number of heat pumps [millions]']
    index_level_1 = ['Air-sourced', 
                     'Ground', 
                     'Total', 
                     'Maximum (830 W) [1]', 
                     'Minimum (6900 W) [1]']
    # Create a MultiIndex
    multi_index = pd.MultiIndex.from_arrays([index_level_0, index_level_1])

    # Define column levels
    col_level_0 = ["2030"]*5 + ["2040"]*4 + ["2050"]*4
    col_level_1 = ["Optimal Renovation and Heating", "Optimal Renovation and Green Heating", 
                   "Limited Renovation and Optimal Heating", "No Renovation and Optimal Heating", 
                   "EU action plan (Announced Pledges Scenario) [2]"] + \
                   ["Optimal Renovation and Heating", "Optimal Renovation and Green Heating", 
                   "Limited Renovation and Optimal Heating", "No Renovation and Optimal Heating"]*2

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
    co2l_limits = {"2030":"0.45", "2040":"0.1", "2050":"0.0"}
    line_limits = {"2030":"v1.15", "2040":"v1.3", "2050":"v1.5"}
    clusters = config["plotting"]["clusters"]
    time_resolution = config["plotting"]["time_resolution"]

    # heat pump techs
    heat_pumps = {"air heat pump": "Air-sourced", "ground heat pump": "Ground"}
    # heat pump min and max capacities in W
    MAX_CAPACITY = 6900 # https://www.energysage.com/electricity/house-watts/how-many-watts-does-an-air-source-heat-pump-use/
    MIN_CAPACITY = 830 # https://www.energysage.com/electricity/house-watts/how-many-watts-does-an-air-source-heat-pump-use/

    # define scenario namings
    scenarios = {"flexible": "Optimal Renovation and Heating", 
                 "retro_tes": "Optimal Renovation and Green Heating", 
                 "flexible-moderate": "Limited Renovation and Optimal Heating", 
                 "rigid": "No Renovation and Optimal Heating"}

    # define heat pumps dataframe
    df_heat_pumps = define_heat_pump_dataframe()

    # heat pumps estimation
    for planning_horizon in ["2030", "2040", "2050"]:
        lineex = line_limits[planning_horizon]
        sector_opts = f"Co2L{co2l_limits[planning_horizon]}-{time_resolution}-T-H-B-I"

        for scenario, nice_name in scenarios.items():
            # load networks
            n = load_network(lineex, clusters, sector_opts, planning_horizon, scenario)

            if n is None:
                # Skip further computation for this scenario if network is not loaded
                print(f"Network is not found for scenario '{scenario}', planning year '{planning_horizon}', and time resolution of '{time_resolution}'. Skipping...")
                continue

            # compute heat pump capacities
            for h, h_name in heat_pumps.items():
                p_nom_opt = n.links.filter(like=h, axis=0).p_nom_opt.sum()
                df_heat_pumps.loc[("Optimal Heat Pump Capacity [MW]", h_name), (planning_horizon, nice_name)] = p_nom_opt
            # total heat pump capacity
            total_p_nom_opt = n.links.filter(like="heat pump", axis=0).p_nom_opt.sum()
            df_heat_pumps.loc[("Optimal Heat Pump Capacity [MW]", "Total"), (planning_horizon, nice_name)] = total_p_nom_opt
            
            # estimate heat pump amount
            max_amount = total_p_nom_opt / MIN_CAPACITY
            min_amount = total_p_nom_opt / MAX_CAPACITY
            df_heat_pumps.loc[('Approximate number of heat pumps [millions]', 'Maximum (830 W) [1]'), (planning_horizon, nice_name)] = round(max_amount, 2)
            df_heat_pumps.loc[('Approximate number of heat pumps [millions]', 'Minimum (6900 W) [1]'), (planning_horizon, nice_name)] = round(min_amount, 2)

    # set heat plan amount based on EU action plan
    df_heat_pumps.loc[('Approximate number of heat pumps [millions]', 'Maximum (830 W) [1]'), ("2030", "EU action plan (Announced Pledges Scenario) [2]")] = 45
    df_heat_pumps.loc[('Approximate number of heat pumps [millions]', 'Minimum (6900 W) [1]'), ("2030", "EU action plan (Announced Pledges Scenario) [2]")] = 45
        
    # move to base directory
    change_path_to_base()

    # save the heat pumps data in Excel format
    df_heat_pumps.to_csv(snakemake.output.table)
