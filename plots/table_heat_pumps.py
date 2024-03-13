import os
import sys
sys.path.append("../submodules/pypsa-eur")
import pypsa
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

# get the base working directory
BASE_PATH = os.path.abspath(os.path.join(__file__ ,"../.."))
# relative path to folder where to store plots
PATH_PLOTS = "plots/results/"


def change_path_to_pypsa_eur():
    # path to pypsa-eur
    pypsa_path = "submodules/pypsa-eur/"
    # absolute path to pypsa-eur
    new_path = os.path.join(BASE_PATH, pypsa_path)
    # change path to pypsa-eur
    os.chdir(new_path)


def load_network(lineex, space_resolution, sector_opts, planning, scenario):
    FILE = f"elec_s_{space_resolution}_l{lineex}__{sector_opts}_{planning}.nc"
    DIR = f"results/{scenario}/postnetworks"
    n = pypsa.Network(os.path.join(DIR, FILE))
    return n


def change_path_to_base():
    os.chdir(BASE_PATH)
    # create folder to store images
    os.makedirs(PATH_PLOTS, exist_ok=True)


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
    col_level_0 = [2030]*5 + [2040]*4 + [2050]*4
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
    # move to submodules/pypsa-eur
    change_path_to_pypsa_eur()
    # network parameters
    co2l_limits = {2030:0.45, 2040:0.1, 2050:0.0}
    line_limits = {2030:"v1.15", 2040:"v1.3", 2050:"v1.5"}
    space_resolution = 48

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
    for planning in [2030, 2040, 2050]:
        lineex = line_limits[planning]
        sector_opts = f"Co2L{co2l_limits[planning]}-1H-T-H-B-I"

        for scenario, nice_name in scenarios.items():
            # load networks
            n = load_network(lineex, space_resolution, sector_opts, planning, scenario)
            # compute heat pump capacities
            for h, h_name in heat_pumps.items():
                p_nom_opt = n.links.filter(like=h, axis=0).p_nom_opt.sum()
                df_heat_pumps.loc[("Optimal Heat Pump Capacity [MW]", h_name), (planning, nice_name)] = p_nom_opt
            # total heat pump capacity
            total_p_nom_opt = n.links.filter(like="heat pump", axis=0).p_nom_opt.sum()
            df_heat_pumps.loc[("Optimal Heat Pump Capacity [MW]", "Total"), (planning, nice_name)] = total_p_nom_opt
            
            # estimate heat pump amount
            max_amount = total_p_nom_opt / MIN_CAPACITY
            min_amount = total_p_nom_opt / MAX_CAPACITY
            df_heat_pumps.loc[('Approximate number of heat pumps [millions]', 'Maximum (830 W) [1]'), (planning, nice_name)] = round(max_amount, 2)
            df_heat_pumps.loc[('Approximate number of heat pumps [millions]', 'Minimum (6900 W) [1]'), (planning, nice_name)] = round(min_amount, 2)

    # set heat plan amount based on EU action plan
    df_heat_pumps.loc[('Approximate number of heat pumps [millions]', 'Maximum (830 W) [1]'), (2030, "EU action plan (Announced Pledges Scenario) [2]")] = 45
    df_heat_pumps.loc[('Approximate number of heat pumps [millions]', 'Minimum (6900 W) [1]'), (2030, "EU action plan (Announced Pledges Scenario) [2]")] = 45
        
    # move to base directory
    change_path_to_base()

    # save the heat pumps data in Excel format
    df_heat_pumps.to_excel(PATH_PLOTS+f"table_heat_pumps_{space_resolution}.xlsx")

    # save the heat pumps data in CSV format
    # df_heat_pumps.to_csv(PATH_PLOTS+"table_heat_pumps.csv")
