import os
import sys
sys.path.append("../submodules/pypsa-eur")
import pypsa
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
from _helpers import mock_snakemake, update_config_from_wildcards

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
    try:
        n = pypsa.Network(os.path.join(DIR, FILE))
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return None
    return n


def change_path_to_base():
    os.chdir(BASE_PATH)
    # create folder to store images
    os.makedirs(PATH_PLOTS, exist_ok=True)


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "get_line_congestion", 
            space_resolution="48",
        )
    # update config based on wildcards
    config = update_config_from_wildcards(snakemake.config, snakemake.wildcards)


    # move to submodules/pypsa-eur
    change_path_to_pypsa_eur()
    # network parameters
    co2l_limits = {"2030":"0.45", "2040":"0.1", "2050":"0.0"}
    line_limits = {"2030":"v1.15", "2040":"v1.3", "2050":"v1.5"}
    space_resolution = config["plotting"]["space_resolution"]
    time_resolution = config["plotting"]["time_resolution"]

    # define scenario namings
    scenarios = {"flexible": "Optimal Renovation and Heating", 
                 "retro_tes": "Optimal Renovation and Green Heating", 
                 "flexible-moderate": "Limited Renovation and Optimal Heating", 
                 "rigid": "No Renovation and Optimal Heating"}

    # define dataframe to store grid congestion
    df_congestion = pd.DataFrame(index=list(scenarios.values()), columns=["2030", "2040", "2050"])

    # line congestion estimation
    for planning in ["2030", "2040", "2050"]:
        lineex = line_limits[planning]
        sector_opts = f"Co2L{co2l_limits[planning]}-{time_resolution}-T-H-B-I"

        for scenario, nice_name in scenarios.items():
            # load networks
            n = load_network(lineex, space_resolution, sector_opts, planning, scenario)

            if n is None:
                # Skip further computation for this scenario if network is not loaded
                print(f"Network is not found for scenario '{scenario}', planning year '{planning}', and time resolution of '{time_resolution}'. Skipping...")
                continue

            # estimate upper and lower limits of congestion of grid
            upper = n.lines_t.mu_upper.multiply(n.snapshot_weightings.stores, axis=0).mean().multiply(n.lines.length).mean()
            lower = n.lines_t.mu_lower.multiply(n.snapshot_weightings.stores, axis=0).mean().multiply(n.lines.length).mean()
            congestion = (upper + lower) / 2
            df_congestion.loc[nice_name, planning] = round((congestion / 1e6), 2)

    # add name for columns
    df_congestion.index.name = "Scenario [m. €/MW]"

    # move to base directory
    change_path_to_base()

    # save the heat pumps data in Excel format
    df_congestion.to_excel(snakemake.output.table)

