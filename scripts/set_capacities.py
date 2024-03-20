import os
import sys
sys.path.append(os.path.abspath(os.path.join(__file__ ,"../../")))
import pypsa
import pandas as pd
import logging
import warnings
warnings.filterwarnings("ignore")
print(os.getcwd())
print(sys.path)
from plots._helpers import mock_snakemake, update_config_from_wildcards, load_network, \
                    change_path_to_pypsa_eur, change_path_to_base


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "set_capacities", 
            clusters="48",
            planning_horizon="2030",
            scenario="flexible"
        )
    # update config based on wildcards
    config = update_config_from_wildcards(snakemake.config, snakemake.wildcards)


    # network parameters
    co2l_limits = {"2030":"0.45", "2040":"0.1", "2050":"0.0"}
    line_limits = {"2030":"v1.15", "2040":"v1.3", "2050":"v1.5"}
    clusters = config["set_capacities"]["clusters"]
    planning_horizon = config["set_capacities"]["planning_horizon"]
    time_resolution = config["set_capacities"]["time_resolution"]
    scenario = config["set_capacities"]["scenario"]
    lineex = line_limits[planning_horizon]
    sector_opts = f"Co2L{co2l_limits[planning_horizon]}-{time_resolution}-T-H-B-I"

    # move to pypsa-eur directory
    change_path_to_pypsa_eur()

    # load solved network
    solved_network = load_network(lineex, clusters, sector_opts, planning_horizon, scenario)
    
    


    # move to base directory
    change_path_to_base()

    # write logs
    with open(snakemake.output.logs, 'w') as f:
        f.write("Done")