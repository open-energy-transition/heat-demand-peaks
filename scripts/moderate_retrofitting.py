import os
import sys
sys.path.append(os.path.abspath(os.path.join(__file__ ,"../../")))
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
from plots._helpers import mock_snakemake, update_config_from_wildcards, load_network, \
                    change_path_to_pypsa_eur, change_path_to_base, load_unsolved_network, \
                    save_unsolved_network


def set_moderate_retrofitting(n_solved, n_unsolved):
    # get optimal retrofitting from the solved network
    retro_opt = n_solved.generators.query("carrier in 'retrofitting'")[["p_nom_opt"]]
    retro_data = retro_opt / 2
    # set p_nom as half of p_nom_opt of flexible scenario
    n_unsolved.generators.loc[retro_data.index, "p_nom"] = retro_data["p_nom_opt"]
    # set retrofitting not extendable
    n_unsolved.generators.loc[retro_data.index, "p_nom_extendable"] = False

    return n_unsolved


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "moderate_retrofitting", 
            clusters="48",
            planning_horizon="2030",
        )
    # update config based on wildcards
    config = update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    # network parameters by year
    co2l_limits = {"2030":"0.45", "2040":"0.1", "2050":"0.0"}
    line_limits = {"2030":"v1.15", "2040":"v1.3", "2050":"v1.5"}

    # network parameters of unsolved network
    clusters = config["moderate_retrofitting"]["clusters"]
    planning_horizon = config["moderate_retrofitting"]["planning_horizon"]
    time_resolution = config["moderate_retrofitting"]["time_resolution"]
    lineex = line_limits[planning_horizon]
    sector_opts = f"Co2L{co2l_limits[planning_horizon]}-{time_resolution}-T-H-B-I"

    # move to pypsa-eur directory
    change_path_to_pypsa_eur()

    # load solved network of flexible scenario
    n_solved = load_network(lineex, clusters, sector_opts, planning_horizon, "flexible")

    # load unsolved network of flexible-moderate scenario
    n_unsolved = load_unsolved_network(lineex, clusters, sector_opts, planning_horizon, "flexible-moderate")

    if not n_solved is None and not n_unsolved is None:
        # update the network by setting p_nom_opt/2 of flexible scenario as p_nom for flexible-moderate scenario
        n_updated = set_moderate_retrofitting(n_solved, n_unsolved)

        # save updated network
        try:
            save_unsolved_network(n_updated, lineex, clusters, sector_opts, planning_horizon, "flexible-moderate")
            success = True
        except Exception as e:
            print(f"Error: {e}")
            success = False
    else:
        success = False

    # move to base directory
    change_path_to_base()

    # write logs
    with open(snakemake.output.logs, 'w') as f:
        f.write(f"""Planning horizon: {planning_horizon} 
                \nClusters: {clusters} 
                \nSuccess: {success}""")
