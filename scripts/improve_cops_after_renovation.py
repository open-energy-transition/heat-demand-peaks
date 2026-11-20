# SPDX-FileCopyrightText:  Open Energy Transition gGmbH
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import os
import sys
sys.path.append(os.path.abspath(os.path.join(__file__ ,"../../")))
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
from plots._helpers import mock_snakemake, update_config_from_wildcards, load_network, \
                    change_path_to_pypsa_eur, change_path_to_base, load_unsolved_network, \
                    save_unsolved_network


def improve_cops_after_renovation(n_solved, n_unsolved):
    # get optimal retrofitting from the solved network
    retro_data = n_solved.generators.query("carrier in 'retrofitting'")[["p_nom_opt"]]
    # set p_nom as p_nom_opt of previous run with sink_T are defined 
    n_unsolved.generators.loc[retro_data.index, "p_nom"] = retro_data["p_nom_opt"]
    # set retrofitting not extendable
    n_unsolved.generators.loc[retro_data.index, "p_nom_extendable"] = False

    return n_unsolved


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "improve_cops_after_renovation", 
            clusters="48",
            planning_horizon="2030",
            scenario="flexible"
        )
    # update config based on wildcards
    config = update_config_from_wildcards(snakemake.config, snakemake.wildcards)


    # network parameters of unsolved network
    clusters = config["improve_cops_after_renovation"]["clusters"]
    planning_horizon = config["improve_cops_after_renovation"]["planning_horizon"]
    scenario = config["improve_cops_after_renovation"]["scenario"]

    # move to pypsa-eur directory
    change_path_to_pypsa_eur()

    # load solved network of scenario
    n_solved = load_network(scenario, planning_horizon)

    # load unsolved network of scenario
    n_unsolved = load_unsolved_network(scenario, planning_horizon)

    if n_solved is None or n_unsolved is None:
        change_path_to_base()
        raise FileNotFoundError(
            f"Missing networks for improve_cops: {scenario} {planning_horizon}"
        )

    # improve cops after renovation
    n_updated = improve_cops_after_renovation(n_solved, n_unsolved)
    save_unsolved_network(n_updated, scenario, planning_horizon)
    success = True

    # move to base directory
    change_path_to_base()

    # write logs
    with open(snakemake.output.logs, 'w') as f:
        f.write(f"""Planning horizon: {planning_horizon} 
                \nClusters: {clusters} 
                \nScenario: {scenario}
                \nSuccess: {success}""")
