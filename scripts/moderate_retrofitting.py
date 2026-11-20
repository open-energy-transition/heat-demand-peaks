# -*- coding: utf-8 -*-
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


def set_moderate_retrofitting(n_solved, n_unsolved):
    # get optimal retrofitting from the solved network
    retro_opt = n_solved.generators.query("carrier in 'retrofitting'")[["p_nom_opt"]]
    retro_data = retro_opt / 2
    # set p_nom_max as half of p_nom_opt of flexible scenario
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


    # network parameters of unsolved network
    clusters = config["moderate_retrofitting"]["clusters"]
    planning_horizon = config["moderate_retrofitting"]["planning_horizon"]

    # move to pypsa-eur directory
    change_path_to_pypsa_eur()

    # load solved network of flexible scenario
    n_solved = load_network("flexible", planning_horizon)

    # load unsolved network of flexible-moderate scenario
    n_unsolved = load_unsolved_network("flexible-moderate", planning_horizon)

    if n_solved is None or n_unsolved is None:
        change_path_to_base()
        raise FileNotFoundError(
            f"Missing networks for moderate retrofitting at {planning_horizon}"
        )

    # set moderate retrofitting
    n_updated = set_moderate_retrofitting(n_solved, n_unsolved)
    save_unsolved_network(n_updated, "flexible-moderate", planning_horizon)
    success = True

    # move to base directory
    change_path_to_base()

    # write logs
    with open(snakemake.output.logs, 'w') as f:
        f.write(f"""Planning horizon: {planning_horizon} 
                \nClusters: {clusters} 
                \nSuccess: {success}""")
