import os
import sys
sys.path.append(os.path.abspath(os.path.join(__file__ ,"../../")))
import pypsa
import pandas as pd
import logging
import warnings
warnings.filterwarnings("ignore")
from plots._helpers import mock_snakemake, update_config_from_wildcards, load_network, \
                    change_path_to_pypsa_eur, change_path_to_base, load_unsolved_network, \
                    save_unsolved_network, LINE_LIMITS, CO2L_LIMITS


def get_capacities(network, capacity="opt"):
    # extendable generators and capacities
    ext_gen = network.generators.query("p_nom_extendable==True")
    ext_gen_cap = ext_gen["p_nom_"+capacity]
    
    # extendable stores and capacities
    ext_store = network.stores.query("e_nom_extendable==True")
    ext_store_cap = ext_store["e_nom_"+capacity]

    # extendable links and capacities
    ext_link = network.links.query("p_nom_extendable==True")
    ext_link_cap = ext_link["p_nom_"+capacity]

    return ext_gen_cap, ext_store_cap, ext_link_cap


def get_common_index(list1, list2):
    # get common index
    common_index = set(list1).intersection(set(list2))
    return common_index


def set_optimal_capacities(solved_network, unsolved_network):
    # get optimal capacities from previous horizon
    opt_gen_cap, opt_store_cap, opt_link_cap = get_capacities(solved_network, "opt")

    # get minimum capacities from planned horizon
    min_gen_cap, min_store_cap, min_link_cap = get_capacities(unsolved_network, "min")

    # get common extendable generators and stores
    common_gen = get_common_index(opt_gen_cap.index, min_gen_cap.index)
    common_store = get_common_index(opt_store_cap.index, min_store_cap.index)
    common_link = get_common_index(opt_link_cap.index, min_link_cap.index)

    # remove generators with coal, gas, and oil carriers
    fossil_carriers = ["coal", "oil", "gas", "co2", "co2 stored", "co2 sequestered",
                       "land transport oil", "solid biomass for industry",
                       "gas for industry", "shipping methanol", "shipping oil",
                       "naphtha for industry", "kerosene for aviation",
                       "process emissions", "coal for industry", "nuclear"]
    fossil_gens = unsolved_network.generators.query("carrier in @fossil_carriers").index
    fossil_stores = unsolved_network.stores.query("carrier in @fossil_carriers").index
    fossil_links = unsolved_network.links.query("carrier in @fossil_carriers").index
    target_gens = list(set(common_gen).difference(set(fossil_gens)))
    target_stores = list(set(common_store).difference(set(fossil_stores)))
    target_links = list(set(common_link).difference(set(fossil_links)))

    # set p_nom_min and e_nom_min capacities
    unsolved_network.generators.loc[target_gens, "p_nom_min"] = opt_gen_cap[target_gens]
    unsolved_network.stores.loc[target_stores, "e_nom_min"] = opt_store_cap[target_stores]
    unsolved_network.links.loc[target_links, "p_nom_min"] = opt_link_cap[target_links]

    # set p_nom_max as max of p_nom_min and p_nom_max
    unsolved_network.generators.loc[target_gens, "p_nom_max"] = unsolved_network.generators.loc[target_gens, ["p_nom_min", "p_nom_max"]].max(axis=1)
    unsolved_network.stores.loc[target_stores, "e_nom_min"] = unsolved_network.stores.loc[target_stores, ["e_nom_min", "e_nom_max"]].max(axis=1)
    unsolved_network.links.loc[target_links, "p_nom_min"] = unsolved_network.links.loc[target_links, ["p_nom_min", "p_nom_max"]].max(axis=1)

    # set p_nom_max as max(p_nom_max, p_nom_min) for retrofitting
    retro_idx = unsolved_network.generators.query("carrier == 'retrofitting'").index
    unsolved_network.generators.loc[retro_idx, "p_nom_max"] = unsolved_network.generators.loc[retro_idx, ["p_nom_max", "p_nom_min"]].max(axis=1)

    return unsolved_network


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "set_capacities", 
            clusters="48",
            planning_horizon="2040",
            scenario="flexible"
        )
    # update config based on wildcards
    config = update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    # network parameters by year
    co2l_limits = CO2L_LIMITS
    line_limits = LINE_LIMITS
    previous_horizons = {"2040":"2030", "2050":"2040"} 

    # network parameters of unsolved network
    clusters = config["set_capacities"]["clusters"]
    planning_horizon = config["set_capacities"]["planning_horizon"]
    time_resolution = config["set_capacities"]["time_resolution"]
    scenario = config["set_capacities"]["scenario"]
    opts = config["set_capacities"]["sector_opts"]
    lineex = line_limits[planning_horizon]
    sector_opts = f"Co2L{co2l_limits[planning_horizon]}-{time_resolution}-{opts}"

    # network parameters of solved network
    previous_horizon = previous_horizons[planning_horizon]
    previous_lineex = line_limits[previous_horizon]
    previous_sector_opts = f"Co2L{co2l_limits[previous_horizon]}-{time_resolution}-{opts}"

    # move to pypsa-eur directory
    change_path_to_pypsa_eur()

    # load solved network
    solved_network = load_network(previous_lineex, clusters, previous_sector_opts, previous_horizon, scenario)
    
    # load unsolved network for future
    unsolved_network = load_unsolved_network(lineex, clusters, sector_opts, planning_horizon, scenario)

    # if both network is present, then set optimal capacities
    if not solved_network is None and not unsolved_network is None:
        updated_network = set_optimal_capacities(solved_network, unsolved_network)
        print(f"Updated: {lineex}, {clusters}, {sector_opts}, {planning_horizon}, {scenario}")
        # save updated network
        try:
            save_unsolved_network(updated_network, lineex, clusters, sector_opts, planning_horizon, scenario)
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
        f.write(f"""Scenarios: {scenario} 
                \nPlanning horizon: {planning_horizon} 
                \nClusters: {clusters} 
                \nSuccess: {success}""")