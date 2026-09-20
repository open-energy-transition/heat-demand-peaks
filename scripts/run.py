# SPDX-FileCopyrightText:  Open Energy Transition gGmbH
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import subprocess
import argparse
import logging
import sys
import os
sys.path.append("plots")
from _helpers import (
    change_path_to_pypsa_eur,
    change_path_to_base,
    load_network,
    get_config_path,
    get_n_clusters,
    compose_target,
    solve_target,
    BASE_PATH,
)

# Set up logging configuration
logging.basicConfig(level=logging.INFO)

def get_scenario():
    parser = argparse.ArgumentParser(description="Running the scenario")
    parser.add_argument("-s", "--scenario", help="Specify the scenario.", required=False, 
                        choices=["flexible", "flexible-moderate", "retro_tes", "rigid", "BAU"])
    parser.add_argument("-c", "--continue_horizon", help="Specify the horizon to continue simulations", 
                        choices=["2030", "2040", "2050"])
    parser.add_argument("-y", "--year", help="Specify a single horizon to simulate", 
                        choices=["2030", "2040", "2050"])
    parser.add_argument("-i", "--improved_cop", help="Specify if improved COP calculation is needed",
                        choices=["true", "false"])
    parser.add_argument(
        "--legacy",
        action="store_true",
        help="Include configs/EEE_study/config.legacy_default.yaml (EEE reproduction).",
    )
    args = parser.parse_args()

    # Access the value of the scenario argument
    scenarios = args.scenario
    if scenarios == None:
        scenarios = ["flexible", "flexible-moderate", "retro_tes", "rigid", "BAU"]
    else:
        scenarios = [scenarios]
    # log scenario name
    logging.info(f"Scenarios: {scenarios}")

    # Access the value of horizon from which simulation is continued
    c = args.continue_horizon
    # Access the value of specific horizon to simulate
    y = args.year

    # Access bool for improved COP
    i = args.improved_cop
    improved_cop = False if i == "false" else True
    use_legacy = args.legacy
    logging.info(f"Legacy config overlay: {'enabled' if use_legacy else 'disabled'}")

    # log scenario name
    if y:
        horizons = [int(y)]
        logging.info(f"Start simulating {scenarios} scenario for {y}")
    elif c:
        horizons = get_horizon_list(int(c))
        logging.info(f"Start simulating {scenarios} scenario for {horizons}")
    else:
        c = 2030
        horizons = get_horizon_list(int(c))
        logging.info("No horizon specified. Starting from default horizon (2030)")

    return scenarios, horizons, improved_cop, use_legacy


def get_horizon_list(start_horizon):
    horizons = [2030, 2040, 2050]
    try:
        start_index = horizons.index(start_horizon)
        return horizons[start_index:]
    except ValueError:
        return []


def copy_custom_data(base_network: str = "entsoegridkit"):
    # busmap: new location + name matching clustering.mode: custom_busmap
    # electricity.base_network: entsoegridkit in your EEE configs
    os.makedirs("submodules/pypsa-eur/data/busmaps", exist_ok=True)
    subprocess.run(
        f"cp data/busmaps/simplified_48_{base_network}.csv "
        "submodules/pypsa-eur/data/busmaps/",
        shell=True, check=True,
    )
    subprocess.run(
        "cp data/custom_powerplants.csv submodules/pypsa-eur/data/",
        shell=True, check=True,
    )
    subprocess.run(
        "cp data/custom_costs.csv submodules/pypsa-eur/data/",
        shell=True, check=True,
    )
    # log the success
    logging.info(f"Copied custom data from data/ folder to submodules/pypsa-eur/data/ folder")


def get_configfiles(scenario, horizon, use_legacy=False):
    """Scenario config on top of upstream defaults; optional legacy overlay for EEE reproduction."""
    scenario_cfg = os.path.relpath(get_config_path(scenario, horizon), os.getcwd())
    if not use_legacy:
        return scenario_cfg
    legacy = os.path.relpath(
        os.path.join(BASE_PATH, "configs", "EEE_study", "config.legacy_default.yaml"),
        os.getcwd(),
    )
    return f"{legacy} {scenario_cfg}"


def compose_network(scenario, horizon, use_legacy=False):
    # change path to pypsa-eur
    change_path_to_pypsa_eur()

    # compose target
    target = compose_target(scenario, horizon)

    # compose network (optional legacy defaults, then scenario overrides)
    command = (
        f"snakemake -call {target} "
        f"--configfile {get_configfiles(scenario, horizon, use_legacy=use_legacy)} "
        f"--force --rerun-incomplete"
    )
    result = subprocess.run(command, shell=True)

    # move to base directory
    change_path_to_base()

    if result.returncode != 0:
        raise RuntimeError(f"Compose failed for {scenario} {horizon}")
    logging.info(f"Composed network for {scenario} {horizon}")


def set_capacities(scenario, horizon):
    error = []
    try:
        # get number of clusters
        clusters = get_n_clusters(scenario, horizon)
        command = f"snakemake -call scripts/logs/set_capacities_{clusters}_{horizon}_{scenario}.txt --forceall"
        subprocess.run(command, shell=True, check=True)
        logging.info(f"Capacities are set to {scenario} scenario in {horizon} horizon!")
    except subprocess.CalledProcessError  as e:
        logging.error(f"Error occurred during command execution: {e}")
        error.append(e)
    return error


def moderate_retrofitting(scenario, horizon):
    error = []
    try:
        # get number of clusters
        clusters = get_n_clusters(scenario, horizon)
        command = f"snakemake -call scripts/logs/set_moderate_retrofitting_{clusters}_{horizon}.txt --forceall"
        subprocess.run(command, shell=True, check=True)
        logging.info(f"Moderate retrofitting capacities are set to {scenario} scenario in {horizon} horizon!")
    except subprocess.CalledProcessError  as e:
        logging.error(f"Error occurred during command execution: {e}")
        error.append(e)
    return error


def improve_cops_after_renovation(scenario, horizon):
    error = []
    try:
        # get number of clusters
        clusters = get_n_clusters(scenario, horizon)
        command = f"snakemake -call scripts/logs/improve_cops_after_renovation_{clusters}_{horizon}_{scenario}.txt --forceall"
        subprocess.run(command, shell=True, check=True)
        logging.info(f"Retrofitting capacities are fixed to {scenario} scenario in {horizon} horizon!")
    except subprocess.CalledProcessError  as e:
        logging.error(f"Error occurred during command execution: {e}")
        error.append(e)
    return error


def solve_network(scenario, horizon, use_legacy=False):
    # change path to pypsa-eur
    change_path_to_pypsa_eur()

    # solve target
    target = solve_target(scenario, horizon)

    # solve the network (optional legacy defaults, then scenario overrides)
    command = (
        f"snakemake -call {target} "
        f"--configfile {get_configfiles(scenario, horizon, use_legacy=use_legacy)} "
        f"--force --rerun-incomplete"
    )
    result = subprocess.run(command, shell=True)

    # move to base directory
    change_path_to_base()

    if result.returncode != 0:
        raise RuntimeError(f"Solve failed for {scenario} {horizon}")
    logging.info(f"Solved network for {scenario} {horizon}")


def get_heat_saved(scenario, horizon):
    # change path to pypsa-eur
    change_path_to_pypsa_eur()

    # load the network
    n = load_network(scenario, horizon)

    # move to base directory
    change_path_to_base()

    if n is None:
        raise FileNotFoundError(f"No solved network for {scenario} {horizon}")

    # calculate saved heat ratio
    retrofitting = (
        n.generators_t.p.filter(like="retrofitting")
        .multiply(n.snapshot_weightings.objective, axis=0).sum().sum()
    )
    heat_demand = (
        n.loads_t.p_set.filter(like="heat")
        .multiply(n.snapshot_weightings.objective, axis=0).sum().sum()
    )
    heat_saved_ratio = retrofitting / heat_demand
    return heat_saved_ratio


def calculate_sink_T(heat_saved_ratio):
    # calculate heat_pump_sink_T
    sink_T = (55-21)*(1-heat_saved_ratio) + 21
    return sink_T


def delete_config_yaml():
    # change path to pypsa-eur
    change_path_to_pypsa_eur()

    # Define the path to the config.yaml file
    config_file_path = os.path.join('config', 'config.yaml')

    # Check if the file exists before attempting to delete it
    if os.path.exists(config_file_path):
        os.remove(config_file_path)
        print(f"Deleted {config_file_path}")
    else:
        print(f"{config_file_path} does not exist") 

    # log changes
    logging.info("config/config.yaml was deleted.")
    
    # move to base directory
    change_path_to_base()
    

def update_sink_T(scenario, horizon, sink_T):
    # Define file path
    key = "heat_pump_sink_T_individual_heating"
    config_path = get_config_path(scenario, horizon)  # absolute path

    # Define the line to be set
    new_line = f'  {key}: {sink_T:.1f}\n'

    # Read the contents of the file
    with open(config_path, 'r') as file:
        lines = file.readlines()

    # Find the index of the line containing the specified text
    index = next((i for i, line in enumerate(lines) if f'  {key}:' in line), None)

    # Insert the new line after the specified line
    if index is not None:
        lines[index] = new_line

    # Write the modified contents back to the file
    with open(config_path, 'w') as file:
        file.writelines(lines)

    # log changes
    logging.info(f"Changed {key} to {sink_T}")


def read_sink_T(scenario, horizon):
    # Define file path
    key = "heat_pump_sink_T_individual_heating"
    config_path = get_config_path(scenario, horizon)  # absolute path

    # Read the contents of the file
    with open(config_path, 'r') as file:
        lines = file.readlines()

    # Find the index of the line containing the specified text
    index = next((i for i, line in enumerate(lines) if f'  {key}:' in line), None)

    # Insert the new line after the specified line
    if index is not None:
        # Extract the line
        line = lines[index]
        
        # Split the line by the specified text and take the part after it
        value_str = line.split(f'  {key}:')[-1].strip()

        # Convert the extracted part to float
        value = float(value_str)
        
        logging.info(f"The extracted {key} value is: {value:.1f}")
    else:
        logging.info(f"The specified {key} was not found in any line.")

    return value


def run_workflow(scenario, horizon, improved_cop=False, use_legacy=False):
    compose_network(scenario=scenario, horizon=horizon, use_legacy=use_legacy)

    # initialize error_capacities and error_moderate
    error_capacities, error_moderate, error_improve_cop = [], [], []
    
    # set capacities if 2040 or 2050
    if not horizon == 2030:
        error_capacities = set_capacities(scenario=scenario, horizon=horizon)

    # set moderate retrofitting
    if scenario == "flexible-moderate":
        error_moderate = moderate_retrofitting(scenario=scenario, horizon=horizon)

    # set retrofitting p_nom for improved COP
    if improved_cop:
        error_improve_cop = improve_cops_after_renovation(scenario=scenario, horizon=horizon)

    # break if error happens
    if error_capacities or error_moderate or error_improve_cop:
        return None

    # solve the network
    solve_network(scenario, horizon, use_legacy=use_legacy)

    return True


if __name__ == "__main__":
    # get scenario from argument
    scenarios, horizons, improved_cop, use_legacy = get_scenario()

    # copy custom data into pypsa-eur/data folder
    copy_custom_data()

    # remove BAU scenario from scenarios list if present (BAU is simulated separately)
    scenario_BAU = "BAU" in scenarios
    if scenario_BAU:
        scenarios.remove("BAU")

    # run model for given horizon
    for horizon in horizons:
        for scenario in scenarios:
            # calculate heat_pump_sink_T for flexible-moderate, set to 55.0 at the beginning of first run for other scenarios
            if scenario == "flexible-moderate" and improved_cop:
                # read heat saved from flexible scenario
                heat_saved_ratio = get_heat_saved("flexible", horizon)
                # calculate sink_T for half of heat saved
                sink_T = calculate_sink_T(heat_saved_ratio/2)
                # update sink_T
                update_sink_T(scenario, horizon, sink_T)
            elif scenario in ["flexible", "flexible-moderate", "retro_tes"]:
                update_sink_T(scenario, horizon, 55.0)

            # run full network preparation and solving workflow 
            run_status = run_workflow(scenario, horizon, use_legacy=use_legacy)

            # stop further execution if workflow did not succeed
            if run_status is None:
                logging.error(f"Workflow of {scenario} scenario for {horizon} broke!")
                continue

            # for Optimal and Limited retrofitting proceed with improved COP
            if scenario in ["flexible", "retro_tes"] and improved_cop:
                # read heat saved
                heat_saved_ratio = get_heat_saved(scenario, horizon)

                # calculate heat_pump_sink_T
                sink_T = calculate_sink_T(heat_saved_ratio)

                # delete config.yaml
                delete_config_yaml()

                # update heat_pump_sink_T
                update_sink_T(scenario, horizon, sink_T)

                # run full network preparation and solving workflow
                run_status = run_workflow(
                    scenario, horizon, improved_cop=improved_cop, use_legacy=use_legacy
                )


    # run BAU scenario
    if scenario_BAU:
        compose_network("BAU", 2020, use_legacy=use_legacy)
        solve_network("BAU", 2020, use_legacy=use_legacy)
