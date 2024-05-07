import subprocess
import argparse
import logging
import yaml
import sys
sys.path.append("plots")
from _helpers import change_path_to_pypsa_eur, change_path_to_base

# Set up logging configuration
logging.basicConfig(level=logging.INFO)

def get_scenario():
    parser = argparse.ArgumentParser(description="Running the scenario")
    parser.add_argument("-s", "--scenario", help="Specify the scenario.", required=True, 
                        choices=["flexible", "flexible-moderate", "retro_tes", "rigid"])
    parser.add_argument("-c", "--continue_horizon", help="Specify the horizon to continue simulations", 
                        choices=["2030", "2040", "2050"])
    args = parser.parse_args()

    # Access the value of the scenario argument
    scenario = args.scenario
    # log scenario name
    logging.info(f"Scenario: {scenario}")

    # Access the value of the horizon argument
    c = args.continue_horizon
    # log scenario name
    if c:
        logging.info(f"Start simulating from {c}")
    else:
        c = 2030
        logging.info("No horizon specified. Starting from default horizon (2030)")

    return scenario, int(c)


def get_horizon_list(start_horizon):
    horizons = [2030, 2040, 2050]
    try:
        start_index = horizons.index(start_horizon)
        return horizons[start_index:]
    except ValueError:
        return []


def get_clusters(scenario, horizon):
    config_path = get_config_path(scenario, horizon)
    # read config file
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    return config["scenario"]["clusters"][0]


def get_config_path(scenario, horizon):
    configname_dict = {"flexible": "flexible-industry",
                       "flexible-moderate": "flexible-moderate",
                       "retro_tes": "retro_tes-industry",
                       "rigid": "rigid-industry"}
    configname = f"config.{configname_dict[scenario]}_{horizon}.yaml"
    configpath = "configs/EEE_study/"
    return configpath + configname


def get_network_name(scenario, horizon):
    config_path = get_config_path(scenario, horizon)
    config_path = "../../" + config_path
    # read config file
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    # read config.default.yaml file
    with open("config/config.default.yaml", 'r') as file:
        config_default = yaml.safe_load(file)
    # scenario params
    params = config["scenario"]
    params_default = config_default["scenario"]
    ll, clusters, sector_opts, planning_horizons = params['ll'][0], params['clusters'][0], params['sector_opts'][0], params['planning_horizons'][0]
    simpl, opts = params_default['simpl'][0], params_default['opts'][0]
    filename = f"elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc"
    return filename


def increase_biomass_potential(factor=1.2):
    # change path to pypsa-eur
    change_path_to_pypsa_eur()

    # Define the file path
    file_path = 'scripts/prepare_sector_network.py'

    # Define the line to be added
    new_line = f'    biomass_potentials = {factor} * biomass_potentials\n'

    # Read the contents of the file
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Find the index of the line containing the specified text
    index = next((i for i, line in enumerate(lines) if 'biomass_potentials = pd.read_csv(snakemake.input.biomass_potentials, index_col=0)' in line), None)

    # Insert the new line after the specified line
    if index is not None:
        lines.insert(index + 1, new_line)

    # Write the modified contents back to the file
    with open(file_path, 'w') as file:
        file.writelines(lines)

    # log changes
    logging.info(f"Increase biomass potentials by {factor} factor")

    # move to base directory
    change_path_to_base()


def revert_biomass_potential():
    # Change path to pypsa-eur
    change_path_to_pypsa_eur()

    # Define the file path
    file_path = 'scripts/prepare_sector_network.py'

    # Read the contents of the file
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Find the index of the line containing the specified text
    index = next((i for i, line in enumerate(lines) if 'biomass_potentials = 1.2 * biomass_potentials' in line), None)

    # Remove the line if found
    if index is not None:
        del lines[index]

    # Write the modified contents back to the file
    with open(file_path, 'w') as file:
        file.writelines(lines)

    # log changes
    logging.info(f"Revert biomass potentials back")

    # Move to base directory
    change_path_to_base()


def prepare_prenetwork(scenario, horizon):
    # get config path
    config_path = get_config_path(scenario, horizon)
    # config path relative to pypsa-eur folder
    config_path = "../../" + config_path

    # change path to pypsa-eur
    change_path_to_pypsa_eur()

    # get .nc filename
    filename = get_network_name(scenario, horizon)
    # run prenetwork
    command = f"snakemake -call results/{scenario}/prenetworks/{filename} --configfile {config_path} --force"
    subprocess.run(command, shell=True)
    logging.info(f"Prenetwork was prepared for {scenario} scenario in {horizon} horizon!")

    # move to base directory
    change_path_to_base()


def set_capacities(scenario, horizon):
    error = []
    try:
        # get number of clusters
        clusters = get_clusters(scenario, horizon)
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
        clusters = get_clusters(scenario, horizon)
        command = f"snakemake -call scripts/logs/set_moderate_retrofitting_{clusters}_{horizon}.txt --forceall"
        subprocess.run(command, shell=True, check=True)
        logging.info(f"Moderate retrofitting capacities are set to {scenario} scenario in {horizon} horizon!")
    except subprocess.CalledProcessError  as e:
        logging.error(f"Error occurred during command execution: {e}")
        error.append(e)
    return error


def solve_network(scenario, horizon):
    # get config path
    config_path = get_config_path(scenario, horizon)
    # config path relative to pypsa-eur folder
    config_path = "../../" + config_path

    # change path to pypsa-eur
    change_path_to_pypsa_eur()

    # solve the network
    command = f"snakemake -call solve_sector_networks --configfile {config_path}"
    subprocess.run(command, shell=True)
    logging.info(f"Network was solved for {scenario} scenario in {horizon} horizon!")

    # move to base directory
    change_path_to_base()


if __name__ == "__main__":
    # get scenario from argument
    scenario, start_horizon = get_scenario()

    # horizons to be simulated
    horizons = get_horizon_list(start_horizon=start_horizon)

    # run model for given horizon
    for horizon in horizons:
        if horizon == 2050:
            increase_biomass_potential()
        # run prenetwork
        prepare_prenetwork(scenario=scenario, horizon=horizon)
        
        # set capacities if 2040 or 2050
        if not horizon == 2030:
            error_capacities = set_capacities(scenario=scenario, horizon=horizon)

        # set moderate retrofitting
        if scenario == "flexible-moderate":
            error_moderate = moderate_retrofitting(scenario=scenario, horizon=horizon)

        # break if error happens
        if error_capacities or error_moderate:
            if horizon == 2050:
                revert_biomass_potential()
            break

        # solve the network
        solve_network(scenario, horizon)

        # revert biomass potential
        if horizon == 2050:
            revert_biomass_potential()
