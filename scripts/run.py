import subprocess
import argparse
import logging
import yaml
import sys
import os
import pypsa
sys.path.append("plots")
from _helpers import change_path_to_pypsa_eur, change_path_to_base, load_network

# Set up logging configuration
logging.basicConfig(level=logging.INFO)

def get_scenario():
    parser = argparse.ArgumentParser(description="Running the scenario")
    parser.add_argument("-s", "--scenario", help="Specify the scenario.", required=True, 
                        choices=["flexible", "flexible-moderate", "retro_tes", "rigid"])
    parser.add_argument("-c", "--continue_horizon", help="Specify the horizon to continue simulations", 
                        choices=["2030", "2040", "2050"])
    parser.add_argument("-y", "--year", help="Specify a single horizon to simulate", 
                        choices=["2030", "2040", "2050"])
    args = parser.parse_args()

    # Access the value of the scenario argument
    scenario = args.scenario
    # log scenario name
    logging.info(f"Scenario: {scenario}")

    # Access the value of horizon from which simulation is continued
    c = args.continue_horizon
    # Access the value of specific horizon to simulate
    y = args.year

    # log scenario name
    if y:
        horizons = [int(y)]
        logging.info(f"Start simulating {scenario} scenario for {y}")
    elif c:
        horizons = get_horizon_list(int(c))
        logging.info(f"Start simulating {scenario} scenario for {horizons}")
    else:
        c = 2030
        horizons = get_horizon_list(int(c))
        logging.info("No horizon specified. Starting from default horizon (2030)")

    return scenario, horizons


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
    command = f"snakemake -call results/{scenario}/prenetworks/{filename} --configfile {config_path} --force --ri"
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


def get_heat_saved(scenario, horizon):
    # change path to pypsa-eur
    change_path_to_pypsa_eur()

    # get .nc filename
    filename = get_network_name(scenario, horizon)

    # load solved network
    n = None
    try:
        n = pypsa.Network(os.path.join(f"results/{scenario}/postnetworks", filename))
        logging.info(f"Loading {filename} for {scenario}")
    except FileNotFoundError as e:
        print(f"Error: {e}")

    # move to base directory
    change_path_to_base()

    # calculate saved heat ratio
    retrofitting = n.generators_t.p.filter(like="retrofitting").multiply(n.snapshot_weightings.objective, axis=0).sum().sum()
    heat_demand = n.loads_t.p_set.filter(like="heat").multiply(n.snapshot_weightings.objective, axis=0).sum().sum()
    heat_saved_ratio = retrofitting / heat_demand
    return heat_saved_ratio


def calculate_sink_T(heat_saved_ratio):
    # calculate heat_pump_sink_T
    sink_T = (55-21)*heat_saved_ratio + 21
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
    # change path to pypsa-eur
    change_path_to_pypsa_eur()

    # Define the file path
    config_path = get_config_path(scenario, horizon)
    config_path = "../../" + config_path

    # Define the line to be set
    new_line = f'  heat_pump_sink_T: {sink_T:.2f}\n'

    # Read the contents of the file
    with open(config_path, 'r') as file:
        lines = file.readlines()

    # Find the index of the line containing the specified text
    index = next((i for i, line in enumerate(lines) if '  heat_pump_sink_T:' in line), None)

    # Insert the new line after the specified line
    if index is not None:
        lines[index] = new_line

    # Write the modified contents back to the file
    with open(config_path, 'w') as file:
        file.writelines(lines)

    # log changes
    logging.info(f"Changed heat_pump_sink_T to {sink_T}")

    # move to base directory
    change_path_to_base()


def run_workflow(scenario, horizon):
    # increase biomass potential for 2050 by 1.2
    if horizon == 2050:
        increase_biomass_potential()
    # run prenetwork
    prepare_prenetwork(scenario=scenario, horizon=horizon)

    # initialize error_capacities and error_moderate
    error_capacities, error_moderate = [], []
    
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
        return None # return None is error happens

    # solve the network
    solve_network(scenario, horizon)

    # revert biomass potential
    if horizon == 2050:
        revert_biomass_potential()

    return True # return True if success


if __name__ == "__main__":
    # get scenario from argument
    scenario, horizons = get_scenario()

    # run model for given horizon
    for horizon in horizons:
        # run full network preparation and solving workflow 
        run_status = run_workflow(scenario, horizon)

        # stop further execution if workflow did not succeed
        if run_workflow is None:
            logging.error("Workflow broke!")
            break

        # for Optimal and Limited retrofitting proceed with improved COP
        if scenario in ["flexible", "flexible-moderate"]:
            # read heat saved
            heat_saved_ratio = get_heat_saved(scenario, horizon)

            # calculate heat_pump_sink_T
            sink_T = calculate_sink_T(heat_saved_ratio)

            # delete config.yaml
            delete_config_yaml()

            # update heat_pump_sink_T
            update_sink_T(scenario, horizon, sink_T)

            # run full network preparation and solving workflow
            run_status = run_workflow(scenario, horizon)

            # revert heat_pump_sink_T to 55.0
            # update_sink_T(scenario, horizon, 55.0)
