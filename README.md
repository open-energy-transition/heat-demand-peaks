<!--
SPDX-FileCopyrightText:  Open Energy Transition gGmbH

SPDX-License-Identifier: AGPL-3.0-or-later
-->

# heat-demand-peaks


<img src="https://github.com/open-energy-transition/heat-demand-peaks/assets/53824825/d305974e-030f-46f9-8b7c-9d7d276a81a3" alt="Open Energy Transition Logo" width="260" height="100" align="right">
<img src="https://github.com/open-energy-transition/heat-demand-peaks/assets/53824825/4470a2e7-42eb-4f97-81e6-7d58de1042e5" alt="Eurima Logo" width="240" height="90">
<img src="https://github.com/open-energy-transition/heat-demand-peaks/assets/53824825/bcf1713b-39e7-463e-a4c5-99f7310e140a" alt="ECF Logo" width="200" height="100">
<img src="https://github.com/open-energy-transition/heat-demand-peaks/assets/53824825/b2a6782b-28db-4ae2-be26-22e97e058b14" alt="ECI Logo" width="260" height="85">
<br>
<br>

[![REUSE status](https://api.reuse.software/badge/git.fsfe.org/reuse/api)](https://api.reuse.software/info/git.fsfe.org/reuse/api)
[![REUSE Compliance Check](https://github.com/yerbol-akhmetov/heat-demand-peaks/actions/workflows/reuse-compliance.yml/badge.svg)](https://github.com/yerbol-akhmetov/heat-demand-peaks/actions/workflows/reuse-compliance.yml)

The project, commissioned by the EEE consortium  and supported by Open Energy Transition (OET), aims to assess the impact of various energy efficiency measures on the European energy system. This study focuses on energy affordability, social and household impacts, and industry considerations in the context of the EU's target of 90-95% emissions reduction by 2040. Utilizing the PyPSA-Eur integrated energy system planning tool, the project evaluates the effectiveness of different renovation scenarios, energy management measures, and demand-side flexibility measures in reducing energy generation needs, flattening the peak demand curve, and influencing energy prices. The results provide insights into the benefits of isolated and combined efficiency measures, contributing to data-driven decision-making in energy policy and planning. The project commenced in January 2024 and is expected to conclude by July 2024.

The study is developed as an extension to the PyPSA-Eur model in this repository.

# Repository Structure

- `submodules` contains the relevant submodules `PyPSA-Eur` and `technology-data`.
- `configs` contains all relevant configuration files. in `EEE_study`, there are all relevant config files for this study, while `Zeyen_etal` contains config files to reproduce [previously published](https://doi.org/10.1016/j.energy.2021.120784) results; `config.plot.yaml` contains network parameters for plotting and network modifications.
- `plots` contains some scripts to generate a few plots for model evaluation purposes
- `scripts` contains scripts that modify existing results in order to create new scenarios

# Installation and Usage

## 1. Installation

Clone the repository including its submodules:

    git clone --recurse-submodules https://github.com/open-energy-transition/heat-demand-peaks

Install the necessary dependencies using `conda` or `mamba`:

    mamba env create -f submodules/pypsa-eur/envs/environment.yaml

Activate `pypsa-eur` environment:

    conda activate pypsa-eur

Navigate into the main Snakemake workflow directory of `PyPSA-Eur`:

    cd submodules/pypsa-eur

## 2. Running scenarios

Before running the scenarios, it is important to know their short working names used in the code. The table below provides the mapping between official scenario names and correspoding short names utilized in the code. To run the scenarios, refer to the coding names.

|Scenario name                                |Coding name      |
|---------------------------------------------|-----------------|
|Optimal Retrofitting & Optimal Heating (OROH)|flexible         |
|Optimal Retrofitting & Green Heating (ORGH)  |retro_tes        |
|Limited Retrofitting & Optimal Heating (LROH)|flexible-moderate|
|No Retrofitting & Green Heating (NRGH)       |rigid            |

**Note!** Running the scenarios requires a high-performance computing environment, as well as a [Gurobi license](https://www.gurobi.com/downloads/gurobi-software/).

### A. Running scenarios using *automated workflow* (the easy way)

To run the simulations for specific scenario for several horizons at once (e.g. *Optimal Renovation & Optimal Heating (OROH)* scenario), run:

    python scripts/run.py -s flexible -c 2040

where `-s` is a mandatory flag used for the scenario selection. Use the coding name of corresponding scenario from the table to trigger the execution. `-c` flag (optional) speficies from which planning horizon the simulations should start. So if `-c 2040` is given, then 2040 and 2050 horizons will be simulated consequitively; if not specified, then 2030 is used as a starting year by default. `-c` flag is useful when simulation is interupted and it needs to be re-run from certain horizon. 

To simulate a single horizon for the scenario (e.g. *Limited Retrofitting & Optimal Heating* scenario for 2050), use `-y` flag as follows:

    python scripts/run.py -s flexible-moderate -y 2050

* **Note!** To run *Limited Renovation & Optimal Heating* (`flexible-moderate`) scenario, the simulation results for *Optimal Renovation & Optimal Heating* (`flexible`) scenario must be present.

The available flags and their details are presented in table below:

|Flag                  |Default   |Description        |Status    |
|----------------------|----------|-------------------|----------|
|`-s`, `--scenario`    |n/a       |Selects scenario   |required  |
|`-c`, `--continue`    |2030      |Selects horizon from which simulation starts (used in simulating multiple horizons). If 2030 is selected, then 2030, 2040, and 2050 is simulated consecutively. If 2040 is selected, then 2040 and 2050 is simulated.|optional  |
|`-y`, `--year`        |n/a       |Selects a single horizon to be simulated. If both `-c` and `-y` are prodived occasionally, then priority is given to `-y`.|optional  |
|`-i`, `--improved_cop`|true      |Enables/disables improved COP workflow. By default, improved COP for heat pumps is used.|optional  |

### 2.2. Running scenarios *manually* using `snakemake` (the hard way)

#### A. Run the scenario

To run the scenario of a particular configuration file (e.g. `configs/EEE_study/config.flexible-industry.yaml`), navigate to `pypsa-eur` directory using:

    cd submodules/pypsa-eur

Then, run the following comamnd to prepare the pre-network:

    snakemake -call prepare_sector_networks --configfile ../../configs/EEE_study/config.flexible-industry_2030.yaml

To solve the network, run:

    snakemake -call solve_sector_networks --configfile ../../configs/EEE_study/config.flexible-industry_2030.yaml 

Please follow the documentation of PyPSA-Eur for more details.

#### B. Setting nominal capacities of retrofitting for *Limited Renovation & Optimal Heating (LROH)* scenario

The nominal capacities of retrofitting for *Limited Renovation & Optimal Heating (LROH)* (flexible-moderate) scenario is set by running:

    snakemake -call set_moderate_retrofitting

* **Note!** This and the following `snakemake` commands must be run in `heat-demands-peak` base directory (not `pypsa-eur` submodule). The command needs to be run after preparation of pre-network.

This command will set `p_nom` for moderate retrofitting network as a half of `p_nom_opt` of solved *Optimal Renovation and Heating* (flexible) scenario. The network parameters, such as `clusters` and `planning_horizon`, are defined in `moderate_retrofitting` section of `configs/config.plot.yaml`. 

As an alternative, the nominal capacities for moderate retrofitting scenario can be set by running the command with wildcards (e.g. scenario for 2030 with 48 clusters):

    snakemake -call scripts/logs/set_moderate_retrofitting_48_2030.txt --force

The resultant file in `scripts/logs/set_moderate_retrofitting_48_2030.txt` contains `Success` parameter which indicates the success status of executed rule. Here, `--force` flag is used to forcefully re-execute the rule.

#### C. Transfering optimal capacities to future horizons

After optimizing scenarios for one horizon, it is important to transfer optimal generation and store capacities into future horizons. To do so, configure `planning_horizon` in `set_capacities` section of `configs/config.plot.yaml` to horizon of interest. By default, `2040` is set as `planning_horizon` in `set_capacities`, which helps to transfer `p_nom_opt` values from solved networks of 2030 into `p_nom_min` of corresponding generators and stores of unsolved network of 2040. The optimal capacities are transfered to corresponding scenarios. To set `p_nom_min` for all scenarios of 2040, run:

    snakemake -call set_capacities

* **Note!** The command needs to be run after preparation of pre-network, but before solving it.

To set minimum capacities for specific scenario (e.g. flexible scenario of 2050 with 48 clusters), you can run run:

    snakemake -call scripts/logs/set_capacities_48_2050_flexible.txt

The resultant file in `scripts/logs/set_capacities_48_2050_flexible.txt` contains `Success` parameter which indicates the success status of executed rule.

#### D. Utilizing improved COP for heat pumps in scenarios with building retrofitting

To use the workflow of improved COP, solve the network regularly and determine the heat saved ratio by building retrofitting from the solved network. Use the ratio to compute sink temperature:

    heat_pump_sink_T = (55 - 21) * (1 - heat_saved_ratio) + 21

Update the `heat_pump_sink_T` parameter in corresponding configuration file. Then, prepare the pre-network and set capacities from previous horizon. To set retrofitting capacities from the first run to current run, execute:

    snakemake -call improve_cops_after_renovation

Finally, the updated pre-network can be solved.

### 3. Plotting

To visualize the simulation results, the following list of `snakemake` rules can be used. 

|Rule name                    |Description        |
|-----------------------------|-------------------|
|`plot_total_costs`           |Plots total costs and optimal capacities of technologies and provides the results in a tabular form|
|`plot_electricity_bills`     |Plots electricity bills per household and energy price per MWh|
|`plot_electricity_for_heats` |Plots electricity and gas consumption profiles to cover heat demands|
|`plot_electricity_generations`|Plots electricity generation mix|
|`plot_curtailments`          |Plots total curtailment and provides detailed curtailment for renewables in a tabular form|
|`plot_hydrogen_productions`  |Plots hydrogen generation and electricity used for its production and provides the results in a tabular form|
|`plot_heat_tech_ratios`      |Plots optimal capacities (GW<sub>el</sub>) of heat generation technologies, such as gas boilers, heat pumps, and resistive heaters|
|`plot_COP`                   |Plots coefficient of performance (COP) for heat pumps|
|`plot_co2_levels`            |Plots CO<sub>2</sub> balances and provides the results in a tabular form|
|`plot_historic_generation`   |Plots historic electricity generation mix for 2022 from [1]|
|`plot_transmission_congestions`|Plots transmission lines congestion|
|`plot_all`                   |Plots all abovementioned figures|
|`get_heat_pumps`             |Provodes estimated heat pumps amount in a tabualr form|
|`get_infra_savings`          |Provides infrastructure savings interms of costs and capacities in a tabular form|

[1] How is EU electricity produced and sold? URL: https://www.consilium.europa.eu/en/infographics/how-is-eu-electricity-produced-and-sold/

To demonstrate how to use the plotting rules to visualize the total system cost for various planning horizons (e.g., 2030, 2040, and 2050) as specified in the `config.plot.yaml` file, run the following command:

    snakemake -call plot_total_costs

**Note!** The total costs will be plotted only for scenarios where the solved networks are present. `--forceall` flag can be used to regenerate and overwrite the already existing plots.

As another example, to estimate number of heat pumps and generate table for all horizons and clusters specified in `config.plot.yaml`, run:

    snakemake -call get_heat_pumps

To generate all aforementioned plots and tables at once, run:

    snakemake -call plot_all --forceall

**Note!** `--forceall` flag is used to force a rerun of rule.
