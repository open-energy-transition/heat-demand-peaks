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

The project, commissioned by the EEE consortium  and supported by Open Energy Transition (OET), aims to assess the impact of various energy efficiency measures on the European energy system. This study focuses on energy affordability, social and household impacts, and industry considerations in the context of the EU's target of 90-95% emissions reduction by 2040. Utilizing the PyPSA-Eur integrated energy system planning tool, the project evaluates the effectiveness of different renovation scenarios, energy management measures, and demand-side flexibility measures in reducing energy generation needs, flattening the peak demand curve, and influencing energy prices. The results provide insights into the benefits of isolated and combined efficiency measures, contributing to data-driven decision-making in energy policy and planning. The project commenced in January 2024 and is expected to conclude by July 2024.

The study is developed as an extension to the PyPSA-Eur model in this repository.

# Repository Structure

- `submodules` contains the relevant submodules `PyPSA-Eur` and `technology-data`.
- `configs` contains all relevant configuration files. in `EEE_study`, there are all relevant config files for this study, while `Zeyen_etal` contains config files to reproduce [previously published](https://doi.org/10.1016/j.energy.2021.120784) results; `config.plot.yaml` contains network parameters for plotting and network modifications.
- `plots` contains some scripts to generate a few plots for model evaluation purposes
- `scripts` contains scripts that modify existing results in order to create new scenarios

# Installation and Usage

### 1. Installation

Clone the repository including its submodules:

    git clone --recurse-submodules https://github.com/open-energy-transition/heat-demand-peaks

Install the necessary dependencies using `conda` or `mamba`:

    mamba env create -f submodules/pypsa-eur/envs/environment.yaml

Navigate into the main Snakemake workflow directory of `PyPSA-Eur`:

    cd submodules/pypsa-eur

### 2. Running scenarios

To run the scenarios of a particular configuration file (e.g. `configs/EEE_study/config.flexible-industry.yaml`), run:

    snakemake -call solve_sector_networks --configfile ../../configs/EEE_study/config.flexible-industry_2030.yaml 

This call requires a high-performance computing environment, as well as a [Gurobi license](https://www.gurobi.com/downloads/gurobi-software/).

Please follow the documentation of PyPSA-Eur for more details.

### 3. Setting nominal capacities of retrofitting for `Limited Renovation & Optimal Heating` scenario

The nominal capacities of retrofitting for `Limited Renovation & Optimal Heating` (moderate retrofitting) scenario is set by running:

    snakemake -call set_moderate_retrofitting

* Note! This and the following `snakemake` commands must be run in `heat-demands-peak` base directory (not `pypsa-eur` submodule).

This command will set `p_nom` for moderate retrofitting network as a half of `p_nom_opt` of solved `Optimal Renovation and Heating` (flexible) scenario. The network parameters, such as `clusters`, `planning_horizon`, and `time_resolution`, are defined in `moderate_retrofitting` section of `configs/config.plot.yaml`. 

As an alternative, the nominal capacities for moderate retrofitting scenario can be set by running the command with wildcards (e.g. scenario for 2030 with 48 clusters):

    snakemake -call scripts/logs/set_moderate_retrofitting_48_2030.txt --force

The resultant file in `scripts/logs/set_moderate_retrofitting_48_2030.txt` contains `Success` parameter which indicates the success status of executed rule. Here, `--force` flag is used to forcefully re-execute the rule.

### 4. Transfering optimal capacities to future horizons

After optimizing scenarios for one horizon, it is important to transfer optimal generation and store capacities into future horizons. To do so, configure `planning_horizon` in `set_capacities` section of `configs/config.plot.yaml`. By default, `2040` is set as `planning_horizon` in `set_capacities`, which helps to transfer `p_nom_opt` values from solved networks of 2030 into `p_nom_min` of corresponding generators and stores of unsolved network of 2040. The optimal capacities are transfered to corresponding scenarios. To set `p_nom_min` for all scenarios of 2040, run:

    snakemake -call set_capacities

To set minimum capacities for specific scenario (e.g. flexible scenario of 2050 with 48 clusters), you can run run:

    snakemake -call scripts/logs/set_capacities_48_2050_flexible.txt

The resultant file in `scripts/logs/set_capacities_48_2050_flexible.txt` contains `Success` parameter which indicates the success status of executed rule.

### 5. Plotting

To plot the total system cost for all planning horizons (i.e. 2030, 2040, and 2050) specified in `config.plot.yaml`, run:

    snakemake -call plot_total_costs

**Note!** The total costs will be plotted only for scenarios where the solved networks are present.

To plot total costs for specific horizon (e.g. 2040) and number of clusters (e.g. 48), run:

    snakemake -call plots/results/plot_total_costs_48_2040.png

To plot electricity bills per household and electricity prices per country for all horizons (i.e. 2030, 2040, and 2050), run:

    snakemake -call plot_electricity_bills

**Note!** The electricity bills and prices will be computed only for scenarios with solved networks.

To plot electricity bills separately for specific horizon (e.g. 2040) and number of clusters (e.g. 48), run:

    snakemake -call plots/results/plot_bill_per_household_48_2040.png

To estimate number of heat pumps and generate table for all horizons and clusters specified in `config.plot.yaml`, run:

    snakemake -call get_heat_pumps

To estimate average grid congestion and generate table for all horizons and clusters specified in `config.plot.yaml`, run:

    snakemake -call get_line_congestions

To generate all aforementioned plots and tables at once, run:

    snakemake -call plot_all --forceall

**Note!** --forceall flag is used to force a rerun of rule.
