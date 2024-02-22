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

- `workflows` contains the relevant submodules `PyPSA-Eur` and `technology-data`.
- `configs` contains all relevant configuration files. in `EEE_study`, there are all relevant config files for this study, while `Zeyen_etal` contains config files to reproduce [previously published](https://doi.org/10.1016/j.energy.2021.120784) results
- `plots` contains some scripts to generate a few plots for model evaluation purposes
- `scripts` contains scripts that modify existing results in order to create new scenarios

# Installation and Usage

Clone the repository including its submodules:

`git clone --recurse-submodules https://github.com/open-energy-transition/heat-demand-peaks`

Install the necessary dependencies using `conda` or `mamba`:

`mamba env create -f workflows/pypsa-eur/envs/environment.yaml`

Navigate into the main Snakemake workflow directory of `PyPSA-Eur`:

`cd workflows/pypsa-eur`

To run the scenarios of a particular configuration file (e.g. `configs/EEE_study/config.flexible-industry.yaml`), run:

`snakemake -call --configfile ../../configs/EEE_study/config.flexible-industry.yaml solve_sector_networks`

This call requires a high-performance computing environment, as well as a [Gurobi license](https://www.gurobi.com/downloads/gurobi-software/).

Please follow the documentation of PyPSA-Eur for more details.
