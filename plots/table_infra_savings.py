import os
import sys
sys.path.append("../submodules/pypsa-eur")
import pypsa
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
from _helpers import mock_snakemake, update_config_from_wildcards, load_network, \
                     change_path_to_pypsa_eur, change_path_to_base


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "get_infra_savings", 
            clusters="48",
        )
    # update config based on wildcards
    config = update_config_from_wildcards(snakemake.config, snakemake.wildcards)


    # move to submodules/pypsa-eur
    change_path_to_pypsa_eur()
    # network parameters
    co2l_limits = {"2030":"0.45", "2040":"0.1", "2050":"0.0"}
    line_limits = {"2030":"v1.15", "2040":"v1.3", "2050":"v1.5"}
    clusters = config["plotting"]["clusters"]
    time_resolution = config["plotting"]["time_resolution"]

    # define scenario namings
    scenarios = {"flexible": "Optimal Renovation and Heating", 
                 "retro_tes": "Optimal Renovation and Green Heating", 
                 "flexible-moderate": "Limited Renovation and Optimal Heating", 
                 }

    # define dataframe to store infra savings
    df_savings = pd.DataFrame(
        index=list(scenarios.values()),
        columns=[
            ("2030", "wind"), ("2030", "solar"), ("2030", "gas"),
            ("2040", "wind"), ("2040", "solar"), ("2040", "gas"),
            ("2050", "wind"), ("2050", "solar"), ("2050", "gas")
        ]
    )
    df_savings.columns = pd.MultiIndex.from_tuples(df_savings.columns, names=['horizon','tech'])

    for planning_horizon in ["2030", "2040", "2050"]:
        lineex = line_limits[planning_horizon]
        sector_opts = f"Co2L{co2l_limits[planning_horizon]}-{time_resolution}-T-H-B-I"

        # benchmark network
        b = load_network(lineex, clusters, sector_opts, planning_horizon, "rigid")
        for scenario, nice_name in scenarios.items():
            # load networks
            n = load_network(lineex, clusters, sector_opts, planning_horizon, scenario)

            if n is None:
                # Skip further computation for this scenario if network is not loaded
                print(f"Network is not found for scenario '{scenario}', planning year '{planning_horizon}', and time resolution of '{time_resolution}'. Skipping...")
                continue

            # estimate upper and lower limits of congestion of grid
            solar_carriers = ["solar", "solar rooftop"]
            solar = (
                n.generators.query("carrier in @solar_carriers").p_nom_opt.sum() -
                b.generators.query("carrier in @solar_carriers").p_nom_opt.sum()
            )/1e3
            wind_carriers = ["onwind", "offwind-ac", "offwind-dc"]
            wind = (
                n.generators.query("carrier in @wind_carriers").p_nom_opt.sum() -
                b.generators.query("carrier in @wind_carriers").p_nom_opt.sum()
            )/1e3
            OCGT_carriers = ["OCGT"]
            gas = (
                n.links.query("carrier in @OCGT_carriers").p_nom_opt.multiply(n.links.efficiency).sum() -
                b.links.query("carrier in @OCGT_carriers").p_nom_opt.multiply(b.links.efficiency).sum()
            )/1e3

            df_savings.loc[nice_name, (planning_horizon, "solar")] = solar
            df_savings.loc[nice_name, (planning_horizon, "wind")] = wind
            df_savings.loc[nice_name, (planning_horizon, "gas")] = gas

    # add name for columns
    df_savings.index.name = "Scenario [GW]"

    # move to base directory
    change_path_to_base()

    # save the heat pumps data in Excel format
    df_savings.to_csv(snakemake.output.table)

