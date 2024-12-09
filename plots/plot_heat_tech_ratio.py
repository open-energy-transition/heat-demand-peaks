# SPDX-FileCopyrightText:  Open Energy Transition gGmbH
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import sys
sys.path.append("../submodules/pypsa-eur")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
import colors as c
from _helpers import mock_snakemake, update_config_from_wildcards, load_network, \
                     change_path_to_pypsa_eur, change_path_to_base, \
                     LINE_LIMITS, CO2L_LIMITS, BAU_HORIZON, PATH_PLOTS, replace_multiindex_values


def get_heat_capacities(n, nice_name):
    # geat heat pump capacity [MW_el]
    heat_pump_carriers = [x for x in n.links.carrier.unique() if "heat pump" in x]
    heat_pump_cap = n.links.query("carrier in @heat_pump_carriers").p_nom_opt.sum()
    # get resistive heater capacity [MW_el]
    resistive_heater_carriers = [x for x in n.links.carrier.unique() if "resistive heater" in x]
    resistive_heater_cap = n.links.query("carrier in @resistive_heater_carriers").p_nom_opt.sum()
    # get gas boiler capacity [MW_el]
    gas_boiler_carriers = [x for x in n.links.carrier.unique() if "gas boiler" in x]
    gas_boiler_cap = n.links.query("carrier in @gas_boiler_carriers").p_nom_opt.sum()

    heat_tech_dict = {"heat pump": heat_pump_cap,
                      "resistive heater": resistive_heater_cap,
                      "gas boiler": gas_boiler_cap}
    heat_tech_caps = pd.Series(heat_tech_dict)
    heat_tech_caps.name = nice_name

    return heat_tech_caps


def plot_capacities(capacities_df, clusters, planning_horizon, plot_width=7):
    # MW to GW
    df = capacities_df / 1e3

    # narrow plot if BAU
    if planning_horizon == BAU_HORIZON:
        plot_width = 1.5

    fig, ax = plt.subplots(figsize=(plot_width, 9))

    df.T.plot(
        kind="bar",
        ax=ax,
        stacked=True,
        color=[c.tech_colors[i] for i in df.index],
        zorder=1,
    )

    handles, labels = ax.get_legend_handles_labels()
    handles.reverse()
    labels.reverse()

    plt.xticks(rotation=0, fontsize=10)
    ax.set_ylabel("Capacity [$\mathrm{GW_{el}}$]")
    ax.set_xlabel("")
    ax.set_ylim([0, 2000])
    x_ticks = list(df.columns)
    if planning_horizon in ["2040", "2050"] and "Limited\nRenovation" in x_ticks:
        # replace name for Limited Renovation scenario for 2030 to be LROH
        x_ticks[x_ticks.index("Limited\nRenovation")] = "Limited\nRenovation &\nElectrification"

    ax.set_xticklabels(x_ticks)

    # Turn off both horizontal and vertical grid lines
    ax.grid(False, which='both')
    ax.legend(
        handles, labels, ncol=1, loc="upper left", bbox_to_anchor=[1, 1], frameon=False
    )
    ax.set_facecolor('white')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    ax.grid(axis='y', linestyle='--', linewidth=0.5, color='gray', zorder=0)
    plt.savefig(f"{PATH_PLOTS}/plot_heat_tech_ratio_{clusters}_{planning_horizon}.png", dpi=600, bbox_inches = 'tight')


def define_table_df(scenarios):
    # Define column levels
    col_level_0 = ["2030"]*4 + ["2040"]*4 + ["2050"]*4
    col_level_1 = list(scenarios.values()) * 3
    # Create a MultiColumns
    multi_cols = pd.MultiIndex.from_arrays([col_level_0, col_level_1], names=['Year', 'Scenario'])
    df = pd.DataFrame(columns=multi_cols)
    return df


def fill_table_df(df, planning_horizon, scenarios, values):
    for scenario in scenarios.values():
        for tech_name, _ in values.iterrows():
            if scenario in values.columns:
                df.loc[tech_name, (planning_horizon, scenario)] = values.loc[tech_name, scenario]
    return df


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "plot_heat_tech_ratio", 
            clusters="48",
            planning_horizon=["2030"],
        )
    # update config based on wildcards
    config = update_config_from_wildcards(snakemake.config, snakemake.wildcards)


    # network parameters
    co2l_limits = CO2L_LIMITS
    line_limits = LINE_LIMITS
    clusters = config["plotting"]["clusters"]
    planning_horizons = config["plotting"]["planning_horizon"]
    planning_horizons = [str(x) for x in planning_horizons]
    opts = config["plotting"]["sector_opts"]

    # define scenario namings
    scenarios = {"flexible": "Widespread\nRenovation",
                 "retro_tes": "Widespread\nRenovation &\nElectrification",
                 "flexible-moderate": "Limited\nRenovation",
                 "rigid": "Business\nas Usual &\nElectrification"}

    # initialize df for storing table information
    table_cap_df = define_table_df(scenarios)
    
    for planning_horizon in planning_horizons:
        lineex = line_limits[planning_horizon]
        sector_opts = f"Co2L{co2l_limits[planning_horizon]}-{opts}"

        # if planning_horizon is 2020
        if planning_horizon == BAU_HORIZON:
            scenarios = {"BAU": "BASE 2023"}
        else:
            scenarios = {"flexible": "Widespread\nRenovation",
                        "retro_tes": "Widespread\nRenovation &\nElectrification",
                        "flexible-moderate": "Limited\nRenovation",
                        "rigid": "Business\nas Usual &\nElectrification"}

        # move to submodules/pypsa-eur
        change_path_to_pypsa_eur()
        
        # define dataframe to store heat tech capacities
        capacities_df = pd.DataFrame()

        # compute heat tech capacities
        for scenario, nice_name in scenarios.items():
            # load the network
            n = load_network(lineex, clusters, sector_opts, planning_horizon, scenario)

            if n is None:
                # Skip further computation for this scenario if network is not loaded
                print(f"Network is not found for scenario '{scenario}', planning year '{planning_horizon}'. Skipping...")
                continue

            # get heat tech capacities
            heat_tech_caps = get_heat_capacities(n, nice_name)
            capacities_df = pd.concat([capacities_df, heat_tech_caps], axis=1)

        # move to base directory
        change_path_to_base()

        # plot capacities
        if not capacities_df.empty:
            # plot figures
            plot_capacities(capacities_df, clusters, planning_horizon)
            # store to df
            table_cap_df = fill_table_df(table_cap_df, planning_horizon, scenarios, capacities_df)
            
    if not table_cap_df.empty:
        # save to csv
        table_cap_df.index.name = "Capacity [MW_el]"
        table_cap_df.columns = replace_multiindex_values(table_cap_df.columns, 
                                                         ("2040", "Limited\nRenovation"),
                                                         ("2040","Limited\nRenovation &\nElectrification"))
        table_cap_df.columns = replace_multiindex_values(table_cap_df.columns, 
                                                         ("2050", "Limited\nRenovation"),
                                                         ("2050","Limited\nRenovation &\nElectrification"))
        table_cap_df.to_csv(snakemake.output.table)
        
