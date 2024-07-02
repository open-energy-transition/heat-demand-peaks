# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  Open Energy Transition gGmbH
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import os
import sys
sys.path.append("../submodules/pypsa-eur")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import logging
import colors as c
import warnings
warnings.filterwarnings("ignore")
from _helpers import mock_snakemake, update_config_from_wildcards, load_network, \
                     change_path_to_pypsa_eur, change_path_to_base, \
                     LINE_LIMITS, CO2L_LIMITS, BAU_HORIZON, replace_multiindex_values

logger = logging.getLogger(__name__)

RESULTS_DIR = "plots/results"

GAS_BOILERS = [
    "Link:residential rural gas boiler",
    "Link:residential urban decentral gas boiler",
    "Link:services rural gas boiler",
    "Link:services urban decentral gas boiler",
    "Link:urban central gas boiler",
]

PREFIX_TO_REMOVE = [
    "residential ",
    "services ",
    "urban ",
    "rural ",
    "central ",
    "decentral ",
]

RENAME_IF_CONTAINS = [
    "solid biomass CHP",
    "gas CHP",
    "gas boiler",
    "biogas",
    "solar thermal",
    "air heat pump",
    "ground heat pump",
    "resistive heater",
    "Fischer-Tropsch",
]

RENAME_IF_CONTAINS_DICT = {
    "water tanks": "TES",
    "retrofitting": "building retrofitting",
    # "H2 Electrolysis": "hydrogen storage",
    # "H2 Fuel Cell": "hydrogen storage",
    # "H2 pipeline": "hydrogen storage",
    "battery": "battery storage",
    # "CC": "CC"
}

RENAME = {
    "Solar": "solar PV",
    "solar": "solar PV",
    "Sabatier": "methanation",
    "helmeth" : "methanation",
    "Offshore Wind (AC)": "offshore wind",
    "Offshore Wind (DC)": "offshore wind",
    "Onshore Wind": "onshore wind",
    "offwind-ac": "offshore wind",
    "offwind-dc": "offshore wind",
    "Run of River": "hydroelectricity",
    "Run of river": "hydroelectricity",
    "Reservoir & Dam": "hydroelectricity",
    "Pumped Hydro Storage": "hydroelectricity",
    "PHS": "hydroelectricity",
    "NH3": "ammonia",
    "co2 Store": "DAC",
    "co2 stored": "CO2 sequestration",
    "AC": "transmission lines",
    "DC": "transmission lines",
    "B2B": "transmission lines",
    "solid biomass for industry": "solid biomass",
    "solid biomass for industry CC": "solid biomass",
    "electricity distribution grid": "distribution lines",
    "Open-Cycle Gas":"OCGT",
    "Combined-Cycle Gas":"CCGT",
    "gas": "gas storage",
    'gas pipeline new': 'gas pipeline',
    "gas for industry CC": "gas for industry",
    "SMR CC": "SMR",
    "process emissions CC": "process emissions",
    "Battery Storage": "battery storage",
    'H2 Store': "H2 storage",
    'Hydrogen Storage': "H2 storage",
    'co2 sequestered': "CO2 sequestration",
    "solid biomass transport": "solid biomass",
    "uranium": "nuclear",
}

PREFERRED_ORDER = pd.Index(
    [
        "uranium",
        "nuclear",
        "solid biomass",
        "biogas",
        "gas for industry",
        "coal for industry",
        "methanol",
        "oil",
        "lignite",
        "coal",
        "shipping oil",
        "shipping methanol",
        "naphtha for industry",
        "land transport oil",
        "kerosene for aviation",
        
        "transmission lines",
        "distribution lines",
        "gas pipeline",
        "H2 pipeline",
        
        "H2 Electrolysis",
        "H2 Fuel Cell",
        "DAC",
        "Fischer-Tropsch",
        "methanation",
        "BEV charger",
        "V2G",
        "SMR",
        "methanolisation",
        
        "battery storage",
        "gas storage",
        "H2 storage",
        "TES",
        
        "hydroelectricity",
        "OCGT",
        "CCGT",
        "onshore wind",
        "offshore wind",
        "solar PV",
        "solar thermal",
        "solar rooftop",

        "co2",
        "CO2 sequestration",
        "process emissions",

        "gas CHP",
        "solid biomass CHP",
        "resistive heater",
        "air heat pump",
        "ground heat pump",
        "gas boiler",
        "biomass boiler",
        "WWHRS",
        "building retrofitting",
        "WWHRS",
     ]
)


def rename_techs(label):

    for ptr in PREFIX_TO_REMOVE:
        if label[: len(ptr)] == ptr:
            label = label[len(ptr) :]

    for rif in RENAME_IF_CONTAINS:
        if rif in label:
            label = rif

    for old, new in RENAME_IF_CONTAINS_DICT.items():
        if old in label:
            label = new

    for old, new in RENAME.items():
        if old == label:
            label = new
    return label


def compute_capacity_expansion(n, nice_name):
    opt_capacities = n.statistics()["Optimal Capacity"]
    nom_capacities = n.statistics()["Installed Capacity"]
    diff_capacities = opt_capacities - nom_capacities
    diff_capacities = diff_capacities[diff_capacities.index.get_level_values(0).isin(["Generator","Link","StorageUnit"])]
    diff_caps = diff_capacities.droplevel(0).to_frame()
    diff_caps = diff_caps.groupby(diff_caps.index).sum()
    diff_caps.columns = [nice_name]
    return diff_caps


def plot_capacities(caps_df, clusters, planning_horizon, plot_width=7):
    df = caps_df.groupby(caps_df.index).sum()

    # drop solid biomass transport
    df = df[df.index != 'solid biomass transport']

    # convert to GW
    df = df / 1e3
    df = df.groupby(df.index.map(rename_techs)).sum()

    caps_threshold = 10
    to_drop = df.index[df.max(axis=1) < caps_threshold]  #df <

    logger.info(
        f"Dropping technology with capacity below {caps_threshold} GW"
    )
    logger.debug(df.loc[to_drop])

    df = df.drop(to_drop)

    logger.info(f"Total optimal capacity is {round(df.sum())} GW")

    new_index = PREFERRED_ORDER.intersection(df.index).append(
        df.index.difference(PREFERRED_ORDER)
    )

    _, ax = plt.subplots(figsize=(plot_width, 10))

    df.loc[new_index].T.plot(
        kind="bar",
        ax=ax,
        stacked=True,
        color=[c.tech_colors[i] for i in new_index],
        zorder=1,
    )

    handles, labels = ax.get_legend_handles_labels()

    handles.reverse()
    labels.reverse()

    plt.xticks(rotation=0, fontsize=10)

    ax.set_ylabel("Capacity expansion [GW]")

    ax.set_xlabel("")
    ax.set_ylim([0,16000])
    ax.set_yticks(np.arange(0, 17000, 2000))

    x_ticks = list(df.columns)
    if planning_horizon in ["2040", "2050"] and "Limited \nRenovation &\nCost-Optimal Heating" in x_ticks:
        # replace name for Limited Renovation scenario for 2030 to be LROH
        x_ticks[x_ticks.index("Limited \nRenovation &\nCost-Optimal Heating")] = "Limited \nRenovation &\nElectric Heating"

    ax.set_xticklabels(x_ticks)

    # Turn off both horizontal and vertical grid lines
    ax.grid(False, which='both')

    ax.legend(
        handles, labels, ncol=1, loc="upper left", bbox_to_anchor=[1, 1], frameon=False
    )
    if planning_horizon == BAU_HORIZON:
        ax.set_title("BAU", fontsize=12)
    else:
        ax.set_title(planning_horizon, fontsize=12)
    
    ax.set_facecolor('white')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    ax.grid(axis='y', linestyle='--', linewidth=0.5, color='gray', zorder=0)
    plt.savefig(f"{RESULTS_DIR}/plot_capacity_expansion_{clusters}_{planning_horizon}.png", dpi=600, bbox_inches = 'tight')
    return df.loc[new_index[::-1]]


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
            "plot_capacity_expansion", 
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
    planning_horizons = [str(x) for x in planning_horizons if not str(x) == BAU_HORIZON]
    opts = config["plotting"]["sector_opts"]

    # define scenario namings
    scenarios = {"flexible": "Optimal \nRenovation &\nCost-Optimal Heating", 
                "retro_tes": "Optimal \nRenovation &\nElectric Heating", 
                "flexible-moderate": "Limited \nRenovation &\nCost-Optimal Heating", 
                "rigid": "No \nRenovation &\nElectric Heating"}

    # initialize df for storing capacity expansion in table form
    table_cap_df = define_table_df(scenarios)

    for planning_horizon in planning_horizons:
        lineex = line_limits[planning_horizon]
        sector_opts = f"Co2L{co2l_limits[planning_horizon]}-{opts}"

        # move to submodules/pypsa-eur
        change_path_to_pypsa_eur()
    
        # load networks
        networks = {}
        capacities_df = pd.DataFrame()
        p_nom_exp_df = pd.DataFrame()
        for scenario, nice_name in scenarios.items():
            n = load_network(lineex, clusters, sector_opts, planning_horizon, scenario)

            if n is None:
                # Skip further computation for this scenario if network is not loaded
                print(f"Network is not found for scenario '{scenario}', planning year '{planning_horizon}'. Skipping...")
                continue

            # calculate capacity expansion for scenario
            capacities = compute_capacity_expansion(n, nice_name)
            capacities_df = capacities_df.join(capacities, how="outer").fillna(0)


        # move to base directory
        change_path_to_base()

        # plot capacities
        if not capacities_df.empty:
            processed_capacities_df = plot_capacities(capacities_df, clusters, planning_horizon)
            table_cap_df = fill_table_df(table_cap_df, planning_horizon, scenarios, processed_capacities_df)


    # add BAU
    BAU_horizon = BAU_HORIZON
    scenario = "BAU"
    lineex = line_limits[BAU_horizon]
    sector_opts = f"Co2L{co2l_limits[BAU_horizon]}-{opts}"
    
    # move to submodules/pypsa-eur
    change_path_to_pypsa_eur()

    n = load_network(lineex, clusters, sector_opts, BAU_horizon, scenario)

    # move to base directory
    change_path_to_base()

    if n is None:
        # Skip further computation for this scenario if network is not loaded
        print(f"Network is not found for scenario '{scenario}', BAU year '{BAU_horizon}'. Skipping...")
    else:
        capacities_BAU = compute_capacity_expansion(n, "BAU")
        if not table_cap_df.empty and not capacities_BAU.empty:
            processed_capacities_BAU = plot_capacities(capacities_BAU, clusters, BAU_horizon, plot_width=1.6)
            table_cap_df = fill_table_df(table_cap_df, BAU_horizon, {"BAU":"BAU"}, processed_capacities_BAU)


    # save all capacities to csv
    if not table_cap_df.empty:
        table_cap_df.index.name = "Capacity expansion [GW]"
        table_cap_df.columns = replace_multiindex_values(table_cap_df.columns, 
                                                         ("2040", "Limited \nRenovation &\nCost-Optimal Heating"),
                                                         ("2040","Limited \nRenovation &\nElectric Heating"))
        table_cap_df.columns = replace_multiindex_values(table_cap_df.columns, 
                                                         ("2050", "Limited \nRenovation &\nCost-Optimal Heating"),
                                                         ("2050","Limited \nRenovation &\nElectric Heating"))

        table_cap_df.to_csv(snakemake.output.capacities) 
