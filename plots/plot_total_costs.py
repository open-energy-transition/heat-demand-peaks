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


def compute_costs(n, nice_name, cost_type):
    assert cost_type in ["Operational", "Capital"], "Type variable must be 'Operational' or 'Capital'"
    # fix gas storage manually
    n.stores.loc["EU gas Store", "e_nom_opt"] = (
        n.stores_t.e["EU gas Store"].max() - n.stores_t.e["EU gas Store"].min()
    )
    costs = n.statistics()[[f"{cost_type} Expenditure"]]
    new_index = [':'.join(idx) for idx in costs.index]
    costs.index = new_index
    costs.columns = [nice_name]
    return costs


def plot_costs(cost_df, clusters, planning_horizon, plot_width=7):
    df = cost_df.groupby(cost_df.index).sum()

    # convert to billions
    df = df / 1e9
    df = df.groupby(df.index.map(rename_techs)).sum()

    costs_threshold = 0.5
    to_drop = df.index[df.max(axis=1) < costs_threshold]  #df <

    logger.info(
        f"Dropping technology with costs below {costs_threshold} EUR billion per year"
    )
    logger.debug(df.loc[to_drop])

    df = df.drop(to_drop)

    logger.info(f"Total system cost of {round(df.sum())} EUR billion per year")

    new_index = PREFERRED_ORDER.intersection(df.index).append(
        df.index.difference(PREFERRED_ORDER)
    )

    new_columns = df.sum().sort_values().index  


    fig, ax = plt.subplots(figsize=(plot_width, 9))

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

    ax.set_ylabel("System Cost [EUR billion per year]")

    ax.set_xlabel("")
    ax.set_ylim([0,1100])
    ax.set_yticks(np.arange(0, 1200, 100))

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
    
    ax.set_facecolor('white')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    ax.grid(axis='y', linestyle='--', linewidth=0.5, color='gray', zorder=0)
    plt.savefig(f"{RESULTS_DIR}/plot_total_costs_{clusters}_{planning_horizon}.png", dpi=600, bbox_inches = 'tight')
    return df.loc[new_index[::-1]]
    

def compute_capacities(n, nice_name):
    capacities = n.statistics()[["Optimal Capacity"]]
    capacities = capacities[capacities.index.get_level_values(0).isin(["Generator","Link","StorageUnit"])]
    full_caps = capacities.sum(axis=1).droplevel(0).to_frame()
    full_caps = full_caps.groupby(full_caps.index).sum()
    full_caps.columns = [nice_name]
    return full_caps


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

    _, ax = plt.subplots(figsize=(plot_width, 9))

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

    ax.set_ylabel("Installed capacities [GW]")

    ax.set_xlabel("")
    ax.set_ylim([0,20000])
    ax.set_yticks(np.arange(0, 21000, 2000))

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
    
    ax.set_facecolor('white')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    ax.grid(axis='y', linestyle='--', linewidth=0.5, color='gray', zorder=0)
    plt.savefig(f"{RESULTS_DIR}/plot_total_capacities_{clusters}_{planning_horizon}.png", dpi=600, bbox_inches = 'tight')
    return df.loc[new_index[::-1]]


def get_p_nom_opt(n, nice_name):
    capacities = n.statistics()[["Optimal Capacity"]]
    new_index = [':'.join(idx) for idx in capacities.index]
    capacities.index = new_index
    capacities.columns = [nice_name]
    return capacities


def sum_costs(cap_cost_df, op_cost_df):
    total_cost = cap_cost_df + op_cost_df
    new_index = [x.split(":")[1] for x in total_cost.index]
    total_cost.index = new_index
    return total_cost


def get_common_index(list1, list2):
    # get common index
    common_index = set(list1).intersection(set(list2))
    return common_index


def update_capital_cost(cap_costs_dict, p_nom_opt_dict, planning_horizon):
    delta_horizon = 10
    planning_horizon_init = int(planning_horizon)
    planning_horizon = int(planning_horizon)
    full_cost = pd.DataFrame()
    missing_scenarios = []

    while planning_horizon > 2030:
        # previous horizon's year
        previous_horizon = planning_horizon - delta_horizon
        # delta p_nom_opt
        delta_p_nom_opt = p_nom_opt_dict[str(planning_horizon)] - p_nom_opt_dict[str(previous_horizon)]
        delta_p_nom_opt = delta_p_nom_opt.clip(lower=0)
        # built_fraction = delta(p_nom_opt) / p_nom_opt
        built_fraction = (delta_p_nom_opt/p_nom_opt_dict[str(planning_horizon)]).fillna(0)
        # cost amount in planning horizon
        cost_horizon = cap_costs_dict[str(planning_horizon)] * built_fraction
        # sum cost
        full_cost = full_cost.add(cost_horizon, fill_value=0)

        # missing scenario
        missing_scenario = cost_horizon.columns[cost_horizon.isna().all()]
        missing_scenarios.extend(missing_scenario)

        # go to previous horizon
        planning_horizon = previous_horizon
    else:
        cost_horizon = cap_costs_dict[str(planning_horizon)]
        full_cost = full_cost.add(cost_horizon, fill_value=0)
        # remove missing scenario
        full_cost = full_cost.loc[:, ~full_cost.columns.isin(missing_scenarios)]

    # set gas boilers manually after finished
    full_cost.loc[GAS_BOILERS] = cap_costs_dict[str(planning_horizon_init)].loc[GAS_BOILERS]
    return full_cost


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
            "plot_total_costs", 
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

    # initialize capital cost and p_nom_opt storing dictionary for different horizons
    cap_costs_dict = {}
    p_nom_opt_dict = {}
    # initialize df for storing table information
    table_cost_df = define_table_df(scenarios)
    table_cap_df = define_table_df(scenarios)

    for planning_horizon in planning_horizons:
        lineex = line_limits[planning_horizon]
        sector_opts = f"Co2L{co2l_limits[planning_horizon]}-{opts}"

        # move to submodules/pypsa-eur
        change_path_to_pypsa_eur()
    
        # load networks
        networks = {}
        cap_cost_df = pd.DataFrame()
        op_cost_df = pd.DataFrame()
        capacities_df = pd.DataFrame()
        p_nom_opt_df = pd.DataFrame()
        cost_df = pd.DataFrame()
        for scenario, nice_name in scenarios.items():
            n = load_network(lineex, clusters, sector_opts, planning_horizon, scenario)

            if n is None:
                # Skip further computation for this scenario if network is not loaded
                print(f"Network is not found for scenario '{scenario}', planning year '{planning_horizon}'. Skipping...")
                continue

            # calculate capital costs for scenario
            cap_costs = compute_costs(n, nice_name, "Capital")
            cap_cost_df = cap_cost_df.join(cap_costs, how="outer").fillna(0)

            # calculate operational costs for scenario
            op_costs = compute_costs(n, nice_name, "Operational")
            op_cost_df = op_cost_df.join(op_costs, how="outer").fillna(0)

            # calculate capacities for scenario
            capacities = compute_capacities(n, nice_name)
            capacities_df = capacities_df.join(capacities, how="outer").fillna(0)

            # get p_nom_opt for scenario
            p_nom_opt = get_p_nom_opt(n, nice_name)
            p_nom_opt_df = p_nom_opt_df.join(p_nom_opt, how="outer").fillna(0)

        # update capital costs based on previous horizons
        cap_costs_dict[planning_horizon] = cap_cost_df
        p_nom_opt_dict[planning_horizon] = p_nom_opt_df

        # if capital costs data is present (not empty)
        if not any([x.empty for x in cap_costs_dict.values()]):
            updated_caps_df = update_capital_cost(cap_costs_dict, p_nom_opt_dict, planning_horizon)

            # add capital and operational costs
            cost_df = sum_costs(updated_caps_df, op_cost_df)
            # reorder scenarios
            reorder_columns = [s for s in scenarios.values() if s in cost_df.columns]
            cost_df = cost_df[reorder_columns]

        # move to base directory
        change_path_to_base()

        # plot costs
        if not cost_df.empty:
            processed_cost_df = plot_costs(cost_df, clusters, planning_horizon)
            table_cost_df = fill_table_df(table_cost_df, planning_horizon, scenarios, processed_cost_df)

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
        cap_costs = compute_costs(n, "BAU", "Capital")
        op_costs = compute_costs(n, "BAU", "Operational")
        cost_BAU = sum_costs(cap_costs, op_costs)
        capacities_BAU = compute_capacities(n, "BAU")
        if not table_cost_df.empty and not cost_BAU.empty:
            processed_cost_BAU = plot_costs(cost_BAU, clusters, BAU_horizon, plot_width=1.6)
            table_cost_df = fill_table_df(table_cost_df, BAU_horizon, {"BAU":"BAU"}, processed_cost_BAU)

        if not table_cap_df.empty and not capacities_BAU.empty:
            processed_capacities_BAU = plot_capacities(capacities_BAU, clusters, BAU_horizon, plot_width=1.6)
            table_cap_df = fill_table_df(table_cap_df, BAU_horizon, {"BAU":"BAU"}, processed_capacities_BAU)


    # save all costs to csv
    if not table_cost_df.empty:
        table_cost_df.index.name = "System cost [EUR billion per year]"
        table_cost_df.columns = replace_multiindex_values(table_cost_df.columns, 
                                                          ("2040", "Limited \nRenovation &\nCost-Optimal Heating"),
                                                          ("2040","Limited \nRenovation &\nElectric Heating"))
        table_cost_df.columns = replace_multiindex_values(table_cost_df.columns, 
                                                          ("2050", "Limited \nRenovation &\nCost-Optimal Heating"),
                                                          ("2050","Limited \nRenovation &\nElectric Heating"))

        table_cost_df.to_csv(snakemake.output.costs)

    # save all capacities to csv
    if not table_cap_df.empty:
        table_cap_df.index.name = "Installed capacity [GW]"
        table_cap_df.columns = replace_multiindex_values(table_cap_df.columns, 
                                                         ("2040", "Limited \nRenovation &\nCost-Optimal Heating"),
                                                         ("2040","Limited \nRenovation &\nElectric Heating"))
        table_cap_df.columns = replace_multiindex_values(table_cap_df.columns, 
                                                         ("2050", "Limited \nRenovation &\nCost-Optimal Heating"),
                                                         ("2050","Limited \nRenovation &\nElectric Heating"))

        table_cap_df.to_csv(snakemake.output.capacities) 
