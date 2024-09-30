# SPDX-FileCopyrightText:  Open Energy Transition gGmbH
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import os
import sys
sys.path.append("../submodules/pypsa-eur")
import pypsa
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import geopandas as gpd
import cartopy.crs as ccrs
import logging
import colors as c
import warnings
warnings.filterwarnings("ignore")
from _helpers import mock_snakemake, update_config_from_wildcards, load_network, \
                     change_path_to_pypsa_eur, change_path_to_base, \
                     LINE_LIMITS, CO2L_LIMITS, BAU_HORIZON, replace_multiindex_values

logger = logging.getLogger(__name__)

RESULTS_DIR = "plots/results"

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
    "solid biomass for industry": "solid biomass for industry",
    "solid biomass for industry CC": "solid biomass for industry",
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
        "co2",
        "DAC",
        "biogas",
        "solid biomass for industry",
        "solid biomass CHP",

        "uranium",
        "nuclear",
        "solid biomass",
        "biogas",
        "solid biomass for industry",
        "gas for industry",
        "coal for industry",
        "shipping oil",
        "shipping methanol",
        "naphtha for industry",
        "land transport oil",
        "kerosene for aviation",
        "process emissions",
        "SMR",
        "coal",
        "lignite",
        "methanol",
        "oil",

        "gas pipeline",
        "H2 pipeline",
        
        "H2 Electrolysis",
        "H2 Fuel Cell",
        "DAC",
        "Fischer-Tropsch",
        "methanation",
        "BEV charger",
        "V2G",
        "methanolisation",
        
        "battery storage",
        "gas storage",
        "H2 storage",
        "TES",
        "OCGT",
        "CCGT",

        "co2",
        "CO2 sequestration",

        "gas CHP",
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

def get_co2_balance(n, nice_name):
    co2_carrier = 'co2'
    balance = n.statistics.energy_balance()
    # get co2 balance
    df = balance.loc[:, :, co2_carrier].droplevel(0).groupby(level=0).sum().to_frame()
    df.columns = [nice_name]
    return df


def plot_co2_balance(co2_df, clusters, planning_horizon, plot_width=7):
    # convert to Billion tCO2_eq
    co2_df = co2_df / 1e9
    # filter out technologies with very small emission
    max_emissions = co2_df.abs().sum().max() / 2
    co2_threshold = max_emissions / 100 # 1% of max
    to_drop = co2_df.index[co2_df.abs().max(axis=1) < co2_threshold]
    logger.info(
        f"Dropping technology with co2 balance below {co2_threshold} Bton CO2_eq per year"
    )
    logger.debug(co2_df.loc[to_drop])
    co2_df = co2_df.drop(to_drop)

    # reorder index
    new_index = PREFERRED_ORDER.intersection(co2_df.index).append(
        co2_df.index.difference(PREFERRED_ORDER)
    )
    co2_reordered = co2_df.T[new_index]

    # get colors
    colors = c.tech_colors
    # plot co2 balance
    fig, ax = plt.subplots(figsize=(plot_width,9))
    co2_reordered.plot.bar(ax=ax, 
                           stacked=True, 
                           color=[colors[x] for x in co2_reordered.columns])
    # configure plot
    handles, labels = ax.get_legend_handles_labels()
    handles.reverse()
    labels.reverse()
    # reverse negatives co2 technologies
    negative_co2_names = co2_reordered.loc[:, co2_reordered.mean()<0].columns
    neg_co2_length = len(negative_co2_names)
    handles = handles[:-neg_co2_length] + handles[-neg_co2_length:][::-1]
    labels = labels[:-neg_co2_length] + labels[-neg_co2_length:][::-1]

    plt.xticks(rotation=0, fontsize=10)
    ax.set_ylabel("CO$_2$ emissions [BtCO$_{2-eq}$]")
    ax.set_xlabel("")
    ax.set_yticks(np.arange(-4.0, 4.0, 0.5))
    ax.set_ylim([-3.6,3.6])

    x_ticks = list(co2_df.columns)
    if planning_horizon in ["2040", "2050"] and "Limited \nRenovation" in x_ticks:
        # replace name for Limited Renovation scenario for 2030 to be LROH
        x_ticks[x_ticks.index("Limited \nRenovation")] = "Limited \nRenovation &\nElectrification"

    ax.set_xticklabels(x_ticks)

    # Turn off both horizontal and vertical grid lines
    ax.grid(False, which='both')
    ax.legend(
        handles, labels, ncol=1, loc="upper left", bbox_to_anchor=[1, 1], frameon=False
    )

    if planning_horizon == BAU_HORIZON:
        ax.set_title("2023", fontsize=12)
    else:
        ax.set_title(planning_horizon, fontsize=12)

    # Percentage drop for renovation scenarios
    if scenarios["rigid"] in co2_df.columns:
        total_co2_rigid = co2_df[co2_df>0].sum(axis=0)[scenarios["rigid"]]
        percentage_lower = 100 * (co2_df[co2_df>0].sum(axis=0) - total_co2_rigid) / total_co2_rigid
        
        # Calculate x-coordinates for the groups
        unique_x_coords = sorted(list(set([bar.get_x() + bar.get_width() / 2 for bar in ax.patches])))
        
        # Add arrows and percentage texts with corrected positions
        arrowprops = dict(facecolor='red', shrink=0.05, width=1, headwidth=8)
        
        # Annotate each group
        for i, x in enumerate(unique_x_coords[:-1]):
            if percentage_lower[i] > 0:
                color_percent = 'red'
            else:
                color_percent = 'green'
            plt.annotate(
                f'{percentage_lower[i]:.2f}%', 
                xy=(x, co2_df[co2_df>0].iloc[:,i].sum()), 
                xytext=(x, co2_df[co2_df>0].iloc[:,i].sum()+0.15), 
                fontsize=12, 
                color=color_percent, 
                ha='center'
            )

        # add horizontal line
        bar_width = unique_x_coords[1] - unique_x_coords[0]
        line_start = unique_x_coords[0] - bar_width * 0.25
        line_end = unique_x_coords[-1] + bar_width * 0.25
        ax.plot([line_start, line_end], [total_co2_rigid, total_co2_rigid], color='red', linestyle='--', linewidth=2)

    ax.set_facecolor('white')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    ax.grid(axis='y', linestyle='--', linewidth=0.5, color='gray', zorder=0)
    plt.savefig(f"{RESULTS_DIR}/plot_co2_level_{clusters}_{planning_horizon}.png", dpi=600, bbox_inches='tight')


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
            "plot_co2_level", 
            clusters="48",
        )
    # update config based on wildcards
    config = update_config_from_wildcards(snakemake.config, snakemake.wildcards)


    # network parameters
    co2l_limits = CO2L_LIMITS
    line_limits = LINE_LIMITS
    clusters = config["plotting"]["clusters"]
    opts = config["plotting"]["sector_opts"]
    planning_horizons = config["plotting"]["planning_horizon"]
    planning_horizons = [str(x) for x in planning_horizons if not str(x) == BAU_HORIZON]

    # define scenario namings
    scenarios = {"flexible": "Widespread \nRenovation",
                "retro_tes": "Widespread \nRenovation &\nElectrification",
                "flexible-moderate": "Limited \nRenovation",
                "rigid": "Business\nas Usual &\nElectrification"}


    # initialize df for storing co2 balance information
    table_co2_df = define_table_df(scenarios)

    # initialize df for storing CO2 emissions
    table_emissions_df = define_table_df(scenarios)

    for planning_horizon in planning_horizons:
        lineex = line_limits[planning_horizon]
        sector_opts = f"Co2L{co2l_limits[planning_horizon]}-{opts}"

        # move to submodules/pypsa-eur
        change_path_to_pypsa_eur()
    
        # load networks
        co2_df = pd.DataFrame()
        for scenario, nice_name in scenarios.items():
            n = load_network(lineex, clusters, sector_opts, planning_horizon, scenario)

            if n is None:
                # Skip further computation for this scenario if network is not loaded
                print(f"Network is not found for scenario '{scenario}', planning year '{planning_horizon}'. Skipping...")
                continue
            
            co2_balance = get_co2_balance(n, nice_name)
            co2_df = co2_df.join(co2_balance, how="outer").fillna(0)
        # grouping and renaming technologies
        co2_df = co2_df.groupby(co2_df.index.map(rename_techs)).sum()

        # move to base directory
        change_path_to_base()

        # plot co2 balance
        if not co2_df.empty:
            plot_co2_balance(co2_df, clusters, planning_horizon)
            table_co2_df = fill_table_df(table_co2_df, planning_horizon, scenarios, co2_df)
            co2_emission = -co2_df.loc["co2", :].to_frame().T
            table_emissions_df = fill_table_df(table_emissions_df, planning_horizon, scenarios, co2_emission)
    
    # Add BAU scenario
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
        print(f"Network is not found for scenario '{scenario}', planning year '{BAU_horizon}'. Skipping...")
    else:
        # get co2 balance for BAU and group technologies
        co2_BAU = get_co2_balance(n, "Baseline 2023")
        co2_BAU = co2_BAU.groupby(co2_BAU.index.map(rename_techs)).sum()
        if not table_co2_df.empty and not co2_BAU.empty:
            plot_co2_balance(co2_BAU, clusters, BAU_horizon, plot_width=1.5)
            table_co2_df = fill_table_df(table_co2_df, BAU_horizon, {"BAU":"Baseline 2023"}, co2_BAU)
            table_co2_df = table_co2_df.fillna(0)
            co2_emission_BAU = -co2_BAU.loc["co2"].values[0]

            # CO2 savings compared to BAU
            co2_savings_df = co2_emission_BAU - table_emissions_df
            co2_savings_df.index = ["CO2 emission savings [tCO2_eq]"]
            # saving in cubic meter of gas (1 m3 natural gas = 1,9 kg CO2) Source: https://www.eeagrants.gov.pt/media/2776/conversion-guidelines.pdf
            co2_savings_df.loc["Equivalent gas savings [trillion m3]"] = co2_savings_df.loc["CO2 emission savings [tCO2_eq]", :] * 1000 / 1.9 / 1e12


    # move to base directory
    change_path_to_base()

    # store to csv
    if not table_co2_df.empty:
        # save to csv
        table_co2_df.index.name = "CO2 emissions [tCO2_eq]"
        table_co2_df.columns = replace_multiindex_values(table_co2_df.columns, 
                                                         ("2040", "Limited \nRenovation"),
                                                         ("2040","Limited \nRenovation &\nElectrification"))
        table_co2_df.columns = replace_multiindex_values(table_co2_df.columns, 
                                                         ("2050", "Limited \nRenovation"),
                                                         ("2050","Limited \nRenovation &\nElectrification"))
        table_co2_df.to_csv(snakemake.output.table)

    if not co2_savings_df.empty:
        co2_savings_df.columns = replace_multiindex_values(co2_savings_df.columns, 
                                                           ("2040", "Limited \nRenovation"),
                                                           ("2040","Limited \nRenovation &\nElectrification"))
        co2_savings_df.columns = replace_multiindex_values(co2_savings_df.columns, 
                                                           ("2050", "Limited \nRenovation"),
                                                           ("2050","Limited \nRenovation &\nElectrification"))
        co2_savings_df.to_csv(snakemake.output.table_savings)
