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
from _helpers import mock_snakemake, update_config_from_wildcards
logger = logging.getLogger(__name__)

# get the current working directory
BASE_PATH = os.path.abspath(os.path.join(__file__ ,"../.."))
# relative path to folder where to store plots
PATH_PLOTS = "plots/results/"

DONT_PLOT = ["gas storage"]

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
    "gas": "gas storage",
    'gas pipeline new': 'gas pipeline',
    "gas for industry CC": "gas for industry",
    "SMR CC": "SMR",
    "process emissions CC": "process emissions",
    "Battery Storage": "battery storage",
    'H2 Store': "H2 storage",
    'Hydrogen Storage': "H2 storage",
    'co2 sequestered': "CO2 sequestration",
    "solid biomass transport": "solid biomass"
}

PREFERRED_ORDER = pd.Index(
    [
        "nuclear",
        "solid biomass",
        "biogas",
        "gas for industry",
        "methanol",
        "oil",
        "coal",
        
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
        "building retrofitting",
     ]
)


def change_path_to_pypsa_eur():
    # path to pypsa-eur
    pypsa_path = "submodules/pypsa-eur/"
    # absolute path to pypsa-eur
    new_path = os.path.join(BASE_PATH, pypsa_path)
    # change path to pypsa-eur
    os.chdir(new_path)


def change_path_to_base():
    os.chdir(BASE_PATH)
    # create folder to store images
    os.makedirs(PATH_PLOTS, exist_ok=True)


def load_network(lineex, space_resolution, sector_opts, planning, scenario):
    FILE = f"elec_s_{space_resolution}_l{lineex}__{sector_opts}_{planning}.nc"
    DIR = f"results/{scenario}/postnetworks"
    try:
        n = pypsa.Network(os.path.join(DIR, FILE))
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return None
    return n


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


def compute_costs(n, nice_name):
    costs = n.statistics()[["Capital Expenditure", "Operational Expenditure"]].dropna()
    full_costs = costs.sum(axis=1).droplevel(0).to_frame()
    full_costs.columns = [nice_name]
    return full_costs


def plot_costs(cost_df):
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

    for remove_tech in DONT_PLOT:
        new_index = new_index.drop(remove_tech)

    new_columns = df.sum().sort_values().index  


    fig, ax = plt.subplots(figsize=(6, 8))

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

    costs_max = cost_df.sum().max() / 1e9
    ax.set_ylim([0, costs_max])
    plt.xticks(rotation=10, fontsize=12)

    ax.set_ylabel("System Cost [EUR billion per year]")

    ax.set_xlabel("")

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
    plt.savefig(snakemake.output.figure, dpi=600, bbox_inches = 'tight')
    

if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "plot_total_costs", 
            space_resolution="48",
            planning="2030",
        )
    # update config based on wildcards
    config = update_config_from_wildcards(snakemake.config, snakemake.wildcards)


    # move to submodules/pypsa-eur
    change_path_to_pypsa_eur()

    # network parameters
    co2l_limits = {"2030":"0.45", "2040":"0.1", "2050":"0.0"}
    line_limits = {"2030":"v1.15", "2040":"v1.3", "2050":"v1.5"}
    space_resolution = config["plotting"]["space_resolution"]
    planning = config["plotting"]["planning"]
    time_resolution = config["plotting"]["time_resolution"]
    lineex = line_limits[planning]
    sector_opts = f"Co2L{co2l_limits[planning]}-{time_resolution}-T-H-B-I"

    # define scenario namings
    scenarios = {"flexible": "Efficient Heating", 
                 "retro_tes": "Efficient Green Heating", 
                 "flexible-moderate": "Semi-Efficient Heating", 
                 "rigid": "Non-efficient Heating"}

    # load networks
    networks = {}
    cost_df = pd.DataFrame()
    for scenario, nice_name in scenarios.items():
        n = load_network(lineex, space_resolution, sector_opts, planning, scenario)

        if n is None:
            # Skip further computation for this scenario if network is not loaded
            print(f"Network is not found for scenario '{scenario}', planning year '{planning}', and time resolution of '{time_resolution}'. Skipping...")
            continue

        # calculate costs for scenario
        full_costs = compute_costs(n, nice_name)
        cost_df = cost_df.join(full_costs, how="outer").fillna(0)

    # drop oil from plot
    if "oil" in cost_df.index:
        cost_df = cost_df.drop("oil", axis=0)

    # move to base directory
    change_path_to_base()

    # plot costs
    if not cost_df.empty:
        plot_costs(cost_df)
