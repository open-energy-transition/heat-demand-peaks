# %%
"""
# 0. Importing libraries
"""

# %%
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
logger = logging.getLogger(__name__)
plt.style.use("ggplot")

# get the current working directory
base_path = os.path.abspath(os.path.join(__file__ ,"../../submodules"))
# path to pypsa-eur
pypsa_path = "pypsa-eur/"
# absolute path to pypsa-eur

new_path = os.path.join(base_path, pypsa_path)
# change path to pypsa-eur
os.chdir(new_path)

# %%
"""
# 1. Define functions
"""

dont_plot = ["gas storage"]

prefix_to_remove = [
    "residential ",
    "services ",
    "urban ",
    "rural ",
    "central ",
    "decentral ",
]

rename_if_contains = [
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

rename_if_contains_dict = {
    "water tanks": "TES",
    "retrofitting": "building retrofitting",
    # "H2 Electrolysis": "hydrogen storage",
    # "H2 Fuel Cell": "hydrogen storage",
    # "H2 pipeline": "hydrogen storage",
    "battery": "battery storage",
    # "CC": "CC"
}

rename = {
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

# %%
# renaming function
def rename_techs(label):

    for ptr in prefix_to_remove:
        if label[: len(ptr)] == ptr:
            label = label[len(ptr) :]

    for rif in rename_if_contains:
        if rif in label:
            label = rif

    for old, new in rename_if_contains_dict.items():
        if old in label:
            label = new

    for old, new in rename.items():
        if old == label:
            label = new
    return label


preferred_order = pd.Index(
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

# %%
"""
# 1. Loading the networks
"""

co2l_limits = {2030:0.45, 2040:0.1, 2050:0.0}
line_limits = {2030:"v1.15", 2040:"v1.3", 2050:"v1.5"}
space_resolution = 48
planning = 2030
lineex = line_limits[planning]
sector_opts = f"Co2L{co2l_limits[planning]}-1H-T-H-B-I"




# %%
FILE = f"elec_s_{space_resolution}_l{lineex}__{sector_opts}_{planning}.nc"
DIR = "results/rigid/postnetworks"
n_rigid = pypsa.Network(os.path.join(DIR, FILE))

# %%
FILE = f"elec_s_{space_resolution}_l{lineex}__{sector_opts}_{planning}.nc"
DIR = "results/flexible/postnetworks"
n_flex = pypsa.Network(os.path.join(DIR, FILE))

# %%
FILE = f"elec_s_{space_resolution}_l{lineex}__{sector_opts}_{planning}.nc"
DIR = "results/retro_tes/postnetworks"
n_igas_tes = pypsa.Network(os.path.join(DIR, FILE))

# %%
FILE = f"elec_s_{space_resolution}_l{lineex}__{sector_opts}_{planning}.nc"
DIR = "results/flexible-moderate/postnetworks"
n_mod = pypsa.Network(os.path.join(DIR, FILE))

# %%
network = {"rigid":n_rigid, "igas+tes":n_igas_tes, "flexible":n_flex, "mod":n_mod}

# change directory back to original
os.chdir(base_path)
# relative path to folder where to store plots
PATH_PLOTS = "../plots/results/"
# create folder to store images
os.makedirs(PATH_PLOTS, exist_ok=True)
# %%
"""
# 2. Total System Cost estimation (Fig. 5)
"""

# %%
"""
## 2A. Cost estimation
"""

# %%
costs_rigid = n_rigid.statistics()[["Capital Expenditure", "Operational Expenditure"]].dropna()
full_costs_rigid = costs_rigid.sum(axis=1).droplevel(0).to_frame()
full_costs_rigid.columns = ["Supplied Heating"]

# %%
costs_flex = n_flex.statistics()[["Capital Expenditure", "Operational Expenditure"]].dropna()
full_costs_flex = costs_flex.sum(axis=1).droplevel(0).to_frame()
full_costs_flex.columns = ["Efficient Heating"]

# %%
costs_igas_tes = n_igas_tes.statistics()[["Capital Expenditure", "Operational Expenditure"]].dropna()
full_costs_igas_tes = costs_igas_tes.sum(axis=1).droplevel(0).to_frame()
full_costs_igas_tes.columns = ["Efficient Green Heating"]

# %%
costs_moderate = n_mod.statistics()[["Capital Expenditure", "Operational Expenditure"]].dropna()
full_costs_moderate = costs_moderate.sum(axis=1).droplevel(0).to_frame()
full_costs_moderate.columns = ["Semi-efficient Heating"]

# %%
cost_df = full_costs_flex.join(full_costs_igas_tes, how="outer").join(full_costs_moderate, how="outer").join(full_costs_rigid, how="outer").fillna(0)
cost_df = cost_df.drop("oil", axis=0)

# %%
"""
## 2B. Plot 
"""

# %%
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

    new_index = preferred_order.intersection(df.index).append(
        df.index.difference(preferred_order)
    )

    for remove_tech in dont_plot:
        new_index = new_index.drop(remove_tech)

    new_columns = df.sum().sort_values().index  


    fig, ax = plt.subplots(figsize=(6, 8))

    df.loc[new_index].T.plot(
        kind="bar",
        ax=ax,
        stacked=True,
        color=[c.tech_colors[i] for i in new_index],
    )

    handles, labels = ax.get_legend_handles_labels()

    handles.reverse()
    labels.reverse()

    costs_max = cost_df.sum().max() / 1e9
    ax.set_ylim([0, costs_max])
    plt.xticks(rotation=10, fontsize=12)

    ax.set_ylabel("System Cost [EUR billion per year]")

    ax.set_xlabel("")

    ax.grid(axis="x")

    ax.legend(
        handles, labels, ncol=1, loc="upper left", bbox_to_anchor=[1, 1], frameon=False
    )
    
    ax.set_facecolor('white')
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    ax.grid(axis='y', linestyle='--', linewidth=0.5, color='gray')
    plt.savefig(PATH_PLOTS+f'plot_total_costs_{space_resolution}_{planning}.png', dpi=600, bbox_inches = 'tight')
    
    
plot_costs(cost_df)
