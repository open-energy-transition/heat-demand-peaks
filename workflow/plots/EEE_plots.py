# %%
"""
# 0. Importing libraries
"""

# %%
import os
import sys
sys.path.append("../pypsa-eur")
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
current_path = os.getcwd()
# path to pypsa-eur
pypsa_path = "workflow/pypsa-eur"
# absolute path to pypsa-eur
new_path = os.path.join(current_path, pypsa_path)
# change path to pypsa-eur
os.chdir(new_path)

# %%
"""
# 1. Define functions
"""

# %%
# renaming function
def rename_techs(label):
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
        
    }

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
        "solid biomass",
        "biogas",
        "gas for industry",
        "methanol",
        "oil",
        
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

        "gas CHP",
        "solid biomass CHP",
        "resistive heater",
        "air heat pump",
        "ground heat pump",
        "gas boiler",
        "biomass boiler",
        "building retrofitting",
        
        "co2",
        "CO2 sequestration",
        "process emissions"
     ]
)

# %%
# renaming function
def rename_techs2(label):
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
        "water tanks": "water tanks discharger",
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
        
    }

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
        "solid biomass",
        "biogas",
        "gas for industry",
        "methanol",
        "oil",
        
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
        "water tanks discharger",
        
        "hydroelectricity",
        "OCGT",
        "onshore wind",
        "offshore wind",
        "solar PV",
        "solar thermal",
        "solar rooftop",
        
        "building retrofitting",
        "solid biomass CHP",
        "gas CHP",
        "biomass boiler",
        "gas boiler",
        "resistive heater",
        "air heat pump",
        "ground heat pump",
        
        "co2",
        "CO2 sequestration",
        "process emissions"
     ]
)

# %%
"""
# 1. Loading the networks
"""

# %%
FILE = "elec_s_48_lcopt__Co2L0-2H-T-H-B_2030.nc"
DIR = "results/rigid/postnetworks"
n_rigid = pypsa.Network(os.path.join(DIR, FILE))

# %%
FILE = "elec_s_48_lcopt__Co2L0-2H-T-H-B_2030.nc"
DIR = "results/flexible/postnetworks"
n_flex = pypsa.Network(os.path.join(DIR, FILE))

# %%
FILE = "elec_s_48_lcopt__Co2L0-2H-T-H-B_2030.nc"
DIR = "results/igas_tes/postnetworks"
n_igas_tes = pypsa.Network(os.path.join(DIR, FILE))

# %%
network = {"rigid":n_rigid, "igas+tes":n_igas_tes, "flexible":n_flex}

# change directory back to original
os.chdir(current_path)
# relative path to folder where to store plots
PATH_PLOTS = "workflow/plots/results/"
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
full_costs_rigid.columns = ["rigid"]

# %%
costs_flex = n_flex.statistics()[["Capital Expenditure", "Operational Expenditure"]].dropna()
full_costs_flex = costs_flex.sum(axis=1).droplevel(0).to_frame()
full_costs_flex.columns = ["flexible"]

# %%
costs_igas_tes = n_igas_tes.statistics()[["Capital Expenditure", "Operational Expenditure"]].dropna()
full_costs_igas_tes = costs_igas_tes.sum(axis=1).droplevel(0).to_frame()
full_costs_igas_tes.columns = ["igas+tes"]

# %%
cost_df = full_costs_flex.join(full_costs_igas_tes, how="outer").join(full_costs_rigid, how="outer").fillna(0)
cost_df

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
    plt.xticks(rotation=0, fontsize=12)

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
    plt.savefig(PATH_PLOTS+'costs.png', dpi=600, bbox_inches = 'tight')
    
    
plot_costs(cost_df)

# %%
"""
# 3. Installed peak capacituies (Fig. 6)
"""

# %%
"""
## 3A. Calculation
"""

# %%
common_rural_techs = ["solar thermal", "ground heat pump", "resistive heater", 
                      "gas boiler", "water tanks discharger" ]#, "biomass boiler" ]
common_urban_techs = ["solar thermal", "air heat pump", "resistive heater", 
                      "gas boiler", "water tanks discharger"] #, "biomass boiler" ]
urban_central_techs = ["solar thermal", "air heat pump", "resistive heater", 
                       "gas boiler", "gas CHP", "gas CHP CC", 
                       "solid biomass CHP","solid biomass CHP CC", "water tanks discharger"]

# %%
techs = {"residential rural":["residential rural "+ x for x in common_rural_techs],
         "services rural":["services rural "+ x for x in common_rural_techs],
         "residential urban decentral":["residential urban decentral "+ x for x in common_urban_techs],
         "services urban decentral":["services urban decentral "+ x for x in common_urban_techs],
         "urban central":["urban central "+ x for x in urban_central_techs]+["H2 Fuel Cell"]#, "Fischer-Tropsch"]
        }

# %%
techs

# %%
cap_full = pd.DataFrame(columns=["rigid", "tes", "igas", "igas+tes", "retro", "retro+tes", "retro+igas", "flexible"])
gen_full = pd.DataFrame(columns=["rigid", "tes", "igas", "igas+tes", "retro", "retro+tes", "retro+igas", "flexible"])
cap_full

# %%
# Calculation of installed capacities and generation
for name, n in network.items():
    print("Case: ", name)
    for r in techs:
        for i in techs[r]:
            gens_cap = n.generators.filter(like=i, axis=0).p_nom_opt
            links_cap = n.links.filter(like=i, axis=0).p_nom_opt
            if not gens_cap.empty:
                cap_full.loc[i, name] = gens_cap.sum()
            if not links_cap.empty:
                cap_full.loc[i, name] = links_cap.sum()
                
            gens_gen = n.generators_t.p.filter(like=i, axis=1).multiply(n_rigid.snapshot_weightings.generators,axis=0)
            if "CHP" in i or "Fuel Cell" in i:
                links_gen = -n.links_t.p2.filter(like=i, axis=1).multiply(n_rigid.snapshot_weightings.stores,axis=0)
            elif "Fischer-Tropsch" in i:
                links_gen = -n.links_t.p3.filter(like=i, axis=1).multiply(n_rigid.snapshot_weightings.stores,axis=0)
            else:
                links_gen = -n.links_t.p1.filter(like=i, axis=1).multiply(n_rigid.snapshot_weightings.stores,axis=0)  # for all links excep Fischer-Tropsch, H2 Fuel Cell, and CHPs
            # gas CHP bus2 
            if not gens_gen.empty:
                gen_full.loc[i, name] = gens_gen.sum().sum()
            if not links_gen.empty:
                gen_full.loc[i, name] = links_gen.sum().sum()
print(cap_full)
print(gen_full)

# %%
"""
## 3B. Plot
"""

# %%
fig, axs = plt.subplots(2, 3, figsize=(12, 6))
plt.subplots_adjust(wspace=0.2, hspace=0.6)
zones = {"rural":"rural", "urban individual":"urban decentral", "district heating":"urban central"}
k = 0
for nice_name, tech_name in zones.items():
    cur_zones = [x for x in techs.keys() if tech_name in x]  # eg. residential rural, services rural
    list_tech = [techs[x] for x in cur_zones]  # list of technology lists
    cur_techs = [t for sublist in list_tech for t in sublist]  # list of technologies
    
    cf = cap_full.loc[cur_techs,:] / 1e3
    cf = cf.groupby(cf.index.map(rename_techs2)).sum()
    logger.info(f"Total installed capacity of {cf.sum().astype(int).values} GW")
    
    gf = gen_full.loc[cur_techs, :] / 1e6
    gf = gf.groupby(gf.index.map(rename_techs2)).sum()
    logger.info(f"Total Generation of {gf.sum().astype(int).values} TWh")
    
    
    new_index = preferred_order.intersection(cf.index).append(
        cf.index.difference(preferred_order)
    )
    
    cf.loc[new_index].T.plot(
        kind="bar",
        ax=axs[0,k],
        stacked=True,
        color=[c.tech_colors[i] for i in new_index],
        legend=False,
        title = nice_name
    )
    gf.loc[new_index].T.plot(
        kind="bar",
        ax=axs[1,k],
        stacked=True,
        color=[c.tech_colors[i] for i in new_index],
        legend=False,
        title = nice_name
    )
    
    axs[0,k].set_facecolor('white')
    axs[0,k].spines['left'].set_color('black')
    axs[0,k].spines['bottom'].set_color('black')
    axs[0,k].grid(axis='y', linestyle='--', linewidth=0.5, color='gray')
    axs[1,k].set_facecolor('white')
    axs[1,k].spines['left'].set_color('black')
    axs[1,k].spines['bottom'].set_color('black')
    axs[1,k].grid(axis='y', linestyle='--', linewidth=0.5, color='gray')
    k += 1


handles1, labels1 = axs[0,0].get_legend_handles_labels()
handles2, labels2 = axs[0,2].get_legend_handles_labels()

handles = handles1[2:3] + handles2 + handles1[-1:]
labels = labels1[2:3] + labels2 + labels1[-1:]

handles.reverse()
labels.reverse()

axs[0,2].legend(
    handles, labels, ncol=1, loc="upper left", bbox_to_anchor=[1, 1], frameon=False
)
axs[0, 0].set_ylabel("capacity [GW]")
axs[1, 0].set_ylabel("generation [TWh]")
plt.savefig(PATH_PLOTS+'capacities.png', dpi=600, bbox_inches = 'tight')


# %%
"""
# 4.Heat supply profile (Fig. 7) 
"""

# %%
"""
## 4A. Estimation and plot
"""

# %%
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
import matplotlib.dates as mdates

# %%
first_date = "2013-01-10"
second_date = "2013-01-17"
bus = "RO1 0"

# %%
techs

# %%
network_small = {"rigid":n_rigid, "igas+tes":n_igas_tes, "flexible":n_flex}
fig, axs = plt.subplots(5, 3, figsize=(12, 6))
i = 0

for name, n in network_small.items():
    print("Case: ", name)
    j = 0
    for r in techs:
        current_technologies = techs[r]
        
        time = n.generators_t.p.multiply(n.snapshot_weightings.generators, 
                                              axis=0).loc[first_date:second_date].index
        
        solar_gen = pd.Series(index=time).fillna(0)
        fischer_gen = pd.Series(index=time).fillna(0)
        gas_CHP_gen = pd.Series(index=time).fillna(0)
        biomass_CHP_gen = pd.Series(index=time).fillna(0)
        fuel_cell_gen = pd.Series(index=time).fillna(0)
        ground_heat_gen = pd.Series(index=time).fillna(0)
        air_heat_gen = pd.Series(index=time).fillna(0)
        res_heater_gen = pd.Series(index=time).fillna(0)
        gas_boiler_gen = pd.Series(index=time).fillna(0)
        biomass_boiler_gen = pd.Series(index=time).fillna(0)
        water_dischar_gen = pd.Series(index=time).fillna(0)
        
        
#         idx_gen = n.generators.query("carrier in @curent_technologies").filter(like=bus, axis=0).index
#         idx_link = n.links.query("carrier in @curent_technologies").filter(like=bus, axis=0).index
        
        
    
        # generation of solar thermal collectors
        solar_gen = n.generators_t.p.multiply(n.snapshot_weightings.generators, 
                                              axis=0).loc[first_date:second_date].filter(like=bus).filter(like="solar thermal").filter(like=r).squeeze().fillna(0)
        if solar_gen.empty:
            solar_gen = pd.Series(index=time).fillna(0)
            
        if r == "urban central":
            # generation of Fischer-Tropsch
            fischer_gen = -n.links_t.p3.multiply(n.snapshot_weightings.stores, 
                                                 axis=0).loc[first_date:second_date].filter(like=bus).filter(like="Fischer").squeeze().fillna(0)
            if fischer_gen.empty:
                fischer_gen = pd.Series(index=time).fillna(0)
            
            # generation of gas CHP
            gas_CHP_gen = -n.links_t.p2.multiply(n.snapshot_weightings.stores,
                                                 axis=0).loc[first_date:second_date].filter(like=bus).filter(like="gas CHP").sum(axis=1).squeeze().fillna(0)
            if gas_CHP_gen.empty:
                gas_CHP_gen = pd.Series(index=time).fillna(0)
                
            # generation of solid biomass CHP
            biomass_CHP_gen = -n.links_t.p2.multiply(n.snapshot_weightings.stores, 
                                                     axis=0).loc[first_date:second_date].filter(like=bus).filter(like="biomass CHP").sum(axis=1).squeeze().fillna(0)
            if biomass_CHP_gen.empty:
                biomass_CHP_gen = pd.Series(index=time).fillna(0)
                
            # generation of Fuel Cell
            fuel_cell_gen = -n.links_t.p2.multiply(n.snapshot_weightings.stores,
                                                   axis=0).loc[first_date:second_date].filter(like=bus).filter(like="Fuel Cell").squeeze().fillna(0)
            if fuel_cell_gen.empty:
                fuel_cell_gen = pd.Series(index=time).fillna(0)
        
        # generation of ground heat pump
        ground_heat_gen = -n.links_t.p1.multiply(n.snapshot_weightings.stores,
                                                 axis=0).loc[first_date:second_date].filter(like=bus).filter(like="ground heat pump").filter(like=r).squeeze()
        if ground_heat_gen.empty:
            ground_heat_gen = pd.Series(index=time).fillna(0)
        
        # generation of air heat pump
        air_heat_gen = -n.links_t.p1.multiply(n.snapshot_weightings.stores,
                                              axis=0).loc[first_date:second_date].filter(like=bus).filter(like="air heat pump").filter(like=r).squeeze()
        if air_heat_gen.empty:
            air_heat_gen = pd.Series(index=time).fillna(0)
        
        # generation of resistive heater
        res_heater_gen = -n.links_t.p1.multiply(n.snapshot_weightings.stores,
                                                axis=0).loc[first_date:second_date].filter(like=bus).filter(like="resistive heater").filter(like=r).squeeze().fillna(0)
        if res_heater_gen.empty:
            res_heater_gen = pd.Series(index=time).fillna(0)
        
        # generation of gas boiler
        gas_boiler_gen = -n.links_t.p1.multiply(n.snapshot_weightings.stores,
                                                axis=0).loc[first_date:second_date].filter(like=bus).filter(like="gas boiler").filter(like=r).squeeze().fillna(0)
        if gas_boiler_gen.empty:
            gas_boiler_gen = pd.Series(index=time).fillna(0)

        # generation of biomass boiler
        biomass_boiler_gen = -n.links_t.p1.multiply(n.snapshot_weightings.stores,
                                                    axis=0).loc[first_date:second_date].filter(like=bus).filter(like="biomass boiler").filter(like=r).squeeze().fillna(0)
        if biomass_boiler_gen.empty:
            biomass_boiler_gen = pd.Series(index=time).fillna(0)

        # generation of water tanks discharger
        water_dischar_gen = -n.links_t.p1.multiply(n.snapshot_weightings.stores,
                                                   axis=0).loc[first_date:second_date].filter(like=bus).filter(like="water tanks discharger").filter(like=r).squeeze().fillna(0)
        if water_dischar_gen.empty:
            water_dischar_gen = pd.Series(index=time).fillna(0)
            
        # storage of water tanks charger
        water_char_store = n.links_t.p1.multiply(n.snapshot_weightings.stores,
                                                 axis=0).loc[first_date:second_date].filter(like=bus).filter(like="water tanks charger").filter(like=r).squeeze().fillna(0)
        if water_char_store.empty:
            water_char_store = pd.Series(index=time).fillna(0)
            
        # heat to DAC
        dac_need = n.links_t.p1.multiply(n.snapshot_weightings.stores,
                                         axis=0).loc[first_date:second_date].filter(like=bus).filter(like="DAC").filter(like=r).squeeze().fillna(0)
        if dac_need.empty:
            dac_need = pd.Series(index=time).fillna(0)
 

        # heat demand
        heat_demand = n.loads_t.p_set.multiply(n.snapshot_weightings.generators,
                                                   axis=0).loc[first_date:second_date].filter(like=bus).filter(like=r).squeeze().fillna(0)
        retrofit = n.generators_t.p.multiply(n.snapshot_weightings.generators,
                                             axis=0).loc[first_date:second_date].filter(like=bus).filter(like="retrofitting").filter(like=r).sum(axis=1).squeeze().fillna(0)
        if retrofit.empty:
            retrofit = pd.Series(index=time).fillna(0)
        total_heat_demand = heat_demand - retrofit
        
        # technologies
        items = ["ground heat pump", "air heat pump", "resistive heater",
                 "gas boiler", "biomass boiler", "gas CHP", "solid biomass CHP",
                 "Fischer-Tropsch", "H2 Fuel Cell", "water tanks discharger", "solar thermal"]
        # stacked plot for supply
        axs[j,i].stackplot(time, ground_heat_gen/1e3, air_heat_gen/1e3, res_heater_gen/1e3, 
                          gas_boiler_gen/1e3, biomass_boiler_gen/1e3, gas_CHP_gen/1e3, biomass_CHP_gen/1e3, 
                          fischer_gen/1e3, fuel_cell_gen/1e3, water_dischar_gen/1e3, solar_gen/1e3,
                          colors=[c.tech_colors[i] for i in items], zorder = 3)
        # stacked plot for demand
        axs[j,i].stackplot(time, water_char_store/1e3, dac_need/1e3,
                          colors=[c.tech_colors[i] for i in ["water tanks charger", "DAC"]], zorder = 3)
        # heat demand - retrofit
        line1, = axs[j,i].plot(time, total_heat_demand/1e3, '--', color='#000000', linewidth=2, label="Heat demand", zorder=10)
        
        
        # heat demand
        line2, = axs[j,i].plot(time, heat_demand/1e3, '--', color='#000000', linewidth=2, label="Heat demand", zorder=11)
        
        j += 1  # plot next graph
        
        # configuration of axis
        if j<5:
            axs[j-1,i].set_xticklabels([])
        axs[4,i].set_xticklabels(labels=time.day)
        axs[4,i].set_xlabel("Jan")
        axs[0,i].set_yticks([-10,0,10,20,30]) # limits and ticks for row 1
        axs[1,i].set_yticks([-2,0,4,8]) # limits and ticks for row 2
        axs[2,i].set_yticks([-10,0,10,20,30]) # limits and ticks for row 3
        axs[3,i].set_yticks([-3,0,3,6]) # limits and ticks for row 4
        axs[4,i].set_yticks([-5,0,5,10,15,20]) # limits and ticks for row 1
        axs[0,i].set_title(name)
        axs[j-1,i].set_facecolor('white')
        axs[j-1,i].spines['left'].set_color('black')
        axs[j-1,i].spines['bottom'].set_color('black')
        axs[j-1,i].grid(axis='y', linestyle='--', linewidth=0.5, color='gray')
        
        if i==0:
            yloc = [0.8, 0.65, 0.49, 0.33, 0.17]
            plt.gcf().text(0.92, yloc[j-1], f'{r}', fontsize=10, va='center')
         
    i += 1

handles = []
for i in items + ["water tanks charger", "DAC"]:
    current_patch = mpatches.Patch(color=c.tech_colors[i], label=i)
    handles.append(current_patch)
handles.append(line1)
lgd = axs[4,1].legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.5), ncol=5, facecolor="white")
ylabel = axs[2,0].set_ylabel("Energy [GWh]")
plt.savefig(PATH_PLOTS+'gen_profiles.png', dpi=600, bbox_inches = 'tight')


    



# %%
"""
# 5. Heat demand saved per country (Fig. 8)
"""

# %%
"""
## A. Importing shapes
"""

# %%
# set a custom path to `resources` folder
resources_path = "workflow/pypsa-eur/resources/flexible"
# onshore and offshore shapes
onshore_cl_fl = os.path.join(resources_path, "regions_onshore_elec_s_48.geojson")
# read onshore shapes
onshore_cl_df = gpd.read_file(onshore_cl_fl)

# %%
"""
## B. Prepare dataframes and plot
"""

# %%
heat_demand = n_flex.loads_t.p_set.multiply(n.snapshot_weightings.generators,
                                             axis=0).filter(like="heat").groupby(n.buses.location, 
                                                                                 axis=1).sum().sum().iloc[1:].to_frame()
heat_demand.columns = ["heat_demand"]
heat_demand.reset_index(inplace=True)
heat_demand.rename(columns={"location":"name"}, inplace=True)
heat_demand.head()

# %%
data_gdf = pd.concat([onshore_cl_df, heat_demand], axis=1, join="inner")
data_gdf = data_gdf.loc[:,~data_gdf.columns.duplicated()]
data_gdf.head()

# %%
# produces an empty subplot due to current mathplotlib logic
# sorry for that and see https://github.com/matplotlib/matplotlib/issues/18138
fig, ax = plt.subplots(1, 1, figsize=(15, 10))
# font = {'family' : 'normal',
#         'size'   : 22}
# plt.rc('font', **font)

# %%
column_name = "heat_demand"

data_gdf.plot(
    column=column_name, 
    legend=True,
    categorical=False,
    cmap='OrRd',
    edgecolor='#36454F',
    linewidth=0.2,
    figsize=(20, 10)
)
plt.title(column_name)
plt.axis("off")
plt.savefig(PATH_PLOTS+'heat_demand.png', dpi=600, bbox_inches = 'tight')

# %%
"""
# 6. Marginal prices for electricity (Fig. 9)
"""

# %%
"""
## 6A. Plot
"""

# %%
network_two = {"rigid":n_rigid, "flexible":n_flex}
all_prices = pd.DataFrame()
all_annual_cost = pd.DataFrame()

fig, axs = plt.subplots(2, 1, figsize=(12, 7))
plt.subplots_adjust(wspace=0.4, hspace=0.4)

# Plot prices per country
for name, n in network_two.items():
#     prices = n.buses_t.marginal_price.groupby(n.buses.country,axis=1).mean().mean().iloc[1:].sort_values()
    prices = n.loads_t.p_set.mul(n.buses_t.marginal_price).sum().groupby(n.buses.country).sum().div(n.loads_t.p_set.groupby(n.buses.country, axis=1).sum().sum()).iloc[1:].sort_values()
    prices.name = name
    all_prices = pd.concat([all_prices, prices], axis=1)

all_prices.plot.bar(ax=axs[0], width=0.6) # plotting bar chart

axs[0].set_yticks([0,20,40,60,80])
axs[0].set_ylim([0, 90])
axs[0].legend(loc="upper left", facecolor="white")
ylabel = axs[0].set_ylabel("EUR/MWh")
xlabel = axs[0].set_xlabel("countries")
axs[0].set_title("Electricity price per country")
axs[0].set_facecolor('white')
axs[0].spines['left'].set_color('black')
axs[0].spines['bottom'].set_color('black')
axs[0].grid(axis='y', linestyle='--', linewidth=0.5, color='gray')

# Plot prices during year
for name, n in network_two.items():
    annual_cost = n.buses_t.marginal_price.mean(axis=1)
    annual_cost.name = name
    all_annual_cost = pd.concat([all_annual_cost, annual_cost], axis=1)
    
A = all_annual_cost.plot(ax=axs[1]) # plotting annual cost variation

axs[1].set_facecolor('white')
axs[1].spines['left'].set_color('black')
axs[1].spines['bottom'].set_color('black')
axs[1].grid(axis='y', linestyle='--', linewidth=0.5, color='gray')
yticks = axs[1].set_yticks([0,200,400,650])
axs[1].set_ylim([0, 700])
lgd = axs[1].legend(loc="upper right", facecolor="white")
date_rng = pd.date_range(start='2013-01-01', periods=13, freq='MS')
axs[1].set_xticks(date_rng)
a = axs[1].set_xticklabels(labels=date_rng.strftime('%B'))
ylabel = axs[1].set_ylabel("EUR/MWh (avg of countries)")
axs[1].set_title("Electricity price during year")
plt.savefig(PATH_PLOTS+'marginal_costs.png', dpi=600, bbox_inches = 'tight')

# %%


# %%
