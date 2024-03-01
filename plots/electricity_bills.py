import os
import sys
sys.path.append("../submodules/pypsa-eur")
import pypsa
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

# get the base working directory
BASE_PATH = os.path.abspath(os.path.join(__file__ ,"../.."))
# relative path to folder where to store plots
PATH_PLOTS = "plots/results/"

def change_path_to_pypsa_eur():
    # path to pypsa-eur
    pypsa_path = "submodules/pypsa-eur/"
    # absolute path to pypsa-eur
    new_path = os.path.join(BASE_PATH, pypsa_path)
    # change path to pypsa-eur
    os.chdir(new_path)


def load_networks(lineex, space_resolution, sector_opts, planning):
    FILE = f"elec_s_{space_resolution}_l{lineex}__{sector_opts}_{planning}_flexible.nc"
    DIR = "results/flexible/postnetworks"
    n_flex = pypsa.Network(os.path.join(DIR, FILE))

    FILE = f"elec_s_{space_resolution}_l{lineex}__{sector_opts}_{planning}_retro_tes.nc"
    DIR = "results/retro_tes/postnetworks"
    n_retro_tes = pypsa.Network(os.path.join(DIR, FILE))

    FILE = f"elec_s_{space_resolution}_l{lineex}__{sector_opts}_{planning}_flexible-moderate.nc"
    DIR = "results/flexible-moderate/postnetworks"
    n_flex_mod = pypsa.Network(os.path.join(DIR, FILE))

    FILE = f"elec_s_{space_resolution}_l{lineex}__{sector_opts}_{planning}_rigid.nc"
    DIR = "results/rigid/postnetworks"
    n_rigid = pypsa.Network(os.path.join(DIR, FILE))

    return n_flex, n_flex_mod, n_retro_tes, n_rigid


def change_path_to_base():
    os.chdir(BASE_PATH)
    # create folder to store images
    os.makedirs(PATH_PLOTS, exist_ok=True)


def get_households():
    household_data = pd.read_excel("plots/Households.xlsx",skiprows=1,usecols="B:D")
    households = household_data.set_index("Country")[ "Households (thousands)"]
    households.index.name=""
    return households


def electricity_bills(network, households):
    n = network.copy()
    
    rh_techs_elec = ['residential rural ground heat pump',
                     'residential rural resistive heater', 
                     'residential urban decentral air heat pump',
                     'residential urban decentral resistive heater']
    rh_techs_gas = ['residential rural gas boiler', 'residential urban decentral gas boiler']
    rh_techs_mCHP = ['residential rural micro gas CHP', 'residential urban decentral micro gas CHP']
    
    ev_tech_charge = ['BEV charger']
    ev_tech_discharge = ["V2G"]
    
    # map loads to corresponding buses (AL1 0 in load to AL1 0 low voltage)
    n.loads_t.p_set =  n.loads_t.p_set.rename(columns=n.loads.bus.to_dict())
    
    # low voltage load
    lv_buses = n.buses.query("carrier == 'low voltage'").index
    lv_load = n.loads_t.p_set[lv_buses]
    
    # electricity consumption of residential heat techs in low_voltage bus of links
    rh_elec_links = n.links.query("carrier in @rh_techs_elec").index
    rh_techs_mapping = n.links.loc[rh_elec_links, "bus0"].to_dict()
    rh_techs_consume = n.links_t.p0[rh_elec_links]
    rh_techs_consume = rh_techs_consume.rename(columns=rh_techs_mapping)
    rh_techs_consume = rh_techs_consume.groupby(level=0, axis=1).sum()
    
    # electricity consumption of land transport EV in low_voltage bus in links
    ev_charge_links = n.links.query("carrier in @ev_tech_charge").index
    ev_charge_mapping = n.links.loc[ev_charge_links, "bus0"].to_dict()
    ev_tech_consume = n.links_t.p0[ev_charge_links]
    ev_tech_consume = ev_tech_consume.rename(columns=ev_charge_mapping)
    ev_tech_consume = ev_tech_consume.groupby(level=0, axis=1).sum()
    
    # electricity prosumption of land transport EV to low_voltage bus in links
    ev_discharge_links = n.links.query("carrier in @ev_tech_discharge").index
    ev_discharge_mapping = n.links.loc[ev_discharge_links, "bus1"].to_dict()
    ev_tech_prosume = n.links_t.p1[ev_discharge_links]
    ev_tech_prosume = ev_tech_prosume.rename(columns=ev_discharge_mapping)
    ev_tech_prosume = ev_tech_prosume.groupby(level=0, axis=1).sum()
    
    # electricity prosumption of micro CHP to low_voltage bus in links
    rh_mCHP_links = n.links.query("carrier in @rh_techs_mCHP").index
    rh_mCHP_lv_mapping = n.links.loc[rh_mCHP_links, "bus1"].to_dict()
    rh_mCHP_prosume = n.links_t.p1[rh_mCHP_links]
    rh_mCHP_prosume = rh_mCHP_prosume.rename(columns=rh_mCHP_lv_mapping)
    rh_mCHP_prosume = rh_mCHP_prosume.groupby(level=0, axis=1).sum()
    
    # gas consumption of micro CHP in EU_gas bus in links
    rh_mCHP_gas_mapping = n.links.loc[rh_mCHP_links, "bus0"].to_dict()
    rh_mCHP_consume = n.links_t.p0[rh_mCHP_links]
    rh_mCHP_consume = rh_mCHP_consume.rename(columns=rh_mCHP_gas_mapping)
    rh_mCHP_consume = rh_mCHP_consume.groupby(level=0, axis=1).sum()
    
    # gas consumption of gas boilers in EU_gas bus in links
    rh_gas_links = n.links.query("carrier in @rh_techs_gas").index
    if not rh_gas_links.empty:   
        rh_gas_mapping = n.links.loc[rh_gas_links, "bus0"].to_dict()
        rh_gas_consume = n.links_t.p0[rh_gas_links]
        rh_gas_consume = rh_gas_consume.rename(columns=rh_gas_mapping)
        rh_gas_consume = rh_gas_consume.groupby(level=0, axis=1).sum()
    else:
        rh_gas_consume = pd.DataFrame()
    
    # total electricity consumption in low_voltage nodes
    total_elec_load = (lv_load + rh_techs_consume + ev_tech_consume + ev_tech_prosume + rh_mCHP_prosume).multiply(n.snapshot_weightings.stores, axis=0)
    
    # total gas consumption in EU_gas node
    if not rh_gas_links.empty: 
        total_gas_load = (rh_mCHP_consume + rh_gas_consume).multiply(n.snapshot_weightings.stores, axis=0)
    else:
        total_gas_load = rh_mCHP_consume.multiply(n.snapshot_weightings.stores, axis=0)
    
    # total electricity cost
    total_cost = total_elec_load.multiply(n.buses_t.marginal_price[total_elec_load.columns]).sum()
    total_cost.index = [x[:2] for x in total_cost.index]
    total_cost_country = total_cost.groupby(level=0).sum()
        
    # electricity bill per household [EUR/household] (households given in thousands)
    elec_bills_household = total_cost_country / (households*1e3)
    
    return elec_bills_household


def plot_electricity_bills(total_elec_bills):
    # sort countries by flexible bills' cost
    sorted_bills_household = total_elec_bills.sort_values(by="Efficient Heating", axis=1)

    # color codes for legend
    color_codes = {"Efficient Heating":"purple", "Efficient Green Heating":"limegreen", "Semi-Efficient Heating":"royalblue", "Non-efficient Heating":"#f4b609"}

    # plot as bar plot
    fig, ax = plt.subplots(figsize=(12,4))
    sorted_bills_household.T.plot.bar(ax=ax, width=0.7, color=color_codes)
    # define plot parameters
    ax.set_facecolor("white")
    ax.legend(loc="upper left", facecolor="white")
    ylabel = ax.set_ylabel("EUR/household")
    xlabel = ax.set_xlabel("countries")
    ax.set_title("Electricity price per country")
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    ax.grid(axis='y', linestyle='--', linewidth=0.5, color='gray')
    # save figure
    plt.savefig(PATH_PLOTS+"bill_per_household.png", bbox_inches='tight', dpi=600)

if __name__ == "__main__":
    # network parameters
    lineex = "v1.0"
    space_resolution = 48
    sector_opts = "Co2L0-1H-T-H-B-I"
    planning = 2030

    # move to submodules/pypsa-eur
    change_path_to_pypsa_eur()

    # load networks
    n_flex, n_flex_mod, n_retro_tes, n_rigid = load_networks(lineex, space_resolution, sector_opts, planning)
    networks = {"Efficient Heating": n_flex, 
                "Efficient Green Heating": n_retro_tes, 
                "Semi-Efficient Heating": n_flex_mod, 
                "Non-efficient Heating": n_rigid}

    # move to base directory
    change_path_to_base()

    # load country-wise households
    households = get_households()
    
    # calculate electricity bills per household for each network
    total_elec_bills = pd.DataFrame()
    for name, network in networks.items():
        # get electricity bills
        elec_bills_household = electricity_bills(network, households)
        # rename series name to scenario name
        elec_bills_household.name = name
        # concatenate current results
        total_elec_bills = pd.concat([total_elec_bills, elec_bills_household.to_frame().T], axis=0)

    # plot and store electricity bills
    plot_electricity_bills(total_elec_bills)


        







