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


def get_households():
    household_data = pd.read_csv("data/Households.csv")
    households = household_data.set_index("Country")[ "Households (thousands)"]
    households.index.name=""
    return households


def electricity_bills(network, households):
    n = network
    
    rh_techs_elec = ['residential rural ground heat pump',
                     'residential rural resistive heater', 
                     'residential urban decentral air heat pump',
                     'residential urban decentral resistive heater',
                     'urban central resistive heater',
                     'urban central air heat pump',
                     'residential rural air heat pump']
    rh_techs_gas = ['residential rural gas boiler', 'residential urban decentral gas boiler', 'urban central gas boiler']
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
    rh_techs_consume = rh_techs_consume.T.groupby(level=0).sum().T
    
    # electricity consumption of land transport EV in low_voltage bus in links
    ev_charge_links = n.links.query("carrier in @ev_tech_charge").index
    ev_charge_mapping = n.links.loc[ev_charge_links, "bus0"].to_dict()
    ev_tech_consume = n.links_t.p0[ev_charge_links]
    ev_tech_consume = ev_tech_consume.rename(columns=ev_charge_mapping)
    ev_tech_consume = ev_tech_consume.T.groupby(level=0).sum().T
    
    # electricity prosumption of land transport EV to low_voltage bus in links
    ev_discharge_links = n.links.query("carrier in @ev_tech_discharge").index
    ev_discharge_mapping = n.links.loc[ev_discharge_links, "bus1"].to_dict()
    ev_tech_prosume = n.links_t.p1[ev_discharge_links]
    ev_tech_prosume = ev_tech_prosume.rename(columns=ev_discharge_mapping)
    ev_tech_prosume = ev_tech_prosume.T.groupby(level=0).sum().T
    
    # electricity prosumption of micro CHP to low_voltage bus in links
    rh_mCHP_links = n.links.query("carrier in @rh_techs_mCHP").index
    rh_mCHP_lv_mapping = n.links.loc[rh_mCHP_links, "bus1"].to_dict()
    rh_mCHP_prosume = n.links_t.p1[rh_mCHP_links]
    rh_mCHP_prosume = rh_mCHP_prosume.rename(columns=rh_mCHP_lv_mapping)
    rh_mCHP_prosume = rh_mCHP_prosume.T.groupby(level=0).sum().T
    
    # gas consumption of micro CHP in EU_gas bus in links
    rh_mCHP_gas_mapping = n.links.loc[rh_mCHP_links, "bus0"].to_dict()
    rh_mCHP_consume = n.links_t.p0[rh_mCHP_links]
    rh_mCHP_consume = rh_mCHP_consume.rename(columns=rh_mCHP_gas_mapping)
    rh_mCHP_consume = rh_mCHP_consume.T.groupby(level=0).sum().T
    
    # gas consumption of gas boilers in EU_gas bus in links
    rh_gas_links = n.links.query("carrier in @rh_techs_gas").index
    if not rh_gas_links.empty:   
        rh_gas_mapping = n.links.loc[rh_gas_links, "bus0"].to_dict()
        rh_gas_consume = n.links_t.p0[rh_gas_links]
        rh_gas_consume = rh_gas_consume.rename(columns=rh_gas_mapping)
        rh_gas_consume = rh_gas_consume.T.groupby(level=0).sum().T
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


def plot_electricity_cost(df_prices, name):
    # Check if name is one of the allowed values
    allowed_names = ["bills", "prices"]
    if name not in allowed_names:
        raise ValueError("name must be one of {}".format(allowed_names))

    # sort countries by flexible's electricity cost
    sorted_df_prices = df_prices.sort_values(by=df_prices.index[0], axis=1)

    # color codes for legend
    color_codes = {"Optimal Renovation and Heating":"purple", 
                   "Optimal Renovation and Green Heating":"limegreen", 
                   "Limited Renovation and Optimal Heating":"royalblue", 
                   "No Renovation and Optimal Heating":"#f4b609"}

    # plot as bar plot
    fig, ax = plt.subplots(figsize=(7,3))
    sorted_df_prices.T.plot.bar(ax=ax, width=0.7, color=color_codes)
    # define plot parameters
    ax.set_facecolor("white")
    ax.legend(loc="upper left", facecolor="white", fontsize='x-small')
    xlabel = ax.set_xlabel("countries")
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    ax.grid(axis='y', linestyle='--', linewidth=0.5, color='gray')

    max_price = sorted_df_prices.max().max()
    ax.set_ylim(0, max_price * 1.2) 

    if name == "bills":
        ax.set_title("Electricity bills")
        ylabel = ax.set_ylabel("EUR/household")
        plt.savefig(snakemake.output.figure_bills, bbox_inches='tight', dpi=600)
    elif name == "prices":
        ax.set_title("Energy price per country")
        ylabel = ax.set_ylabel("EUR/MWh")
        plt.savefig(snakemake.output.figure_price, bbox_inches='tight', dpi=600)


def electricity_prices(network, households):
    n = network
    
    # rename AL1 0 load to AL1 0 low voltage
    n.loads_t.p_set =  n.loads_t.p_set.rename(columns=n.loads.bus.to_dict())
    load_cols = n.loads_t.p_set.columns
    
    # sum total electricity price per country in EUR using loads and marginal prices ar buses
    prices = n.loads_t.p_set.multiply(n.snapshot_weightings.stores, axis=0).multiply(n.buses_t.marginal_price[load_cols]).sum()
    
    # rename indexes to 2-letter country code
    prices.index = [x[:2] for x in prices.index]
    
    # group by country
    prices = prices.groupby(level=0).sum()

    # total load
    total_load = n.loads_t.p_set.multiply(n.snapshot_weightings.stores, axis=0).sum()
    total_load.index = [x[:2] for x in total_load.index]
    total_load_country = total_load.groupby(level=0).sum()

    # electricity bill per MWh [EUR/MWh]
    energy_price_MWh = prices / total_load_country

    # drop EU
    if "EU" in energy_price_MWh.index:
        energy_price_MWh.drop("EU", axis=0, inplace=True)
    
    return energy_price_MWh



if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "plot_electricity_bill", 
            clusters="48",
            planning_horizon="2030",
        )
    # update config based on wildcards
    config = update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    # network parameters
    co2l_limits = {"2030":"0.45", "2040":"0.1", "2050":"0.0"}
    line_limits = {"2030":"v1.15", "2040":"v1.3", "2050":"v1.5"}
    clusters = config["plotting"]["clusters"]
    planning_horizon = config["plotting"]["planning_horizon"]
    time_resolution = config["plotting"]["time_resolution"]
    lineex = line_limits[planning_horizon]
    sector_opts = f"Co2L{co2l_limits[planning_horizon]}-{time_resolution}-T-H-B-I"

    # move to submodules/pypsa-eur
    change_path_to_pypsa_eur()

    # define scenario namings
    scenarios = {"flexible": "Optimal Renovation and Heating", 
                 "retro_tes": "Optimal Renovation and Green Heating", 
                 "flexible-moderate": "Limited Renovation and Optimal Heating", 
                 "rigid": "No Renovation and Optimal Heating"}

    # load networks
    networks = {}
    for scenario, nice_name in scenarios.items():
        n = load_network(lineex, clusters, sector_opts, planning_horizon, scenario)
        networks[nice_name] = n

    # move to base directory
    change_path_to_base()

    # load country-wise households
    households = get_households()
    
    # calculate electricity bills per household for each network and electricity prices
    total_elec_bills = pd.DataFrame()
    total_elec_prices = pd.DataFrame()
    for name, network in networks.items():
        if network is None:
            # Skip further computation for this scenario if network is not loaded
            print(f"Network is not found for scenario '{scenario}', planning year '{planning_horizon}', and time resolution of '{time_resolution}'. Skipping...")
            continue
        # get electricity bills
        elec_bills_household = electricity_bills(network, households)
        # get electricity prices
        elec_bills_MWh = electricity_prices(network, households)
        # rename series name to scenario name
        elec_bills_household.name = name
        elec_bills_MWh.name = name
        # concatenate current results
        total_elec_bills = pd.concat([total_elec_bills, elec_bills_household.to_frame().T], axis=0)
        total_elec_prices = pd.concat([total_elec_prices, elec_bills_MWh.to_frame().T], axis=0)

    # plot and store electricity bills
    if not total_elec_bills.empty:
        plot_electricity_cost(total_elec_bills, "bills")
    
    # plot and store electricity prices
    if not total_elec_prices.empty:
        plot_electricity_cost(total_elec_prices, "prices")

