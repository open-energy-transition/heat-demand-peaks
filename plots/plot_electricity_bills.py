# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  Open Energy Transition gGmbH
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import sys
sys.path.append("../submodules/pypsa-eur")
import matplotlib.pyplot as plt
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
from _helpers import mock_snakemake, update_config_from_wildcards, load_network, \
                     change_path_to_pypsa_eur, change_path_to_base, \
                     CO2L_LIMITS, LINE_LIMITS, BAU_HORIZON


def get_households():
    household_filepath = "submodules/pypsa-eur/data/retro/households.csv"
    household_data = pd.read_csv(household_filepath)
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
    rh_techs_gasCHP = ['urban central gas CHP', 'urban central gas CHP CC']

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
    if ev_tech_consume.empty:
        ev_tech_consume = pd.DataFrame(0, index=lv_load.index, columns=lv_load.columns)
    
    # electricity prosumption of land transport EV to low_voltage bus in links
    ev_discharge_links = n.links.query("carrier in @ev_tech_discharge").index
    ev_discharge_mapping = n.links.loc[ev_discharge_links, "bus1"].to_dict()
    ev_tech_prosume = n.links_t.p1[ev_discharge_links]
    ev_tech_prosume = ev_tech_prosume.rename(columns=ev_discharge_mapping)
    ev_tech_prosume = ev_tech_prosume.T.groupby(level=0).sum().T
    if ev_tech_prosume.empty:
        ev_tech_prosume = pd.DataFrame(0, index=lv_load.index, columns=lv_load.columns)
    
    # electricity prosumption of micro CHP to low_voltage bus in links
    rh_mCHP_links = n.links.query("carrier in @rh_techs_mCHP").index
    rh_mCHP_lv_mapping = n.links.loc[rh_mCHP_links, "bus1"].to_dict()
    rh_mCHP_prosume = n.links_t.p1[rh_mCHP_links]
    rh_mCHP_prosume = rh_mCHP_prosume.rename(columns=rh_mCHP_lv_mapping)
    rh_mCHP_prosume = rh_mCHP_prosume.T.groupby(level=0).sum().T
    
    # gas consumption of micro CHP in EU_gas bus in links
    rh_mCHP_consume = n.links_t.p2[rh_mCHP_links].divide(n.links.loc[rh_mCHP_links, "efficiency2"], axis=1)
    
    # gas consumption by gas CHPs
    rh_gasCHP_links = n.links.query("carrier in @rh_techs_gasCHP").index
    rh_gasCHP_consume = n.links_t.p2[rh_gasCHP_links].divide(n.links.loc[rh_gasCHP_links, "efficiency2"], axis=1)

    # gas consumption of gas boilers in EU_gas bus in links
    rh_gas_links = n.links.query("carrier in @rh_techs_gas").index
    if not rh_gas_links.empty:
        rh_gas_consume = n.links_t.p1[rh_gas_links].divide(n.links.loc[rh_gas_links, "efficiency"], axis=1)
    else:
        rh_gas_consume = pd.DataFrame()
    
    # total electricity consumption in low_voltage nodes
    total_elec_load = (lv_load + rh_techs_consume + ev_tech_consume + ev_tech_prosume + rh_mCHP_prosume).multiply(n.snapshot_weightings.stores, axis=0)
    
    # total gas consumption in EU_gas node
    if not rh_gas_links.empty:
        rh_mCHP_consume.columns = [" ".join(x.split(" ")[0:4]) for x in rh_mCHP_consume.columns]
        rh_gas_consume.columns = [" ".join(x.split(" ")[0:4]) for x in rh_gas_consume.columns]
        rh_gasCHP_consume.columns = [" ".join(x.split(" ")[0:4]) for x in rh_gasCHP_consume.columns]
        # group gas consumption of CHP for each country (because we have gas CHP and gas CHP CC)
        rh_gasCHP_consume = rh_gasCHP_consume.groupby(rh_gasCHP_consume.columns, axis=1).sum()
        total_gas_load = rh_mCHP_consume.add(rh_gas_consume, fill_value=0).add(rh_gasCHP_consume, fill_value=0).multiply(n.snapshot_weightings.stores, axis=0)
    else:
        total_gas_load = rh_mCHP_consume.add(rh_gasCHP_consume, fill_value=0).multiply(n.snapshot_weightings.stores, axis=0)
    
     # total gas cost
    total_gas_cost = (-n.generators.loc["EU gas"].marginal_cost * total_gas_load).sum()
    total_gas_cost.index = [x[:2] for x in total_gas_cost.index]
    total_gas_cost_country = total_gas_cost.groupby(level=0).sum()

    # total electricity cost
    total_cost = total_elec_load.multiply(n.buses_t.marginal_price[total_elec_load.columns]).sum()
    total_cost.index = [x[:2] for x in total_cost.index]
    total_cost_country = total_cost.groupby(level=0).sum()

    # gas plus electricity costs
    total_cost_country += total_gas_cost_country
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
    color_codes = {"Optimal Renovation and Cost-Optimal Heating":"purple", 
                   "Optimal Renovation and Electric Heating":"limegreen", 
                   "Limited Renovation and Cost-Optimal Heating":"royalblue", 
                   "No Renovation and Electric Heating":"#f4b609",
                   "BAU": "grey"}

    # plot as bar plot
    fig, ax = plt.subplots(figsize=(7,3))
    sorted_df_prices.T.plot.bar(ax=ax, width=0.7, color=color_codes)
    # define plot parameters
    ax.set_facecolor("white")
    # modify the name for LR in legend for 2040 and 2050
    handles, labels = ax.get_legend_handles_labels()
    if planning_horizon in ["2040", "2050"]:
        labels = ["Limited Renovation and Electric Heating" if label == "Limited Renovation and Cost-Optimal Heating" else label for label in labels]
    ax.legend(handles, labels, loc="upper left", facecolor="white", fontsize='x-small')
    xlabel = ax.set_xlabel("countries")
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    ax.grid(axis='y', linestyle='--', linewidth=0.5, color='gray')

    max_price = sorted_df_prices.max().max()
    ax.set_ylim(0, max_price * 1.2) 

    if name == "bills":
        ax.set_title("Electricity bills")
        ylabel = ax.set_ylabel("EUR/household")
        ax.set_ylim([0,5000])
        plt.savefig(snakemake.output.figure_bills, bbox_inches='tight', dpi=600)
    elif name == "prices":
        ax.set_title("Energy price per country")
        ylabel = ax.set_ylabel("EUR/MWh")
        ax.set_ylim([0, 300])
        plt.savefig(snakemake.output.figure_price, bbox_inches='tight', dpi=600)


def plot_industry_opex(df):

    # color codes f
    color_codes = {
        "Coal":"#545454",
        "Methanol":"#468c8b",
        "Biomass":"#baa741",
        "Hydrogen":"#f073da",
        "Oil products": "#aaaaaa",
        "Gas": "#e05b09",
        "Heat": "#cc1f1f",
        "Electricity": "#110d63",
    }

    fig, ax = plt.subplots(figsize=(3,3))
    df.T.plot.bar(ax=ax, width=0.7, color=color_codes, stacked=True)
    ax.set_facecolor("white")
    handles, labels = ax.get_legend_handles_labels()
    if planning_horizon in ["2040", "2050"]:
        labels = ["LREH" if label == "LROH" else label for label in labels]
    ax.legend(handles[::-1], labels[::-1], loc=[1.05,0], ncol=1, facecolor="white", fontsize='x-small')
    ax.set_title("")
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')

    ylabel = ax.set_ylabel("operating expenses [bn EUR]")
    ax.set_ylim([0, 550])
    plt.savefig(snakemake.output.figure_opex, bbox_inches='tight', dpi=600)

def electricity_prices(network):
    n = network
    
    rh_techs_elec = ['residential rural ground heat pump',
                     'residential rural resistive heater', 
                     'residential urban decentral air heat pump',
                     'residential urban decentral resistive heater',
                     'urban central resistive heater',
                     'urban central air heat pump',
                     'residential rural air heat pump']
    rh_techs_mCHP = ['residential rural micro gas CHP', 'residential urban decentral micro gas CHP']
    rh_techs_gasCHP = ['urban central gas CHP', 'urban central gas CHP CC']

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
    if ev_tech_consume.empty:
        ev_tech_consume = pd.DataFrame(0, index=lv_load.index, columns=lv_load.columns)
    
    # electricity prosumption of land transport EV to low_voltage bus in links
    ev_discharge_links = n.links.query("carrier in @ev_tech_discharge").index
    ev_discharge_mapping = n.links.loc[ev_discharge_links, "bus1"].to_dict()
    ev_tech_prosume = n.links_t.p1[ev_discharge_links]
    ev_tech_prosume = ev_tech_prosume.rename(columns=ev_discharge_mapping)
    ev_tech_prosume = ev_tech_prosume.T.groupby(level=0).sum().T
    if ev_tech_prosume.empty:
        ev_tech_prosume = pd.DataFrame(0, index=lv_load.index, columns=lv_load.columns)
    
    # electricity prosumption of micro CHP to low_voltage bus in links
    rh_mCHP_links = n.links.query("carrier in @rh_techs_mCHP").index
    rh_mCHP_lv_mapping = n.links.loc[rh_mCHP_links, "bus1"].to_dict()
    rh_mCHP_prosume = n.links_t.p1[rh_mCHP_links]
    rh_mCHP_prosume = rh_mCHP_prosume.rename(columns=rh_mCHP_lv_mapping)
    rh_mCHP_prosume = rh_mCHP_prosume.T.groupby(level=0).sum().T
    
    # total electricity consumption in low_voltage nodes
    total_elec_load = (lv_load + rh_techs_consume + ev_tech_consume + ev_tech_prosume + rh_mCHP_prosume).multiply(n.snapshot_weightings.stores, axis=0)

    # total electricity cost per country
    total_cost = total_elec_load.multiply(n.buses_t.marginal_price[total_elec_load.columns]).sum()
    total_cost.index = [x[:2] for x in total_cost.index]
    prices = total_cost.groupby(level=0).sum()

    # total electricity consumption by countries
    total_elec_load.columns = [x[:2] for x in total_elec_load.columns]
    total_load_country = total_elec_load.groupby(level=0, axis=1).sum().sum()

    # electricity bill per MWh [EUR/MWh]
    energy_price_MWh = prices / total_load_country

    # drop EU
    if "EU" in energy_price_MWh.index:
        energy_price_MWh.drop("EU", axis=0, inplace=True)
    
    return energy_price_MWh

def calc_opex(network, scenarioname):
    mapping = network.loads.bus.map(network.buses.country)

    query = "carrier == 'industry electricity'"
    load = network.loads.query(query).groupby(mapping).sum().p_set
    per_MWh = network.buses_t.marginal_price[network.buses.query("carrier == 'AC'").index].T.groupby(network.buses.country).mean().T
    elec = (load * per_MWh).sum().sum()

    query = "carrier == 'solid biomass for industry'"
    load = network.loads.query(query).groupby(mapping).sum().p_set
    per_MWh = network.buses_t.marginal_price[network.buses.query("carrier == 'solid biomass'").index].T.groupby(network.buses.country).mean().T
    solid_biomass = (load * per_MWh).sum().sum()

    query = "carrier == 'gas for industry'"
    load = network.loads.query(query).groupby(mapping).sum().p_set
    per_MWh = network.buses_t.marginal_price[network.buses.query("carrier == 'gas for industry'").index].T.groupby(network.buses.country).mean().T
    gas = (load * per_MWh).sum().values[0]

    query = "carrier == 'coal for industry'"
    load = network.loads.query(query).groupby(mapping).sum().p_set.sum()
    per_MWh = network.buses_t.marginal_price[network.buses.query("carrier == 'coal'").index].T.groupby(network.buses.country).mean().T
    coal = (load * per_MWh).sum().values[0]

    query = "carrier == 'H2 for industry'"
    load = network.loads.query(query).groupby(mapping).sum().p_set
    per_MWh = network.buses_t.marginal_price[network.buses.query("carrier == 'H2'").index].T.groupby(network.buses.country).mean().T
    H2 = (load * per_MWh).sum().sum()

    query = "carrier == 'H2 for industry'"
    load = network.loads.query(query).groupby(mapping).sum().p_set
    per_MWh = network.buses_t.marginal_price[network.buses.query("carrier == 'H2'").index].T.groupby(network.buses.country).mean().T
    H2 = (load * per_MWh).sum().sum()

    query = "carrier in ['industry methanol', 'shipping methanol']"
    load = network.loads.query(query).p_set
    per_MWh = network.buses_t.marginal_price[network.buses.query("carrier in ['industry methanol', 'shipping methanol']").index]
    methanol = (load * per_MWh).sum().sum()

    query = "carrier in ['shipping oil', 'naphtha for industry']"
    load = network.loads.query(query).p_set
    per_MWh = network.buses_t.marginal_price[network.buses.query("carrier in ['shipping oil', 'naphtha for industry']").index]
    oil = (load * per_MWh).sum().sum()

    query = "carrier == 'low-temperature heat for industry'"
    load = network.loads.query(query).set_index("bus").p_set
    per_MWh = network.buses_t.marginal_price[network.buses.query("carrier in ['urban central heat', 'services urban decentral heat']").index].clip(lower=0)
    heat = (load * per_MWh).fillna(0).sum().sum()

    opex = pd.DataFrame(
        index=["Electricity", "Biomass", "Gas", "Coal", "Hydrogen", "Methanol", "Oil products", "Heat"],
        columns=[scenarioname],
        data=[elec, solid_biomass, gas, coal, H2, methanol, oil, heat]
    )/1e9
    return (
        opex.loc[["Coal", "Methanol", "Hydrogen", "Oil products", "Biomass", "Heat", "Gas", "Electricity"]]
    )


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
    co2l_limits = CO2L_LIMITS
    line_limits = LINE_LIMITS
    clusters = config["plotting"]["clusters"]
    planning_horizon = config["plotting"]["planning_horizon"]
    opts = config["plotting"]["sector_opts"]
    lineex = line_limits[planning_horizon]
    sector_opts = f"Co2L{co2l_limits[planning_horizon]}-{opts}"

    # move to submodules/pypsa-eur
    change_path_to_pypsa_eur()

    # define scenario namings
    if planning_horizon == BAU_HORIZON:
        scenarios = {"BAU": "BAU"}
        scenario_abbrev = scenarios
    else:
        scenarios = {"flexible": "Optimal Renovation and Cost-Optimal Heating",
                     "retro_tes": "Optimal Renovation and Electric Heating",
                     "flexible-moderate": "Limited Renovation and Cost-Optimal Heating",
                     "rigid": "No Renovation and Electric Heating"}
        scenario_abbrev = {"Optimal Renovation and Cost-Optimal Heating": "OROH",
                     "Optimal Renovation and Electric Heating": "OREH",
                     "Limited Renovation and Cost-Optimal Heating": "LROH",
                     "No Renovation and Electric Heating": "NREH"}

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
    industry_opex = pd.DataFrame()
    for name, network in networks.items():
        if network is None:
            # Skip further computation for this scenario if network is not loaded
            print(f"Network is not found for scenario '{scenario}', planning year '{planning_horizon}'. Skipping...")
            continue
        # get electricity bills
        elec_bills_household = electricity_bills(network, households).rename(name)
        # get electricity prices
        elec_bills_MWh = electricity_prices(network).rename(name)
        # OPEX for industry
        opex = calc_opex(network, scenario_abbrev[name])
        # rename series name to scenario name
        elec_bills_household.name = name
        elec_bills_MWh.name = name
        # concatenate current results
        total_elec_bills = pd.concat([total_elec_bills, elec_bills_household.to_frame().T], axis=0)
        total_elec_prices = pd.concat([total_elec_prices, elec_bills_MWh.to_frame().T], axis=0)
        industry_opex = pd.concat([industry_opex, opex], axis=1)

    # plot and store electricity bills
    if not total_elec_bills.empty:
        plot_electricity_cost(total_elec_bills, "bills")
    
    # plot and store electricity prices
    if not total_elec_prices.empty:
        plot_electricity_cost(total_elec_prices, "prices")

    # plot and store industry opex
    if not industry_opex.empty:
        plot_industry_opex(industry_opex)

