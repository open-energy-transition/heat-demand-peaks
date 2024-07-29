# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  Open Energy Transition gGmbH
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import os
import sys
sys.path.append("../submodules/pypsa-eur")
import pypsa
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
from _helpers import mock_snakemake, update_config_from_wildcards, load_network, \
                     change_path_to_pypsa_eur, change_path_to_base, \
                     LINE_LIMITS, CO2L_LIMITS, BAU_HORIZON


def add_new_carriers(n):
    n.generators['new_carrier'] = (
        n.generators.apply(
            lambda b: "solar thermal" if "solar thermal" in b.carrier else b.carrier, axis=1
        )
    )

    n.links['new_carrier'] = (
        n.links.apply(
            lambda b: "ground heat pump" if "ground heat pump" in b.carrier else (
                "water tank" if "water tank" in b.carrier else (
                    "resistive heater" if "resistive heater" in b.carrier else (
                        "gas boiler" if "gas boiler" in b.carrier else (
                            "micro gas CHP" if "micro gas CHP" in b.carrier else (
                                "air heat pump" if "air heat pump" in b.carrier else b.carrier
                            )
                        )
                    )
                )
            ), axis=1
        )
    )


def get_elec_consumption_for_heat(n):
    supp_techs = ["ground heat pump", "water tank", "resistive heater", "gas boiler", "micro gas CHP", "air heat pump", "solar thermal"]
    carriers = ['residential rural heat', 'services rural heat', 'residential urban decentral heat', 'services urban decentral heat', 'urban central heat']
    elec_buses = n.buses.query("carrier == 'low voltage'").index
    heat_buses = n.buses.query("carrier in @carriers").index


    supply = pd.DataFrame(index=n.snapshots, columns=["gas boiler"], data=0)
    if "gas boiler" in n.links.new_carrier.unique():
        supply[["gas boiler"]] -= (
            n.links_t.p1.T.groupby(n.links.query("bus1 in @heat_buses").new_carrier).sum().T[["gas boiler"]]
        )
    else:
        supply["gas boiler"] = 0

    elec_demand_f_heating = pd.DataFrame(index=n.snapshots, columns=["micro gas CHP", "air heat pump", "ground heat pump", "resistive heater"], data=0)
    elec_demand_f_heating[["air heat pump", "ground heat pump", "resistive heater"]] = n.links_t.p0.T.groupby(n.links.query("new_carrier in @supp_techs and bus0 in @elec_buses").new_carrier).sum().T
    elec_demand_f_heating[["micro gas CHP"]] = n.links_t.p1.T.groupby(n.links.query("new_carrier in @supp_techs and bus1 in @elec_buses").new_carrier).sum().T.clip(0)
    elec_demand_f_heating["gas boiler"] = supply["gas boiler"]

    return elec_demand_f_heating


def plot_elec_consumption_for_heat(dict_elec, horizon):
    # set heights for each subplots
    if "BAU" in dict_elec.keys():
        heights = [1.4]
    else:
        heights = [1.4] * 4
    fig = plt.figure(figsize=(6.4, sum(heights)))
    # fig, axes = plt.subplots(len(dict_elec.keys()),1) #, figsize=(11,5))
    # axes = [ax1, ax2, ax3, ax4]
    gs = gridspec.GridSpec(len(heights), 1, height_ratios=heights)
    axes = [fig.add_subplot(gs[i]) for i in range(len(heights))]
    if "BAU" in dict_elec.keys():
        axes = axes * 4
    i=0
    for name, elec_demand_f_heating in dict_elec.items():
        ax = axes[i]
        ax.set_facecolor("whitesmoke")

        print("Hard coded coordinates on x-axis, selected for 3H granularity")
        where = [359, 695]

        elec_demand_f_heating = elec_demand_f_heating.reset_index()[["ground heat pump", "air heat pump", "resistive heater", "gas boiler"]]
        colors = ["#2fb537", "#48f74f", "#d8f9b8", "#db6a25"]

        (elec_demand_f_heating/1e3).iloc[where[0] : where[1]].plot(
            kind='area', stacked=True, color=colors, legend=False,
            linewidth=0, ax=ax
        )

        cumulative = (elec_demand_f_heating / 1e3).cumsum(axis=1)
        ax.plot(cumulative.iloc[where[0]:where[1]], color='black', linewidth=0.5)

        ax.set_xlabel("", fontsize=12)
        # get in y limits for all horizons
        y_limits = {'2020':1600, '2030':1100, '2040':900, '2050': 900}

        ax.set_ylim([1, y_limits[horizon]])

        if i < 3:
            ax.set_xticks([])
            ax.set_xlabel("")
        if i == 3 or "BAU" in dict_elec.keys():
            ticks = [i for i in range(where[0], where[1], 48)]
            ax.set_xticks(ticks)  # Set the tick positions
            ticks = [
                str(n.snapshots[i]).split(" ")[0].split("-")[2]+"."+str(n.snapshots[i]).split(" ")[0].split("-")[1]
                for i in ticks
            ]
            ax.set_xticklabels(ticks, fontsize=10)  # Set the tick labels
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)

        ax.set_xlim([where[0], where[1]-1])
        # change name to LR for 2040 and 2050
        name = "Limited Renovation and Electric Heating" if name == "Limited Renovation and Cost-Optimal Heating" and planning_horizon in ["2040", "2050"] else name
        ax.set_title(name, fontsize=10)
        i+= 1

    handles1, _ = axes[0].get_legend_handles_labels()
    axes[0].legend(
        reversed(handles1[0:7]),
        [
            "Gas Boilers (gas)",
            "Resistive Heaters (electricity)",
            "Air-Sourced Heat Pumps (electricity)",
            "Ground-Sourced Heat Pumps (electricity)"
        ],
        loc=[1.02, -.2], fontsize=10
    )
    if len(dict_elec.keys()) == 1:
        ylabel = "Electricity and Gas \nConsumption to Cover \nHeating Demands [GW]"
        axes[0].set_ylabel(ylabel, fontsize=10)
    else:
        ylabel = "Electricity and Gas Consumption to \n Cover Heating Demands [GW]"
        axes[2].set_ylabel(ylabel, fontsize=10)
        axes[2].yaxis.set_label_coords(-0.1, 1)
    plt.subplots_adjust(hspace=0.3)
    plt.savefig(snakemake.output.figure, bbox_inches='tight', dpi=600)



if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "plot_electricity_for_heat", 
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
    else:
        scenarios = {"flexible": "Optimal Renovation and Cost-Optimal Heating", 
                     "retro_tes": "Optimal Renovation and Electric Heating", 
                     "flexible-moderate": "Limited Renovation and Cost-Optimal Heating", 
                     "rigid": "No Renovation and Electric Heating"}

    # load networks
    networks = {}
    for scenario, nice_name in scenarios.items():
        n = load_network(lineex, clusters, sector_opts, planning_horizon, scenario)
        if n is None:
            # Skip further computation for this scenario if network is not loaded
            print(f"Network is not found for scenario '{scenario}', planning year '{planning_horizon}'. Skipping...")
            continue
        
        add_new_carriers(n)
        networks[nice_name] = n

    # move to base directory
    change_path_to_base()
   
    total_elec_consumption_for_heat = {}
    for name, network in networks.items():
        if network is None:
            # Skip further computation for this scenario if network is not loaded
            print(f"Network is not found for scenario '{scenario}', planning year '{planning_horizon}'. Skipping...")
            continue
        elec_consumption_for_heat = get_elec_consumption_for_heat(network)
        total_elec_consumption_for_heat[name] = elec_consumption_for_heat
    
    # plot and store electricity prices
    plot_elec_consumption_for_heat(total_elec_consumption_for_heat, planning_horizon)

