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


def get_heat_data(n, with_flexibility=False):
    # get heat load
    heat_load = n.loads_t.p_set.filter(like="heat").multiply(n.snapshot_weightings.objective, axis=0).sum(axis=1)
    # get heat savings by retrofitting
    heat_savings = n.generators_t.p.filter(like="retrofitting").multiply(n.snapshot_weightings.objective, axis=0).sum(axis=1)
    # get heat flexibility
    if with_flexibility:
        heat_flex_carr = ['residential rural heat', 'residential urban decentral heat', 'urban central heat']
        heat_flexibility_techs = n.stores.query("carrier in @heat_flex_carr").index
        heat_flexibility = n.stores_t.p[heat_flexibility_techs].multiply(n.snapshot_weightings.objective, axis=0).sum(axis=1)
        
        # compute heat load after renovation and flexbility
        net_heat_load = heat_load - heat_savings - heat_flexibility
        # concatenate net heat load and heat_savings
        heat_data = pd.concat([net_heat_load, heat_flexibility + heat_savings], axis=1)
        heat_data.columns = ["Net heat demand", "Heat savings by renovation"]
    else:
        # compute heat load after renovation
        net_heat_load = heat_load - heat_savings
        # concatenate net heat load and heat_savings
        heat_data = pd.concat([net_heat_load, heat_savings], axis=1)
        heat_data.columns = ["Net heat demand", "Heat savings by renovation"]
    return heat_data


def get_flexibility(network):
    # heat flexibility
    heat_flex_carr = ['residential rural heat', 'residential urban decentral heat', 'urban central heat']
    heat_flexibility_techs = n.stores.query("carrier in @heat_flex_carr").index
    heat_flexibility = n.stores_t.p[heat_flexibility_techs].multiply(n.snapshot_weightings.objective, axis=0).sum(axis=1)
    # EV flexibility
    ev_flex_carr = ['Li ion']
    ev_flexibility_techs = n.stores.query("carrier in @ev_flex_carr").index
    ev_flexibility = n.stores_t.p[ev_flexibility_techs].multiply(n.snapshot_weightings.objective, axis=0).sum(axis=1)
    # concatenate flexibility
    flexibility = pd.concat([ev_flexibility, heat_flexibility], axis=1)
    flexibility.columns = ["EV", "heat"]
    return flexibility

def plot_elec_consumption_for_heat(dict_elec, full_year=False):
    # set heights for each subplots
    if "Baseline 2023" in dict_elec.keys():
        heights = [1.3]
    else:
        heights = [1.4] * 3
    fig = plt.figure(figsize=(6.4, sum(heights)))
    gs = gridspec.GridSpec(len(heights), 1, height_ratios=heights)
    axes = [fig.add_subplot(gs[i]) for i in range(len(heights))]
    if "Baseline 2023" in dict_elec.keys():
        axes = axes * 3
    i=0
    for name, heat_demand in dict_elec.items():
        ax = axes[i]
        ax.set_facecolor("whitesmoke")

        # heat saved ratio
        heat_saved_ratio = heat_demand.sum()["Heat savings by renovation"] / heat_demand.sum().sum()

        print("Hard coded coordinates on x-axis, selected for 3H granularity")
        if not full_year:
            where = [359, 695]
        else:
            where = [0, 365]
            heat_demand = heat_demand.resample('24H').max()
        
        heat_demand = heat_demand.reset_index(drop=True)
        colors = ["#a63f3f", "#e39191"]

        (heat_demand/1e3).iloc[where[0] : where[1]].plot(
            kind='area', stacked=True, color=colors, legend=False,
            linewidth=0, ax=ax
        )

        cumulative = (heat_demand / 1e3).cumsum(axis=1)
        ax.plot(cumulative.iloc[where[0]:where[1]], color='black', linewidth=0.5)

        # x position for arrow
        if not full_year:
            x_loc = 436
            y_shift_arrow = 0
            x_shift_arrow = -24
            x_shift_after = 33
            x_shift_before = -65
            x_shift_mid = -12
            y_after = 1550
            y_before = 2000
            y_mid = 1800
        else:
            x_loc = 90
            y_shift_arrow = -50
            x_shift_arrow = 0
            x_shift_after = 33
            x_shift_before = -15
            x_shift_mid = 20
            y_after = 746
            y_before = 1600
            y_mid = 1150

        if name == "Widespread Renovation":
            # Heating demand before renovation
            ax.annotate('heat demand before renovation', xy=(x_loc+x_shift_arrow, cumulative['Heat savings by renovation'][x_loc]+y_shift_arrow), xytext=(x_loc+x_shift_before, y_before),
                arrowprops=dict(facecolor='black', edgecolor='black', arrowstyle='->', linewidth=0.75), fontsize=9, color="black")
            # Heating demand after renovation
            ax.annotate('heat demand after renovation', xy=(x_loc, cumulative['Net heat demand'][x_loc]+y_shift_arrow), xytext=(x_loc+x_shift_after, y_after),
                arrowprops=dict(facecolor='black', edgecolor='black', arrowstyle='->', linewidth=0.75), fontsize=9, color="black")
        if not name == "Baseline 2023":
            # Heating savings
            mid_point = cumulative.sum(axis=1) / 2
            ax.annotate(f'saved heating demand ({100*heat_saved_ratio:.1f}%)', xy=(x_loc, mid_point[x_loc]+y_shift_arrow), xytext=(x_loc+x_shift_mid, y_mid),
                arrowprops=dict(facecolor='black', edgecolor='black', arrowstyle='->', linewidth=0.75), fontsize=9, color="#a63f3f")
            
        ax.set_xlabel("", fontsize=12)
        if full_year:
            ax.set_ylim([1,2000])
        else:
            ax.set_ylim([1,2200])

        if i < 2:
            ax.set_xticks([])
            ax.set_xlabel("")
        if i == 2 or "Baseline 2023" in dict_elec.keys():
            if not full_year:
                ticks = [i for i in range(where[0], where[1], 48)]
                ax.set_xticks(ticks)  # Set the tick positions
                ticks = [
                    str(n.snapshots[i]).split(" ")[0].split("-")[2]+"."+str(n.snapshots[i]).split(" ")[0].split("-")[1]
                    for i in ticks
                ]
                ax.set_xticklabels(ticks, fontsize=10)  # Set the tick labels
            else:
                snapshots_series = pd.Series(1, index=n.snapshots)
                first_day_each_month = snapshots_series.resample('MS').first().index
                days_diff = first_day_each_month.to_series().diff().dt.days
                days_diff.iloc[0] = 0
                cumulative_days = days_diff.cumsum()
                cumulative_days_list = cumulative_days.tolist()
                ax.set_xticks([int(x) for x in cumulative_days_list])
                ticks = [str(i).split(" ")[0].replace("-",".")[5:] for i in first_day_each_month]
                ax.set_xticklabels(ticks, fontsize=10, rotation=15)  # Set the tick labels
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)

        ax.set_xlim([where[0], where[1]-1])
        # change name to LR for 2040 and 2050
        name = "Limited Renovation & Electrification" if name == "Limited Renovation" and planning_horizon in ["2040", "2050"] else name
        ax.set_title(name, fontsize=10)
        i+= 1

    # handles1, _ = axes[0].get_legend_handles_labels()
    # axes[0].legend(
    #     reversed(handles1[0:7]),
    #     [
    #         "Heat savings by renovation",
    #         "Net heat demand",
    #     ],
    #     loc=[1.02, -.2], fontsize=10
    # )
    if len(dict_elec.keys()) == 1:
        ylabel = "Heat demand [$GW_{th}$]"
        axes[0].set_ylabel(ylabel, fontsize=10)
    else:
        ylabel = "Heat demand [$GW_{th}$]"
        axes[1].set_ylabel(ylabel, fontsize=10)
        # axes[1].yaxis.set_label_coords(-0.1, 1)
    plt.subplots_adjust(hspace=0.3)
    if not full_year:
        plt.savefig(snakemake.output.figure, bbox_inches='tight', dpi=600)
    else:
        plt.savefig(snakemake.output.figure_full, bbox_inches='tight', dpi=600)


def plot_flexibility(dict_elec):
    # set heights for each subplots
    if "Baseline 2023" in dict_elec.keys():
        heights = [1.2]
    else:
        heights = [1.4] * 3
    fig = plt.figure(figsize=(6.4, sum(heights)))
    gs = gridspec.GridSpec(len(heights), 1, height_ratios=heights)
    axes = [fig.add_subplot(gs[i]) for i in range(len(heights))]
    if "Baseline 2023" in dict_elec.keys():
        axes = axes * 3
    i=0
    colors = ["red", "black"]
    for name, heat_demand in dict_elec.items():
        ax = axes[i]
        ax.set_facecolor("whitesmoke")

        print("Hard coded coordinates on x-axis, selected for 3H granularity")
        where = [0, 8759]
        heat_demand = heat_demand.reset_index(drop=True)

        ax.plot((heat_demand["EV"]/1e3).iloc[where[0]:where[1]], color="red", linewidth=0.5, label="EV flexibility")
        ax.plot((heat_demand["heat"]/1e3).iloc[where[0]:where[1]], color="black", linewidth=0.5, label="Heat flexibility")

        ax.set_xlabel("", fontsize=12)
        ax.set_ylim([-550,550])
        ax.set_yticks(np.arange(-500, 600, 250))

        if i < 2:
            ax.set_xticks([])
            ax.set_xlabel("")
        if i == 2 or "Baseline 2023" in dict_elec.keys():
            snapshots_series = pd.Series(1, index=n.snapshots)
            first_day_each_month = snapshots_series.resample('MS').first().index
            hours_diff = first_day_each_month.to_series().diff().dt.total_seconds() / 3600
            hours_diff.iloc[0] = 0
            cumulative_hours = hours_diff.cumsum()
            cumulative_hours_list = cumulative_hours.tolist()
            ax.set_xticks([int(x) for x in cumulative_hours_list])
            ticks = [str(i).split(" ")[0].replace("-",".")[5:] for i in first_day_each_month]
            ax.set_xticklabels(ticks, fontsize=10, rotation=15)  # Set the tick labels
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)

        ax.set_xlim([where[0], where[1]-1])
        # change name to LR for 2040 and 2050
        name = "Limited Renovation & Electrification" if name == "Limited Renovation" and planning_horizon in ["2040", "2050"] else name
        ax.set_title(name, fontsize=10)
        i+= 1

    axes[0].legend(loc=[1.02, -.2], fontsize=10)
    if len(dict_elec.keys()) == 1:
        ylabel = "Flexibility [GW]"
        axes[0].set_ylabel(ylabel, fontsize=10)
    else:
        ylabel = "Flexibility [GW]"
        axes[1].set_ylabel(ylabel, fontsize=10)
    plt.subplots_adjust(hspace=0.3)
    plt.savefig(snakemake.output.figure_flexibility, bbox_inches='tight', dpi=600)



if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "plot_heat_saving", 
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
        scenarios = {"BAU": "Baseline 2023"}
    else:
        scenarios = {"flexible": "Widespread Renovation",
                     "retro_tes": "Widespread Renovation & Electrification",
                     "flexible-moderate": "Limited Renovation"} #, 
                     #"rigid": "No Renovation and Electric Heating"}

    # load networks
    networks = {}
    for scenario, nice_name in scenarios.items():
        n = load_network(lineex, clusters, sector_opts, planning_horizon, scenario)
        if n is None:
            # Skip further computation for this scenario if network is not loaded
            print(f"Network is not found for scenario '{scenario}', planning year '{planning_horizon}'. Skipping...")
            continue
        
        networks[nice_name] = n

    # move to base directory
    change_path_to_base()
   
    total_heat_data = {}
    total_heat_data_with_flex = {}
    total_flexibility = {}
    for name, network in networks.items():
        if network is None:
            # Skip further computation for this scenario if network is not loaded
            print(f"Network is not found for scenario '{scenario}', planning year '{planning_horizon}'. Skipping...")
            continue
        heat_data= get_heat_data(network)
        total_heat_data[name] = heat_data
        heat_data_with_flex = get_heat_data(network, with_flexibility=True)
        total_heat_data_with_flex[name] = heat_data_with_flex
        # get flexibility data
        flexibility = get_flexibility(network)
        total_flexibility[name] = flexibility
    
    # plot heat demand data for short period
    plot_elec_consumption_for_heat(total_heat_data_with_flex, full_year=False)
    # plot heat demand data for full year
    plot_elec_consumption_for_heat(total_heat_data_with_flex, full_year=True)
    # plot heat flexibility
    plot_flexibility(total_flexibility)