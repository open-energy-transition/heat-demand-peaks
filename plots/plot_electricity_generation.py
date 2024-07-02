# SPDX-FileCopyrightText:  Open Energy Transition gGmbH
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import sys
sys.path.append("../submodules/pypsa-eur")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
from _helpers import mock_snakemake, update_config_from_wildcards, load_network, \
                     change_path_to_pypsa_eur, change_path_to_base, \
                     LINE_LIMITS, CO2L_LIMITS, BAU_HORIZON


def autopct_format_inner(value, threshold=1):
    return f'{value:.1f}%' if value >= threshold else ''


def autopct_format_outer(value, threshold=3):
    return f'{value:.1f}%' if value >= threshold else ''


def filter_labels(labels, autopct_values):
    return [label if autopct_values[i] else None for i, label in enumerate(labels)]


def plot_pies(ax, elec_mix_array):
    size = 0.4
    vals = np.array(elec_mix_array)

    cmap = plt.colormaps["tab20c"]
    outer_colors = ["#ff8c00", "#0fa101", "#A18181"]
    inner_colors = cmap([1, 2, 5, 6, 9, 10])
    inner_colors = [
        "#ff8c00", "#ff8c00", "#ff8c00", "#ff8c00", "#ff8c00", # all Nuclear
        "#6895dd", "#235ebc", "#3dbfb0", "#f9d002", '#ffea80', # different VRES
        "#545454", "#db6a25", "#db6a25", "#db6a25", "#db6a25" # coal, rest gas
    ]

    outer_labels = ["Nuclear", "VRES", "Fossil"]

    # threshold to show outer label
    threshold = 2.0
    # get percentage of outer layer  
    autopct_values = [100*value/vals.sum(axis=1).sum() for value in vals.sum(axis=1)]
    # filter out percentages lower than threshold
    autopct_values = [value if value >= threshold else None for value in autopct_values]
    # outer labels that has value >= threshold 
    valid_outer_labels = filter_labels(outer_labels, autopct_values)

    _, texts, autotexts = ax.pie(vals.sum(axis=1), radius=1, colors=outer_colors,
        wedgeprops=dict(width=size, edgecolor='w', linewidth=0.2), 
        autopct=autopct_format_outer, textprops={'fontsize': 4}, pctdistance=1.5,
        labels=valid_outer_labels)
    
    # Adjust the position of autopct labels
    for autotext, label in zip(autotexts, texts):
        x, y = label.get_position()  # Get position of corresponding wedge label
        autotext.set_position((x , y - 0.11))  # Set position of autopct label below the wedge label
        autotext.set_fontsize(3)
        align = "right" if x < 0 else "left"
        autotext.set_horizontalalignment(align)

    all_vals = vals.flatten()
    ax.pie(
        all_vals[all_vals!=0], radius=1.15-size,
        colors=[inner_colors[i] for i in range(len(inner_colors)) if all_vals[i] != 0],
        wedgeprops=dict(width=size, edgecolor='w', linewidth=0.2),
        autopct=autopct_format_inner, textprops={'fontsize': 3}, pctdistance=0.7
    )
    
    # Calculate total generated electricity
    total_electricity = np.sum(vals)
    # Add total generated electricity to the center of the pie chart
    ax.text(0, 0, f"{total_electricity/1e6:.2f}"+"\nTWh", ha='center', va='center', fontsize=4)
    ax.set(aspect="equal")


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "plot_electricity_generation", 
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
        scenarios = {"flexible": "OROH", 
                     "retro_tes": "OREH", 
                     "flexible-moderate": "LROH", 
                     "rigid": "NREH"}

    # define figure
    count_scenarios = len(scenarios.keys())
    fig, axes = plt.subplots(1, count_scenarios, figsize=(1.6*count_scenarios, 1.9)) #each scenario one axis

    # load networks
    i = 0
    for scenario, nice_name in scenarios.items():
        ax = axes[i] if isinstance(axes, np.ndarray) else axes
        i += 1

        n = load_network(lineex, clusters, sector_opts, planning_horizon, scenario)

        if n is None:
            # Skip further computation for this scenario if network is not loaded
            print(f"Network is not found for scenario '{scenario}', planning year '{planning_horizon}'. Skipping...")
            continue

        # definde elec buses
        elec = ["AC", "low voltage"]
        elec = n.buses.query("carrier in @elec").index

        # elec mix from generators
        gens = n.generators.query("bus in @elec").index
        elec_mix = n.generators_t.p[gens].multiply(n.snapshot_weightings.objective,axis=0).T.groupby(n.generators.carrier).sum().T.sum()

        # elec mix storage units
        elec_mix_hydro = n.storage_units_t.p.multiply(n.snapshot_weightings.objective,axis=0).T.groupby(n.storage_units.carrier).sum().T.sum()

        # elec mix links (gas)
        buses = n.stores.query("carrier == 'gas'").bus
        links = n.links.query("bus0 in @buses and bus1 in @elec").index
        elec_mix_gas = -1*n.links_t.p1[links].multiply(n.snapshot_weightings.objective,axis=0).T.groupby(n.links.carrier).sum().T.sum()
        # elec mix links (nuclear)
        buses = n.stores.query("carrier == 'uranium'").bus
        links = n.links.query("bus0 in @buses and bus1 in @elec").index
        elec_mix_nuc = -1*n.links_t.p1[links].multiply(n.snapshot_weightings.objective,axis=0).T.groupby(n.links.carrier).sum().T.sum()
        # elec mix links (coal/lignite)
        buses = n.stores.query("carrier in ['coal', 'lignite']").bus
        links = n.links.query("bus0 in @buses and bus1 in @elec").index
        elec_mix_coal = -1*n.links_t.p1[links].multiply(n.snapshot_weightings.objective,axis=0).T.groupby(n.links.carrier).sum().T.sum()

        elec_mix_array = [
            [elec_mix_nuc.sum(), 0, 0, 0, 0],
            [
                elec_mix.loc[["offwind-ac", "offwind-dc"]].sum(),
                elec_mix.loc[["onwind"]].sum(),
                (elec_mix.loc[["ror"]] + elec_mix_hydro.loc[["hydro"]].sum()).sum(),
                elec_mix.loc[["solar"]].sum(), 
                elec_mix.loc[["solar rooftop"]].sum(),
            ],
            [elec_mix_coal.sum().sum(), elec_mix_gas.sum(), 0, 0, 0]
        ]

        plot_pies(ax, elec_mix_array)


        ax.set(aspect="equal")
        # set nice_name to LRGH for LR in 2040 and 2050
        nice_name = "LREH" if nice_name == "LROH" and planning_horizon in ["2040", "2050"] else nice_name
        ax.set_title(nice_name, fontsize=6)
    
    nuclear_patch = mpatches.Patch(color='#ff8c00', label='Nuclear')
    onwind_patch = mpatches.Patch(color='#235ebc', label='Onshore Wind')
    offwind_patch = mpatches.Patch(color='#6895dd', label='Offshore Wind')
    ror_patch = mpatches.Patch(color='#3dbfb0', label='Hydropower')
    solar_patch = mpatches.Patch(color='#f9d002', label='Solar utility')
    solar_rooftop_patch = mpatches.Patch(color='#ffea80', label='Solar rooftop')
    gas_patch = mpatches.Patch(color='#db6a25', label='Gas')
    vres_patch = mpatches.Patch(color='#0fa101', label='VRES')
    coal_patch = mpatches.Patch(color='#545454', label='Hard Coal')
    fossil_patch = mpatches.Patch(color='#A18181', label='Fossil Fuel')

    if isinstance(axes, np.ndarray):
        axes[1].legend(
        handles=[nuclear_patch, vres_patch, gas_patch, coal_patch, fossil_patch, onwind_patch, offwind_patch, ror_patch, solar_patch, solar_rooftop_patch],
        loc="lower center", ncol=5, fontsize=4, bbox_to_anchor=(1.1, -0.15)
        )
    else:
        axes.legend(
        handles=[nuclear_patch, vres_patch, gas_patch, coal_patch, fossil_patch, onwind_patch, offwind_patch, ror_patch, solar_patch, solar_rooftop_patch],
        loc="lower center", ncol=5, fontsize=4, bbox_to_anchor=(0.5, -0.15)
        )
    
    # move to base directory
    change_path_to_base()

    plt.savefig(snakemake.output.figure, dpi=600, bbox_inches="tight")

