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
import logging
from _helpers import PATH_PLOTS

logger = logging.getLogger(__name__)


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
        "#235ebc", "#3dbfb0", "#f9d002", '#baa741', "#baa741", # different VRES
        "#545454", "#db6a25", "#A18181", "#A18181", "#A18181" # coal, gas, and others
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
    total_electricity = total_elec
    # Add total generated electricity to the center of the pie chart
    ax.text(0, 0, f"{total_electricity:.2f}"+"\nTWh", ha='center', va='center', fontsize=4)
    ax.set(aspect="equal")


if __name__ == "__main__":
    # define figure
    fig, ax = plt.subplots(1, 1, figsize=(1.6, 1.9)) #each scenario one axis

    # electricity generation mix for cited source https://www.consilium.europa.eu/en/infographics/how-is-eu-electricity-produced-and-sold/
    total_elec = 2641
    mix_nuclear = 0.219
    mix_coal = 0.158
    mix_gas = 0.196
    mix_other = 0.033
    mix_wind = 0.159
    mix_solar = 0.076
    mix_hydro = 0.113
    mix_biomass = 0.044
    # form electricity mix array
    elec_mix_array = [
        [total_elec*mix_nuclear, 0, 0, 0, 0],
        [total_elec*mix_wind, total_elec*mix_hydro, 
         total_elec*mix_solar, total_elec*mix_biomass, 0],
        [total_elec*mix_coal, total_elec*mix_gas, total_elec*mix_other, 0, 0]
    ]

    plot_pies(ax, elec_mix_array)

    ax.set(aspect="equal")
    ax.set_title("2022", fontsize=6)
    
    nuclear_patch = mpatches.Patch(color='#ff8c00', label='Nuclear')
    wind_patch = mpatches.Patch(color='#235ebc', label='Wind')
    hydro_patch = mpatches.Patch(color='#3dbfb0', label='Hydropower')
    solar_patch = mpatches.Patch(color='#f9d002', label='Solar')
    gas_patch = mpatches.Patch(color='#db6a25', label='Gas')
    biomass_patch = mpatches.Patch(color='#baa741', label='Biomass')
    vres_patch = mpatches.Patch(color='#0fa101', label='VRES')
    coal_patch = mpatches.Patch(color='#545454', label='Hard Coal')
    fossil_patch = mpatches.Patch(color='#A18181', label='Fossil Fuel')
    other_patch = mpatches.Patch(color='#A18181', label='Other')

    ax.legend(
    handles=[nuclear_patch, vres_patch, gas_patch, coal_patch, fossil_patch, wind_patch, hydro_patch, solar_patch, biomass_patch, other_patch],
    loc="lower center", ncol=5, fontsize=4, bbox_to_anchor=(0.5, -0.15)
    )

    plt.savefig(PATH_PLOTS+"plot_historic_generation.png", dpi=600, bbox_inches="tight")
    logger.info("Saved plot of historic electricity generation mix for 2022!")
