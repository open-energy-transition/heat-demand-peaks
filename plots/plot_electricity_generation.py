import sys
sys.path.append("../submodules/pypsa-eur")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
from _helpers import mock_snakemake, update_config_from_wildcards, load_network, \
                     change_path_to_pypsa_eur, change_path_to_base


def plot_pies(ax, elec_mix_array):
    size = 0.4
    vals = np.array(elec_mix_array)

    cmap = plt.colormaps["tab20c"]
    outer_colors = ["#ff8c00", "#0fa101", "#db6a25"]
    inner_colors = cmap([1, 2, 5, 6, 9, 10])
    inner_colors = ["#ff8c00", "#6895dd", "#235ebc", "#3dbfb0", "#f9d002", '#ffea80', "#db6a25"]

    outer_labels = ["Nuclear", "VRES", "Gas"]

    wedges, texts, autotexts = ax.pie(vals.sum(axis=1), radius=1, colors=outer_colors,
        wedgeprops=dict(width=size, edgecolor='w', linewidth=0.2), 
        autopct='%.1f%%', textprops={'fontsize': 4}, pctdistance=1.5,
        labels=outer_labels)
    
    # Adjust the position of autopct labels
    for autotext, label in zip(autotexts, texts):
        x, y = label.get_position()  # Get position of corresponding wedge label
        autotext.set_position((x , y - 0.11))  # Set position of autopct label below the wedge label
        autotext.set_fontsize(3)
        align = "right" if x < 0 else "left"
        autotext.set_horizontalalignment(align)

    all_vals = vals.flatten()
    ax.pie(all_vals[all_vals!=0], radius=1.15-size, colors=inner_colors,
       wedgeprops=dict(width=size, edgecolor='w', linewidth=0.2),
       autopct='%.1f%%', textprops={'fontsize': 3}, pctdistance=0.7)

    ax.set(aspect="equal")


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "plot_electricity_generation_pies", 
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
    scenarios = {"flexible": "OROH", 
                 "retro_tes": "ORGH", 
                 "flexible-moderate": "LROH", 
                 "rigid": "NROH"}

    # define figure
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4) #each scenario one axis
    axes = [ax1, ax2, ax3, ax4]

    # load networks
    i = 0
    for scenario, nice_name in scenarios.items():
        ax = axes[i]
        i += 1
        
        n = load_network(lineex, clusters, sector_opts, planning_horizon, scenario)

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

        elec_mix_array = [
            [elec_mix.loc[["nuclear"]].sum(), 0, 0 , 0, 0],
            [
                elec_mix.loc[["offwind-ac", "offwind-dc"]].sum(),
                elec_mix.loc[["onwind"]].sum(),
                (elec_mix.loc[["ror"]] + elec_mix_hydro.loc[["hydro"]].sum()).sum(),
                elec_mix.loc[["solar"]].sum(), 
                elec_mix.loc[["solar rooftop"]].sum(),
            ],
            [elec_mix_gas.sum(), 0, 0, 0, 0]
        ]

        plot_pies(ax, elec_mix_array)


        ax.set(aspect="equal")
        ax.set_title(nice_name, fontsize=6)
    
    nuclear_patch = mpatches.Patch(color='#ff8c00', label='Nuclear')
    onwind_patch = mpatches.Patch(color='#235ebc', label='Onshore Wind')
    offwind_patch = mpatches.Patch(color='#6895dd', label='Offshore Wind')
    ror_patch = mpatches.Patch(color='#3dbfb0', label='Run of River')
    solar_patch = mpatches.Patch(color='#f9d002', label='Solar utility')
    solar_rooftop_patch = mpatches.Patch(color='#ffea80', label='Solar rooftop')
    gas_patch = mpatches.Patch(color='#db6a25', label='Gas')
    vres_patch = mpatches.Patch(color='#0fa101', label='VRES')

    axes[1].legend(
    handles=[nuclear_patch, vres_patch, gas_patch, onwind_patch, offwind_patch, ror_patch, solar_patch, solar_rooftop_patch],
    loc="lower center", ncol=4, fontsize=4, bbox_to_anchor=(1.1, -0.15)
    )
    
    # move to base directory
    change_path_to_base()

    plt.savefig(snakemake.output.figure, dpi=300, bbox_inches="tight")

