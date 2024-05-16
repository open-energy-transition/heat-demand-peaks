import sys
sys.path.append("../submodules/pypsa-eur")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import pandas as pd
import numpy as np
import cartopy.crs as ccrs
import logging
import warnings
import math
warnings.filterwarnings("ignore")
from _helpers import mock_snakemake, update_config_from_wildcards, load_network, \
                     change_path_to_pypsa_eur, change_path_to_base, \
                     LINE_LIMITS, CO2L_LIMITS, BAU_HORIZON

logger = logging.getLogger(__name__)

"""
Shadow (dual) price of the line volume constraint in multiplies of underground cable costs.
"""


def get_congestion_spatial(n, scaling_factor = 0.4):

    drop = n.buses[n.buses.x.isna()].index
    n.mremove("Bus", drop)

    drop = n.buses[n.buses.x == -5.5].index
    n.mremove("Bus", drop)

    test = n.links.query("carrier == 'DC' and length > 0").length
    if len(test) == len(n.links.query("carrier == 'DC' and length == 0").index):
        test.index = n.links.query("carrier == 'DC' and length == 0").index
        n.links.loc[n.links.query("carrier == 'DC' and length == 0").index, "length"] = test

    # https://pdf.sciencedirectassets.com/318494/1-s2.0-S2589004221X00120/1-s2.0-S2589004221014668/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEJj%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJGMEQCIH76%2F5RRp%2FsyvzTic70TZIlp86HT4kbfRopf%2BPfCiKlIAiBqjv49cp5HNULujFZfLq19GyetA%2F809nhKn69hEosbdyq8BQjw%2F%2F%2F%2F%2F%2F%2F%2F%2F%2F8BEAUaDDA1OTAwMzU0Njg2NSIMMCqRnIBKGAzBgNaHKpAFjqw6qsazYTtiYyRP17t2V2FEz8Kvij96D6Io2g6uSMPQMK5%2FCHjGAAfz1lVHhlsAjR1PjMJDzUtkEc%2BNsxNb7Dlfl80MXXYSwrkBXbn9JPbwuAFSdcAxQ46zMZ9a83oHcBIT27rV%2BKCBNqyu6Y39mAV6yUx6xh82OCFidkFybcvPs59aNR4NrN6mg6OvD%2FYcB4FKebXtllUJ6s7F07CXdT%2Fk%2BDaRrGQu2MStEUjUiBAMXUdf%2BWCzSsEvfz9HJJNmZE5B3vUrRSKI6%2B5rbD%2Bu9RkuT28RFfb%2FHvQTAAIIhFq91XePPtsJIDTSFpw0R1U7YMBofS4AH0MaDXETHM2ZdOot4yQZD%2BkskD0PZW2LRzqp6g3OJ2TV8E%2FViHw7vYF8ad1e2tWU%2FGIs%2Bdt20l%2FKu%2BKtJ6pr1kOCkRGmU3oqxc2DI0opTwJLaUQ2LvRs3%2FmIL%2FRmvRzg%2F6XUS0IZhyjUysqcyI%2FkVb%2B52uTGM5GhhAYpQcPDGbeyXp1T3Pd0fqKxlDxn3uKaIuAbB8OspJfBBpPVRF6Juye4Z11l6qhi5cWV1uDNj1FzUnbxQhq27Upqf3RRzqdSy085BPhsqhrfem06LygmcaLbVxjqj%2FRMG2tUVbSj%2Fbbl4sfAMCW0CHZAmzoN3ZKT2auVJ07jGdRbhKud14oKlzGDtZRSbGGxdsacuAav65srDgPQx2jdKn5oysvWUbkeWb7Ers1KwsiGa%2FPUstkxV%2BdZh3j%2BW2Jw9YJz4V3XwXzT3yTGHAi9mAHlOuNDtD2pqaDQXNTyXuYLvKeLDmg44CACOD72UdxdB3qS81HLQGwwKzTKxMfT2GfxoixIHUo2qqtACKTPuf%2FuZuIMmyDAbqHzz8uyfJhupdAwu6%2FusQY6sgEkY6yDtQ0mlyGtyXl%2BwEU1JQzkZDWO7kedvhtEUzH%2BS2tYGu5TDRnwHXIEKZ3py%2BJ5hHlHuO5qNhzlu7aDkREluSr5DdGiBvINuAGGlfKLcApKATslo43d1xYupSmS2lN8GfrmWR00C%2F9%2F40HcRBAdWENO%2BIjO2vldofIdGM%2Bm7%2Bn1NeQzwnKOJgDgam%2BLEd8PHjyppJZPYsqJ4lJ6h4rxe3iO6Rw%2F%2BtSXK2akqcJkDnkn&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20240508T155938Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTY7JI3FMEE%2F20240508%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=dbbc05d918718a07e6d93d5e80a132a7a580c594d44027197c1a1b07e868f0eb&hash=17515ee88de6288f8910103d9e5e5400e91b9e9644facac0bc3fd21b4bbd2597&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S2589004221014668&tid=spdf-77f5d991-24e5-49bb-82f7-b5bffd495306&sid=55b3e9e9166f494fcd2bfde0146533710f85gxrqa&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=01025a56530f005b595b&rr=880aa879ac193a8a&cc=de

    link_widths = pd.Series(index=n.links.index, data=0)
    link_widths = n.links_t.mu_lower.multiply(n.snapshot_weightings.objective, axis=0).sum(axis=0)
    link_widths += -1*n.links_t.mu_upper.multiply(n.snapshot_weightings.objective, axis=0).sum(axis=0)
    link_widths /= n.links.length # shadow price is EUR / MVA
    link_widths /= 100 # costs DC: 100 EUR / MWkm
    link_widths[n.links.query("carrier != 'DC'").index] = 0
    link_widths = link_widths.apply(lambda b: math.floor(b))*scaling_factor

    line_widths = n.lines_t.mu_lower.multiply(n.snapshot_weightings.objective, axis=0).sum(axis=0)
    line_widths += -1*n.lines_t.mu_upper.multiply(n.snapshot_weightings.objective, axis=0).sum(axis=0)
    line_widths /= n.lines.length
    line_widths /= 200 # costs AC: 200 EUR / MWkm
    line_widths = line_widths.apply(lambda b: math.floor(b))*scaling_factor

    line_color = pd.Series(index=n.lines.index, data="black")
    line_color[line_widths==0] = "green"

    link_color = pd.Series(index=n.links.index, data="black")
    link_color[link_widths==0] = "green"

    return line_widths, link_widths, line_color, link_color


def add_legend(axes, scaling_factor):
    # Define line styles for the legend
    green_line = mlines.Line2D([], [], color='green', linewidth=1, label="below today's costs")
    black_line_2 = mlines.Line2D([], [], color='black', linewidth=2*scaling_factor, label="2x today's cost")
    black_line_3 = mlines.Line2D([], [], color='black', linewidth=3*scaling_factor, label="3x today's cost")
    black_line_5 = mlines.Line2D([], [], color='black', linewidth=5*scaling_factor, label="5x today's cost")

    # Add legend to the plot
    if isinstance(axes, np.ndarray):
        legend = axes[1,1].legend(handles=[green_line, black_line_2, black_line_3, black_line_5],
                                  loc='lower center', bbox_to_anchor=(-0.2, -0.35), 
                                  ncol=2, fontsize='x-small', title="Marginal costs of lines",
                                  title_fontsize=8)
    else:
        legend = axes.legend(handles=[green_line, black_line_2, black_line_3, black_line_5],
                             loc='lower center', bbox_to_anchor=(0.5, -0.18), 
                             ncol=2, fontsize='small', title="Marginal costs of lines",
                             title_fontsize=10)
    return legend


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "plot_transmission_congestion", 
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
    time_resolution = config["plotting"]["time_resolution"]
    opts = config["plotting"]["sector_opts"]

    # define scenario namings
    scenarios = {"flexible": "OROH", 
                "retro_tes": "ORGH", 
                "flexible-moderate": "LROH", 
                "rigid": "NRGH"}

    table = pd.Series(index=scenarios.values(), data=0)
    if planning_horizon != BAU_HORIZON:
        lineex = line_limits[planning_horizon]
        sector_opts = f"Co2L{co2l_limits[planning_horizon]}-{time_resolution}-{opts}"

        # move to submodules/pypsa-eur
        change_path_to_pypsa_eur()
    
        # load networks
        networks = {}
        _, axes = plt.subplots(2, 2, subplot_kw={"projection":ccrs.EqualEarth()})

        for scenario, short_name in scenarios.items():
            n = load_network(lineex, clusters, sector_opts, planning_horizon, scenario)
            if scenario == "flexible": ax = axes[0,0]
            if scenario == "retro_tes": ax = axes[0,1]
            if scenario == "flexible-moderate": ax = axes[1,0]
            if scenario == "rigid": ax = axes[1,1]

            if n is None:
                # Skip further computation for this scenario if network is not loaded
                print(f"Network is not found for scenario '{scenario}', planning year '{planning_horizon}', and time resolution of '{time_resolution}'. Skipping...")
                continue

            scaling_factor = 4e-1
            line_widths, link_widths, line_color, link_color = get_congestion_spatial(n, scaling_factor)

            table.loc[short_name] = (line_widths/scaling_factor).sum()+(link_widths/scaling_factor).sum()
            line_widths[line_widths==0] = 1
            link_widths[link_widths==0] = (
                link_widths.loc[n.links.query("carrier == 'DC'").index].apply(lambda b: b if b > 0 else 1*scaling_factor)
            )

            n.plot(
                ax=ax, color_geomap={"land": "ghostwhite"},
                line_colors=line_color, link_colors=link_color,
                link_widths=link_widths, line_widths=line_widths
            )
            ax.set_title(short_name)

        add_legend(axes, scaling_factor)    
        
        # move to base directory
        change_path_to_base()
        plt.savefig(snakemake.output.plot, bbox_inches="tight", dpi=200)
        table.to_csv(snakemake.output.table)

    # add BAU
    BAU_horizon = BAU_HORIZON
    if BAU_horizon in config["plotting"]["planning_horizon"]:
        scenario, short_name = "BAU", "BAU"
        lineex = line_limits[BAU_horizon]
        sector_opts = f"Co2L{co2l_limits[BAU_horizon]}-{time_resolution}-{opts}"
        
        # move to submodules/pypsa-eur
        change_path_to_pypsa_eur()

        n = load_network(lineex, clusters, sector_opts, BAU_horizon, scenario)

        # move to base directory
        change_path_to_base()

        if n is None:
            # Skip further computation for this scenario if network is not loaded
            print(f"Network is not found for scenario '{scenario}', BAU year '{BAU_horizon}', and time resolution of '{time_resolution}'. Skipping...")
        else:
            _, ax = plt.subplots(subplot_kw={"projection":ccrs.EqualEarth()})
            plt.savefig("test.png", bbox_inches="tight", dpi=200)
            scaling_factor = 4e-1
            line_widths, link_widths, line_color, link_color = get_congestion_spatial(n, scaling_factor)

            table.loc[short_name] = (line_widths/scaling_factor).sum()+(link_widths/scaling_factor).sum()
            line_widths[line_widths==0] = 1
            link_widths[link_widths==0] = (
                link_widths.loc[n.links.query("carrier == 'DC'").index].apply(lambda b: b if b > 0 else 1)
            )

            n.plot(
                ax=ax, color_geomap={"land": "ghostwhite"},
                line_colors=line_color, link_colors=link_color,
                link_widths=link_widths, line_widths=line_widths
            )

            add_legend(ax, scaling_factor)

            plt.savefig(snakemake.output.plot, bbox_inches="tight", dpi=200)