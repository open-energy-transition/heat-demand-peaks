import os
import sys
sys.path.append("../submodules/pypsa-eur")
import matplotlib.pyplot as plt
import pandas as pd
import logging
import yaml
import warnings
warnings.filterwarnings("ignore")
from _helpers import mock_snakemake, update_config_from_wildcards, load_network, \
                     change_path_to_pypsa_eur, change_path_to_base


logger = logging.getLogger(__name__)

RESULTS_DIR = "plots/results"
DEFAULT_CONFIG_DIR = "config/config.default.yaml"
TEMP_LOW = -20
TEMP_HIGH = 40


def coefficient_of_performance(delta_T, source="air"):
    if source == "air":
        return 6.81 - 0.121 * delta_T + 0.000630 * delta_T**2
    elif source == "soil":
        return 8.77 - 0.150 * delta_T + 0.000734 * delta_T**2
    else:
        raise NotImplementedError("'source' must be one of  ['air', 'soil']")


def get_heat_pump_sink_T():
    # read heat_pump_sink_T from default config
    with open(DEFAULT_CONFIG_DIR, 'r') as file:
        config_default = yaml.safe_load(file)

    heat_pump_sink_T = config_default["sector"]["heat_pump_sink_T"]
    return heat_pump_sink_T


def plot_COP(t, heat_pump_sink_T):
    # define sources of heat pumps
    media = {"air": "air heat pumps", 
             "soil": "ground heat pumps"}
    
    # define df for storing COP vs temperatures
    cop_df = pd.DataFrame(index=t, columns=media.values())

    # delta_T
    delta_T = heat_pump_sink_T - t

    # compute COP for each sources and temperatures
    for source, tech_name in media.items():
        cop = coefficient_of_performance(delta_T, source=source)
        cop_df.loc[:, tech_name] = cop

    # plot COP
    fig, ax = plt.subplots(figsize=(7, 3))
    cop_df.plot(ax=ax)
    ax.set_ylabel("COP")
    ax.set_xlabel("Source temperature (°C)")
    ax.legend(loc="upper left", facecolor="white", fontsize='x-small')
    plt.savefig(snakemake.output.figure, dpi=600, bbox_inches = 'tight')


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "plot_COP", 
        )
    # update config based on wildcards
    config = update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    # move to submodules/pypsa-eur
    change_path_to_pypsa_eur()

    # get_heat_pump_sink_T
    heat_pump_sink_T = get_heat_pump_sink_T()

    # move to base directory
    change_path_to_base()
    
    # temperature range for plotting COP
    t = pd.Series(range(TEMP_LOW, TEMP_HIGH))

    # plot COP for air and ground heat pumps
    plot_COP(t, heat_pump_sink_T)
