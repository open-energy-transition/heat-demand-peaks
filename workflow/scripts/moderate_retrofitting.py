# %%
"""
# 0. Import libraties
"""

# %%
import os
import sys
sys.path.append("../pypsa-eur")
import pypsa
import pandas as pd
import logging
import warnings
warnings.filterwarnings("ignore")
logger = logging.getLogger(__name__)

# get the current working directory
base_path = os.path.abspath(os.path.join(__file__ ,"../.."))
# path to pypsa-eur
pypsa_path = "pypsa-eur/"
# absolute path to pypsa-eur

new_path = os.path.join(base_path, pypsa_path)
# change path to pypsa-eur
os.chdir(new_path)

# %%
"""
# 1. Import network
"""

# %%
"""
### 1a. Import solved flexible network
"""

# %%
FILE_flex = "elec_s_48_lvopt__Co2L0-1H-T-H-B-I_2030.nc"
DIR_flex = "results/flexible/postnetworks"
n_flex = pypsa.Network(os.path.join(DIR_flex, FILE_flex))

# %%
"""
### 1b. Import unsolved flexible network
"""

# %%
FILE = "elec_s_48_lvopt__Co2L0-1H-T-H-B-I_2030.nc"
#FILE = "chpCC_elec_s_48_lcopt__Co2L0-2H-T-H-B_2030.nc"
DIR = "results/flexible/prenetworks"
n = pypsa.Network(os.path.join(DIR, FILE))

# %%
"""
# 2. Set moderate retrofitting
"""

# %%
# get optimal retrofitting from the solved network
retro_opt = n_flex.generators.query("carrier in 'retrofitting'")[["p_nom_opt"]]
retro_data = retro_opt / 2

# %%
# set p_nom as half of p_nom_opt of flexible scenario
n.generators.loc[retro_data.index, "p_nom"] = retro_data["p_nom_opt"]
# set retrofitting not extendable
n.generators.loc[retro_data.index, "p_nom_extendable"] = False

# %%
# save changed network
output_filename = "results/flexible-moderate/prenetworks/" + FILE
n.export_to_netcdf(output_filename)

# %%
