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
simpl = ""
opts = ""
clusters = "48"
ll = "vopt"
sector_opts = "Co2L0-1H-T-H-B-I"
planning_horizons = "2030"
FILE = f"elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc"
DIR_flex = "results/flexible/postnetworks"
n_flex = pypsa.Network(os.path.join(DIR_flex, FILE))

# %%
"""
### 1b. Import unsolved flexible network
"""

# %%
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
output_path = "results/flexible-moderate/prenetworks/"
if not os.path.exists(output_path):
    os.makedirs(output_path)
    print(f"Folder created at: {output_path}")
output_filename = output_path + FILE
n.export_to_netcdf(output_filename)

# %%
