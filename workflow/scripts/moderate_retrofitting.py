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
FILE = "elec_s_48_lvopt__Co2L0-1H-T-H-B_2030.nc"
#FILE = "chpCC_elec_s_48_lcopt__Co2L0-2H-T-H-B_2030.nc"
DIR = "results/flexible/prenetworks"
n = pypsa.Network(os.path.join(DIR, FILE))

# %%
"""
# 2. Set moderate retrofitting
"""

# %%
# import data for moderate retrofitting
filename = "../data/flexible_retro_half.csv"
retro_data = pd.read_csv(os.path.join(os.getcwd(), filename), index_col=0)

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
