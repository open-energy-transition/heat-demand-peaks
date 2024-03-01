import pypsa
import matplotlib.pyplot as plt
import geopandas as gpd
import numpy as np
import pandas as pd
from pathlib import Path
import seaborn as sns
from datetime import datetime
from cartopy import crs as ccrs

network_path = "data/"
network_file = "elec_s_48_lvopt__Co2L0-1H-T-H-B-I_2030_flexible.nc"

building_data_dir = "data/"
building_stock_fl = "data_building_stock.csv"
# extracted from build_retro.py of PyPSA-Eur
thermal_quality_fl = "heat_transfer.csv"

building_stock_df = pd.read_csv(building_data_dir + building_stock_fl)
thermal_quality_df = pd.read_csv(building_data_dir + thermal_quality_fl)


def extract_E_ann(n, country):
    Q_heat_ann = (
        n.loads_t.p_set
        .filter(like=country)
        .filter(like="heat")
        .sum(axis=1)
        .mul(n.snapshot_weightings.generators)
        .sum()/1e6
    )
    return Q_heat_ann

def extract_relative_dE_ann(retro_scheme, n, country):

    dE_retro_ambit_ann = (
        n.generators_t.p
        .filter(like=country)
        .filter(like="retrofitting")
        .filter(like=retro_scheme)
        .sum(axis=1)
        .mul(n.snapshot_weightings.generators)
        .sum()/1e6
    )

    Q_heat_ann = (
        n.loads_t.p_set
        .filter(like=country)
        .filter(like="heat")
        .sum(axis=1)
        .mul(n.snapshot_weightings.generators)
        .sum()/1e6
    )

    return(dE_retro_ambit_ann/Q_heat_ann)

# TODO Address data loss for the missed countries codes
def build_national_stock_stats(
    country,
    thermal_df,
    buildings_df,
):

    thermal_quality_in_country = (
        thermal_df
        .query("country_code == @country")
    )
    
    n_buildings_in_country = (
        buildings_df
        .query("country_code == @country")
        .query("topic == 'BUILDING'")
        .query("feature == 'Area'")
        .query("subsector != 'Total'")
        .query("type == 'Number of buildings [Mil.]'")
    )
    n_buildings_in_country["assumed_subsector"] = n_buildings_in_country.subsector
    n_buildings_in_country.subsector = n_buildings_in_country.assumed_subsector
    
    n_U_original_in_country = (
        buildings_df
        .query("country_code == @country")
        .query("topic == 'BUILDING'")
        .query("feature == 'Construction features (U-values)' | feature == 'Construction features (U-value)'")
        .query("subsector != 'Total'")
    )
    n_U_original_in_country["assumed_subsector"] = n_U_original_in_country.subsector
    n_U_original_in_country.subsector = n_U_original_in_country.assumed_subsector

    n_U_original_in_country_aggreg = (
        n_U_original_in_country[["country_code", "subsector", "bage", "value"]]
        .groupby(["country_code", "subsector", "bage"])
        .mean()
    )
    
    national_buildings_df = (
        thermal_quality_in_country
        .set_index(["country_code", "subsector", "bage"])
        .join(
            n_buildings_in_country[["country_code", "subsector", "bage", "value"]]
            .set_index(["country_code", "subsector", "bage"]),
            rsuffix="_mln_buildings"
            )
        )
    national_buildings_df = (
        national_buildings_df
        .join(
            n_U_original_in_country[["country_code", "subsector", "bage", "value"]]
            .groupby(["country_code", "subsector", "bage"])
            .mean(),
            #.set_index(["country_code", "subsector", "bage"]),
            rsuffix="_U_before_retrofit"
        )
    )

    return national_buildings_df

def prepare_retro_df(
    national_buildings_df, 
    retrofitted_sectors=["AB", "MFH", "SFH"]
    ):
    bg = ['Before 1945', '1945 - 1969', '1970 - 1979', '1980 - 1989',
        '1990 - 1999', '2000 - 2010', 'Post 2010']
    # TODO 1 "value_mln_buildings" column may be not available
    # TODO 2 The heat losses values should account for a mixture of ambitious and moderate schemes
    national_buildings_df["demand_weight"] = national_buildings_df.apply(
    # heat_transfer is still [W/K], no need to multiply on A_envelope
        lambda row: row["new_U_0.07"]*row["value_mln_buildings"],
        axis=1
    )
    
    national_buildings_df["demand_share"] = (
        national_buildings_df["demand_weight"] / 
        national_buildings_df["demand_weight"].sum()
    )
    
    residential_categories = [
        i for i in retrofitted_sectors if i in national_buildings_df.index.get_level_values(1).unique()
    ]
    other_categories = [
        i for i in national_buildings_df.index.get_level_values(1).unique() if not i in retrofitted_sectors
    ]
    sk = residential_categories + other_categories
    
    #new_index = national_buildings_df.index.droplevel(level=0)
    national_buildings_df2 = national_buildings_df.copy()
    new_index_clean = pd.MultiIndex.from_product(
        [sk, bg],
        names=["subsector", "bage"]
    )
    national_buildings_df2.reset_index(level=0, drop=True, inplace=True)
    national_buildings_df2 = (
        national_buildings_df2
        .reindex(new_index_clean)
    )
    return national_buildings_df2

def eval_shell_renovation_rate(national_buildings_df, demand_decrease):
    national_buildings_df["demand_share_cumsum"] = (
        national_buildings_df.demand_share.cumsum(axis=0)
    )
    national_buildings_df["retrofitted"] = (
        national_buildings_df["demand_share_cumsum"] <= demand_decrease
    )
    
    shell_renovation_rate = (
        national_buildings_df.query("retrofitted")
        .eval("A_envelope * value_mln_buildings").sum()/
        national_buildings_df.eval("A_envelope * value_mln_buildings").sum()
    )
    
    demand_share_bulk_covered_by_retro = (
        national_buildings_df.query("retrofitted")["demand_share"].sum()
    )

    # the next buildings cohort may be too numerous to cover the remained demand part
    remained_demand_share_to_cover = (
       demand_decrease_overall - demand_share_bulk_covered_by_retro
    )
    if remained_demand_share_to_cover:
    
        k_retro_share_cohort = (
            national_buildings_df.loc[national_buildings_df[["retrofitted"]].idxmin()]["demand_share"] /
            remained_demand_share_to_cover
        )
    
        shell_renovation_rate_cohort = (
            national_buildings_df.loc[national_buildings_df[["retrofitted"]].idxmin()]
            .eval("A_envelope * value_mln_buildings").sum()/
            national_buildings_df.eval("A_envelope * value_mln_buildings").sum()
        )
    
        shell_renovation_to_add = k_retro_share_cohort[0] * shell_renovation_rate_cohort
    else:
        shell_renovation_to_add = 0    
    
    shell_renovation_rate = (
        demand_share_bulk_covered_by_retro + shell_renovation_to_add
    )

    #print("demand_share covered")
    # print(
    #     overall_building_stock_df2.query("retrofitted")["demand_share"].sum() + 
    #     remained_demand_share_to_cover
    # )

    return shell_renovation_rate



# energy modeling datasets ----------------------------------------------------

n = pypsa.Network(network_path + network_file)

heat_savings = pd.DataFrame(
   {
       # drop empty country codes
       "country": n.generators.country.unique()[n.generators.country.unique().astype(bool)],
        "E_ann": -99999
    }
)

heat_savings.E_ann = heat_savings.country.apply(
    lambda x: (extract_E_ann(n=n, country=x))
)
heat_savings["ambit_retro_rel_dE_ann"] = heat_savings.country.apply(
    lambda x: (extract_relative_dE_ann(retro_scheme="ambitious", n=n, country=x))
)
heat_savings["moderate_retro_rel_dE_ann"] = heat_savings.country.apply(
    lambda x: (extract_relative_dE_ann(retro_scheme="moderate", n=n, country=x))
)
heat_savings["overall_dE_ann"] = heat_savings["moderate_retro_rel_dE_ann"] + heat_savings["ambit_retro_rel_dE_ann"]

print("heat_savings")
print(heat_savings)

# EU-wide datasets ------------------------------------------------------------

building_stock_clean = building_stock_df.copy()
rename_sectors = {
    "Single family- Terraced houses": "SFH",
    "Multifamily houses": "MFH",
    "Appartment blocks": "AB",
    "Hotels and restaurants": "Hotels and Restaurants"
}

building_stock_clean.subsector.replace(rename_sectors, inplace=True)
building_stock_clean.btype.replace(rename_sectors, inplace=True)
# for missing weighting of surfaces of building types assume MFH
building_stock_clean["assumed_subsector"] = building_stock_clean.subsector
building_stock_clean.loc[
    ~building_stock_clean.subsector.isin(rename_sectors.values()), "assumed_subsector"
] = "MFH"
building_stock_clean.country_code.replace({"UK": "GB"}, inplace=True)
building_stock_clean.bage.replace({"Berfore 1945": "Before 1945"}, inplace=True)
building_stock_clean = building_stock_clean[~building_stock_clean.bage.isna()]
building_stock_clean.country_code = building_stock_clean.country_code.str.upper()

# country-specific data -------------------------------------------------------
country_in_focus = "PL"

national_thermal_quality_df = build_national_stock_stats(
    buildings_df=building_stock_clean,
    country=country_in_focus,
    thermal_df=thermal_quality_df,
)

overall_building_stock_df2 = prepare_retro_df(
    national_buildings_df=national_thermal_quality_df, 
    retrofitted_sectors=["AB", "MFH", "SFH"]
    )

demand_decrease_overall = (
    heat_savings.loc[
        heat_savings.country == country_in_focus, 
        ["ambit_retro_rel_dE_ann", "moderate_retro_rel_dE_ann"]
    ].sum()
    .sum()
)

shell_renovation_rate = eval_shell_renovation_rate(
    national_buildings_df=overall_building_stock_df2,
    demand_decrease=demand_decrease_overall
)

print("country_in_focus")
print(country_in_focus)
print("demand_decrease_overall")
print(demand_decrease_overall)
print("shell_renovation_rate")
print(shell_renovation_rate)