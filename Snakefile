# SPDX-FileCopyrightText:  Open Energy Transition gGmbH
#
# SPDX-License-Identifier: AGPL-3.0-or-later

configfile: "configs/config.plot.yaml"

RESULTS = "plots/results/"


wildcard_constraints:
    clusters="[0-9]+",
    planning_horizon="[0-9]{4}",


localrules:
    all,


rule plot_total_cost:
    params:
        clusters=config["plotting"]["clusters"],
        planning_horizon=config["plotting"]["planning_horizon"],
    output:
        costs=RESULTS+"table_total_costs_{clusters}.csv",
        capacities=RESULTS+"table_total_capacities_{clusters}.csv",
    resources:
        mem_mb=20000,
    script:
        "plots/plot_total_costs.py"


rule plot_total_costs:
    input:
        expand(RESULTS
            + "table_total_costs_{clusters}.csv",
            **config["plotting"],
        ),
        expand(RESULTS
            + "table_total_capacities_{clusters}.csv",
            **config["plotting"],
        ),


rule plot_capacity_expansion:
    params:
        clusters=config["plotting"]["clusters"],
        planning_horizon=config["plotting"]["planning_horizon"],
    output:
        capacities=RESULTS+"table_capacity_expansions_{clusters}.csv",
    resources:
        mem_mb=20000,
    script:
        "plots/plot_capacity_expansion.py"


rule plot_capacity_expansions:
    input:
        expand(RESULTS
            + "table_capacity_expansions_{clusters}.csv",
            **config["plotting"],
        ),


rule plot_electricity_bill:
    params:
        clusters=config["plotting"]["clusters"],
        planning_horizon=config["plotting"]["planning_horizon"],
    output:
        figure_bills=RESULTS+"plot_bill_per_household_{clusters}_{planning_horizon}.png",
        figure_price=RESULTS+"plot_prices_per_MWh_{clusters}_{planning_horizon}.png",
    resources:
        mem_mb=20000,
    script:
        "plots/plot_electricity_bills.py"


rule plot_electricity_for_heat:
    params:
        clusters=config["plotting"]["clusters"],
        planning_horizon=config["plotting"]["planning_horizon"],
    output:
        figure=RESULTS+"plot_electricity_for_heat_{clusters}_{planning_horizon}.png",
    resources:
        mem_mb=20000,
    script:
        "plots/plot_electricity_for_heat.py"


rule plot_electricity_bills:
    input:
        expand(
            RESULTS
            + "plot_bill_per_household_{clusters}_{planning_horizon}.png",
            **config["plotting"],
        ),
        expand(
            RESULTS
            + "plot_prices_per_MWh_{clusters}_{planning_horizon}.png",
            **config["plotting"],
        ),


rule plot_electricity_for_heats:
    input:
        expand(
            RESULTS
            + "plot_electricity_for_heat_{clusters}_{planning_horizon}.png",
            **config["plotting"],
        ),


rule plot_heat_saving:
    params:
        clusters=config["plotting"]["clusters"],
        planning_horizon=config["plotting"]["planning_horizon"],
    output:
        figure=RESULTS+"plot_heat_savings_{clusters}_{planning_horizon}.png",
        figure_full=RESULTS+"plot_heat_savings_full_year_{clusters}_{planning_horizon}.png",
    resources:
        mem_mb=10000,
    script:
        "plots/plot_heat_savings.py"


rule plot_heat_savings:
    input:
        expand(
            RESULTS
            + "plot_heat_savings_{clusters}_{planning_horizon}.png",
            **config["plotting"],
        ),
        expand(
            RESULTS
            + "plot_heat_savings_full_year_{clusters}_{planning_horizon}.png",
            **config["plotting"],
        ),


rule plot_electricity_generation:
    params:
        clusters=config["plotting"]["clusters"],
        planning_horizon=config["plotting"]["planning_horizon"],
    output:
        figure=RESULTS+"plot_elec_generation_{clusters}_{planning_horizon}.png",
    resources:
        mem_mb=20000,
    script:
        "plots/plot_electricity_generation.py"


rule plot_electricity_generations:
    input:
        expand(
            RESULTS
            + "plot_elec_generation_{clusters}_{planning_horizon}.png",
            **config["plotting"],
        ),


rule plot_curtailment:
    params:
        clusters=config["plotting"]["clusters"],
    output:
        figure=RESULTS+"plot_curtailment_{clusters}.png",
        table=RESULTS+"table_curtailment_{clusters}.csv",
    resources:
        mem_mb=20000,
    script:
        "plots/plot_curtailment.py"


rule plot_curtailments:
    input:
        expand(
            RESULTS
            + "plot_curtailment_{clusters}.png",
            **config["plotting"],
        ),
        expand(
            RESULTS
            + "table_curtailment_{clusters}.csv",
            **config["plotting"],
        ),


rule plot_hydrogen_production:
    params:
        clusters=config["plotting"]["clusters"],
    output:
        figure_elec=RESULTS+"plot_H2_prod_elec_use_{clusters}.png",
        figure_prod=RESULTS+"plot_H2_prod_{clusters}.png",
        table=RESULTS+"table_H2_prod_{clusters}.csv",
    resources:
        mem_mb=20000,
    script:
        "plots/plot_hydrogen_production.py"


rule plot_hydrogen_productions:
    input:
        expand(
            RESULTS
            + "plot_H2_prod_elec_use_{clusters}.png",
            **config["plotting"],
        ),
        expand(
            RESULTS
            + "plot_H2_prod_{clusters}.png",
            **config["plotting"],
        ),
        expand(
            RESULTS
            + "table_H2_prod_{clusters}.csv",
            **config["plotting"],
        ),


rule plot_heat_tech_ratio:
    params:
        clusters=config["plotting"]["clusters"],
        planning_horizon=config["plotting"]["planning_horizon"],
    output:
        table=RESULTS+"table_heat_tech_ratio_{clusters}.csv",
    resources:
        mem_mb=20000,
    script:
        "plots/plot_heat_tech_ratio.py"


rule plot_heat_tech_ratios:
    input:
        expand(
            RESULTS
            + "table_heat_tech_ratio_{clusters}.csv",
            **config["plotting"],
        ),


rule plot_COP:
    output:
        figure=RESULTS+"plot_COP.png",
    resources:
        mem_mb=20000,
    script:
        "plots/plot_COP.py"


rule plot_co2_level:
    params:
        clusters=config["plotting"]["clusters"],
    output:
        table=RESULTS+"table_co2_level_{clusters}.csv",
        table_savings=RESULTS+"table_co2_savings_{clusters}.csv",
    resources:
        mem_mb=20000,
    script:
        "plots/plot_co2_level.py"


rule plot_co2_levels:
    input:
        expand(
            RESULTS
            + "table_co2_level_{clusters}.csv",
            **config["plotting"],
        ),
        expand(
            RESULTS
            + "table_co2_savings_{clusters}.csv",
            **config["plotting"],
        ),


rule plot_historic_generation:
    shell:
        "python plots/plot_historic_generation.py"


rule get_heat_pump:
    params:
        clusters=config["plotting"]["clusters"],
    output:
        table=RESULTS+"table_heat_pumps_{clusters}.csv",
    resources:
        mem_mb=20000,
    script:
        "plots/table_heat_pumps.py"


rule get_heat_pumps:
    input:
        expand(
            RESULTS
            + "table_heat_pumps_{clusters}.csv",
            **config["plotting"],
        ),


rule get_infra_saving:
    params:
        clusters=config["plotting"]["clusters"],
    output:
        table_cap=RESULTS+"table_infra_savings_caps_{clusters}.csv",
        table_costs=RESULTS+"table_infra_savings_costs_{clusters}.csv",
    resources:
        mem_mb=20000,
    script:
        "plots/table_infra_savings.py"


rule get_infra_savings:
    input:
        expand(
            RESULTS
            + "table_infra_savings_caps_{clusters}.csv",
            **config["plotting"],
        ),
        expand(
            RESULTS
            + "table_infra_savings_costs_{clusters}.csv",
            **config["plotting"],
        ),


rule plot_transmission_congestion:
    params:
        clusters=config["plotting"]["clusters"],
    output:
        plot=RESULTS+"plot_line_congestion_{clusters}_{planning_horizon}.png",
        table=RESULTS+"table_line_congestion_{clusters}_{planning_horizon}.csv",
    script:
        "plots/plot_transmission_congestion.py"


rule plot_transmission_congestions:
    input:
        expand(RESULTS
            + "plot_line_congestion_{clusters}_{planning_horizon}.png",
            **config["plotting"],
        ),
        expand(RESULTS
            + "table_line_congestion_{clusters}_{planning_horizon}.csv",
            **config["plotting"],
        ),


rule get_res_share:
    params:
        clusters=config["plotting"]["clusters"],
    output:
        table=RESULTS+"table_res_share_{clusters}.csv",
    resources:
        mem_mb=20000,
    script:
        "plots/table_res_share.py"


rule get_res_shares:
    input:
        expand(
            RESULTS
            + "table_res_share_{clusters}.csv",
            **config["plotting"],
        ),


rule plot_all:
    input:
        expand(RESULTS
            + "table_total_costs_{clusters}.csv",
            **config["plotting"],
        ),
        expand(RESULTS
            + "table_total_capacities_{clusters}.csv",
            **config["plotting"],
        ),
        expand(RESULTS
            + "table_capacity_expansions_{clusters}.csv",
            **config["plotting"],
        ),
        expand(
            RESULTS
            + "plot_bill_per_household_{clusters}_{planning_horizon}.png",
            **config["plotting"],
        ),
        expand(
            RESULTS
            + "plot_prices_per_MWh_{clusters}_{planning_horizon}.png",
            **config["plotting"],
        ),
        expand(
            RESULTS
            + "plot_electricity_for_heat_{clusters}_{planning_horizon}.png",
            **config["plotting"],
        ),
        expand(
            RESULTS
            + "plot_heat_savings_{clusters}_{planning_horizon}.png",
            **config["plotting"],
        ),
        expand(
            RESULTS
            + "plot_heat_savings_full_year_{clusters}_{planning_horizon}.png",
            **config["plotting"],
        ),
        expand(
            RESULTS
            + "table_heat_pumps_{clusters}.csv",
            **config["plotting"],
        ),
        expand(
            RESULTS
            + "plot_elec_generation_{clusters}_{planning_horizon}.png",
            **config["plotting"],
        ),
        expand(
            RESULTS
            + "plot_curtailment_{clusters}.png",
            **config["plotting"],
        ),
        expand(
            RESULTS
            + "table_curtailment_{clusters}.csv",
            **config["plotting"],
        ),
        expand(
            RESULTS
            + "plot_H2_prod_elec_use_{clusters}.png",
            **config["plotting"],
        ),
        expand(
            RESULTS
            + "plot_H2_prod_{clusters}.png",
            **config["plotting"],
        ),
        expand(
            RESULTS
            + "plot_COP.png",
            **config["plotting"],
        ),
        expand(
            RESULTS
            + "table_co2_level_{clusters}.csv",
            **config["plotting"],
        ),
        expand(
            RESULTS
            + "table_co2_savings_{clusters}.csv",
            **config["plotting"],
        ),
        expand(
            RESULTS
            + "table_H2_prod_{clusters}.csv",
            **config["plotting"],
        ),
        expand(
            RESULTS
            + "plot_line_congestion_{clusters}_{planning_horizon}.png",
            **config["plotting"],
        ),
        expand(
            RESULTS
            + "table_heat_tech_ratio_{clusters}.csv",
            **config["plotting"],
        ),
        expand(
            RESULTS
            + "table_infra_savings_caps_{clusters}.csv",
            **config["plotting"],
        ),
        expand(
            RESULTS
            + "table_infra_savings_costs_{clusters}.csv",
            **config["plotting"],
        ),
        expand(
            RESULTS
            + "table_res_share_{clusters}.csv",
            **config["plotting"],
        ),


rule set_capacity:
    params:
        clusters=config["set_capacities"]["clusters"],
        planning_horizon=config["set_capacities"]["planning_horizon"],
        scenario=config["set_capacities"]["scenario"]
    output:
        logs="scripts/logs/set_capacities_{clusters}_{planning_horizon}_{scenario}.txt",
    resources:
        mem_mb=20000,
    script:
        "scripts/set_capacities.py"


rule set_capacities:
    input:
        expand(
            "scripts/logs/set_capacities_{clusters}_{planning_horizon}_{scenario}.txt",
            **config["set_capacities"],
        ),


rule moderate_retrofitting:
    params:
        clusters=config["moderate_retrofitting"]["clusters"],
        planning_horizon=config["moderate_retrofitting"]["planning_horizon"],
    output:
        logs="scripts/logs/set_moderate_retrofitting_{clusters}_{planning_horizon}.txt",
    resources:
        mem_mb=20000,
    script:
        "scripts/moderate_retrofitting.py"


rule set_moderate_retrofitting:
    input:
        expand(
            "scripts/logs/set_moderate_retrofitting_{clusters}_{planning_horizon}.txt",
            **config["moderate_retrofitting"],
        ),


rule improve_cops_after_renovation:
    params:
        clusters=config["improve_cops_after_renovation"]["clusters"],
        planning_horizon=config["improve_cops_after_renovation"]["planning_horizon"],
        scenario=config["improve_cops_after_renovation"]["scenario"],
    output:
        logs="scripts/logs/improve_cops_after_renovation_{clusters}_{planning_horizon}_{scenario}.txt",
    resources:
        mem_mb=20000,
    script:
        "scripts/improve_cops_after_renovation.py"


rule improve_cops_after_renovations:
    input:
        expand(
            "scripts/logs/improve_cops_after_renovation_{clusters}_{planning_horizon}_{scenario}.txt",
            **config["improve_cops_after_renovation"],
        ),
