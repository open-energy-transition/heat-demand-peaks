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
        figure=RESULTS+"plot_total_costs_{clusters}_{planning_horizon}.png",
        capacities=RESULTS+"plot_total_capacities_{clusters}_{planning_horizon}.png",
    resources:
        mem_mb=20000,
    script:
        "plots/plot_total_costs.py"


rule plot_total_costs:
    input:
        expand(
            RESULTS
            + "plot_total_costs_{clusters}_{planning_horizon}.png",
            **config["plotting"],
        ),
        expand(
            RESULTS
            + "plot_total_capacities_{clusters}_{planning_horizon}.png",
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


rule get_line_congestion:
    params:
        clusters=config["plotting"]["clusters"],
    output:
        table=RESULTS+"table_line_congestion_{clusters}.csv",
    resources:
        mem_mb=20000,
    script:
        "plots/table_line_congestion.py"


rule get_line_congestions:
    input:
        expand(
            RESULTS
            + "table_line_congestion_{clusters}.csv",
            **config["plotting"],
        ),


rule plot_all:
    input:
        expand(
            RESULTS
            + "plot_total_costs_{clusters}_{planning_horizon}.png",
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
            + "table_heat_pumps_{clusters}.csv",
            **config["plotting"],
        ),
        expand(
            RESULTS
            + "table_line_congestion_{clusters}.csv",
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
