
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


rule get_infra_savings:
    params:
        clusters=config["plotting"]["clusters"],
    output:
        table_cap=RESULTS+"table_infra_savings_caps_{clusters}.csv",
        table_costs=RESULTS+"table_infra_savings_costs_{clusters}.csv",
    resources:
        mem_mb=20000,
    script:
        "plots/table_infra_savings.py"


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
        expand(RESULTS
            + "table_total_costs_{clusters}.csv",
            **config["plotting"],
        ),
        expand(RESULTS
            + "table_total_capacities_{clusters}.csv",
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
