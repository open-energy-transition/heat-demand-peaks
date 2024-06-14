# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  Open Energy Transition gGmbH
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import os
from pathlib import Path
import pypsa
import logging
import yaml
import pandas as pd

# get the base working directory
BASE_PATH = os.path.abspath(os.path.join(__file__ ,"../.."))
# relative path to folder where to store plots
PATH_PLOTS = "plots/results/"

# Co2L limits
CO2L_LIMITS = {"2020": "0.725",
               "2030": "0.45", 
               "2040": "0.1", 
               "2050": "0.0"}
# Line limits
LINE_LIMITS = {"2020": "v1.0",
               "2030": "v1.15",
               "2040": "v1.3",
               "2050": "v1.5"}
# BAU year
BAU_HORIZON = "2020"


def mock_snakemake(
    rulename,
    root_dir=None,
    configfiles=None,
    submodule_dir="workflow/submodules/pypsa-eur",
    **wildcards,
):
    """
    This function is expected to be executed from the 'scripts'-directory of '
    the snakemake project. It returns a snakemake.script.Snakemake object,
    based on the Snakefile.

    If a rule has wildcards, you have to specify them in **wildcards.

    Parameters
    ----------
    rulename: str
        name of the rule for which the snakemake object should be generated
    root_dir: str/path-like
        path to the root directory of the snakemake project
    configfiles: list, str
        list of configfiles to be used to update the config
    submodule_dir: str, Path
        in case PyPSA-Eur is used as a submodule, submodule_dir is
        the path of pypsa-eur relative to the project directory.
    **wildcards:
        keyword arguments fixing the wildcards. Only necessary if wildcards are
        needed.
    """
    import os

    import snakemake as sm
    from pypsa.descriptors import Dict
    from snakemake.api import Workflow
    from snakemake.common import SNAKEFILE_CHOICES
    from snakemake.script import Snakemake
    from snakemake.settings import (
        ConfigSettings,
        DAGSettings,
        ResourceSettings,
        StorageSettings,
        WorkflowSettings,
    )

    script_dir = Path(__file__).parent.resolve()
    if root_dir is None:
        root_dir = script_dir.parent
    else:
        root_dir = Path(root_dir).resolve()

    user_in_script_dir = Path.cwd().resolve() == script_dir
    if str(submodule_dir) in __file__:
        # the submodule_dir path is only need to locate the project dir
        os.chdir(Path(__file__[: __file__.find(str(submodule_dir))]))
    elif user_in_script_dir:
        os.chdir(root_dir)
    elif Path.cwd().resolve() != root_dir:
        raise RuntimeError(
            "mock_snakemake has to be run from the repository root"
            f" {root_dir} or scripts directory {script_dir}"
        )
    try:
        for p in SNAKEFILE_CHOICES:
            if os.path.exists(p):
                snakefile = p
                break
        if configfiles is None:
            configfiles = []
        elif isinstance(configfiles, str):
            configfiles = [configfiles]

        resource_settings = ResourceSettings()
        config_settings = ConfigSettings(configfiles=map(Path, configfiles))
        workflow_settings = WorkflowSettings()
        storage_settings = StorageSettings()
        dag_settings = DAGSettings(rerun_triggers=[])
        workflow = Workflow(
            config_settings,
            resource_settings,
            workflow_settings,
            storage_settings,
            dag_settings,
            storage_provider_settings=dict(),
        )
        workflow.include(snakefile)

        if configfiles:
            for f in configfiles:
                if not os.path.exists(f):
                    raise FileNotFoundError(f"Config file {f} does not exist.")
                workflow.configfile(f)

        workflow.global_resources = {}
        rule = workflow.get_rule(rulename)
        dag = sm.dag.DAG(workflow, rules=[rule])
        wc = Dict(wildcards)
        job = sm.jobs.Job(rule, dag, wc)

        def make_accessable(*ios):
            for io in ios:
                for i, _ in enumerate(io):
                    io[i] = os.path.abspath(io[i])

        make_accessable(job.input, job.output, job.log)
        snakemake = Snakemake(
            job.input,
            job.output,
            job.params,
            job.wildcards,
            job.threads,
            job.resources,
            job.log,
            job.dag.workflow.config,
            job.rule.name,
            None,
        )
        # create log and output dir if not existent
        for path in list(snakemake.log) + list(snakemake.output):
            Path(path).parent.mkdir(parents=True, exist_ok=True)

    finally:
        if user_in_script_dir:
            os.chdir(script_dir)
    return snakemake


def update_config_from_wildcards(config, w):
    if w.get("planning_horizon"):
        planning_horizon = w.planning_horizon
        config["plotting"]["planning_horizon"] = planning_horizon
        config["set_capacities"]["planning_horizon"] = planning_horizon
        config["moderate_retrofitting"]["planning_horizon"] = planning_horizon
        config["improve_cops_after_renovation"]["planning_horizon"] = planning_horizon
    if w.get("clusters"):
        clusters = w.clusters
        config["plotting"]["clusters"] = clusters
        config["set_capacities"]["clusters"] = clusters
        config["moderate_retrofitting"]["clusters"] = clusters
        config["improve_cops_after_renovation"]["clusters"] = clusters
    if w.get("scenario"):
        scenario = w.scenario
        config["set_capacities"]["scenario"] = scenario
        config["improve_cops_after_renovation"]["scenario"] = scenario
    return config


def load_network(lineex, clusters, sector_opts, planning_horizon, scenario):
    FILE = f"elec_s_{clusters}_l{lineex}__{sector_opts}_{planning_horizon}.nc"
    DIR = f"results/{scenario}/postnetworks"
    try:
        n = pypsa.Network(os.path.join(DIR, FILE))
        logging.info(f"Loading {FILE} in {DIR}")
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return None
    return n


def load_unsolved_network(lineex, clusters, sector_opts, planning_horizon, scenario):
    FILE = f"elec_s_{clusters}_l{lineex}__{sector_opts}_{planning_horizon}.nc"
    DIR = f"results/{scenario}/prenetworks"
    try:
        n = pypsa.Network(os.path.join(DIR, FILE))
        logging.info(f"Loading {FILE} in {DIR}")
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return None
    return n


def save_unsolved_network(network, lineex, clusters, sector_opts, planning_horizon, scenario):
    FILE = f"elec_s_{clusters}_l{lineex}__{sector_opts}_{planning_horizon}.nc"
    DIR = f"results/{scenario}/prenetworks/"
    network.export_to_netcdf(DIR+FILE)
    logging.info(f"Saving {FILE} to {DIR}")


def change_path_to_pypsa_eur():
    # path to pypsa-eur
    pypsa_path = "submodules/pypsa-eur/"
    # absolute path to pypsa-eur
    new_path = os.path.join(BASE_PATH, pypsa_path)
    # change path to pypsa-eur
    os.chdir(new_path)


def change_path_to_base():
    os.chdir(BASE_PATH)
    # create folder to store images
    os.makedirs(PATH_PLOTS, exist_ok=True)


def get_config_path(scenario, horizon):
    configname_dict = {"flexible": "flexible-industry",
                       "flexible-moderate": "flexible-moderate",
                       "retro_tes": "retro_tes-industry",
                       "rigid": "rigid-industry"}
    configname = f"config.{configname_dict[scenario]}_{horizon}.yaml"
    configpath = "configs/EEE_study/"
    return configpath + configname


def get_config(scenario, horizon):
    config_path = get_config_path(scenario, horizon)
    # read config file
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    return config


def replace_multiindex_values(multiindex, old_value, new_value):
    # Create a new MultiIndex with replaced values
    new_tuples = [new_value if item == old_value else item for item in multiindex]
    return pd.MultiIndex.from_tuples(new_tuples, names=multiindex.names)
