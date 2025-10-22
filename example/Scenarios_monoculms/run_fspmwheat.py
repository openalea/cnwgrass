# -*- coding: latin-1 -*-

from __future__ import print_function

import os
import sys
import getopt

import pandas as pd
from math import exp

from openalea.fspmwheat import fspmwheat_postprocessing

from example.Scenarios_monoculms import main
from example.Scenarios_monoculms import tools
from example.Scenarios_monoculms import additional_graphs


def exponential_fertilization_rate(V0, K, t, dt, plant_density):
    ferti_per_plant = V0 / K * (exp(-K * (t + dt) / 168) - exp(-K * t / 168))  # g N per plant
    return ferti_per_plant * plant_density * (10 ** 6) / 14  # µmol N m-2


def run_fspmwheat(scenario_id=1, inputs_dir_path=None, outputs_dir_path='outputs'):
    """
    Run the main.py of fspmwheat using data from a specific scenario

    :param int scenario_id: the index of the scenario to be read in the CSV file containing the list of scenarios
    :param str inputs_dir_path: the path directory of inputs
    :param str outputs_dir_path: the path to save outputs
    """

    # Path of the directory which contains the inputs of the model
    if inputs_dir_path:
        INPUTS_DIRPATH = inputs_dir_path
    else:
        INPUTS_DIRPATH = 'inputs'

    # Scenario to be run
    scenarios_df = pd.read_csv(os.path.join(INPUTS_DIRPATH, 'scenarios_list.csv'), index_col='Scenario')
    scenario = scenarios_df.loc[scenario_id].to_dict()
    scenario_name = 'Scenario_%.4d' % scenario_id

    # -- SIMULATION PARAMETERS --

    # Create dict of parameters for the scenario
    scenario_parameters = tools.buildDic(scenario)

    # Do run the simulation?
    RUN_SIMU = scenario_parameters.get('Run_Simulation', True)

    SIMULATION_LENGTH = scenario_parameters.get('Simulation_Length', 3000)

    # Do run the simulation from the output files ?
    RUN_FROM_OUTPUTS = scenario_parameters.get('Run_From_Outputs', False)

    # Do run the postprocessing?
    RUN_POSTPROCESSING = scenario_parameters.get('Run_Postprocessing', True)  #: TODO separate postprocessings coming from other models

    # Do generate the graphs?
    GENERATE_GRAPHS = scenario_parameters.get('Generate_Graphs', False)  #: TODO separate postprocessings coming from other models

    if RUN_SIMU or RUN_POSTPROCESSING or GENERATE_GRAPHS:

        # -- SIMULATION DIRECTORIES --

        # Path of the directory which contains the outputs of the model
        if not outputs_dir_path:
            outputs_dir_path = 'outputs'
        scenario_dirpath = os.path.join(outputs_dir_path, scenario_name)
        if not os.path.exists(outputs_dir_path):
            os.mkdir(outputs_dir_path)

        # Create the directory of the Scenario where results will be stored
        if not os.path.exists(scenario_dirpath):
            os.mkdir(scenario_dirpath)

        # Create directory paths for graphs, outputs and postprocessings of this scenario
        scenario_graphs_dirpath = os.path.join(scenario_dirpath, 'graphs')
        if not os.path.exists(scenario_graphs_dirpath):
            os.mkdir(scenario_graphs_dirpath)
        # Outputs
        scenario_outputs_dirpath = os.path.join(scenario_dirpath, 'outputs')
        if not os.path.exists(scenario_outputs_dirpath):
            os.mkdir(scenario_outputs_dirpath)
        # Postprocessings
        scenario_postprocessing_dirpath = os.path.join(scenario_dirpath, 'postprocessing')
        if not os.path.exists(scenario_postprocessing_dirpath):
            os.mkdir(scenario_postprocessing_dirpath)

        # -- SIMULATION CONDITIONS

        # Scenario input dirpath
        SCENARIO_INPUTS_PLANTSOIL_DIRPATH = scenario_parameters.get('Inputs_PlantSoil_Dirpath', INPUTS_DIRPATH)
        SCENARIO_INPUT_METEO_DIRPATH = scenario_parameters.get('Input_Meteo_Dirpath', INPUTS_DIRPATH)

        # Plant density and inter-row
        PLANT_DENSITY = {1: scenario_parameters.get('Plant_Density', 250.)}
        INTER_ROW = scenario_parameters.get('Inter_Row', 0.15)

        # Build N Fertilizations dict
        N_FERTILIZATIONS = {}
        if 'constant_Conc_Nitrates' in scenario_parameters:
            N_FERTILIZATIONS = {'constant_Conc_Nitrates': scenario_parameters.get('constant_Conc_Nitrates')}
        # Setup N_fertilizations time if time interval is given:
        if 'fertilization_interval' in scenario_parameters and 'fertilization_quantity' in scenario_parameters:
            fertilization_times = range(0, SIMULATION_LENGTH, scenario_parameters['fertilization_interval'])
            N_FERTILIZATIONS = {t: scenario_parameters['fertilization_quantity'] for t in fertilization_times}
        if 'fertilization_interval' in scenario_parameters and 'fertilization_expo_rate' in scenario_parameters:
            fertilization_times = range(0, SIMULATION_LENGTH, scenario_parameters['fertilization_interval'])
            K = scenario_parameters['fertilization_expo_rate']['K']
            V0 = scenario_parameters['fertilization_expo_rate']['V0']
            dt = scenario_parameters['fertilization_interval']
            N_FERTILIZATIONS = {t: exponential_fertilization_rate(V0=V0, K=K, t=t, dt=dt, plant_density=PLANT_DENSITY[1]) for t in fertilization_times}
            N0 = scenario_parameters['fertilization_expo_rate'].get('N0_U', 0) * (10 ** 5) / 14  # Initial nitrates in the soil (conversion from kg N ha-1 to µmol m-2)
            N_FERTILIZATIONS[0] += N0

        # -- RUN main fspmwheat --
        try:
            main.main(simulation_length=SIMULATION_LENGTH, forced_start_time=0,
                      run_simu=RUN_SIMU, run_postprocessing=RUN_POSTPROCESSING,
                      generate_graphs=GENERATE_GRAPHS, run_from_outputs=RUN_FROM_OUTPUTS,
                      show_3Dplant=False, heterogeneous_canopy=True,
                      N_fertilizations=N_FERTILIZATIONS,
                      PLANT_DENSITY=PLANT_DENSITY,
                      INTER_ROW=INTER_ROW,
                      INPUTS_PLANTSOIL_DIRPATH=SCENARIO_INPUTS_PLANTSOIL_DIRPATH,
                      INPUT_METEO_DIRPATH=SCENARIO_INPUT_METEO_DIRPATH,
                      METEO_FILENAME=scenario.get('METEO_FILENAME'),
                      GRAPHS_DIRPATH=scenario_graphs_dirpath,
                      OUTPUTS_DIRPATH=scenario_outputs_dirpath,
                      POSTPROCESSING_DIRPATH=scenario_postprocessing_dirpath,
                      update_parameters_all_models=scenario_parameters)
            if GENERATE_GRAPHS:
                additional_graphs.graph_summary(scenario_id, scenario_graphs_dirpath,
                                                graph_list=['LAI', 'sum_dry_mass_axis', 'shoot_roots_ratio_axis', 'N_content_shoot_axis', 'Conc_Amino_acids_phloem', 'Conc_Sucrose_phloem', 'leaf_Lmax',
                                                            'green_area_blade'])
            if RUN_POSTPROCESSING:
                fspmwheat_postprocessing.leaf_traits(scenario_outputs_dirpath, scenario_postprocessing_dirpath)
                fspmwheat_postprocessing.table_C_usages(scenario_postprocessing_dirpath)
                fspmwheat_postprocessing.calculate_performance_indices(scenario_outputs_dirpath, scenario_postprocessing_dirpath, os.path.join(INPUTS_DIRPATH, scenario.get('METEO_FILENAME')),
                                                                       scenario.get('Plant_Density', 250.))
                fspmwheat_postprocessing.canopy_dynamics(scenario_postprocessing_dirpath, os.path.join(INPUTS_DIRPATH, scenario.get('METEO_FILENAME')),
                                                         scenario.get('Plant_Density', 250.))

        except Exception as ex:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            template = "An exception of type {0} occurred. Arguments:\n{1!r}"
            message = template.format(type(ex).__name__, ex.args)
            print(message, fname, exc_tb.tb_lineno)


if __name__ == '__main__':
    inputs = None
    outputs = None
    scenario = 2

    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:o:s:d", ["inputs=", "outputs=", "scenario="])
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-i", "--inputs"):
            inputs = arg
        elif opt in ("-o", "--outputs"):
            outputs = arg
        elif opt in ("-s", "--scenario"):
            scenario = int(arg)

    run_fspmwheat(inputs_dir_path=inputs, outputs_dir_path=outputs, scenario_id=scenario)
