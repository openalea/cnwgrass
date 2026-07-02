# -*- coding: latin-1 -*-

from __future__ import print_function

import os
import sys
import getopt

import pandas as pd

from openalea.cnwgrass.integration.runner import run as runner

from example.Scenarios_monoculms import tools
from example.Scenarios_monoculms import additional_postprocessing
from example.Scenarios_monoculms import additional_graphs


"""
    run_scenario
    ~~~~~~~~~~~~

    Script adapted from the paper Gauthier et al. 2021 (https://doi.org/10.1093/insilicoplants/diab034), simulating the vegetative response of wheat monoculms grown
    in contrasting densities, N fertilisation and light conditions.
    This example also illustrates how the model can be run from a computing server.


"""

def ADEL_MTG(adel_wheat, INPUTS_DIRPATH, nff):
    """
    Builds a new ADEL mtg

    :param alinea.adel.adel_dynamic.AdelDyn adel_wheat: AdelDyn object
    :param str INPUTS_DIRPATH: dirpath of axeTable.csv file
    :param int nff: final leaf number

    :return: g
    """

    axeT_df = pd.read_csv(os.path.join(INPUTS_DIRPATH, 'axeTable.csv'))
    # Final leaf number
    axeT_df.HS_final = nff
    axeT_df.nf = nff
    axeT_df.nf_end = nff

    # mtg
    g = adel_wheat.build_stand(axeT_df)

    return g


def run_scenario(scenario_data, inputs_dir_path, outputs_dir_path='outputs'):
    """
    Run the simulation using data from a specific scenario

    :param pandas.Series scenario_data: a Series with the input data to be used for a given scenario
    :param str or None inputs_dir_path: the path directory of inputs
    :param str or None outputs_dir_path: the path to save outputs
    """

    # Path of the directory which contains the inputs of the model
    if inputs_dir_path:
        INPUTS_DIRPATH = inputs_dir_path
    else:
        INPUTS_DIRPATH = scenario_data.get('INPUTS_DIRPATH', 'inputs')


    # -- SIMULATION PARAMETERS --
    scenario_name = str(scenario_data.name)

    # Create dict of parameters for the scenario
    scenario_parameters = tools.buildDic(scenario_data.to_dict())

    # Do run the simulation?
    RUN_SIMU = scenario_parameters.get('Run_Simulation', True)

    SIMULATION_LENGTH = scenario_parameters.get('Simulation_Length', 3000)

    # Do run the simulation from the output files ?
    RUN_FROM_OUTPUTS = scenario_parameters.get('Run_From_Outputs', False)

    # Do run the postprocessing?
    RUN_POSTPROCESSING = scenario_parameters.get('Run_Postprocessing', True)

    # Do generate the graphs?
    GENERATE_GRAPHS = scenario_parameters.get('Generate_Graphs', False)

    if RUN_SIMU or RUN_POSTPROCESSING or GENERATE_GRAPHS:

        # -- SIMULATION DIRECTORIES --

        # Path of the directory which contains the outputs of the model
        if not outputs_dir_path:
            OUTPUTS_DIRPATH = 'outputs'
        else:
            OUTPUTS_DIRPATH = outputs_dir_path
        scenario_dirpath = os.path.join(OUTPUTS_DIRPATH, scenario_name)

        if not os.path.exists(OUTPUTS_DIRPATH):
            os.mkdir(OUTPUTS_DIRPATH)

        # Create the directory of the Scenario where results will be stored
        if not os.path.exists(scenario_dirpath):
            os.mkdir(scenario_dirpath)

        # Create directory paths for graphs, outputs and postprocessing of this scenario
        scenario_graphs_dirpath = os.path.join(scenario_dirpath, 'graphs')
        if not os.path.exists(scenario_graphs_dirpath):
            os.mkdir(scenario_graphs_dirpath)
        # Outputs
        scenario_outputs_dirpath = os.path.join(scenario_dirpath, 'outputs')
        if not os.path.exists(scenario_outputs_dirpath):
            os.mkdir(scenario_outputs_dirpath)
        # Postprocessing
        scenario_postprocessing_dirpath = os.path.join(scenario_dirpath, 'postprocessing')
        if not os.path.exists(scenario_postprocessing_dirpath):
            os.mkdir(scenario_postprocessing_dirpath)


        # -- RUN main integration --
        try:
            runner(simulation_length=SIMULATION_LENGTH, forced_start_time=0,
                      run_simu=RUN_SIMU, run_postprocessing=RUN_POSTPROCESSING,generate_graphs=GENERATE_GRAPHS, run_from_outputs=RUN_FROM_OUTPUTS,
                      show_3Dplant=False, heterogeneous_canopy=True,
                      INPUTS_DIRPATH=INPUTS_DIRPATH,
                      METEO_FILENAME=scenario_data.get('METEO_FILENAME'), MANAGEMENT_FILENAME=scenario_data.get('Management_filename'),
                      GRAPHS_DIRPATH=scenario_graphs_dirpath,
                      OUTPUTS_DIRPATH=scenario_outputs_dirpath,
                      POSTPROCESSING_DIRPATH=scenario_postprocessing_dirpath,
                      update_parameters_all_models=scenario_parameters,
                      step_callback={'ADEL_mtg': ADEL_MTG})
            if GENERATE_GRAPHS:
                additional_graphs.graph_summary(scenario_id, scenario_graphs_dirpath,
                                                graph_list=['LAI', 'sum_dry_mass_axis', 'shoot_roots_ratio_axis', 'N_content_shoot_axis', 'Conc_Amino_acids_phloem', 'Conc_Sucrose_phloem', 'leaf_Lmax',
                                                            'green_area_blade'])
            if RUN_POSTPROCESSING:
                additional_postprocessing.leaf_traits(scenario_outputs_dirpath, scenario_postprocessing_dirpath)
                additional_postprocessing.table_C_usages(scenario_postprocessing_dirpath)
                additional_postprocessing.calculate_performance_indices(scenario_outputs_dirpath, scenario_postprocessing_dirpath, os.path.join(INPUTS_DIRPATH, scenario_data.get('METEO_FILENAME')),
                                                                       scenario_data.get('Plant_Density', 250.))
                additional_postprocessing.canopy_dynamics(scenario_postprocessing_dirpath, os.path.join(INPUTS_DIRPATH, scenario_data.get('METEO_FILENAME')),
                                                         scenario_data.get('Plant_Density', 250.))

        except Exception as ex:
            print(f"Scenario {scenario_name} :")
            raise ex


if __name__ == '__main__':
    # Default argument
    inputs = None
    outputs = None
    scenario_id_argument = 1#None  # None for multiprocessing or a single scenario ID

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

    scenarios_list_path = os.path.join('inputs', 'scenarios_list.csv')
    scenarios_df = pd.read_csv(scenarios_list_path, index_col='Scenario')

    # ==========================================================================================
    # CASE 1 : A unique scenario is requested (local test or run from a computing server batch)#
    # ==========================================================================================
    if scenario_id_argument is not None:
        if scenario_id_argument not in scenarios_df.index:
            raise ValueError(f"Scenario {scenario_id_argument} does not exist in the CSV file.")
        scenario_df = scenarios_df.loc[scenario_id_argument]

        run_scenario(inputs_dir_path=inputs, outputs_dir_path=outputs, scenario_data=scenario_df)

    # =========================================
    # CASE 2 : For a local multiprocessing run#
    # =========================================
    else:
        import multiprocessing as mp
        scenarios_to_run = [1, 6]

        pool_args = []
        for scenario_id in scenarios_to_run:
            if scenario_id in scenarios_df.index:
                scenario_df = scenarios_df.loc[scenario_id]
                pool_args.append((scenario_df, inputs, outputs))

        num_processes = min(mp.cpu_count() - 1, len(pool_args))
        with mp.Pool(num_processes) as p:
            p.starmap(run_scenario, pool_args)
