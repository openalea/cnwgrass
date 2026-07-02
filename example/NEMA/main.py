# -*- coding: latin-1 -*-

from __future__ import print_function
import datetime
import logging
import os
import random
import time
import sys

import numpy as np
import pandas as pd

from openalea.adel.adel_dynamic import AdelDyn
from openalea.cnwgrass.cnmetabolism import tools as cnmetabolism_tools
from openalea.cnwgrass.integration import cnmetabolism_facade, gasexchange_facade, senescence_facade, growth_facade, caribu_facade, morphogenesis_facade

"""
    main
    ~~~~

    Script adapted from the paper Barillot et al. 2016 (https://doi.org/10.1093/aob/mcw144), 
    simulating postflowering stages of wheat plants grown at different N fertilisation regimes.
    This example uses the format MTG to exchange data between the models, while the original study only involved the current CN-Metabolism model.


"""

random.seed(1234)
np.random.seed(1234)

HOUR_TO_SECOND_CONVERSION_FACTOR = 3600

AXES_INDEX_COLUMNS = ['t', 'plant', 'axis']
ELEMENTS_INDEX_COLUMNS = ['t', 'plant', 'axis', 'metamer', 'organ', 'element']
ORGANS_INDEX_COLUMNS = ['t', 'plant', 'axis', 'organ']
SOILS_INDEX_COLUMNS = ['t', 'plant', 'axis']

# Define culm density (culm m-2)
DENSITY = 410.
NPLANTS = 1
CULM_DENSITY = {i: DENSITY / NPLANTS for i in range(1, NPLANTS + 1)}

INPUTS_OUTPUTS_PRECISION = 5  # 10

LOGGING_CONFIG_FILEPATH = os.path.join('logging.json')

LOGGING_LEVEL = logging.WARNING  # can be one of: DEBUG, INFO, WARNING, ERROR, CRITICAL

cnmetabolism_tools.setup_logging(LOGGING_CONFIG_FILEPATH, LOGGING_LEVEL, log_model=False, log_compartments=False, log_derivatives=False)


def calculate_PARa_from_df(g, Eabs_df, PARi, multiple_sources=False, ratio_diffus_PAR=None):
    """
    Compute PARa from an input dataframe having Eabs values.

    :param openalea.mtg.mtg.MTG g: the MTG of the simulation
    :param pandas.DataFrame Eabs_df: the dataframe with Eabs values of each element calculated a priori
    :param float PARi: the incident PAR at time t (µmol m-2 s-1)
    :param bool multiple_sources: whether to use direct + diffuse sources
    :param float or None ratio_diffus_PAR: the ratio of diffuse PAR to incident PAR
    """

    Eabs_df_grouped = Eabs_df.groupby(['plant', 'metamer', 'organ'])

    #: the name of the elements modeled by Gas-Exchange
    CARIBU_ELEMENTS_NAMES = {'StemElement', 'LeafElement1'}

    PARa_element_data_dict = {}
    # traverse the MTG recursively from top ...
    for mtg_plant_vid in g.components_iter(g.root):
        mtg_plant_index = int(g.index(mtg_plant_vid))
        for mtg_axis_vid in g.components_iter(mtg_plant_vid):
            for mtg_metamer_vid in g.components_iter(mtg_axis_vid):
                mtg_metamer_index = int(g.index(mtg_metamer_vid))
                for mtg_organ_vid in g.components_iter(mtg_metamer_vid):
                    mtg_organ_label = g.label(mtg_organ_vid)
                    for mtg_element_vid in g.components_iter(mtg_organ_vid):
                        mtg_element_label = g.label(mtg_element_vid)
                        if mtg_element_label not in CARIBU_ELEMENTS_NAMES:
                            continue
                        element_id = (mtg_plant_index, mtg_metamer_index, mtg_organ_label)
                        if element_id in Eabs_df_grouped.groups.keys():
                            if PARi == 0:
                                PARa_element_data_dict[mtg_element_vid] = 0
                            elif multiple_sources:
                                PARa_diffuse = Eabs_df_grouped.get_group(element_id)['Eabs_diffuse'].iloc[0] * PARi * ratio_diffus_PAR
                                PARa_direct = Eabs_df_grouped.get_group(element_id)['Eabs_direct'].iloc[0] * PARi * (1 - ratio_diffus_PAR)
                                PARa_element_data_dict[mtg_element_vid] = PARa_diffuse + PARa_direct
                            else:
                                PARa_element_data_dict[mtg_element_vid] = Eabs_df_grouped.get_group(element_id)['Eabs'].iloc[0] * PARi

    return PARa_element_data_dict


def main(simulation_length, NEMA_scenario, run_simu=True, make_graphs=True):
    """
    Run the simulation for a given NEMA scenario

    :param int simulation_length: duration of the simulation (hours)
    :param str NEMA_scenario: name of NEMA scenario to be run (must be "NEMA_HO", "NEMA_H3" or "NEMA_H15")
    :param bool run_simu: whether to run the simulation
    :param bool make_graphs: whether to make the graphs

    """

    INPUTS_DIRPATH = os.path.join(NEMA_scenario, 'inputs')
    GRAPHS_DIRPATH = os.path.join(NEMA_scenario, 'graphs')

    # adelwheat inputs at t0
    ADELWHEAT_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH,
                                            'adelwheat')  # the directory adelwheat must contain files 'adel0000.pckl' and 'scene0000.bgeom'

    # cnmetabolism inputs at t0
    CNMETABOLISM_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH, 'cnmetabolism')
    CNMETABOLISM_AXES_INPUTS_FILEPATH = os.path.join(CNMETABOLISM_INPUTS_DIRPATH, 'SAM_inputs.csv')
    CNMETABOLISM_ORGANS_INPUTS_FILEPATH = os.path.join(CNMETABOLISM_INPUTS_DIRPATH, 'organs_inputs.csv')
    CNMETABOLISM_HIDDENZONE_INPUTS_FILEPATH = os.path.join(CNMETABOLISM_INPUTS_DIRPATH, 'hiddenzones_inputs.csv')
    CNMETABOLISM_ELEMENTS_INPUTS_FILEPATH = os.path.join(CNMETABOLISM_INPUTS_DIRPATH, 'elements_inputs.csv')
    CNMETABOLISM_SOILS_INPUTS_FILEPATH = os.path.join(CNMETABOLISM_INPUTS_DIRPATH, 'soils_inputs.csv')

    # gasexchange inputs at t0
    GASEXCHANGE_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH, 'gasexchange')
    GASEXCHANGE_INPUTS_FILEPATH = os.path.join(GASEXCHANGE_INPUTS_DIRPATH, 'inputs.csv')
    GASEXCHANGE_AXES_INPUTS_FILEPATH = os.path.join(GASEXCHANGE_INPUTS_DIRPATH, 'SAM_inputs.csv')
    METEO_FILEPATH = os.path.join(GASEXCHANGE_INPUTS_DIRPATH, 'meteo_Clermont_rebuild.csv')
    CARIBU_FILEPATH = os.path.join(GASEXCHANGE_INPUTS_DIRPATH, 'inputs_eabs.csv')

    # morphogenesis inputs at t0
    MORPHOGENESIS_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH, 'morphogenesis')
    MORPHOGENESIS_HZ_INPUTS_FILEPATH = os.path.join(MORPHOGENESIS_INPUTS_DIRPATH, 'hiddenzones_inputs.csv')
    MORPHOGENESIS_ELEMENTS_INPUTS_FILEPATH = os.path.join(MORPHOGENESIS_INPUTS_DIRPATH, 'elements_inputs.csv')
    MORPHOGENESIS_AXES_INPUTS_FILEPATH = os.path.join(MORPHOGENESIS_INPUTS_DIRPATH, 'SAM_inputs.csv')

    # senescence inputs at t0
    SENESCENCE_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH, 'senescence')
    SENESCENCE_ROOTS_INPUTS_FILEPATH = os.path.join(SENESCENCE_INPUTS_DIRPATH, 'roots_inputs.csv')
    SENESCENCE_AXES_INPUTS_FILEPATH = os.path.join(SENESCENCE_INPUTS_DIRPATH, 'SAM_inputs.csv')
    SENESCENCE_ELEMENTS_INPUTS_FILEPATH = os.path.join(SENESCENCE_INPUTS_DIRPATH, 'elements_inputs.csv')

    # growth inputs at t0
    GROWTH_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH, 'growth')
    GROWTH_HIDDENZONE_INPUTS_FILEPATH = os.path.join(GROWTH_INPUTS_DIRPATH, 'hiddenzones_inputs.csv')
    GROWTH_ORGANS_INPUTS_FILEPATH = os.path.join(GROWTH_INPUTS_DIRPATH, 'organs_inputs.csv')
    GROWTH_ROOTS_INPUTS_FILEPATH = os.path.join(GROWTH_INPUTS_DIRPATH, 'roots_inputs.csv')
    GROWTH_AXES_INPUTS_FILEPATH = os.path.join(GROWTH_INPUTS_DIRPATH, 'SAM_inputs.csv')

    # the path of the CSV files where to save the states of the modeled system at each step
    OUTPUTS_DIRPATH = os.path.join(NEMA_scenario, 'outputs')
    AXES_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'axes_states.csv')
    ORGANS_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'organs_states.csv')
    ELEMENTS_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'elements_states.csv')
    SOILS_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'soils_states.csv')

    # post-processing directory path
    POSTPROCESSING_DIRPATH = os.path.join(NEMA_scenario, 'postprocessing')
    AXES_POSTPROCESSING_FILEPATH = os.path.join(POSTPROCESSING_DIRPATH, 'axes_postprocessing.csv')
    ORGANS_POSTPROCESSING_FILEPATH = os.path.join(POSTPROCESSING_DIRPATH, 'organs_postprocessing.csv')
    ELEMENTS_POSTPROCESSING_FILEPATH = os.path.join(POSTPROCESSING_DIRPATH, 'elements_postprocessing.csv')
    SOILS_POSTPROCESSING_FILEPATH = os.path.join(POSTPROCESSING_DIRPATH, 'soils_postprocessing.csv')

    if run_simu:
        print(run_simu)
        meteo = pd.read_csv(METEO_FILEPATH, index_col='t')
        Eabs_df = pd.read_csv(CARIBU_FILEPATH)

        # define the time step in hours for each simulator
        senescence_ts = 2
        growth_ts = 2
        gasexchange_ts = 2
        morphogenesis_ts = 2
        cnmetabolism_ts = 1

        hour_to_second_conversion_factor = 3600

        # read adelwheat inputs at t0
        adel_wheat = AdelDyn(seed=1234, scene_unit='m')
        g = adel_wheat.load(directory=ADELWHEAT_INPUTS_DIRPATH)

        # adel_wheat.plot(g)

        # create empty dataframes to shared data between the models
        shared_axes_inputs_outputs_df = pd.DataFrame()
        shared_organs_inputs_outputs_df = pd.DataFrame()
        shared_hiddenzones_inputs_outputs_df = pd.DataFrame()
        shared_elements_inputs_outputs_df = pd.DataFrame()
        shared_soils_inputs_outputs_df = pd.DataFrame()

        # read the inputs at t0 and create the facades

        # caribu
        caribu_facade_ = caribu_facade.CaribuFacade(g,
                                                    shared_elements_inputs_outputs_df,
                                                    adel_wheat)

        # senescence
        senescence_roots_inputs_t0 = pd.read_csv(SENESCENCE_ROOTS_INPUTS_FILEPATH)
        senescence_axes_inputs_t0 = pd.read_csv(SENESCENCE_AXES_INPUTS_FILEPATH)
        senescence_elements_inputs_t0 = pd.read_csv(SENESCENCE_ELEMENTS_INPUTS_FILEPATH)
        senescence_facade_ = senescence_facade.SENESCENCEFacade(g,
                                                                senescence_ts * hour_to_second_conversion_factor,
                                                                senescence_roots_inputs_t0,
                                                                senescence_axes_inputs_t0,
                                                                senescence_elements_inputs_t0,
                                                                shared_organs_inputs_outputs_df,
                                                                shared_axes_inputs_outputs_df,
                                                                shared_elements_inputs_outputs_df)
        # growth
        growth_hiddenzones_inputs_t0 = pd.read_csv(GROWTH_HIDDENZONE_INPUTS_FILEPATH)
        growth_organ_inputs_t0 = pd.read_csv(GROWTH_ORGANS_INPUTS_FILEPATH)
        growth_root_inputs_t0 = pd.read_csv(GROWTH_ROOTS_INPUTS_FILEPATH)
        growth_axes_inputs_t0 = pd.read_csv(GROWTH_AXES_INPUTS_FILEPATH)
        growth_facade_ = growth_facade.GrowthFacade(g,
                                                    growth_ts * hour_to_second_conversion_factor,
                                                    growth_hiddenzones_inputs_t0,
                                                    growth_organ_inputs_t0,
                                                    growth_root_inputs_t0,
                                                    growth_axes_inputs_t0,
                                                    shared_organs_inputs_outputs_df,
                                                    shared_hiddenzones_inputs_outputs_df,
                                                    shared_elements_inputs_outputs_df,
                                                    shared_axes_inputs_outputs_df)

        # gasexchange
        gasexchange_elements_inputs_t0 = pd.read_csv(GASEXCHANGE_INPUTS_FILEPATH)
        gasexchange_axes_inputs_t0 = pd.read_csv(GASEXCHANGE_AXES_INPUTS_FILEPATH)
        # Use the initial version of the photosynthesis sub-model (as in Barillot et al. 2016, and in Gauthier et al. 2020)
        update_parameters_gasexchange = {'SurfacicProteins': False, 'NSC_Retroinhibition': False}

        gasexchange_facade_ = gasexchange_facade.GasExchangeFacade(g,
                                                                   gasexchange_elements_inputs_t0,
                                                                   gasexchange_axes_inputs_t0,
                                                                   shared_elements_inputs_outputs_df,
                                                                   update_parameters=update_parameters_gasexchange)

        # morphogenesis # Only for temperature related computations
        morphogenesis_hiddenzones_inputs_t0 = pd.read_csv(MORPHOGENESIS_HZ_INPUTS_FILEPATH)
        morphogenesis_elements_inputs_t0 = pd.read_csv(MORPHOGENESIS_ELEMENTS_INPUTS_FILEPATH)
        morphogenesis_axes_inputs_t0 = pd.read_csv(MORPHOGENESIS_AXES_INPUTS_FILEPATH)

        morphogenesis_facade_ = morphogenesis_facade.MorphogenesisFacade(g,
                                                                         morphogenesis_ts * hour_to_second_conversion_factor,
                                                                         morphogenesis_axes_inputs_t0,
                                                                         morphogenesis_hiddenzones_inputs_t0,
                                                                         morphogenesis_elements_inputs_t0,
                                                                         shared_axes_inputs_outputs_df,
                                                                         shared_hiddenzones_inputs_outputs_df,
                                                                         shared_elements_inputs_outputs_df,
                                                                         adel_wheat, option_static=True)
        # cnmetabolism
        cnmetabolism_axes_inputs_t0 = pd.read_csv(CNMETABOLISM_AXES_INPUTS_FILEPATH)
        cnmetabolism_organs_inputs_t0 = pd.read_csv(CNMETABOLISM_ORGANS_INPUTS_FILEPATH)
        cnmetabolism_hiddenzones_inputs_t0 = pd.read_csv(CNMETABOLISM_HIDDENZONE_INPUTS_FILEPATH)
        cnmetabolism_elements_inputs_t0 = pd.read_csv(CNMETABOLISM_ELEMENTS_INPUTS_FILEPATH)
        cnmetabolism_soils_inputs_t0 = pd.read_csv(CNMETABOLISM_SOILS_INPUTS_FILEPATH)
        update_cnmetabolism_parameters = {'roots': {'K_AMINO_ACIDS_EXPORT': 25*3E-5,
                                               'K_NITRATE_EXPORT': 25*1E-6}}

        cnmetabolism_facade_ = cnmetabolism_facade.CNMetabolismFacade(g,
                                                                      cnmetabolism_ts * hour_to_second_conversion_factor,
                                                                      CULM_DENSITY,
                                                                      update_cnmetabolism_parameters,
                                                                      cnmetabolism_axes_inputs_t0,
                                                                      cnmetabolism_organs_inputs_t0,
                                                                      cnmetabolism_hiddenzones_inputs_t0,
                                                                      cnmetabolism_elements_inputs_t0,
                                                                      cnmetabolism_soils_inputs_t0,
                                                                      shared_axes_inputs_outputs_df,
                                                                      shared_organs_inputs_outputs_df,
                                                                      shared_hiddenzones_inputs_outputs_df,
                                                                      shared_elements_inputs_outputs_df,
                                                                      shared_soils_inputs_outputs_df)

        # adel_wheat.update_geometry(g) # NE FONCTIONNE PAS car MTG non compatible (pas de top et base element)
        # adel_wheat.plot(g)

        # define organs for which the variable 'max_proteins' is fixed
        forced_max_protein_elements = {(1, 'MS', 9, 'blade', 'LeafElement1'), (1, 'MS', 10, 'blade', 'LeafElement1'), (1, 'MS', 11, 'blade', 'LeafElement1'), (2, 'MS', 9, 'blade', 'LeafElement1'),
                                       (2, 'MS', 10, 'blade', 'LeafElement1'), (2, 'MS', 11, 'blade', 'LeafElement1')}

        # define the start and the end of the whole simulation (in hours)
        start_time = 0

        # define lists of dataframes to store the inputs and the outputs of the models at each step.
        axes_all_data_list = []
        organs_all_data_list = []  # organs which belong to axes: roots, phloem, grains
        elements_all_data_list = []
        soils_all_data_list = []

        all_simulation_steps = []  # to store the steps of the simulation

        # run the simulators
        current_time_of_the_system = time.time()

        try:

            for t_morphogenesis in range(start_time, simulation_length, morphogenesis_ts):  # Only to compute temperature related variable

                # run Morphogenesis
                print('t morphogenesis is {}'.format(t_morphogenesis))
                Tair, Tsoil = meteo.loc[t_morphogenesis, ['air_temperature', 'air_temperature']]
                morphogenesis_facade_.run(Tair, Tsoil)

                for t_senescence in range(t_morphogenesis, t_morphogenesis + morphogenesis_ts, senescence_ts):

                    # run Senescence
                    print('t senescence is {}'.format(t_senescence))
                    senescence_facade_.run(forced_max_protein_elements, postflowering_stages=True)

                    # Test for fully senesced shoot tissues  #TODO: Make the model to work even if the whole shoot is dead but the roots are alived
                    if sum(senescence_facade_._shared_elements_inputs_outputs_df['green_area']) <= 0.25E-6:
                        break

                    for t_growth in range(t_senescence, t_senescence + senescence_ts, growth_ts):
                        # run Growth
                        print('t growth is {}'.format(t_growth))
                        growth_facade_.run(postflowering_stages=True)

                        for t_gasexchange in range(t_growth, t_growth + growth_ts, gasexchange_ts):
                            # get the meteo of the current step
                            Tair, ambient_CO2, RH, Ur, PARi = meteo.loc[t_gasexchange, ['air_temperature', 'ambient_CO2', 'humidity', 'Wind', 'PARi']]
                            # get PARa for current step
                            aggregated_PARa = calculate_PARa_from_df(g, Eabs_df, PARi, multiple_sources=False)
                            print('t caribu is {}'.format(t_gasexchange))
                            # caribu_facade_.run(energy=PARi,sun_sky_option='sky')
                            caribu_facade_.update_shared_MTG({'PARa': aggregated_PARa})
                            caribu_facade_.update_shared_dataframes({'PARa':aggregated_PARa})
                            # run Gas-Exchange
                            print('t farquhar is {}'.format(t_gasexchange))
                            gasexchange_facade_.run(Tair, ambient_CO2, RH, Ur)

                            for t_cnmetabolism in range(t_gasexchange, t_gasexchange + senescence_ts, cnmetabolism_ts):
                                Tair, Tsoil = meteo.loc[t_cnmetabolism, ['air_temperature', 'air_temperature']]
                                # run CN-Metabolism
                                print('t cnmetabolism is {}'.format(t_cnmetabolism))
                                cnmetabolism_facade_.run(Tair=Tair, Tsoil=Tsoil)

                                # append the inputs and outputs at current step to global lists
                                all_simulation_steps.append(t_cnmetabolism)
                                axes_all_data_list.append(shared_axes_inputs_outputs_df.copy())
                                organs_all_data_list.append(shared_organs_inputs_outputs_df.copy())
                                elements_all_data_list.append(shared_elements_inputs_outputs_df.copy())
                                soils_all_data_list.append(shared_soils_inputs_outputs_df.copy())
                else:
                    continue
                break

            execution_time = int(time.time() - current_time_of_the_system)
            print('\n', 'Simulation run in ', str(datetime.timedelta(seconds=execution_time)))

        except Exception as ex:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            template = "An exception of type {0} occurred. Arguments:\n{1!r}"
            message = template.format(type(ex).__name__, ex.args)
            print(message, fname, exc_tb.tb_lineno)

        finally:
            # write all inputs and outputs to CSV files
            all_axes_inputs_outputs = pd.concat(axes_all_data_list, keys=all_simulation_steps)
            all_axes_inputs_outputs.reset_index(0, inplace=True)
            all_axes_inputs_outputs.rename({'level_0': 't'}, axis=1, inplace=True)
            all_axes_inputs_outputs.to_csv(AXES_STATES_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION))

            all_organs_inputs_outputs = pd.concat(organs_all_data_list, keys=all_simulation_steps)
            all_organs_inputs_outputs.reset_index(0, inplace=True)
            all_organs_inputs_outputs.rename({'level_0': 't'}, axis=1, inplace=True)
            all_organs_inputs_outputs.to_csv(ORGANS_STATES_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION))

            all_elements_inputs_outputs = pd.concat(elements_all_data_list, keys=all_simulation_steps)
            all_elements_inputs_outputs.reset_index(0, inplace=True)
            all_elements_inputs_outputs.rename({'level_0': 't'}, axis=1, inplace=True)
            all_elements_inputs_outputs.to_csv(ELEMENTS_STATES_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION))

            all_soils_inputs_outputs = pd.concat(soils_all_data_list, keys=all_simulation_steps)
            all_soils_inputs_outputs.reset_index(0, inplace=True)
            all_soils_inputs_outputs.rename({'level_0': 't'}, axis=1, inplace=True)
            all_soils_inputs_outputs.to_csv(SOILS_STATES_FILEPATH, na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION))

    if make_graphs:
        graphs = os.listdir(GRAPHS_DIRPATH)
        for graph in graphs:
            if graph.endswith(".PNG"):
                os.remove(os.path.join(GRAPHS_DIRPATH, graph))

        # -POST-PROCESSING##
        states_df_dict = {}
        for states_filepath in (AXES_STATES_FILEPATH,
                                ORGANS_STATES_FILEPATH,
                                ELEMENTS_STATES_FILEPATH,
                                SOILS_STATES_FILEPATH):
            # assert states_filepaths were not opened during simulation run meaning that other filenames were saved
            path, filename = os.path.split(states_filepath)
            filename = os.path.splitext(filename)[0]
            newfilename = 'ACTUAL_{}.csv'.format(filename)
            newpath = os.path.join(path, newfilename)
            assert not os.path.isfile(newpath), \
                "File {} was saved because {} was opened during simulation run. Rename it before running postprocessing".format(newfilename, states_filepath)

            # Retrieve outputs dataframes from precedent simulation run
            states_df = pd.read_csv(states_filepath)
            states_file_basename = os.path.basename(states_filepath).split('.')[0]
            states_df_dict[states_file_basename] = states_df
        time_grid = states_df_dict['elements_states']['t']
        delta_t = (time_grid.unique()[1] - time_grid.unique()[0]) * HOUR_TO_SECOND_CONVERSION_FACTOR

        # run the postprocessing
        axes_postprocessing_file_basename = os.path.basename(AXES_POSTPROCESSING_FILEPATH).split('.')[0]
        organs_postprocessing_file_basename = os.path.basename(ORGANS_POSTPROCESSING_FILEPATH).split('.')[0]
        elements_postprocessing_file_basename = os.path.basename(ELEMENTS_POSTPROCESSING_FILEPATH).split('.')[0]
        soils_postprocessing_file_basename = os.path.basename(SOILS_POSTPROCESSING_FILEPATH).split('.')[0]

        postprocessing = cnmetabolism_facade.CNMetabolismFacade.postprocessing(axes_outputs_df=states_df_dict[os.path.basename(AXES_STATES_FILEPATH).split('.')[0]],
                                                                               organs_outputs_df=states_df_dict[os.path.basename(ORGANS_STATES_FILEPATH).split('.')[0]],
                                                                               hiddenzone_outputs_df=None,
                                                                               elements_outputs_df=states_df_dict[os.path.basename(ELEMENTS_STATES_FILEPATH).split('.')[0]],
                                                                               soils_outputs_df=states_df_dict[os.path.basename(SOILS_STATES_FILEPATH).split('.')[0]],
                                                                               delta_t=delta_t)

        # save the postprocessing to disk
        for postprocessing_filepath, postprocessing_filename, index_columns in (('axes', AXES_POSTPROCESSING_FILEPATH, AXES_INDEX_COLUMNS),
                                                                                     ('organs', ORGANS_POSTPROCESSING_FILEPATH, ORGANS_INDEX_COLUMNS),
                                                                                     ('elements', ELEMENTS_POSTPROCESSING_FILEPATH, ELEMENTS_INDEX_COLUMNS),
                                                                                     ('soils', SOILS_POSTPROCESSING_FILEPATH, SOILS_INDEX_COLUMNS)):
            postprocessing_df = postprocessing[postprocessing_filepath]
            postprocessing_df.rename({'level_0': 't'}, axis=1, inplace=True)
            postprocessing_df = postprocessing_df.reindex(index_columns + postprocessing_df.columns.difference(index_columns).tolist(), axis=1, copy=False)
            postprocessing_df.to_csv(postprocessing_filename, na_rep='NA', index=False, float_format='%.{}f'.format(8))

        # - GRAPHS

        # Retrieve last computed post-processing dataframes
        axes_postprocessing_file_basename = os.path.basename(AXES_POSTPROCESSING_FILEPATH).split('.')[0]
        organs_postprocessing_file_basename = os.path.basename(ORGANS_POSTPROCESSING_FILEPATH).split('.')[0]
        elements_postprocessing_file_basename = os.path.basename(ELEMENTS_POSTPROCESSING_FILEPATH).split('.')[0]
        soils_postprocessing_file_basename = os.path.basename(SOILS_POSTPROCESSING_FILEPATH).split('.')[0]
        postprocessing_df_dict = {}
        for (postprocessing_filepath, postprocessing_file_basename) in ((AXES_POSTPROCESSING_FILEPATH, axes_postprocessing_file_basename),
                                                                        (ORGANS_POSTPROCESSING_FILEPATH, organs_postprocessing_file_basename),
                                                                        (ELEMENTS_POSTPROCESSING_FILEPATH, elements_postprocessing_file_basename),
                                                                        (SOILS_POSTPROCESSING_FILEPATH, soils_postprocessing_file_basename)):
            postprocessing_df = pd.read_csv(postprocessing_filepath)
            postprocessing_df_dict[postprocessing_file_basename] = postprocessing_df

        # Generate graphs
        cnmetabolism_facade.CNMetabolismFacade.graphs(axes_postprocessing_df=postprocessing_df_dict[axes_postprocessing_file_basename],
                                                      hiddenzones_postprocessing_df=None,
                                                      organs_postprocessing_df=postprocessing_df_dict[organs_postprocessing_file_basename],
                                                      elements_postprocessing_df=postprocessing_df_dict[elements_postprocessing_file_basename],
                                                      soils_postprocessing_df=postprocessing_df_dict[soils_postprocessing_file_basename],
                                                      graphs_dirpath=GRAPHS_DIRPATH)


if __name__ == '__main__':
    scenario_list = ['NEMA_H0', 'NEMA_H3', 'NEMA_H15']

    sim_length = 1200
    run_simulation = True
    make_graphs = True

    if len(scenario_list) == 1:
        main(sim_length, scenario_list[0], run_simu=run_simulation, make_graphs=make_graphs)
    else:
        import multiprocessing as mp
        num_processes = mp.cpu_count() - 1

        pool_args = [(sim_length, sc, run_simulation, make_graphs) for sc in scenario_list]
        with mp.Pool(num_processes) as p:
            p.starmap(main, pool_args)
