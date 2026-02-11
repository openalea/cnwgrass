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
from openalea.cnwheat import tools as cnwheat_tools
from openalea.fspmwheat import cnwheat_facade, farquharwheat_facade, senescwheat_facade, growthwheat_facade, caribu_facade, elongwheat_facade

"""
    main
    ~~~~

    Script readatpted from example NEMA_H3 used in the paper Barillot et al. (2016).
    This example uses the format MTG to exchange data between the models.

"""

random.seed(1234)
np.random.seed(1234)

HOUR_TO_SECOND_CONVERSION_FACTOR = 3600

INPUTS_DIRPATH = 'inputs'
GRAPHS_DIRPATH = 'graphs'

# adelwheat inputs at t0
ADELWHEAT_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH, 'adelwheat')  #  the directory adelwheat must contain files 'adel0000.pckl' and 'scene0000.bgeom'

# cnwheat inputs at t0
CNWHEAT_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH, 'cnwheat')
CNWHEAT_PLANTS_INPUTS_FILEPATH = os.path.join(CNWHEAT_INPUTS_DIRPATH, 'plants_inputs.csv')
CNWHEAT_AXES_INPUTS_FILEPATH = os.path.join(CNWHEAT_INPUTS_DIRPATH, 'axes_inputs.csv')
CNWHEAT_METAMERS_INPUTS_FILEPATH = os.path.join(CNWHEAT_INPUTS_DIRPATH, 'metamers_inputs.csv')
CNWHEAT_ORGANS_INPUTS_FILEPATH = os.path.join(CNWHEAT_INPUTS_DIRPATH, 'organs_inputs.csv')
CNWHEAT_HIDDENZONE_INPUTS_FILEPATH = os.path.join(CNWHEAT_INPUTS_DIRPATH, 'hiddenzones_inputs.csv')
CNWHEAT_ELEMENTS_INPUTS_FILEPATH = os.path.join(CNWHEAT_INPUTS_DIRPATH, 'elements_inputs.csv')
CNWHEAT_SOILS_INPUTS_FILEPATH = os.path.join(CNWHEAT_INPUTS_DIRPATH, 'soils_inputs.csv')

# farquharwheat inputs at t0
FARQUHARWHEAT_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH, 'farquharwheat')
FARQUHARWHEAT_INPUTS_FILEPATH = os.path.join(FARQUHARWHEAT_INPUTS_DIRPATH, 'inputs.csv')
FARQUHARWHEAT_AXES_INPUTS_FILEPATH = os.path.join(FARQUHARWHEAT_INPUTS_DIRPATH, 'SAM_inputs.csv')
METEO_FILEPATH = os.path.join(FARQUHARWHEAT_INPUTS_DIRPATH, 'meteo_Clermont_rebuild.csv')
CARIBU_FILEPATH = os.path.join(FARQUHARWHEAT_INPUTS_DIRPATH, 'inputs_eabs.csv')

# elongwheat inputs at t0
ELONGWHEAT_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH, 'elongwheat')
ELONGWHEAT_HZ_INPUTS_FILEPATH = os.path.join(ELONGWHEAT_INPUTS_DIRPATH, 'hiddenzones_inputs.csv')
ELONGWHEAT_ELEMENTS_INPUTS_FILEPATH = os.path.join(ELONGWHEAT_INPUTS_DIRPATH, 'elements_inputs.csv')
ELONGWHEAT_AXES_INPUTS_FILEPATH = os.path.join(ELONGWHEAT_INPUTS_DIRPATH, 'SAM_inputs.csv')

# senescwheat inputs at t0
SENESCWHEAT_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH, 'senescwheat')
SENESCWHEAT_ROOTS_INPUTS_FILEPATH = os.path.join(SENESCWHEAT_INPUTS_DIRPATH, 'roots_inputs.csv')
SENESCWHEAT_AXES_INPUTS_FILEPATH = os.path.join(SENESCWHEAT_INPUTS_DIRPATH, 'SAM_inputs.csv')
SENESCWHEAT_ELEMENTS_INPUTS_FILEPATH = os.path.join(SENESCWHEAT_INPUTS_DIRPATH, 'elements_inputs.csv')

# growthwheat inputs at t0
GROWTHWHEAT_INPUTS_DIRPATH = os.path.join(INPUTS_DIRPATH, 'growthwheat')
GROWTHWHEAT_HIDDENZONE_INPUTS_FILEPATH = os.path.join(GROWTHWHEAT_INPUTS_DIRPATH, 'hiddenzones_inputs.csv')
GROWTHWHEAT_ORGANS_INPUTS_FILEPATH = os.path.join(GROWTHWHEAT_INPUTS_DIRPATH, 'organs_inputs.csv')
GROWTHWHEAT_ROOTS_INPUTS_FILEPATH = os.path.join(GROWTHWHEAT_INPUTS_DIRPATH, 'roots_inputs.csv')
GROWTHWHEAT_AXES_INPUTS_FILEPATH = os.path.join(GROWTHWHEAT_INPUTS_DIRPATH, 'SAM_inputs.csv')

# the path of the CSV files where to save the states of the modeled system at each step
OUTPUTS_DIRPATH = 'outputs'
AXES_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'axes_states.csv')
ORGANS_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'organs_states.csv')
ELEMENTS_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'elements_states.csv')
SOILS_STATES_FILEPATH = os.path.join(OUTPUTS_DIRPATH, 'soils_states.csv')

# post-processing directory path
POSTPROCESSING_DIRPATH = 'postprocessing'
AXES_POSTPROCESSING_FILEPATH = os.path.join(POSTPROCESSING_DIRPATH, 'axes_postprocessing.csv')
ORGANS_POSTPROCESSING_FILEPATH = os.path.join(POSTPROCESSING_DIRPATH, 'organs_postprocessing.csv')
HIDDENZONES_POSTPROCESSING_FILEPATH = os.path.join(POSTPROCESSING_DIRPATH, 'hiddenzones_postprocessing.csv')
ELEMENTS_POSTPROCESSING_FILEPATH = os.path.join(POSTPROCESSING_DIRPATH, 'elements_postprocessing.csv')
SOILS_POSTPROCESSING_FILEPATH = os.path.join(POSTPROCESSING_DIRPATH, 'soils_postprocessing.csv')

AXES_INDEX_COLUMNS = ['t', 'plant', 'axis']
ELEMENTS_INDEX_COLUMNS = ['t', 'plant', 'axis', 'metamer', 'organ', 'element']
ORGANS_INDEX_COLUMNS = ['t', 'plant', 'axis', 'organ']
SOILS_INDEX_COLUMNS = ['t', 'plant', 'axis']

# Define culm density (culm m-2)
DENSITY = 410.
NPLANTS = 1
CULM_DENSITY = {i: DENSITY / NPLANTS for i in range(1, NPLANTS + 1)}

INPUTS_OUTPUTS_PRECISION = 5  # 10

# LOGGING_CONFIG_FILEPATH = os.path.join('logging.json')
#
# LOGGING_LEVEL = logging.INFO  # can be one of: DEBUG, INFO, WARNING, ERROR, CRITICAL
#
# cnwheat_tools.setup_logging(LOGGING_CONFIG_FILEPATH, LOGGING_LEVEL, log_model=False, log_compartments=False, log_derivatives=False)


def calculate_PARa_from_df(g, Eabs_df, PARi, multiple_sources=False, ratio_diffus_PAR=None):
    """
    Compute PARa from an input dataframe having Eabs values.
    """

    Eabs_df_grouped = Eabs_df.groupby(['plant', 'metamer', 'organ'])

    #: the name of the elements modeled by FarquharWheat
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


def main(stop_time, run_simu=True, make_graphs=True):
    if run_simu:
        meteo = pd.read_csv(METEO_FILEPATH, index_col='t')
        Eabs_df = pd.read_csv(CARIBU_FILEPATH)

        # define the time step in hours for each simulator
        senescwheat_ts = 2
        growthwheat_ts = 2
        farquharwheat_ts = 2
        elongwheat_ts = 2
        cnwheat_ts = 1

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

        # senescwheat
        senescwheat_roots_inputs_t0 = pd.read_csv(SENESCWHEAT_ROOTS_INPUTS_FILEPATH)
        senescwheat_axes_inputs_t0 = pd.read_csv(SENESCWHEAT_AXES_INPUTS_FILEPATH)
        senescwheat_elements_inputs_t0 = pd.read_csv(SENESCWHEAT_ELEMENTS_INPUTS_FILEPATH)
        senescwheat_facade_ = senescwheat_facade.SenescWheatFacade(g,
                                                                   senescwheat_ts * hour_to_second_conversion_factor,
                                                                   senescwheat_roots_inputs_t0,
                                                                   senescwheat_axes_inputs_t0,
                                                                   senescwheat_elements_inputs_t0,
                                                                   shared_organs_inputs_outputs_df,
                                                                   shared_axes_inputs_outputs_df,
                                                                   shared_elements_inputs_outputs_df)
        # growthwheat
        growthwheat_hiddenzones_inputs_t0 = pd.read_csv(GROWTHWHEAT_HIDDENZONE_INPUTS_FILEPATH)
        growthwheat_organ_inputs_t0 = pd.read_csv(GROWTHWHEAT_ORGANS_INPUTS_FILEPATH)
        growthwheat_root_inputs_t0 = pd.read_csv(GROWTHWHEAT_ROOTS_INPUTS_FILEPATH)
        growthwheat_axes_inputs_t0 = pd.read_csv(GROWTHWHEAT_AXES_INPUTS_FILEPATH)
        growthwheat_facade_ = growthwheat_facade.GrowthWheatFacade(g,
                                                                   growthwheat_ts * hour_to_second_conversion_factor,
                                                                   growthwheat_hiddenzones_inputs_t0,
                                                                   growthwheat_organ_inputs_t0,
                                                                   growthwheat_root_inputs_t0,
                                                                   growthwheat_axes_inputs_t0,
                                                                   shared_organs_inputs_outputs_df,
                                                                   shared_hiddenzones_inputs_outputs_df,
                                                                   shared_elements_inputs_outputs_df,
                                                                   shared_axes_inputs_outputs_df)

        # farquharwheat
        farquharwheat_elements_inputs_t0 = pd.read_csv(FARQUHARWHEAT_INPUTS_FILEPATH)
        farquharwheat_axes_inputs_t0 = pd.read_csv(FARQUHARWHEAT_AXES_INPUTS_FILEPATH)
        # Use the initial version of the photosynthesis sub-model (as in Barillot et al. 2016, and in Gauthier et al. 2020)
        update_parameters_farquharwheat = {'SurfacicProteins': False, 'NSC_Retroinhibition': False}

        farquharwheat_facade_ = farquharwheat_facade.FarquharWheatFacade(g,
                                                                         farquharwheat_elements_inputs_t0,
                                                                         farquharwheat_axes_inputs_t0,
                                                                         shared_elements_inputs_outputs_df,
                                                                         update_parameters=update_parameters_farquharwheat)

        # elongwheat # Only for temperature related computations
        elongwheat_hiddenzones_inputs_t0 = pd.read_csv(ELONGWHEAT_HZ_INPUTS_FILEPATH)
        elongwheat_elements_inputs_t0 = pd.read_csv(ELONGWHEAT_ELEMENTS_INPUTS_FILEPATH)
        elongwheat_axes_inputs_t0 = pd.read_csv(ELONGWHEAT_AXES_INPUTS_FILEPATH)

        elongwheat_facade_ = elongwheat_facade.ElongWheatFacade(g,
                                                                elongwheat_ts * hour_to_second_conversion_factor,
                                                                elongwheat_axes_inputs_t0,
                                                                elongwheat_hiddenzones_inputs_t0,
                                                                elongwheat_elements_inputs_t0,
                                                                shared_axes_inputs_outputs_df,
                                                                shared_hiddenzones_inputs_outputs_df,
                                                                shared_elements_inputs_outputs_df,
                                                                adel_wheat, option_static=True)
        # cnwheat
        cnwheat_organs_inputs_t0 = pd.read_csv(CNWHEAT_ORGANS_INPUTS_FILEPATH)
        cnwheat_hiddenzones_inputs_t0 = pd.read_csv(CNWHEAT_HIDDENZONE_INPUTS_FILEPATH)
        cnwheat_elements_inputs_t0 = pd.read_csv(CNWHEAT_ELEMENTS_INPUTS_FILEPATH)
        cnwheat_soils_inputs_t0 = pd.read_csv(CNWHEAT_SOILS_INPUTS_FILEPATH)
        update_cnwheat_parameters = {'roots': {'K_AMINO_ACIDS_EXPORT': 25*3E-5,
                                               'K_NITRATE_EXPORT': 25*1E-6}}

        cnwheat_facade_ = cnwheat_facade.CNWheatFacade(g,
                                                       cnwheat_ts * hour_to_second_conversion_factor,
                                                       CULM_DENSITY,
                                                       update_cnwheat_parameters,
                                                       cnwheat_organs_inputs_t0,
                                                       cnwheat_hiddenzones_inputs_t0,
                                                       cnwheat_elements_inputs_t0,
                                                       cnwheat_soils_inputs_t0,
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
        stop_time = stop_time

        #  define lists of dataframes to store the inputs and the outputs of the models at each step.
        axes_all_data_list = []
        organs_all_data_list = []  # organs which belong to axes: roots, phloem, grains
        elements_all_data_list = []
        soils_all_data_list = []

        all_simulation_steps = []  # to store the steps of the simulation

        # run the simulators
        current_time_of_the_system = time.time()

        try:

            for t_elongwheat in range(start_time, stop_time, elongwheat_ts):  # Only to compute temperature related variable

                # run ElongWheat
                Tair, Tsoil = meteo.loc[t_elongwheat, ['air_temperature', 'air_temperature']]
                elongwheat_facade_.run(Tair, Tsoil, option_static=True)

                for t_senescwheat in range(t_elongwheat, t_elongwheat + elongwheat_ts, senescwheat_ts):

                    # run SenescWheat
                    senescwheat_facade_.run(forced_max_protein_elements, postflowering_stages=True)

                    # Test for fully senesced shoot tissues  #TODO: Make the model to work even if the whole shoot is dead but the roots are alived
                    if sum(senescwheat_facade_._shared_elements_inputs_outputs_df['green_area']) <= 0.25E-6:
                        break

                    for t_growthwheat in range(t_senescwheat, t_senescwheat + senescwheat_ts, growthwheat_ts):
                        # run GrowthWheat
                        growthwheat_facade_.run(postflowering_stages=True)

                        for t_farquharwheat in range(t_growthwheat, t_growthwheat + growthwheat_ts, farquharwheat_ts):
                            # get the meteo of the current step
                            Tair, ambient_CO2, RH, Ur, PARi = meteo.loc[t_farquharwheat, ['air_temperature', 'ambient_CO2', 'humidity', 'Wind', 'PARi']]
                            # get PARa for current step
                            aggregated_PARa = calculate_PARa_from_df(g, Eabs_df, PARi, multiple_sources=False)
                            # caribu_facade_.run(energy=PARi,sun_sky_option='sky')
                            caribu_facade_.update_shared_MTG({'PARa': aggregated_PARa})
                            caribu_facade_.update_shared_dataframes({'PARa': aggregated_PARa})
                            # run FarquharWheat
                            farquharwheat_facade_.run(Tair, ambient_CO2, RH, Ur)

                            for t_cnwheat in range(t_farquharwheat, t_farquharwheat + senescwheat_ts, cnwheat_ts):
                                Tair, Tsoil = meteo.loc[t_cnwheat, ['air_temperature', 'air_temperature']]
                                # run CNWheat
                                print('t cnwheat is {}'.format(t_cnwheat))
                                cnwheat_facade_.run(Tair=Tair, Tsoil=Tsoil)

                                # append the inputs and outputs at current step to global lists
                                all_simulation_steps.append(t_cnwheat)
                                axes_all_data_list.append(shared_axes_inputs_outputs_df.copy())
                                organs_all_data_list.append(shared_organs_inputs_outputs_df.copy())
                                elements_all_data_list.append(shared_elements_inputs_outputs_df.copy())
                                soils_all_data_list.append(shared_soils_inputs_outputs_df.copy())

                else:
                    continue
                break

        except Exception as ex:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            template = "An exception of type {0} occurred. Arguments:\n{1!r}"
            message = template.format(type(ex).__name__, ex.args)
            print(message, fname, exc_tb.tb_lineno)

        finally:
            execution_time = int(time.time() - current_time_of_the_system)
            print('\n', 'Simulation run in ', str(datetime.timedelta(seconds=execution_time)))

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

    # -POST-PROCESSING##

    if make_graphs:
        generate_graphs()


def generate_graphs():

    # POST PROCESSINGS

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
    postprocessing_df_dict = {}
    (postprocessing_df_dict[axes_postprocessing_file_basename],
     _,
     postprocessing_df_dict[organs_postprocessing_file_basename],
     postprocessing_df_dict[elements_postprocessing_file_basename],
     postprocessing_df_dict[soils_postprocessing_file_basename]) \
        = cnwheat_facade.CNWheatFacade.postprocessing(axes_outputs_df=states_df_dict[os.path.basename(AXES_STATES_FILEPATH).split('.')[0]],
                                                      organs_outputs_df=states_df_dict[os.path.basename(ORGANS_STATES_FILEPATH).split('.')[0]],
                                                      hiddenzone_outputs_df=None,
                                                      elements_outputs_df=states_df_dict[os.path.basename(ELEMENTS_STATES_FILEPATH).split('.')[0]],
                                                      soils_outputs_df=states_df_dict[os.path.basename(SOILS_STATES_FILEPATH).split('.')[0]],
                                                      delta_t=delta_t)

    # save the postprocessing to disk
    for postprocessing_file_basename, postprocessing_filepath in ((axes_postprocessing_file_basename, AXES_POSTPROCESSING_FILEPATH),
                                                                  (organs_postprocessing_file_basename, ORGANS_POSTPROCESSING_FILEPATH),
                                                                  (elements_postprocessing_file_basename, ELEMENTS_POSTPROCESSING_FILEPATH),
                                                                  (soils_postprocessing_file_basename, SOILS_POSTPROCESSING_FILEPATH)):
        postprocessing_df_dict[postprocessing_file_basename].to_csv(postprocessing_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION + 5))

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
    cnwheat_facade.CNWheatFacade.graphs(axes_postprocessing_df=postprocessing_df_dict[axes_postprocessing_file_basename],
                                        hiddenzones_postprocessing_df=None,
                                        organs_postprocessing_df=postprocessing_df_dict[organs_postprocessing_file_basename],
                                        elements_postprocessing_df=postprocessing_df_dict[elements_postprocessing_file_basename],
                                        soils_postprocessing_df=postprocessing_df_dict[soils_postprocessing_file_basename],
                                        graphs_dirpath=GRAPHS_DIRPATH)
    #
    # x_name = 't'
    # x_label='Time (Hour)'
    #
    # # 1) Photosynthetic organs
    # ph_elements_output_df = pd.read_csv(ELEMENTS_STATES_FILEPATH)
    #
    # graph_variables_ph_elements = {'PARa': u'Absorbed PAR (µmol m$^{-2}$ s$^{-1}$)', 'Ag': u'Gross photosynthesis (µmol m$^{-2}$ s$^{-1}$)','An': u'Net photosynthesis (µmol m$^{-2}$ s$^{-1}$)', 'Tr':u'Organ surfacic transpiration rate (mmol H$_{2}$0 m$^{-2}$ s$^{-1}$)', 'Transpiration':u'Organ transpiration rate (mmol H$_{2}$0 s$^{-1}$)', 'Rd': u'Mitochondrial respiration rate of organ in light (µmol C h$^{-1}$)', 'Ts': u'Temperature surface (°C)', 'gs': u'Conductance stomatique (mol m$^{-2}$ s$^{-1}$)',
    #                    'Conc_TriosesP': u'[TriosesP] (µmol g$^{-1}$ mstruct)', 'Conc_Starch':u'[Starch] (µmol g$^{-1}$ mstruct)', 'Conc_Sucrose':u'[Sucrose] (µmol g$^{-1}$ mstruct)', 'Conc_Fructan':u'[Fructan] (µmol g$^{-1}$ mstruct)',
    #                    'Conc_Nitrates': u'[Nitrates] (µmol g$^{-1}$ mstruct)', 'Conc_Amino_Acids': u'[Amino_Acids] (µmol g$^{-1}$ mstruct)', 'Conc_Proteins': u'[Proteins] (g g$^{-1}$ mstruct)',
    #                    'Nitrates_import': u'Total nitrates imported (µmol h$^{-1}$)', 'Amino_Acids_import': u'Total amino acids imported (µmol N h$^{-1}$)',
    #                    'S_Amino_Acids': u'[Rate of amino acids synthesis] (µmol N g$^{-1}$ mstruct h$^{-1}$)', 'S_Proteins': u'Rate of protein synthesis (µmol N g$^{-1}$ mstruct h$^{-1}$)', 'D_Proteins': u'Rate of protein degradation (µmol N g$^{-1}$ mstruct h$^{-1}$)', 'k_proteins': u'Relative rate of protein degradation (s$^{-1}$)',
    #                    'Loading_Sucrose': u'Loading Sucrose (µmol C sucrose h$^{-1}$)', 'Loading_Amino_Acids': u'Loading Amino acids (µmol N amino acids h$^{-1}$)',
    #                    'green_area': u'Green area (m$^{2}$)', 'R_phloem_loading': u'Respiration phloem loading (µmol C h$^{-1}$)', 'R_Nnit_red': u'Respiration nitrate reduction (µmol C h$^{-1}$)', 'R_residual': u'Respiration residual (µmol C h$^{-1}$)', 'R_maintenance': u'Respiration residual (µmol C h$^{-1}$)',
    #                    'mstruct': u'Structural mass (g)', 'Nstruct': u'Structural N mass (g)',
    #                    'Conc_cytokinins':u'[cytokinins] (UA g$^{-1}$ mstruct)', 'D_cytokinins':u'Cytokinin degradation (UA g$^{-1}$ mstruct)', 'cytokinins_import':u'Cytokinin import (UA)'}
    #
    #
    # for org_ph in (['blade'], ['sheath'], ['internode'], ['peduncle', 'ear']):
    #     for variable_name, variable_label in graph_variables_ph_elements.iteritems():
    #         graph_name = variable_name + '_' + '_'.join(org_ph) + '.PNG'
    #         cnwheat_tools.plot_cnwheat_ouputs(ph_elements_output_df,
    #                       x_name = x_name,
    #                       y_name = variable_name,
    #                       x_label=x_label,
    #                       y_label=variable_label,
    #                       filters={'organ': org_ph},
    #                       plot_filepath=os.path.join(GRAPHS_DIRPATH, graph_name),
    #                       explicit_label=False)
    #
    # # 2) Roots, grains and phloem
    # organs_output_df = pd.read_csv(ORGANS_STATES_FILEPATH)
    #
    # graph_variables_organs = {'Conc_Sucrose':u'[Sucrose] (µmol g$^{-1}$ mstruct)', 'Dry_Mass':'Dry mass (g)',
    #                     'Conc_Nitrates': u'[Nitrates] (µmol g$^{-1}$ mstruct)', 'Conc_Amino_Acids':u'[Amino Acids] (µmol g$^{-1}$ mstruct)', 'Proteins_N_Mass': u'[N Proteins] (g)',
    #                     'Uptake_Nitrates':u'Nitrates uptake (µmol h$^{-1}$)', 'Unloading_Sucrose':u'Unloaded sucrose (µmol C g$^{-1}$ mstruct h$^{-1}$)', 'Unloading_Amino_Acids':u'Unloaded Amino Acids (µmol N AA g$^{-1}$ mstruct h$^{-1}$)',
    #                     'S_Amino_Acids': u'Rate of amino acids synthesis (µmol N g$^{-1}$ mstruct h$^{-1}$)', 'S_Proteins': u'Rate of protein synthesis (µmol N h$^{-1}$)', 'Export_Nitrates': u'Total export of nitrates (µmol N h$^{-1}$)', 'Export_Amino_Acids': u'Total export of Amino acids (µmol N h$^{-1}$)',
    #                     'R_Nnit_upt': u'Respiration nitrates uptake (µmol C h$^{-1}$)', 'R_Nnit_red': u'Respiration nitrate reduction (µmol C h$^{-1}$)', 'R_residual': u'Respiration residual (µmol C h$^{-1}$)', 'R_maintenance': u'Respiration residual (µmol C h$^{-1}$)',
    #                     'R_grain_growth_struct': u'Respiration grain structural growth (µmol C h$^{-1}$)', 'R_grain_growth_starch': u'Respiration grain starch growth (µmol C h$^{-1}$)',
    #                     'R_growth': u'Growth respiration of roots (µmol C h$^{-1}$)', 'mstruct': u'Structural mass (g)', 'rate_mstruct_death': u'Rate of structural mass death (g)',
    #                     'C_exudation': u'Carbon lost by root exudation (µmol C g$^{-1}$ mstruct h$^{-1}$', 'N_exudation': u'Nitrogen lost by root exudation (µmol N g$^{-1}$ mstruct h$^{-1}$',
    #                     'Conc_cytokinins':u'[cytokinins] (UA g$^{-1}$ mstruct)', 'S_cytokinins':u'Rate of cytokinins synthesis (UA g$^{-1}$ mstruct)', 'Export_cytokinins': 'Export of cytokinins from roots (UA h$^{-1}$)',
    #                     'HATS_LATS': u'Potential uptake (µmol h$^{-1}$)' , 'regul_transpiration':'Regulating transpiration function'}
    #
    # for org in (['roots'], ['grains'], ['phloem']):
    #     for variable_name, variable_label in graph_variables_organs.iteritems():
    #         graph_name = variable_name + '_' + '_'.join(org) + '.PNG'
    #         cnwheat_tools.plot_cnwheat_ouputs(organs_output_df,
    #                       x_name = x_name,
    #                       y_name = variable_name,
    #                       x_label=x_label,
    #                       y_label=variable_label,
    #                       filters={'organ': org},
    #                       plot_filepath=os.path.join(GRAPHS_DIRPATH, graph_name),
    #                       explicit_label=False)
    #
    # # 3) Soil
    # soil_output_df = pd.read_csv(SOILS_STATES_FILEPATH)
    #
    # fig, (ax1) = plt.subplots(1)
    # conc_nitrates_soil = soil_output_df['Conc_Nitrates_Soil']*14E-6
    # ax1.plot(soil_output_df['t'], conc_nitrates_soil)
    # ax1.set_ylabel(u'[Nitrates] (g m$^{-3}$)')
    # ax1.set_xlabel('Time from flowering (hour)')
    # ax1.set_title = 'Conc Nitrates Soil'
    # plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Conc_Nitrates_Soil.PNG'), format='PNG', bbox_inches='tight')
    # plt.close()


if __name__ == '__main__':
    main(1200, run_simu=True, make_graphs=True)
