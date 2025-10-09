# -*- coding: latin-1 -*-

import datetime
import os
import random
import time
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from openalea.adel.adel_dynamic import AdelDyn
from openalea.adel.echap_leaf import echap_leaves
from openalea.fspmwheat import caribu_facade
from openalea.fspmwheat import cnwheat_facade
from openalea.fspmwheat import elongwheat_facade
from openalea.fspmwheat import farquharwheat_facade
from openalea.fspmwheat import fspmwheat_facade
from openalea.fspmwheat import growthwheat_facade
from openalea.fspmwheat import senescwheat_facade

"""
    main
    ~~~~

    A coupling of CN-Wheat, Farquhar-Wheat, Senesc-Wheat, Elong-Wheat, Growth-Wheat, Adel-Wheat and Caribu.
    This script was used to simulate a field experiment performed in 1998/99 in Grignon (France). 
    Results were published in Gauthier et al. 2020 (https://doi.org/10.1093/jxb/eraa276)
    
    :copyright: Copyright 2014-2016 INRA-ECOSYS, see AUTHORS.
    :license: see LICENSE for details.

"""

random.seed(1234)
np.random.seed(1234)

AXES_INDEX_COLUMNS = ['t', 'plant', 'axis']
ELEMENTS_INDEX_COLUMNS = ['t', 'plant', 'axis', 'metamer', 'organ', 'element']
HIDDENZONES_INDEX_COLUMNS = ['t', 'plant', 'axis', 'metamer']
ORGANS_INDEX_COLUMNS = ['t', 'plant', 'axis', 'organ']
SOILS_INDEX_COLUMNS = ['t', 'plant', 'axis']


def save_df_to_csv(df, outputs_filepath, precision):
    try:
        df.to_csv(outputs_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(precision))
    except IOError as err:
        path, filename = os.path.split(outputs_filepath)
        filename = os.path.splitext(filename)[0]
        newfilename = 'ACTUAL_{}.csv'.format(filename)
        newpath = os.path.join(path, newfilename)
        df.to_csv(newpath, na_rep='NA', index=False, float_format='%.{}f'.format(precision))
        warnings.warn('[{}] {}'.format(err.errno, err.strerror))
        warnings.warn('File will be saved at {}'.format(newpath))


def main(simulation_length, forced_start_time=0, run_simu=True, run_postprocessing=True, generate_graphs=True, run_from_outputs=False, stored_times=None,
         option_static=False, show_3Dplant=True, tillers_replications=None, heterogeneous_canopy=True,
         N_fertilizations=None, PLANT_DENSITY=None, update_parameters_all_models=None,
         INPUTS_DIRPATH='inputs', METEO_FILENAME='meteo.csv',
         OUTPUTS_DIRPATH='outputs', POSTPROCESSING_DIRPATH='postprocessing', GRAPHS_DIRPATH='graphs'):
    """
    Run a simulation of fspmwheat with coupling to several models

    :param int simulation_length: length of the simulation (hours)
    :param int forced_start_time: desired start time (hour)
    :param bool run_simu: whether to run the simulation 
    :param bool run_postprocessing: whether to run the postprocessing
    :param bool generate_graphs: whether to run the generate graphs
    :param bool run_from_outputs: whether to start a simulation from a specific time and initial states as found in previous outputs
    :param str or list stored_times: Time steps when are stored the model outputs. Can be either 'all', a list or an empty list. Default to 'all'
    :param bool option_static: Whether the model should be run for a static plant architecture
    :param bool show_3Dplant: whether to plot the scene in pgl viewer
    :param dict [str, float] tillers_replications: a dictionary with tiller id as key, and weight of replication as value.
    :param bool heterogeneous_canopy: Whether to create a duplicated heterogeneous canopy from the initial mtg.
    :param dict [int, float] or [str, float] N_fertilizations: a dictionary for N fertilisation regime {date: N_input}, with date in hour and N_input in µmol N nitrates
                                               or {'constant_Conc_Nitrates': val} for constant nitrates concentrations
    :param dict [int, int] PLANT_DENSITY: a dict with plant density per plant id (temporary used to account for different cultivars if needed) ; plant m-2
    :param dict update_parameters_all_models: a dict to update model parameters
                                             {'cnwheat': {'organ1': {'param1': 'val1', 'param2': 'val2'},
                                                          'organ2': {'param1': 'val1', 'param2': 'val2'}
                                                         },
                                              'elongwheat': {'param1': 'val1', 'param2': 'val2'}
                                             } 
    :param str or dict INPUTS_DIRPATH: the path directory of inputs, can also be {'adel':str, 'plants':str, 'meteo':str, 'soils':str}
                                                                    #  The directory at path 'adel' must contain files 'adel_pars.RData', 'adel0000.pckl' and 'scene0000.bgeom' for ADELWHEAT
    :param str METEO_FILENAME: the name of the file with meteo data
    :param str OUTPUTS_DIRPATH: the path to save outputs
    :param str POSTPROCESSING_DIRPATH: the path to save postprocessings
    :param str GRAPHS_DIRPATH: the path to save graphs
    
    """
    # ---------------------------------------------
    # ----- CONFIGURATION OF THE SIMULATION -------
    # ---------------------------------------------

    # -- SIMULATION PARAMETERS --

    # Length of the simulation (in hours)
    SIMULATION_LENGTH = simulation_length

    # define the time step in hours for each simulator
    CARIBU_TIMESTEP = 4
    SENESCWHEAT_TIMESTEP = 1
    FARQUHARWHEAT_TIMESTEP = 1
    ELONGWHEAT_TIMESTEP = 1
    GROWTHWHEAT_TIMESTEP = 1
    CNWHEAT_TIMESTEP = 1

    # Define default plant density (culm m-2)
    if PLANT_DENSITY is None:
        PLANT_DENSITY = {1: 250.}

    # precision of floats used to write and format the output CSV files
    OUTPUTS_PRECISION = 8

    # number of seconds in 1 hour
    HOUR_TO_SECOND_CONVERSION_FACTOR = 3600

    # Name of the CSV files which will contain the outputs of the model
    AXES_OUTPUTS_FILENAME = 'axes_outputs.csv'
    ORGANS_OUTPUTS_FILENAME = 'organs_outputs.csv'
    HIDDENZONES_OUTPUTS_FILENAME = 'hiddenzones_outputs.csv'
    ELEMENTS_OUTPUTS_FILENAME = 'elements_outputs.csv'
    SOILS_OUTPUTS_FILENAME = 'soils_outputs.csv'

    # -- INPUTS CONFIGURATION --

    # Path of the directory which contains the inputs of the model
    INPUTS_DIRPATH = INPUTS_DIRPATH

    # Name of the CSV files which describes the initial state of the system
    AXES_INITIAL_STATE_FILENAME = 'axes_initial_state.csv'
    ORGANS_INITIAL_STATE_FILENAME = 'organs_initial_state.csv'
    HIDDENZONES_INITIAL_STATE_FILENAME = 'hiddenzones_initial_state.csv'
    ELEMENTS_INITIAL_STATE_FILENAME = 'elements_initial_state.csv'
    SOILS_INITIAL_STATE_FILENAME = 'soils_initial_state.csv'
    # Read the inputs from CSV files and create inputs dataframes
    inputs_dataframes = {}
    if run_from_outputs:

        previous_outputs_dataframes = {}

        for initial_state_filename, outputs_filename, index_columns in ((AXES_INITIAL_STATE_FILENAME, AXES_OUTPUTS_FILENAME, AXES_INDEX_COLUMNS),
                                                                        (ORGANS_INITIAL_STATE_FILENAME, ORGANS_OUTPUTS_FILENAME, ORGANS_INDEX_COLUMNS),
                                                                        (HIDDENZONES_INITIAL_STATE_FILENAME, HIDDENZONES_OUTPUTS_FILENAME, HIDDENZONES_INDEX_COLUMNS),
                                                                        (ELEMENTS_INITIAL_STATE_FILENAME, ELEMENTS_OUTPUTS_FILENAME, ELEMENTS_INDEX_COLUMNS),
                                                                        (SOILS_INITIAL_STATE_FILENAME, SOILS_OUTPUTS_FILENAME, SOILS_INDEX_COLUMNS)):

            previous_outputs_dataframe = pd.read_csv(os.path.join(OUTPUTS_DIRPATH, outputs_filename))
            # Convert NaN to None
            previous_outputs_dataframes[outputs_filename] = previous_outputs_dataframe.where(previous_outputs_dataframe.notnull(), None)

            assert 't' in previous_outputs_dataframes[outputs_filename].columns
            if forced_start_time > 0:
                new_start_time = forced_start_time + 1
                previous_outputs_dataframes[outputs_filename] = previous_outputs_dataframes[outputs_filename][previous_outputs_dataframes[outputs_filename]['t'] <= forced_start_time]
            else:
                last_t_step = max(previous_outputs_dataframes[outputs_filename]['t'])
                new_start_time = last_t_step + 1

            if initial_state_filename == ELEMENTS_INITIAL_STATE_FILENAME:
                elements_previous_outputs = previous_outputs_dataframes[outputs_filename]
                new_initial_state = elements_previous_outputs[~elements_previous_outputs.is_over.isnull()]
            else:
                new_initial_state = previous_outputs_dataframes[outputs_filename]
            idx = new_initial_state.groupby([col for col in index_columns if col != 't'])['t'].transform('max') == new_initial_state['t']
            inputs_dataframes[initial_state_filename] = new_initial_state[idx].drop(['t'], axis=1)

        # Make sure boolean columns have either type bool or float
        bool_columns = ['is_over', 'is_growing', 'leaf_is_emerged', 'internode_is_visible', 'leaf_is_growing', 'internode_is_growing', 'leaf_is_remobilizing', 'internode_is_remobilizing']
        for df in [inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME], inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME]]:
            for cln in bool_columns:
                if cln in df.keys():
                    df.loc[:, cln] = df.loc[:, cln].replace({'False': 0.0, 'True': 1.0})
                    df.loc[:, cln] = pd.to_numeric(df.loc[:, cln])
    else:
        new_start_time = -1
        for inputs_filename in (AXES_INITIAL_STATE_FILENAME,
                                ORGANS_INITIAL_STATE_FILENAME,
                                HIDDENZONES_INITIAL_STATE_FILENAME,
                                ELEMENTS_INITIAL_STATE_FILENAME,
                                SOILS_INITIAL_STATE_FILENAME):
            inputs_dataframe = pd.read_csv(os.path.join(INPUTS_DIRPATH, inputs_filename))
            inputs_dataframes[inputs_filename] = inputs_dataframe.where(inputs_dataframe.notnull(), None)

    # Start time of the simulation
    START_TIME = max(0, new_start_time)

    # Name of the CSV files which contains the meteo data
    meteo = pd.read_csv(os.path.join(INPUTS_DIRPATH, METEO_FILENAME), index_col='t')

    # -- OUTPUTS CONFIGURATION --

    # Save the outputs with a full scan of the MTG at each time step (or at selected time steps)
    UPDATE_SHARED_DF = False
    if stored_times is None:
        stored_times = 'all'
    if not (stored_times == 'all' or type(stored_times) == list):
        print('stored_times should be either \'all\', a list or an empty list.')
        raise

    # create empty dataframes to shared data between the models
    shared_axes_inputs_outputs_df = pd.DataFrame()
    shared_organs_inputs_outputs_df = pd.DataFrame()
    shared_hiddenzones_inputs_outputs_df = pd.DataFrame()
    shared_elements_inputs_outputs_df = pd.DataFrame()
    shared_soils_inputs_outputs_df = pd.DataFrame()

    # define lists of dataframes to store the inputs and the outputs of the models at each step.
    axes_all_data_list = []
    organs_all_data_list = []  # organs which belong to axes: roots, phloem, grains
    hiddenzones_all_data_list = []
    elements_all_data_list = []
    soils_all_data_list = []

    all_simulation_steps = []  # to store the steps of the simulation

    # -- POSTPROCESSING CONFIGURATION --

    # Name of the CSV files which will contain the postprocessing of the model
    AXES_POSTPROCESSING_FILENAME = 'axes_postprocessing.csv'
    ORGANS_POSTPROCESSING_FILENAME = 'organs_postprocessing.csv'
    HIDDENZONES_POSTPROCESSING_FILENAME = 'hiddenzones_postprocessing.csv'
    ELEMENTS_POSTPROCESSING_FILENAME = 'elements_postprocessing.csv'
    SOILS_POSTPROCESSING_FILENAME = 'soils_postprocessing.csv'

    # -- ADEL and MTG CONFIGURATION --

    # read adelwheat inputs at t0
    adel_wheat = AdelDyn(seed=1, scene_unit='m', leaves=echap_leaves(xy_model='Soissons_byleafclass'))
    g = adel_wheat.load(directory=INPUTS_DIRPATH)

    # ---------------------------------------------
    # ----- CONFIGURATION OF THE FACADES -------
    # ---------------------------------------------

    # -- ELONGWHEAT (created first because it is the only facade to add new metamers) --
    # Initial states
    elongwheat_hiddenzones_initial_state = inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME][
        elongwheat_facade.converter.HIDDENZONE_TOPOLOGY_COLUMNS + [i for i in elongwheat_facade.simulation.HIDDENZONE_INPUTS if i in
                                                                   inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME].columns]].copy()
    elongwheat_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
        elongwheat_facade.converter.ELEMENT_TOPOLOGY_COLUMNS + [i for i in elongwheat_facade.simulation.ELEMENT_INPUTS if i in
                                                                inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns]].copy()
    elongwheat_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
        elongwheat_facade.converter.AXIS_TOPOLOGY_COLUMNS + [i for i in elongwheat_facade.simulation.AXIS_INPUTS if i in inputs_dataframes[AXES_INITIAL_STATE_FILENAME].columns]].copy()

    phytoT_df = pd.read_csv(os.path.join(INPUTS_DIRPATH, 'phytoT.csv'))

    # Update parameters if specified
    if update_parameters_all_models and 'elongwheat' in update_parameters_all_models:
        update_parameters_elongwheat = update_parameters_all_models['elongwheat']
    else:
        update_parameters_elongwheat = None

    # Facade initialisation
    elongwheat_facade_ = elongwheat_facade.ElongWheatFacade(g,
                                                            ELONGWHEAT_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                            elongwheat_axes_initial_state,
                                                            elongwheat_hiddenzones_initial_state,
                                                            elongwheat_elements_initial_state,
                                                            shared_axes_inputs_outputs_df,
                                                            shared_hiddenzones_inputs_outputs_df,
                                                            shared_elements_inputs_outputs_df,
                                                            adel_wheat, phytoT_df,
                                                            update_parameters_elongwheat,
                                                            update_shared_df=UPDATE_SHARED_DF)

    # -- CARIBU --
    caribu_facade_ = caribu_facade.CaribuFacade(g,
                                                shared_elements_inputs_outputs_df,
                                                adel_wheat,
                                                update_shared_df=UPDATE_SHARED_DF)

    # -- SENESCWHEAT --
    # Initial states
    senescwheat_roots_initial_state = inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].loc[inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME]['organ'] == 'roots'][
        senescwheat_facade.converter.ROOTS_TOPOLOGY_COLUMNS +
        [i for i in senescwheat_facade.converter.SENESCWHEAT_ROOTS_INPUTS if i in inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].columns]].copy()

    senescwheat_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
        senescwheat_facade.converter.ELEMENTS_TOPOLOGY_COLUMNS +
        [i for i in senescwheat_facade.converter.SENESCWHEAT_ELEMENTS_INPUTS if i in inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns]].copy()

    senescwheat_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
        senescwheat_facade.converter.AXES_TOPOLOGY_COLUMNS +
        [i for i in senescwheat_facade.converter.SENESCWHEAT_AXES_INPUTS if i in inputs_dataframes[AXES_INITIAL_STATE_FILENAME].columns]].copy()

    # Update parameters if specified
    if update_parameters_all_models and 'senescwheat' in update_parameters_all_models:
        update_parameters_senescwheat = update_parameters_all_models['senescwheat']
    else:
        update_parameters_senescwheat = None

    # Facade initialisation
    senescwheat_facade_ = senescwheat_facade.SenescWheatFacade(g,
                                                               SENESCWHEAT_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                               senescwheat_roots_initial_state,
                                                               senescwheat_axes_initial_state,
                                                               senescwheat_elements_initial_state,
                                                               shared_organs_inputs_outputs_df,
                                                               shared_axes_inputs_outputs_df,
                                                               shared_elements_inputs_outputs_df,
                                                               update_parameters_senescwheat,
                                                               update_shared_df=UPDATE_SHARED_DF)

    # -- FARQUHARWHEAT --
    # Initial states
    farquharwheat_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
        farquharwheat_facade.converter.ELEMENT_TOPOLOGY_COLUMNS +
        [i for i in farquharwheat_facade.converter.FARQUHARWHEAT_ELEMENTS_INPUTS if i in inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns]].copy()

    farquharwheat_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
        farquharwheat_facade.converter.AXIS_TOPOLOGY_COLUMNS +
        [i for i in farquharwheat_facade.converter.FARQUHARWHEAT_AXES_INPUTS if i in inputs_dataframes[AXES_INITIAL_STATE_FILENAME].columns]].copy()

    # Use the initial version of the photosynthesis sub-model (as in Barillot et al. 2016, and in Gauthier et al. 2020)
    update_parameters_farquharwheat = {'SurfacicProteins': False, 'NSC_Retroinhibition': False}

    # Facade initialisation
    farquharwheat_facade_ = farquharwheat_facade.FarquharWheatFacade(g,
                                                                     farquharwheat_elements_initial_state,
                                                                     farquharwheat_axes_initial_state,
                                                                     shared_elements_inputs_outputs_df,
                                                                     update_parameters_farquharwheat,
                                                                     update_shared_df=UPDATE_SHARED_DF)

    # -- GROWTHWHEAT --
    # Initial states
    growthwheat_hiddenzones_initial_state = inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME][
        growthwheat_facade.converter.HIDDENZONE_TOPOLOGY_COLUMNS +
        [i for i in growthwheat_facade.simulation.HIDDENZONE_INPUTS if i in inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME].columns]].copy()

    growthwheat_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
        growthwheat_facade.converter.ELEMENT_TOPOLOGY_COLUMNS +
        [i for i in growthwheat_facade.simulation.ELEMENT_INPUTS if i in inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns]].copy()

    growthwheat_root_initial_state = inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].loc[inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME]['organ'] == 'roots'][
        growthwheat_facade.converter.ROOT_TOPOLOGY_COLUMNS +
        [i for i in growthwheat_facade.simulation.ROOT_INPUTS if i in inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].columns]].copy()

    growthwheat_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
        growthwheat_facade.converter.AXIS_TOPOLOGY_COLUMNS +
        [i for i in growthwheat_facade.simulation.AXIS_INPUTS if i in inputs_dataframes[AXES_INITIAL_STATE_FILENAME].columns]].copy()

    # Update parameters if specified
    if update_parameters_all_models and 'growthwheat' in update_parameters_all_models:
        update_parameters_growthwheat = update_parameters_all_models['growthwheat']
    else:
        update_parameters_growthwheat = None

    # Facade initialisation
    growthwheat_facade_ = growthwheat_facade.GrowthWheatFacade(g,
                                                               GROWTHWHEAT_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                               growthwheat_hiddenzones_initial_state,
                                                               growthwheat_elements_initial_state,
                                                               growthwheat_root_initial_state,
                                                               growthwheat_axes_initial_state,
                                                               shared_organs_inputs_outputs_df,
                                                               shared_hiddenzones_inputs_outputs_df,
                                                               shared_elements_inputs_outputs_df,
                                                               shared_axes_inputs_outputs_df,
                                                               update_parameters_growthwheat,
                                                               update_shared_df=UPDATE_SHARED_DF)

    # -- CNWHEAT --
    # Initial states
    cnwheat_organs_initial_state = inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME][
        [i for i in cnwheat_facade.cnwheat_converter.ORGANS_VARIABLES if i in inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].columns]].copy()

    cnwheat_hiddenzones_initial_state = inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME][
        [i for i in cnwheat_facade.cnwheat_converter.HIDDENZONE_VARIABLES if i in inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME].columns]].copy()

    cnwheat_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
        [i for i in cnwheat_facade.cnwheat_converter.ELEMENTS_VARIABLES if i in inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns]].copy()

    cnwheat_soils_initial_state = inputs_dataframes[SOILS_INITIAL_STATE_FILENAME][
        [i for i in cnwheat_facade.cnwheat_converter.SOILS_VARIABLES if i in inputs_dataframes[SOILS_INITIAL_STATE_FILENAME].columns]].copy()

    # Update parameters if specified
    if update_parameters_all_models and 'cnwheat' in update_parameters_all_models:
        update_parameters_cnwheat = update_parameters_all_models['cnwheat']
    else:
        update_parameters_cnwheat = {}

    # Facade initialisation
    cnwheat_facade_ = cnwheat_facade.CNWheatFacade(g,
                                                   CNWHEAT_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                   PLANT_DENSITY,
                                                   update_parameters_cnwheat,
                                                   cnwheat_organs_initial_state,
                                                   cnwheat_hiddenzones_initial_state,
                                                   cnwheat_elements_initial_state,
                                                   cnwheat_soils_initial_state,
                                                   shared_axes_inputs_outputs_df,
                                                   shared_organs_inputs_outputs_df,
                                                   shared_hiddenzones_inputs_outputs_df,
                                                   shared_elements_inputs_outputs_df,
                                                   shared_soils_inputs_outputs_df,
                                                   update_shared_df=UPDATE_SHARED_DF)

    # Run cnwheat with constant nitrates concentration in the soil if specified
    if N_fertilizations is not None and 'constant_Conc_Nitrates' in N_fertilizations.keys():
        cnwheat_facade_.soils[(1, 'MS')].constant_Conc_Nitrates = True  # TODO: make (1, 'MS') more general
        cnwheat_facade_.soils[(1, 'MS')].nitrates = N_fertilizations['constant_Conc_Nitrates'] * cnwheat_facade_.soils[(1, 'MS')].volume

    # -- FSPMWHEAT --
    # Facade initialisation
    fspmwheat_facade_ = fspmwheat_facade.FSPMWheatFacade(g)

    # Update geometry
    adel_wheat.update_geometry(g)

    # ---------------------------------------------
    # -----      RUN OF THE SIMULATION      -------
    # ---------------------------------------------

    if run_simu:

        try:
            current_time_of_the_system = time.time()
            for t_caribu in range(START_TIME, SIMULATION_LENGTH, SENESCWHEAT_TIMESTEP):
                # run Caribu
                PARi = meteo.loc[t_caribu, ['PARi']].iloc[0]
                DOY = meteo.loc[t_caribu, ['DOY']].iloc[0]
                hour = meteo.loc[t_caribu, ['hour']].iloc[0]
                PARi_next_hours = meteo.loc[range(t_caribu, t_caribu + CARIBU_TIMESTEP), ['PARi']].sum().values[0]

                if (t_caribu % CARIBU_TIMESTEP == 0) and (PARi_next_hours > 0):
                    run_caribu = True
                else:
                    run_caribu = False

                caribu_facade_.run(run_caribu, energy=PARi, DOY=DOY, hourTU=hour, latitude=48.85, sun_sky_option='sky', heterogeneous_canopy=heterogeneous_canopy, plant_density=PLANT_DENSITY[1])

                for t_senescwheat in range(t_caribu, t_caribu + SENESCWHEAT_TIMESTEP, SENESCWHEAT_TIMESTEP):
                    # run SenescWheat
                    senescwheat_facade_.run()

                    # Test for dead plant # TODO: adapt in case of multiple plants
                    if not shared_elements_inputs_outputs_df.empty and \
                            np.nansum(shared_elements_inputs_outputs_df.loc[shared_elements_inputs_outputs_df['element'].isin(['StemElement', 'LeafElement1']), 'green_area']) == 0:
                        # append the inputs and outputs at current step to global lists
                        all_simulation_steps.append(t_senescwheat)
                        axes_all_data_list.append(shared_axes_inputs_outputs_df.copy())
                        organs_all_data_list.append(shared_organs_inputs_outputs_df.copy())
                        hiddenzones_all_data_list.append(shared_hiddenzones_inputs_outputs_df.copy())
                        elements_all_data_list.append(shared_elements_inputs_outputs_df.copy())
                        soils_all_data_list.append(shared_soils_inputs_outputs_df.copy())
                        break

                    # Run the rest of the model if the plant is alive
                    for t_farquharwheat in range(t_senescwheat, t_senescwheat + SENESCWHEAT_TIMESTEP, FARQUHARWHEAT_TIMESTEP):
                        # get the meteo of the current step
                        Ta, ambient_CO2, RH, Ur = meteo.loc[t_farquharwheat, ['air_temperature', 'ambient_CO2', 'humidity', 'Wind']]

                        # run FarquharWheat
                        farquharwheat_facade_.run(Ta, ambient_CO2, RH, Ur)

                        for t_elongwheat in range(t_farquharwheat, t_farquharwheat + FARQUHARWHEAT_TIMESTEP, ELONGWHEAT_TIMESTEP):
                            # run ElongWheat
                            Tair, Tsoil = meteo.loc[t_elongwheat, ['air_temperature', 'soil_temperature']]
                            elongwheat_facade_.run(Tair, Tsoil, option_static=option_static)

                            # Update geometry
                            adel_wheat.update_geometry(g)

                            for t_growthwheat in range(t_elongwheat, t_elongwheat + ELONGWHEAT_TIMESTEP, GROWTHWHEAT_TIMESTEP):
                                # run GrowthWheat
                                growthwheat_facade_.run()

                                for t_cnwheat in range(t_growthwheat, t_growthwheat + GROWTHWHEAT_TIMESTEP, CNWHEAT_TIMESTEP):
                                    print('t cnwheat is {}'.format(t_cnwheat))

                                    # N fertilization if any
                                    if N_fertilizations is not None and len(N_fertilizations) > 0:
                                        if t_cnwheat in N_fertilizations.keys():
                                            cnwheat_facade_.soils[(1, 'MS')].nitrates += N_fertilizations[t_cnwheat]

                                    if t_cnwheat > 0:
                                        # run CNWheat
                                        Tair = meteo.loc[t_elongwheat, 'air_temperature']
                                        Tsoil = meteo.loc[t_elongwheat, 'soil_temperature']
                                        cnwheat_facade_.run(Tair, Tsoil, tillers_replications)

                                    # append outputs at current step to global lists
                                    if (stored_times == 'all') or (t_cnwheat in stored_times):
                                        axes_outputs, elements_outputs, hiddenzones_outputs, organs_outputs, soils_outputs = fspmwheat_facade_.build_outputs_df_from_MTG()

                                        all_simulation_steps.append(t_cnwheat)
                                        axes_all_data_list.append(axes_outputs)
                                        organs_all_data_list.append(organs_outputs)
                                        hiddenzones_all_data_list.append(hiddenzones_outputs)
                                        elements_all_data_list.append(elements_outputs)
                                        soils_all_data_list.append(soils_outputs)

                else:
                    # Continue if SenescWheat loop wasn't broken because of dead plant.
                    continue
                # SenescWheat loop was broken, break the Caribu loop.
                break

            execution_time = int(time.time() - current_time_of_the_system)
            print('\n' 'Simulation run in {}'.format(str(datetime.timedelta(seconds=execution_time))))
            print('Creating postprocessing and graphs...')

        finally:
            # convert list of outputs into dataframes
            outputs_df_dict = {}
            for outputs_df_list, outputs_filename, index_columns in ((axes_all_data_list, AXES_OUTPUTS_FILENAME, AXES_INDEX_COLUMNS),
                                                                     (organs_all_data_list, ORGANS_OUTPUTS_FILENAME, ORGANS_INDEX_COLUMNS),
                                                                     (hiddenzones_all_data_list, HIDDENZONES_OUTPUTS_FILENAME, HIDDENZONES_INDEX_COLUMNS),
                                                                     (elements_all_data_list, ELEMENTS_OUTPUTS_FILENAME, ELEMENTS_INDEX_COLUMNS),
                                                                     (soils_all_data_list, SOILS_OUTPUTS_FILENAME, SOILS_INDEX_COLUMNS)):
                outputs_filepath = os.path.join(OUTPUTS_DIRPATH, outputs_filename)
                outputs_df = pd.concat(outputs_df_list, keys=all_simulation_steps, sort=False)
                outputs_df.reset_index(0, inplace=True)
                outputs_df.rename({'level_0': 't'}, axis=1, inplace=True)
                outputs_df = outputs_df.reindex(index_columns + outputs_df.columns.difference(index_columns).tolist(), axis=1, copy=False)
                if run_from_outputs:
                    outputs_df = pd.concat([previous_outputs_dataframes[outputs_filename], outputs_df], sort=False)
                outputs_df.fillna(value=np.nan, inplace=True)  # Convert back None to NaN
                save_df_to_csv(outputs_df, outputs_filepath, OUTPUTS_PRECISION)
                outputs_file_basename = outputs_filename.split('.')[0]
                outputs_df_dict[outputs_file_basename] = outputs_df.reset_index()

    # ---------------------------------------------
    # -----      POST-PROCESSING      -------
    # ---------------------------------------------

    if run_postprocessing:
        # Retrieve outputs dataframes from precedent simulation run
        delta_t = CNWHEAT_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR

        # run the postprocessing
        axes_postprocessing_file_basename = AXES_POSTPROCESSING_FILENAME.split('.')[0]
        hiddenzones_postprocessing_file_basename = HIDDENZONES_POSTPROCESSING_FILENAME.split('.')[0]
        organs_postprocessing_file_basename = ORGANS_POSTPROCESSING_FILENAME.split('.')[0]
        elements_postprocessing_file_basename = ELEMENTS_POSTPROCESSING_FILENAME.split('.')[0]
        soils_postprocessing_file_basename = SOILS_POSTPROCESSING_FILENAME.split('.')[0]

        postprocessing_df_dict = {}
        (postprocessing_df_dict[axes_postprocessing_file_basename],
         postprocessing_df_dict[hiddenzones_postprocessing_file_basename],
         postprocessing_df_dict[organs_postprocessing_file_basename],
         postprocessing_df_dict[elements_postprocessing_file_basename],
         postprocessing_df_dict[soils_postprocessing_file_basename]) \
            = cnwheat_facade.CNWheatFacade.postprocessing(axes_outputs_df=outputs_df_dict[AXES_OUTPUTS_FILENAME.split('.')[0]],
                                                          hiddenzone_outputs_df=outputs_df_dict[HIDDENZONES_OUTPUTS_FILENAME.split('.')[0]],
                                                          organs_outputs_df=outputs_df_dict[ORGANS_OUTPUTS_FILENAME.split('.')[0]],
                                                          elements_outputs_df=outputs_df_dict[ELEMENTS_OUTPUTS_FILENAME.split('.')[0]],
                                                          soils_outputs_df=outputs_df_dict[SOILS_OUTPUTS_FILENAME.split('.')[0]],
                                                          delta_t=delta_t)

        for postprocessing_file_basename, postprocessing_filename in ((axes_postprocessing_file_basename, AXES_POSTPROCESSING_FILENAME),
                                                                      (hiddenzones_postprocessing_file_basename, HIDDENZONES_POSTPROCESSING_FILENAME),
                                                                      (organs_postprocessing_file_basename, ORGANS_POSTPROCESSING_FILENAME),
                                                                      (elements_postprocessing_file_basename, ELEMENTS_POSTPROCESSING_FILENAME),
                                                                      (soils_postprocessing_file_basename, SOILS_POSTPROCESSING_FILENAME)):
            postprocessing_filepath = os.path.join(POSTPROCESSING_DIRPATH, postprocessing_filename)
            postprocessing_df_dict[postprocessing_file_basename].to_csv(postprocessing_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))

    # ---------------------------------------------
    # -----            GRAPHS               -------
    # ---------------------------------------------

    if generate_graphs:
        # Retrieve last computed post-processing dataframes
        axes_postprocessing_file_basename = AXES_POSTPROCESSING_FILENAME.split('.')[0]
        organs_postprocessing_file_basename = ORGANS_POSTPROCESSING_FILENAME.split('.')[0]
        hiddenzones_postprocessing_file_basename = HIDDENZONES_POSTPROCESSING_FILENAME.split('.')[0]
        elements_postprocessing_file_basename = ELEMENTS_POSTPROCESSING_FILENAME.split('.')[0]
        soils_postprocessing_file_basename = SOILS_POSTPROCESSING_FILENAME.split('.')[0]

        # --- Generate graphs from postprocessing files
        plt.ioff()
        df_elt = postprocessing_df_dict[elements_postprocessing_file_basename]
        df_SAM = pd.read_csv(os.path.join(OUTPUTS_DIRPATH, AXES_OUTPUTS_FILENAME))

        cnwheat_facade.CNWheatFacade.graphs(axes_postprocessing_df=postprocessing_df_dict[axes_postprocessing_file_basename],
                                            hiddenzones_postprocessing_df=postprocessing_df_dict[hiddenzones_postprocessing_file_basename],
                                            organs_postprocessing_df=postprocessing_df_dict[organs_postprocessing_file_basename],
                                            elements_postprocessing_df=postprocessing_df_dict[elements_postprocessing_file_basename],
                                            soils_postprocessing_df=postprocessing_df_dict[soils_postprocessing_file_basename],
                                            graphs_dirpath=GRAPHS_DIRPATH)
    print('End of execution')


if __name__ == '__main__':
    main(2500, forced_start_time=0, run_simu=False, run_postprocessing=False, generate_graphs=False, run_from_outputs=False,
         show_3Dplant=False,
         option_static=False, tillers_replications={'T1': 0.5, 'T2': 0.5, 'T3': 0.5, 'T4': 0.5},
         heterogeneous_canopy=True, N_fertilizations={1440: 357143, 2520: 1000000},
         PLANT_DENSITY={1: 250}, METEO_FILENAME='meteo_Ljutovac2002.csv')
