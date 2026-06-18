# -*- coding: latin-1 -*-

import datetime
import os
import random
import time
import warnings
import ast

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from openalea.adel.adel_dynamic import AdelDyn
from openalea.adel.Stand import AgronomicStand
from openalea.adel.echap_leaf import echap_leaves

from openalea.fspmwheat import caribu_facade
from openalea.fspmwheat import cnwheat_facade
from openalea.fspmwheat import elongwheat_facade
from openalea.fspmwheat import farquharwheat_facade
from openalea.fspmwheat import fspmwheat_facade
from openalea.fspmwheat import growthwheat_facade
from openalea.fspmwheat import senescwheat_facade
from openalea.fspmwheat import turgorgrowth_facade

"""

    Executes the coupling between models CN-Wheat, Farquhar-Wheat, Senesc-Wheat, Elong-Wheat, Growth-Wheat, Adel-Wheat and Caribu.
    This uses the format MTG to exchange data between the models.

"""

random.seed(1234)
np.random.seed(1234)

AXES_INDEX_COLUMNS = ['t', 'plant', 'axis']
ELEMENTS_INDEX_COLUMNS = ['t', 'plant', 'axis', 'metamer', 'organ', 'element']
HIDDENZONES_INDEX_COLUMNS = ['t', 'plant', 'axis', 'metamer']
ORGANS_INDEX_COLUMNS = ['t', 'plant', 'axis', 'organ']
SOILS_INDEX_COLUMNS = ['t', 'plant', 'axis']


def get_management_value(df, variable_name):
    """
    Extracts management variables from the input file.

    :param pandas.DataFrame df: the dataframe with all management variables
    :param str variable_name: name of the current management variable

    return management_input
    :rtype float or dict
    """

    management_input = df.loc[variable_name].iloc[0]

    # In case of missing value in the management file
    if pd.isna(management_input):
        raise ValueError('No input found for {} in the management file.'.format(variable_name))

    management_input = ast.literal_eval(str(management_input))

    return management_input

def save_df_to_csv(df, outputs_filepath, precision):
    """
    Write outputs of the model

    :param pandas.DataFrame df: the current output dataframe
    :param str outputs_filepath: path of the output file
    :param int precision: number of digit

    """

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


def run(simulation_length, forced_start_time=0, run_simu=True, run_postprocessing=True, generate_graphs=True,
         run_from_outputs=False, stored_times=None, show_3Dplant=False,
         hydraulics=False, stomatal_model_name='BWB', drought_trigger=None, rehydration_scenario=None,
         optimal_growth_option=False, option_static=False,
         external_soil_model=False, tillers_replications=None, heterogeneous_canopy=True,
         update_parameters_all_models=None, step_callback=None,
         INPUTS_DIRPATH='inputs', METEO_FILENAME='meteo.csv', MANAGEMENT_FILENAME='management.csv',
         OUTPUTS_DIRPATH='outputs', POSTPROCESSING_DIRPATH='postprocessing', GRAPHS_DIRPATH='graphs'):
    """
    Run a simulation of fspmwheat with coupling to several models

    :param int simulation_length: length of the simulation (hours)
    :param int forced_start_time: desired start time (hour)
    :param bool run_simu: whether to run the simulation
    :param bool run_postprocessing: whether to run the postprocessing
    :param bool generate_graphs: whether to run generate graphs
    :param bool run_from_outputs: whether to start a simulation from a specific time and initial states as found in previous outputs
    :param str or list stored_times: Time steps when are stored the model outputs. Can be either 'all', a list or an empty list. Default to 'all'
    :param bool show_3Dplant: whether to plot the scene in pgl viewer
    :param bool hydraulics: if True the model will assume the coupling to the turgor-driven growth model
    :param str stomatal_model_name: the model of stomatal conductance. Should be one of 'BWB', 'Leuning', 'Tuzet' or 'hydraulics'. 'Tuzet' and 'hydraulics' requires hydraulics to be True
    :param dict or None drought_trigger: a dict for external drought control scenario.
                                            {'trigger_variable': value}. For now, only implemented for 'green_area' variable (value = green area above which the drought will be triggered).
    :param dict or None rehydration_scenario: a dict to specify the rehydration scenario.
                                             {'stop_drought_SRWC': SRWC at which the drought event stops (%),
                                             'SRWC_target': Target SRWC for rehydration (%),
                                             'rehydration_duration': duration of the rehydration period (days)}
    :param bool optimal_growth_option: if True the model will assume optimal growth conditions
    :param bool option_static: Whether the model should be run for a static plant architecture
    :param bool external_soil_model: whether an external soil model is coupled to cnwheat. If True, cnwheat will skip calculations made in soil and uptake N by roots
    :param dict [str, float] tillers_replications: a dictionary with tiller id as key, and weight of replication as value.
    :param bool heterogeneous_canopy: Whether to create a duplicated heterogeneous canopy from the initial mtg.
    :param dict update_parameters_all_models: a dict to update model parameters
                                             {'cnwheat': {'organ1': {'param1': 'val1', 'param2': 'val2'},
                                                          'organ2': {'param1': 'val1', 'param2': 'val2'}
                                                         },
                                              'elongwheat': {'param1': 'val1', 'param2': 'val2'}
                                             }
    :param dict or None step_callback: a dict of functions used to force some external inputs that are natively computed by the model
                                            {'function_name' : function , ...}
    :param str INPUTS_DIRPATH: the path directory of inputs
                                             #  The directory at path 'adel' must contain files 'adel_pars.RData', 'adel0000.pckl' and 'scene0000.bgeom' for ADELWHEAT
    :param str METEO_FILENAME: the name of the file with meteo data
    :param str MANAGEMENT_FILENAME: the name of the file with managment data
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
    ELONGWHEAT_TIMESTEP = 1
    GROWTHWHEAT_TIMESTEP = 1
    CNWHEAT_TIMESTEP = 1
    TURGORGROWTH_TIMESTEP = 1

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
        previous_axes_outputs_dataframe = pd.read_csv(os.path.join(OUTPUTS_DIRPATH, AXES_OUTPUTS_FILENAME))
        assert 't' in previous_axes_outputs_dataframe.columns
        if forced_start_time > 0:
            new_start_time = forced_start_time + 1
        else:
            last_t_step = int(previous_axes_outputs_dataframe['t'].max())
            new_start_time = last_t_step + 1

        previous_outputs_dataframes = {}

        for initial_state_filename, outputs_filename, index_columns in ((AXES_INITIAL_STATE_FILENAME, AXES_OUTPUTS_FILENAME, AXES_INDEX_COLUMNS),
                                                                        (ORGANS_INITIAL_STATE_FILENAME, ORGANS_OUTPUTS_FILENAME, ORGANS_INDEX_COLUMNS),
                                                                        (HIDDENZONES_INITIAL_STATE_FILENAME, HIDDENZONES_OUTPUTS_FILENAME, HIDDENZONES_INDEX_COLUMNS),
                                                                        (ELEMENTS_INITIAL_STATE_FILENAME, ELEMENTS_OUTPUTS_FILENAME, ELEMENTS_INDEX_COLUMNS),
                                                                        (SOILS_INITIAL_STATE_FILENAME, SOILS_OUTPUTS_FILENAME, SOILS_INDEX_COLUMNS)):

            previous_outputs_dataframe = pd.read_csv(os.path.join(OUTPUTS_DIRPATH, outputs_filename))
            # Convert NaN to None
            previous_outputs_dataframes[outputs_filename] = previous_outputs_dataframe.where(previous_outputs_dataframe.notnull(), None)

            # assert 't' in previous_outputs_dataframes[outputs_filename].columns
            if forced_start_time > 0:
                previous_outputs_dataframes[outputs_filename] = previous_outputs_dataframes[outputs_filename][previous_outputs_dataframes[outputs_filename]['t'] <= forced_start_time]

            if initial_state_filename == ELEMENTS_INITIAL_STATE_FILENAME:
                elements_previous_outputs = previous_outputs_dataframes[outputs_filename]
                new_initial_state = elements_previous_outputs[~elements_previous_outputs.is_over.isnull()]
            else:
                new_initial_state = previous_outputs_dataframes[outputs_filename]
            idx = new_initial_state.groupby([col for col in index_columns if col != 't'])['t'].transform(max) == new_initial_state['t']
            inputs_dataframes[initial_state_filename] = new_initial_state[idx].drop(['t'], axis=1)

        # Make sure boolean columns have either type bool or float
        bool_columns = ['is_over', 'is_growing', 'leaf_is_emerged', 'internode_is_visible', 'leaf_is_growing', 'internode_is_growing', 'leaf_is_remobilizing', 'internode_is_remobilizing']
        for df in [inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME], inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME]]:
            for cln in bool_columns:
                if cln in df.keys():
                    df[cln].replace(to_replace='False', value=0.0, inplace=True)
                    df[cln].replace(to_replace='True', value=1.0, inplace=True)
                    df[cln] = pd.to_numeric(df[cln])
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

    # Management data
    management_df = pd.read_csv(os.path.join(INPUTS_DIRPATH, MANAGEMENT_FILENAME), header=0, index_col=0)
    management_variables = {}
    for var_name in management_df.index:
        management_variables[var_name] = get_management_value(management_df, var_name)
    plant_density = management_variables.get('plant_density', {1: 250})
    inter_row = management_variables.get('inter_row', 0.15)
    Zsowing = management_variables.get('sowing_depth', 0.025)
    N_fertilizations = management_variables.get('N_fertilizations', {})


    # -- OUTPUTS CONFIGURATION --

    # Save the outputs with a full scan of the MTG at each time step (or at selected time steps)
    UPDATE_SHARED_DF = False
    if stored_times is None:
        stored_times = 'all'
    if not (stored_times == 'all' or isinstance(stored_times, list)):
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
    organs_all_data_list = []  # organs which belong to the axes: roots, phloem, grains
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
    # Create the stand using density pattern
    stand = AgronomicStand(sowing_density=plant_density[1], plant_density=plant_density[1], inter_row=inter_row, noise=0.) #todo to be adapted if multiple cultivars
    adel_wheat = AdelDyn(seed=1, scene_unit='m', leaves=echap_leaves(xy_model='Soissons_byleafclass'), stand=stand)

    # MTG generation
    if step_callback is not None and 'ADEL_mtg' in step_callback.keys():
        nff = update_parameters_all_models['elongwheat']['max_nb_leaves']
        g = step_callback['ADEL_mtg'](adel_wheat, INPUTS_DIRPATH, nff)  # Create a new MTG
    else:
        g = adel_wheat.load(directory=INPUTS_DIRPATH)  # read adelwheat inputs at t0 from a serialised MTG

    # ---------------------------------------------
    # ----- CONFIGURATION OF THE FACADES -------
    # ---------------------------------------------

    # -- ELONGWHEAT (created first because it is the only facade to add new metamers) --
    # Initial states
    elongwheat_hiddenzones_initial_state = inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME]
    elongwheat_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME]
    elongwheat_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME]

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
                                                            hydraulics=hydraulics,
                                                            optimal_growth_option=optimal_growth_option,
                                                            option_static=option_static,
                                                            update_parameters=update_parameters_elongwheat,
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
                                                               update_parameters=update_parameters_senescwheat,
                                                               update_shared_df=UPDATE_SHARED_DF)

    # -- FARQUHARWHEAT --
    # Initial states
    farquharwheat_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME]
    farquharwheat_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME]

    # Update parameters if specified
    if update_parameters_all_models and 'farquharwheat' in update_parameters_all_models:
        update_parameters_farquharwheat = update_parameters_all_models['farquharwheat']
    else:
        update_parameters_farquharwheat = None

    # Facade initialisation
    farquharwheat_facade_ = farquharwheat_facade.FarquharWheatFacade(g,
                                                                     farquharwheat_elements_initial_state,
                                                                     farquharwheat_axes_initial_state,
                                                                     shared_elements_inputs_outputs_df,
                                                                     stomatal_model_name=stomatal_model_name,
                                                                     hydraulics=hydraulics,
                                                                     update_parameters=update_parameters_farquharwheat,
                                                                     update_shared_df=UPDATE_SHARED_DF)

    # -- GROWTHWHEAT --
    # Initial states
    growthwheat_hiddenzones_initial_state = inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME]
    growthwheat_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME]
    growthwheat_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME]
    growthwheat_root_initial_state = inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME]

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
                                                               hydraulics=hydraulics,
                                                               update_parameters=update_parameters_growthwheat,
                                                               update_shared_df=UPDATE_SHARED_DF)

    # -- CNWHEAT --
    # Initial states
    cnwheat_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
        [i for i in cnwheat_facade.cnwheat_converter.AXES_VARIABLES if i in inputs_dataframes[AXES_INITIAL_STATE_FILENAME].columns]].copy()

    cnwheat_organs_initial_state = inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME][
        [i for i in cnwheat_facade.cnwheat_converter.ORGANS_VARIABLES if i in inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].columns]].copy()

    cnwheat_hiddenzones_initial_state = inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME][
        [i for i in cnwheat_facade.cnwheat_converter.HIDDENZONE_VARIABLES if i in inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME].columns]].copy()

    cnwheat_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
        [i for i in cnwheat_facade.cnwheat_converter.ELEMENTS_VARIABLES if i in inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns]].copy()

    cnwheat_soils_initial_state = inputs_dataframes[SOILS_INITIAL_STATE_FILENAME][
        [i for i in cnwheat_facade.cnwheat_converter.SOILS_VARIABLES if i in inputs_dataframes[SOILS_INITIAL_STATE_FILENAME].columns]].copy()
    if not hydraulics and 'SRWC' not in inputs_dataframes[SOILS_INITIAL_STATE_FILENAME].columns:
        cnwheat_soils_initial_state['SRWC'] = 100
    elif hydraulics and 'SRWC' not in inputs_dataframes[SOILS_INITIAL_STATE_FILENAME].columns:
        raise(ValueError('Hydraulics option is True but SRWC not found in {}.'.format(SOILS_INITIAL_STATE_FILENAME)))

    # Update parameters if specified
    if update_parameters_all_models and 'cnwheat' in update_parameters_all_models:
        update_parameters_cnwheat = update_parameters_all_models['cnwheat']
    else:
        update_parameters_cnwheat = {}

    # Facade initialisation
    cnwheat_facade_ = cnwheat_facade.CNWheatFacade(g,
                                                   CNWHEAT_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                   plant_density,
                                                   update_parameters_cnwheat,
                                                   cnwheat_axes_initial_state,
                                                   cnwheat_organs_initial_state,
                                                   cnwheat_hiddenzones_initial_state,
                                                   cnwheat_elements_initial_state,
                                                   cnwheat_soils_initial_state,
                                                   shared_axes_inputs_outputs_df,
                                                   shared_organs_inputs_outputs_df,
                                                   shared_hiddenzones_inputs_outputs_df,
                                                   shared_elements_inputs_outputs_df,
                                                   shared_soils_inputs_outputs_df,
                                                   tillers_replications=tillers_replications,
                                                   external_soil_model=external_soil_model,
                                                   update_shared_df=UPDATE_SHARED_DF)

    # -- TURGORGROWTH --
    drought_ongoing = False  # Is a drought event ongoing (bool)
    drought_passed = False   # Has a drought event occurred (bool)
    rehydration = False      # Is a rehydration period ongoing (bool)
    turgorgrowth_facade_ = None

    if hydraulics:
        # Initial states
        turgorgrowth_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
            [i for i in turgorgrowth_facade.turgorgrowth_converter.AXES_VARIABLES if
             i in inputs_dataframes[AXES_INITIAL_STATE_FILENAME].columns]].copy()

        turgorgrowth_organs_initial_state = inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME][
            [i for i in turgorgrowth_facade.turgorgrowth_converter.ORGANS_VARIABLES if
             i in inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].columns]].copy()

        turgorgrowth_hiddenzones_initial_state = inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME][
            [i for i in turgorgrowth_facade.turgorgrowth_converter.HIDDENZONE_VARIABLES if
             i in inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME].columns]].copy()

        turgorgrowth_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
            [i for i in turgorgrowth_facade.turgorgrowth_converter.ELEMENTS_VARIABLES if
             i in inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns]].copy()

        # Update parameters if specified
        if update_parameters_all_models and 'turgorgrowth' in update_parameters_all_models:
            update_parameters_turgorgrowth = update_parameters_all_models['turgorgrowth']
        else:
            update_parameters_turgorgrowth = {}

        turgorgrowth_soils_initial_state = inputs_dataframes[SOILS_INITIAL_STATE_FILENAME][
            [i for i in turgorgrowth_facade.turgorgrowth_converter.SOILS_VARIABLES if
             i in inputs_dataframes[SOILS_INITIAL_STATE_FILENAME].columns]].copy()

        # Facade initialisation
        turgorgrowth_facade_ = turgorgrowth_facade.TurgorGrowthFacade(g,
                                                                      TURGORGROWTH_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                                      update_parameters_turgorgrowth,
                                                                      turgorgrowth_axes_initial_state,
                                                                      turgorgrowth_hiddenzones_initial_state,
                                                                      turgorgrowth_elements_initial_state,
                                                                      turgorgrowth_organs_initial_state,
                                                                      turgorgrowth_soils_initial_state,
                                                                      shared_axes_inputs_outputs_df,
                                                                      shared_hiddenzones_inputs_outputs_df,
                                                                      shared_elements_inputs_outputs_df,
                                                                      shared_organs_inputs_outputs_df,
                                                                      shared_soils_inputs_outputs_df,
                                                                      update_shared_df=UPDATE_SHARED_DF)


    # Run cnwheat with constant nitrates concentration in the soil if specified
    if 'constant_Conc_Nitrates' in N_fertilizations.keys():
        cnwheat_facade_.soils[(1, 'MS')].constant_Conc_Nitrates = True  # TODO: make (1, 'MS') more general
        cnwheat_facade_.soils[(1, 'MS')].nitrates = N_fertilizations['constant_Conc_Nitrates'] * cnwheat_facade_.soils[(1, 'MS')].volume

    # Force root nitrate uptake if specified
    if external_soil_model and step_callback is not None:
        try:
            step_callback['nitrate_uptake'](0, cnwheat_facade_.population, g)
        except KeyError:
            print('Function name error in step_callback keys. It should be nitrate_uptake')

    # -- FSPMWHEAT --
    # Facade initialisation
    fspmwheat_facade_ = fspmwheat_facade.FSPMWheatFacade(g, elongwheat_facade_, growthwheat_facade_, farquharwheat_facade_, turgorgrowth_facade_)

    # Update geometry
    adel_wheat.update_geometry(g)
    if show_3Dplant:
        adel_wheat.plot(g)

    # ---------------------------------------------
    # -----      RUN OF THE SIMULATION      -------
    # ---------------------------------------------

    if run_simu:

        try:
            current_time_of_the_system = time.time()
            for t in range(START_TIME, SIMULATION_LENGTH, 1):
                print('t cnwheat is {}'.format(t))
                # if t == 1 or t >= 800:
                #     adel_wheat.scene(g).save('toto_{}.bgeom'.format(t))

                # run Caribu
                PARi = meteo.loc[t, ['PARi']].iloc[0]
                DOY = meteo.loc[t, ['DOY']].iloc[0]
                hour = meteo.loc[t, ['hour']].iloc[0]
                PARi_next_hours = meteo.loc[range(t, t + CARIBU_TIMESTEP), ['PARi']].sum().values[0]

                if (t % CARIBU_TIMESTEP == 0) and (PARi_next_hours > 0) and bool(g.property('geometry')):
                    run_caribu = True
                else:
                    run_caribu = False

                caribu_facade_.run(run_caribu, energy=PARi, DOY=DOY, hourTU=hour, latitude=48.85, sun_sky_option='sky',
                                   heterogeneous_canopy=heterogeneous_canopy, plant_density=plant_density[1], inter_row=inter_row)

                # run SenescWheat
                senescwheat_facade_.run()

                # Test for dead plant # TODO: adapt in case of multiple plants
                if not shared_elements_inputs_outputs_df.empty and \
                        np.nansum(shared_elements_inputs_outputs_df.loc[shared_elements_inputs_outputs_df['element'].isin(['StemElement', 'LeafElement1']), 'green_area']) == 0:
                    # append the inputs and outputs at current step to global lists
                    all_simulation_steps.append(t)
                    axes_all_data_list.append(shared_axes_inputs_outputs_df.copy())
                    organs_all_data_list.append(shared_organs_inputs_outputs_df.copy())
                    hiddenzones_all_data_list.append(shared_hiddenzones_inputs_outputs_df.copy())
                    elements_all_data_list.append(shared_elements_inputs_outputs_df.copy())
                    soils_all_data_list.append(shared_soils_inputs_outputs_df.copy())
                    print('Dead plant')
                    break

                # Run the rest of the model if the plant is alive
                # get the meteo of the current step
                Ta, Tsoil, ambient_CO2, RH, Ur = meteo.loc[t, ['air_temperature', 'soil_temperature', 'ambient_CO2', 'humidity', 'Wind']]

                # run FarquharWheat
                farquharwheat_facade_.run(Ta, ambient_CO2, RH, Ur)

                # run ElongWheat
                Tair, Tsoil = meteo.loc[t, ['air_temperature', 'soil_temperature']]
                elongwheat_facade_.run(Tair, Tsoil, Zsowing)

                # Update geometry
                adel_wheat.update_geometry(g)
                if show_3Dplant:
                    adel_wheat.plot(g)

                # run Turgorgrowth
                if hydraulics and turgorgrowth_facade_ is not None:
                    turgor_soil = turgorgrowth_facade_.soils[(1, 'MS')]
                    # Trigger drought
                    if drought_trigger is not None and 'green_area' in drought_trigger.keys():
                        if  (sum(g.property('green_area').values()) >= drought_trigger['green_area'] or drought_ongoing) and not drought_passed:
                            drought_ongoing = True
                            turgor_soil.constant_water_content = False
                        # Rehydration scenario. Only implemented for a hourly and linear rehydration scenario.
                        if rehydration_scenario is not None:
                            # Maximum of drought, start of rehydration
                            if turgor_soil.SRWC <= rehydration_scenario['stop_drought_SRWC'] and not rehydration:
                                rehydration = True
                                total_irrigation = (rehydration_scenario['SRWC_target'] * turgor_soil.PARAMETERS.AWC) / 100 - turgor_soil.water_content  # Total amount of water to add to the soil in order to reach the target SRWC
                                turgor_soil.hourly_irrigation = total_irrigation / (rehydration_scenario['rehydration_duration'] * 24)  # Amount of water to add each hour to reach the target SRWC at the end of the rehydration period

                            # Ongoing rehydration
                            elif rehydration:
                                # Target SRWC reached after rehydration, end of drought event
                                if turgor_soil.SRWC >= rehydration_scenario['SRWC_target']:
                                    rehydration = False
                                    drought_ongoing = False
                                    drought_passed = True
                                    turgor_soil.water_content = (rehydration_scenario['SRWC_target'] * turgor_soil.PARAMETERS.AWC) / 100
                                    turgor_soil.SRWC = rehydration_scenario['SRWC_target']
                                    turgor_soil.constant_water_content = True
                                    turgor_soil.hourly_rehydration = 0

                    turgorgrowth_facade_.run()

                    # Update geometry
                    adel_wheat.update_geometry(g)
                    if show_3Dplant:
                        adel_wheat.plot(g)
                    # adel_wheat.scene(g).save(r'adel_save\t{}.bgeom'.format(t))

                # run GrowthWheat
                growthwheat_facade_.run()

                # run cnwheat
                # N fertilization if any
                if t in N_fertilizations.keys():
                    cnwheat_facade_.soils[(1, 'MS')].nitrates += N_fertilizations[t]

                # Force root nitrate uptake if specified
                if external_soil_model and step_callback is not None:
                    try:
                        step_callback['nitrate_uptake'](t, cnwheat_facade_.population, g)
                    except KeyError:
                        print(
                            'Function name error in step_callback keys. It should be nitrate_uptake')
                # run CNWheat
                cnwheat_facade_.run(Tair, Tsoil)

                # Adel 3D plant save
                # if t_cnwheat % 24 == 0:
                #     adel_wheat.scene(g).save(os.path.join(OUTPUTS_DIRPATH, 'ADEL', 't{}.bgeom'.format(t_cnwheat)))

                # append outputs at current step to global lists
                if (stored_times == 'all') or (t in stored_times):
                    axes_outputs, elements_outputs, hiddenzones_outputs, organs_outputs, soils_outputs = fspmwheat_facade_.build_outputs_df_from_MTG()

                    all_simulation_steps.append(t)
                    axes_all_data_list.append(axes_outputs)
                    organs_all_data_list.append(organs_outputs)
                    hiddenzones_all_data_list.append(hiddenzones_outputs)
                    elements_all_data_list.append(elements_outputs)
                    soils_all_data_list.append(soils_outputs)

            execution_time = int(time.time() - current_time_of_the_system)
            print('\n' 'Simulation run in {}'.format(str(datetime.timedelta(seconds=execution_time))))

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
        if not run_simu:
            outputs_df_dict = {}

            for outputs_filename in (AXES_OUTPUTS_FILENAME,
                                     ORGANS_OUTPUTS_FILENAME,
                                     HIDDENZONES_OUTPUTS_FILENAME,
                                     ELEMENTS_OUTPUTS_FILENAME,
                                     SOILS_OUTPUTS_FILENAME):
                outputs_filepath = os.path.join(OUTPUTS_DIRPATH, outputs_filename)
                outputs_df = pd.read_csv(outputs_filepath, dtype={'is_over': str, 'is_growing': str})
                outputs_file_basename = outputs_filename.split('.')[0]
                outputs_df_dict[outputs_file_basename] = outputs_df

                # Assert states_filepaths were not opened during simulation run meaning that other filenames were saved
                tmp_filename = 'ACTUAL_{}.csv'.format(outputs_file_basename)
                tmp_path = os.path.join(OUTPUTS_DIRPATH, tmp_filename)
                assert not os.path.isfile(tmp_path), \
                    "File {} was saved because {} was opened during simulation run. Rename it before running postprocessing".format(tmp_filename, outputs_file_basename)

            time_grid = outputs_df_dict['axes_outputs'].t.unique()
            delta_t = (time_grid[1] - time_grid[0]) * HOUR_TO_SECOND_CONVERSION_FACTOR

        else:
            delta_t = CNWHEAT_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR

        # run the postprocessing
        postprocessing = cnwheat_facade.CNWheatFacade.postprocessing(axes_outputs_df=outputs_df_dict[AXES_OUTPUTS_FILENAME.split('.')[0]],
                                                                     hiddenzone_outputs_df=outputs_df_dict[HIDDENZONES_OUTPUTS_FILENAME.split('.')[0]],
                                                                     organs_outputs_df=outputs_df_dict[ORGANS_OUTPUTS_FILENAME.split('.')[0]],
                                                                     elements_outputs_df=outputs_df_dict[ELEMENTS_OUTPUTS_FILENAME.split('.')[0]],
                                                                     soils_outputs_df=outputs_df_dict[SOILS_OUTPUTS_FILENAME.split('.')[0]],
                                                                     delta_t=delta_t)

        if hydraulics:
            turgor_postprocessing = turgorgrowth_facade.TurgorGrowthFacade.postprocessing(axes_outputs_df=outputs_df_dict[AXES_OUTPUTS_FILENAME.split('.')[0]],
                                                                                          hiddenzone_outputs_df=outputs_df_dict[HIDDENZONES_OUTPUTS_FILENAME.split('.')[0]],
                                                                                          elements_outputs_df=outputs_df_dict[ELEMENTS_OUTPUTS_FILENAME.split('.')[0]],
                                                                                          organs_outputs_df=outputs_df_dict[ORGANS_OUTPUTS_FILENAME.split('.')[0]],
                                                                                          soils_outputs_df=outputs_df_dict[SOILS_OUTPUTS_FILENAME.split('.')[0]],
                                                                                          delta_t=delta_t)
            # Merge with cnwheat postprocessing
            mapping_scales = [('axes', AXES_INDEX_COLUMNS),
                ('elements', ELEMENTS_INDEX_COLUMNS),
                ('hiddenzones', HIDDENZONES_INDEX_COLUMNS),
                ('organs', ORGANS_INDEX_COLUMNS),
                ('soils', SOILS_INDEX_COLUMNS)]

            for scale, index_cols in mapping_scales:
                df_cnwheat = postprocessing.get(scale)
                df_turgor = turgor_postprocessing.get(scale)

                turgor_exclusive_cols = index_cols + [col for col in df_turgor.columns if col not in df_cnwheat.columns]
                df_turgor_filtered = df_turgor[turgor_exclusive_cols]

                # Left merge
                postprocessing[scale] = pd.merge(df_cnwheat, df_turgor_filtered, on=index_cols, how='left')



        for postprocessing_file_basename, postprocessing_filename, index_columns in (('axes', AXES_POSTPROCESSING_FILENAME, AXES_INDEX_COLUMNS),
                                                                                     ('hiddenzones', HIDDENZONES_POSTPROCESSING_FILENAME, HIDDENZONES_INDEX_COLUMNS),
                                                                                     ('organs', ORGANS_POSTPROCESSING_FILENAME, ORGANS_INDEX_COLUMNS),
                                                                                     ('elements', ELEMENTS_POSTPROCESSING_FILENAME, ELEMENTS_INDEX_COLUMNS),
                                                                                     ('soils', SOILS_POSTPROCESSING_FILENAME, SOILS_INDEX_COLUMNS)):
            postprocessing_filepath = os.path.join(POSTPROCESSING_DIRPATH, postprocessing_filename)
            postprocessing_df = postprocessing[postprocessing_file_basename]
            postprocessing_df.rename({'level_0': 't'}, axis=1, inplace=True)
            postprocessing_df = postprocessing_df.reindex(index_columns + postprocessing_df.columns.difference(index_columns).tolist(), axis=1, copy=False)
            postprocessing_df.to_csv(postprocessing_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))

    # ---------------------------------------------
    # -----            GRAPHS               -------
    # ---------------------------------------------

    if generate_graphs:
        # Delete previous graphs
        graphs = os.listdir(GRAPHS_DIRPATH)
        for graph in graphs:
            if graph.endswith(".PNG"):
                os.remove(os.path.join(GRAPHS_DIRPATH, graph))

        if not run_postprocessing:
            postprocessing = {}

            for postprocessing_filename in (AXES_POSTPROCESSING_FILENAME,
                                            ORGANS_POSTPROCESSING_FILENAME,
                                            HIDDENZONES_POSTPROCESSING_FILENAME,
                                            ELEMENTS_POSTPROCESSING_FILENAME,
                                            SOILS_POSTPROCESSING_FILENAME):
                postprocessing_filepath = os.path.join(POSTPROCESSING_DIRPATH, postprocessing_filename)
                postprocessing_df = pd.read_csv(postprocessing_filepath)
                postprocessing_file_basename = postprocessing_filename.split('_')[0]
                postprocessing[postprocessing_file_basename] = postprocessing_df

            outputs_df_dict = {}

            for outputs_filename in (AXES_OUTPUTS_FILENAME,
                                     ORGANS_OUTPUTS_FILENAME,
                                     HIDDENZONES_OUTPUTS_FILENAME,
                                     ELEMENTS_OUTPUTS_FILENAME,
                                     SOILS_OUTPUTS_FILENAME):
                outputs_filepath = os.path.join(OUTPUTS_DIRPATH, outputs_filename)
                outputs_df = pd.read_csv(outputs_filepath, dtype={'is_over': str, 'is_growing': str})
                outputs_file_basename = outputs_filename.split('.')[0]
                outputs_df_dict[outputs_file_basename] = outputs_df

                # Assert states_filepaths were not opened during simulation run meaning that other filenames were saved
                tmp_filename = 'ACTUAL_{}.csv'.format(outputs_file_basename)
                tmp_path = os.path.join(OUTPUTS_DIRPATH, tmp_filename)
                assert not os.path.isfile(tmp_path), \
                    "File {} was saved because {} was opened during simulation run. Rename it before running postprocessing".format(
                        tmp_filename, outputs_file_basename)

        # --- Generate graphs from postprocessing files
        plt.ioff()

        cnwheat_facade.CNWheatFacade.graphs(axes_postprocessing_df=postprocessing['axes'],
                                            hiddenzones_postprocessing_df=postprocessing['hiddenzones'],
                                            organs_postprocessing_df=postprocessing['organs'],
                                            elements_postprocessing_df=postprocessing['elements'],
                                            soils_postprocessing_df=postprocessing['soils'],
                                            meteo_data=meteo, graphs_dirpath=GRAPHS_DIRPATH)

        if hydraulics:
            turgorgrowth_facade.TurgorGrowthFacade.graphs(axes_postprocessing_df=postprocessing['axes'],
                                                          hiddenzones_postprocessing_df=postprocessing['hiddenzones'],
                                                          organs_postprocessing_df=postprocessing['organs'],
                                                          elements_postprocessing_df=postprocessing['elements'],
                                                          soils_postprocessing_df=postprocessing['soils'],
                                                          meteo_data=meteo, graphs_dirpath=GRAPHS_DIRPATH)
        # --- Additional graphs
        from openalea.fspmwheat import tools as fspmwheat_tools
        data_obs = pd.read_csv(os.path.join(INPUTS_DIRPATH, 'Ljutovac2002.csv'))
        fspmwheat_tools.additional_graphs(outputs_df_dict['axes_outputs'], outputs_df_dict['hiddenzones_outputs'], outputs_df_dict['elements_outputs'],
                                          postprocessing['axes'], postprocessing['hiddenzones'], postprocessing['elements'], postprocessing['organs'],
                                          plant_density[1], elongwheat_facade_._simulation.model.parameters.RERmax.items(), GRAPHS_DIRPATH, data_obs)
