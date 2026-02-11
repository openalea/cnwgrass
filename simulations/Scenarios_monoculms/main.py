# -*- coding: latin-1 -*-

import os
import random
import warnings
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from openalea.adel.Stand import AgronomicStand
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
    This script was used to simulate monoculms of wheat in different planting densities, light intensities and soil N concentrations. 
    Results were published in Gauthier et al. 2021 (https://doi.org/10.1093/insilicoplants/diab034)

"""

random.seed(1234)
np.random.seed(1234)

AXES_INDEX_COLUMNS = ['t', 'plant', 'axis']
ELEMENTS_INDEX_COLUMNS = ['t', 'plant', 'axis', 'metamer', 'organ', 'element']
HIDDENZONES_INDEX_COLUMNS = ['t', 'plant', 'axis', 'metamer']
ORGANS_INDEX_COLUMNS = ['t', 'plant', 'axis', 'organ']
SOILS_INDEX_COLUMNS = ['t', 'plant', 'axis']


def save_df_to_csv(df, outputs_filepath, precision):
    """
    Save pandas dataframes to csv

    :param pandas.DataFrame df: a pandas dataframe to be saved
    :param str outputs_filepath: the name of the CSV file to be saved
    :param int precision: number of decimals in CSV file
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


def main(simulation_length=2000, forced_start_time=0, run_simu=True, run_postprocessing=True, generate_graphs=True, run_from_outputs=False, stored_times=None,
         option_static=False, show_3Dplant=True, tillers_replications=None, heterogeneous_canopy=True,
         N_fertilizations=None, PLANT_DENSITY=None, INTER_ROW=0.15, update_parameters_all_models=None,
         INPUTS_PLANTSOIL_DIRPATH='inputs', INPUT_METEO_DIRPATH='inputs', METEO_FILENAME='meteo.csv', OUTPUTS_DIRPATH='outputs', POSTPROCESSING_DIRPATH='postprocessing', GRAPHS_DIRPATH='graphs'):
    """
    Run a simulation of fspmwheat with coupling to several models

    :param int simulation_length: length of the simulation (hours)
    :param int forced_start_time: desired start time (hour)
    :param bool run_simu: whether to run the simulation
    :param bool run_postprocessing: whether to run the postprocessing
    :param bool generate_graphs: whether to run the generate graphs
    :param bool run_from_outputs: whether to start a simulation from a specific time and initial states as found in previous outputs
    :param str or list stored_times: can be 'all', a list or an empty list
    :param bool option_static: Whether the model should be run for a static plant architecture
    :param bool show_3Dplant: whether to plot the scene in pgl viewer
    :param dict [str, float] tillers_replications: a dictionary with tiller id as key, and weight of replication as value.
    :param bool heterogeneous_canopy: Whether to create a duplicated heterogeneous canopy from the initial mtg.
    :param dict [int, float] or [str, float] N_fertilizations: a dictionary for N fertilisation regime {date: N_input}, with date in hour and N_input in µmol N nitrates
                                               or {'constant_Conc_Nitrates': val} for constant nitrates concentrations
    :param dict [int, int] PLANT_DENSITY: a dict with plant density per plant id (temporary used to account for different cultivars if needed) ; plant m-2
    :param float INTER_ROW: Inter-row spacing in the stand (m).
    :param dict update_parameters_all_models: a dict to update model parameters
                                             {'cnwheat': {'roots': {'param1': 'val1', 'param2': 'val2'},
                                                          'PhotosyntheticOrgan': {'param1': 'val1', 'param2': 'val2'}
                                                         },
                                              'elongwheat': {'param1': 'val1', 'param2': 'val2'}
                                             }
    :param str INPUTS_PLANTSOIL_DIRPATH: the path directory with plant and soil inputs
    :param str INPUT_METEO_DIRPATH: the path directory of meteo inputs
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
    INPUTS_DIRPATH = INPUTS_PLANTSOIL_DIRPATH

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
                new_initial_state = elements_previous_outputs[~elements_previous_outputs.is_over.isna()]  # TODO : verifier ce que ça donne avec des None
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
    meteo = pd.read_csv(os.path.join(INPUT_METEO_DIRPATH, METEO_FILENAME), index_col='t')

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

    # Create the stand using density pattern
    stand = AgronomicStand(sowing_density=PLANT_DENSITY[1], plant_density=PLANT_DENSITY[1], inter_row=INTER_ROW, noise=0.)

    # Create AdelDyn object and empty mtg
    adel_wheat = AdelDyn(seed=1, scene_unit='m', leaves=echap_leaves(xy_model='Soissons_byleafclass', top_leaves=0), stand=stand)
    axeT_df = pd.read_csv(os.path.join(INPUTS_DIRPATH, 'axeTable.csv'))

    # Final leaf number
    nff = 50
    if update_parameters_all_models and 'elongwheat' in update_parameters_all_models:
        nff = update_parameters_all_models['elongwheat'].get('max_nb_leaves', nff)
    axeT_df.HS_final = nff
    axeT_df.nf = nff
    axeT_df.nf_end = nff

    # mtg
    g = adel_wheat.build_stand(axeT_df)

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
    if show_3Dplant:
        adel_wheat.plot(g)

    # ---------------------------------------------
    # -----      RUN OF THE SIMULATION      -------
    # ---------------------------------------------

    if run_simu:

        try:
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

                caribu_facade_.run(run_caribu, energy=PARi, DOY=DOY, hourTU=hour, latitude=48.85, sun_sky_option='sky', heterogeneous_canopy=heterogeneous_canopy,
                                   plant_density=PLANT_DENSITY[1], inter_row=INTER_ROW)

                for t_senescwheat in range(t_caribu, t_caribu + SENESCWHEAT_TIMESTEP, SENESCWHEAT_TIMESTEP):
                    # run SenescWheat
                    senescwheat_facade_.run()

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
                            if show_3Dplant:
                                adel_wheat.plot(g)

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

                                        # Test for dead plant: if the whole shoot is senesced or if conc_sucrose_phloem < threshold
                                        # TODO: adapt in case of multiple plants
                                        # TODO: create a function with parameters (where ?)
                                        if (np.nansum(elements_outputs.loc[elements_outputs['element'].isin(['StemElement', 'LeafElement1']), 'green_area']) == 0) or \
                                                (organs_outputs.loc[organs_outputs['organ'] == 'phloem', 'sucrose'].values[0] / axes_outputs['mstruct'].values[0] < -50):
                                            break
                                else:
                                    # Continue if CN-Wheat loop wasn't broken because of dead plant.
                                    continue
                                # CN-Wheat loop was broken, break the Caribu loop.
                                break
                            else:
                                # Continue if GrowthWheat loop wasn't broken because of dead plant.
                                continue
                            # GrowthWheat loop was broken, break the Caribu loop.
                            break
                        else:
                            # Continue if ElongWheat loop wasn't broken because of dead plant.
                            continue
                        # ElongWheat loop was broken, break the Caribu loop.
                        break
                    else:
                        # Continue if FarqhuarWheat loop wasn't broken because of dead plant.
                        continue
                    # FarqhuarWheat loop was broken, break the Caribu loop.
                    break

                else:
                    # Continue if SenescWheat loop wasn't broken because of dead plant.
                    continue
                # SenescWheat loop was broken, break the Caribu loop.
                break

        except Exception as ex:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            template = "An exception of type {0} occurred. Arguments:\n{1!r}"
            message = template.format(type(ex).__name__, ex.args)
            print(message, fname, exc_tb.tb_lineno)

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
    # ---------      POST-PROCESSING      ---------
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
                outputs_df = pd.read_csv(outputs_filepath)
                outputs_file_basename = outputs_filename.split('.')[0]
                outputs_df_dict[outputs_file_basename] = outputs_df

                # Assert states_filepaths were not opened during simulation run meaning that other filenames were saved
                tmp_filename = 'ACTUAL_{}.csv'.format(outputs_file_basename)
                tmp_path = os.path.join(OUTPUTS_DIRPATH, tmp_filename)
                assert not os.path.isfile(tmp_path), \
                    "File {} was saved because {} was opened during simulation run. Rename it before running postprocessing".format(tmp_filename, outputs_file_basename)

            time_grid = list(outputs_df_dict.values())[0].t
            delta_t = (time_grid.loc[1] - time_grid.loc[0]) * HOUR_TO_SECOND_CONVERSION_FACTOR

        else:
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
        if not run_postprocessing:
            postprocessing_df_dict = {}

            for postprocessing_filename in (AXES_POSTPROCESSING_FILENAME,
                                            ORGANS_POSTPROCESSING_FILENAME,
                                            HIDDENZONES_POSTPROCESSING_FILENAME,
                                            ELEMENTS_POSTPROCESSING_FILENAME,
                                            SOILS_POSTPROCESSING_FILENAME):
                postprocessing_filepath = os.path.join(POSTPROCESSING_DIRPATH, postprocessing_filename)
                postprocessing_df = pd.read_csv(postprocessing_filepath)
                postprocessing_file_basename = postprocessing_filename.split('.')[0]
                postprocessing_df_dict[postprocessing_file_basename] = postprocessing_df

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

        # --- Additional graphs
        from openalea.cnwheat import tools as cnwheat_tools
        colors = ['blue', 'darkorange', 'green', 'red', 'darkviolet', 'gold', 'magenta', 'brown', 'darkcyan', 'grey', 'lime']
        colors = colors + colors

        # 0) Phyllochron
        df_SAM = df_SAM[df_SAM['axis'] == 'MS']
        df_hz = postprocessing_df_dict[hiddenzones_postprocessing_file_basename]
        grouped_df = df_hz[df_hz['axis'] == 'MS'].groupby(['plant', 'metamer'])[['t', 'leaf_is_emerged']]
        leaf_emergence = {}
        for group_name, data in grouped_df:
            plant, metamer = group_name[0], group_name[1]
            if metamer == 3 or True not in data['leaf_is_emerged'].unique(): continue
            leaf_emergence_t = data[data['leaf_is_emerged'] == True].iloc[0]['t']
            leaf_emergence[(plant, metamer)] = leaf_emergence_t

        phyllochron = {'plant': [], 'metamer': [], 'phyllochron': []}
        for key, leaf_emergence_t in sorted(leaf_emergence.items()):
            plant, metamer = key[0], key[1]
            if metamer == 4:
                continue
            phyllochron['plant'].append(plant)
            phyllochron['metamer'].append(metamer)
            prev_leaf_emergence_t = leaf_emergence[(plant, metamer - 1)]
            if df_SAM[(df_SAM['t'] == leaf_emergence_t) | (df_SAM['t'] == prev_leaf_emergence_t)].sum_TT.count() == 2:
                phyllo_DD = df_SAM[(df_SAM['t'] == leaf_emergence_t)].sum_TT.values[0] - df_SAM[(df_SAM['t'] == prev_leaf_emergence_t)].sum_TT.values[0]
            else:
                phyllo_DD = np.nan
            phyllochron['phyllochron'].append(phyllo_DD)

        if len(phyllochron['metamer']) > 0:
            fig, ax = plt.subplots()
            plt.xlim((int(min(phyllochron['metamer']) - 1), int(max(phyllochron['metamer']) + 1)))
            plt.ylim(ymin=0, ymax=150)
            ax.plot(phyllochron['metamer'], phyllochron['phyllochron'], color='b', marker='o')
            for i, j in zip(phyllochron['metamer'], phyllochron['phyllochron']):
                ax.annotate(str(int(round(j, 0))), xy=(i, j + 2), ha='center')
            ax.set_xlabel('Leaf number')
            ax.set_ylabel('Phyllochron (Degree Day)')
            ax.set_title('phyllochron')
            plt.savefig(os.path.join(GRAPHS_DIRPATH, 'phyllochron' + '.PNG'))
            plt.close()

        # 1) Comparison Dimensions with Ljutovac 2002
        bchmk = pd.read_csv(r'inputs\Ljutovac2002.csv')
        res = pd.read_csv(os.path.join(OUTPUTS_DIRPATH, HIDDENZONES_OUTPUTS_FILENAME))
        res = res[(res['axis'] == 'MS') & (res['plant'] == 1) & ~np.isnan(res.leaf_Lmax)].copy()
        res_IN = res[~ np.isnan(res.internode_Lmax)]
        last_value_idx = res.groupby(['metamer'])['t'].transform('max') == res['t']
        res = res[last_value_idx].copy()
        res['lamina_Wmax'] = res.leaf_Wmax
        res['lamina_W_Lg'] = res.leaf_Wmax / res.lamina_Lmax
        bchmk = bchmk[bchmk.metamer >= min(res.metamer)]
        bchmk['lamina_W_Lg'] = bchmk.lamina_Wmax / bchmk.lamina_Lmax
        last_value_idx = res_IN.groupby(['metamer'])['t'].transform('max') == res_IN['t']
        res_IN = res_IN[last_value_idx].copy()
        res = res[['metamer', 'leaf_Lmax', 'lamina_Lmax', 'sheath_Lmax', 'lamina_Wmax', 'lamina_W_Lg', 'SSLW', 'LSSW']].merge(res_IN[['metamer', 'internode_Lmax']], left_on='metamer',
                                                                                                                              right_on='metamer', how='outer').copy()

        var_list = ['leaf_Lmax', 'lamina_Lmax', 'sheath_Lmax', 'lamina_Wmax', 'internode_Lmax']
        for var in list(var_list):
            fig, ax = plt.subplots()
            plt.xlim((int(min(res.metamer) - 1), int(max(res.metamer) + 1)))
            plt.ylim(ymin=0, ymax=np.nanmax(list(res[var] * 100 * 1.05) + list(bchmk[var] * 1.05)))

            tmp = res[['metamer', var]].drop_duplicates()

            line1 = ax.plot(tmp.metamer, tmp[var] * 100, color='c', marker='o')
            line2 = ax.plot(bchmk.metamer, bchmk[var], color='orange', marker='o')

            ax.set_ylabel(var + ' (cm)')
            ax.set_title(var)
            ax.legend((line1[0], line2[0]), ('Simulation', 'Ljutovac 2002'), loc=2)
            plt.savefig(os.path.join(GRAPHS_DIRPATH, var + '.PNG'))
            plt.close()

        var = 'lamina_W_Lg'
        fig, ax = plt.subplots()
        plt.xlim((int(min(res.metamer) - 1), int(max(res.metamer) + 1)))
        plt.ylim(ymin=0, ymax=np.nanmax(list(res[var] * 1.05) + list(bchmk[var] * 1.05)))
        tmp = res[['metamer', var]].drop_duplicates()
        line1 = ax.plot(tmp.metamer, tmp[var], color='c', marker='o')
        line2 = ax.plot(bchmk.metamer, bchmk[var], color='orange', marker='o')
        ax.set_ylabel(var)
        ax.set_title(var)
        ax.legend((line1[0], line2[0]), ('Simulation', 'Ljutovac 2002'), loc=2)
        plt.savefig(os.path.join(GRAPHS_DIRPATH, var + '.PNG'))
        plt.close()

        # 1bis) Comparison Structural Masses vs. adaptation from Bertheloot 2008

        # SSLW Laminae
        bchmk = pd.DataFrame.from_dict({1: 15, 2: 23, 3: 25, 4: 18, 5: 22, 6: 25, 7: 20, 8: 23, 9: 26, 10: 28, 11: 31}, orient='index').rename(columns={0: 'SSLW'})
        bchmk.index.name = 'metamer'
        bchmk = bchmk.reset_index()
        bchmk = bchmk[bchmk.metamer >= min(res.metamer)]

        fig, ax = plt.subplots()
        plt.xlim((int(min(res.metamer) - 1), int(max(res.metamer) + 1)))
        plt.ylim(ymin=0, ymax=50)

        tmp = res[['metamer', 'SSLW']].drop_duplicates()

        line1 = ax.plot(tmp.metamer, tmp.SSLW, color='c', marker='o')
        line2 = ax.plot(bchmk.metamer, bchmk.SSLW, color='orange', marker='o')

        ax.set_ylabel('Structural Specific Lamina Weight (g.m-2)')
        ax.set_title('Structural Specific Lamina Weight')
        ax.legend((line1[0], line2[0]), ('Simulation', 'adapated from Bertheloot 2008'), loc=3)
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'SSLW.PNG'))
        plt.close()

        # LWS Sheaths
        bchmk = pd.DataFrame.from_dict({1: 0.08, 2: 0.09, 3: 0.11, 4: 0.18, 5: 0.17, 6: 0.21, 7: 0.24, 8: 0.4, 9: 0.5, 10: 0.55, 11: 0.65}, orient='index').rename(columns={0: 'LSSW'})
        bchmk.index.name = 'metamer'
        bchmk = bchmk.reset_index()
        bchmk = bchmk[bchmk.metamer >= min(res.metamer)]

        fig, ax = plt.subplots()
        plt.xlim((int(min(res.metamer) - 1), int(max(res.metamer) + 1)))
        plt.ylim(ymin=0, ymax=0.8)

        tmp = res[['metamer', 'LSSW']].drop_duplicates()

        line1 = ax.plot(tmp.metamer, tmp.LSSW, color='c', marker='o')
        line2 = ax.plot(bchmk.metamer, bchmk.LSSW, color='orange', marker='o')

        ax.set_ylabel('Lineic Structural Sheath Weight (g.m-1)')
        ax.set_title('Lineic Structural Sheath Weight')
        ax.legend((line1[0], line2[0]), ('Simulation', 'adapated from Bertheloot 2008'), loc=2)
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'LSSW.PNG'))
        plt.close()

        # 2) LAI
        df_elt['green_area_rep'] = df_elt.green_area * df_elt.nb_replications
        grouped_df = df_elt[(df_elt.axis == 'MS') & (df_elt.element == 'LeafElement1')].groupby(['t', 'plant'])
        LAI_dict = {'t': [], 'plant': [], 'LAI': []}
        for name, data in grouped_df:
            t, plant = name[0], name[1]
            LAI_dict['t'].append(t)
            LAI_dict['plant'].append(plant)
            LAI_dict['LAI'].append(data['green_area_rep'].sum() * PLANT_DENSITY[plant])

        cnwheat_tools.plot_cnwheat_ouputs(pd.DataFrame(LAI_dict), 't', 'LAI', x_label='Time (Hour)', y_label='LAI', plot_filepath=os.path.join(GRAPHS_DIRPATH, 'LAI.PNG'), explicit_label=False)

        # 3) RUE
        df_elt['day'] = df_elt['t'] // 24 + 1
        days = df_elt['day'].unique()
        df_elt['PARa_MJ'] = df_elt['PARa'] * df_elt['green_area'] * df_elt['nb_replications'] * 3600 / 4.6 * 10 ** -6  # Il faudrait idealement utiliser les calculcs green_area et PARa des talles
        df_elt['RGa_MJ'] = df_elt['PARa'] * df_elt['green_area'] * df_elt['nb_replications'] * 3600 / 2.02 * 10 ** -6  # Il faudrait idealement utiliser les calculcs green_area et PARa des talles
        PARa = df_elt.groupby(['day'])['PARa_MJ'].agg('sum')
        PARa_cum = np.cumsum(PARa)

        df_axe = postprocessing_df_dict[axes_postprocessing_file_basename]
        df_axe['day'] = df_axe['t'] // 24 + 1
        # Total_Photosynthesis = df_axe.groupby(['day'])['Tillers_Photosynthesis'].agg('sum')
        sum_dry_mass_shoot = df_axe.groupby(['day'])['sum_dry_mass_shoot'].agg('max')
        sum_dry_mass = df_axe.groupby(['day'])['sum_dry_mass'].agg('max')

        RUE_shoot = np.polyfit(PARa_cum, sum_dry_mass_shoot, 1)[0]
        RUE_plant = np.polyfit(PARa_cum, sum_dry_mass, 1)[0]

        fig, ax = plt.subplots()
        ax.plot(PARa_cum, sum_dry_mass_shoot, label='Shoot dry mass (g)')
        ax.plot(PARa_cum, sum_dry_mass, label='Plant dry mass (g)')
        ax.legend(prop={'size': 10}, framealpha=0.5, loc='center left', bbox_to_anchor=(1, 0.815), borderaxespad=0.)
        ax.set_xlabel('Cumulative absorbed PAR (MJ)')
        ax.set_ylabel('Dry mass (g)')
        ax.set_title('RUE')
        plt.text(max(PARa_cum) * 0.02, max(sum_dry_mass) * 0.95, 'RUE shoot : {0:.2f} , RUE plant : {1:.2f}'.format(round(RUE_shoot, 2), round(RUE_plant, 2)))
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'RUE.PNG'), dpi=200, format='PNG', bbox_inches='tight')
        plt.close()

        fig, ax = plt.subplots()
        ax.plot(days, sum_dry_mass_shoot, label='Shoot dry mass (g)')
        ax.plot(days, sum_dry_mass, label='Plant dry mass (g)')
        ax.plot(days, PARa_cum, label='Cumulative absorbed PAR (MJ)')
        ax.legend(prop={'size': 10}, framealpha=0.5, loc='center left', bbox_to_anchor=(1, 0.815), borderaxespad=0.)
        ax.set_xlabel('Days')
        ax.set_title('RUE investigations')
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'RUE2.PNG'), dpi=200, format='PNG', bbox_inches='tight')
        plt.close()

        # 3bis) Photosynthetic efficiency of the plant
        df_elt['Photosynthesis_tillers'] = df_elt.Ag * df_elt.green_area * df_elt.nb_replications.fillna(1.)
        df_elt['PARa_tot_tillers'] = df_elt.PARa * df_elt.green_area * df_elt.nb_replications.fillna(1.)
        df_elt['green_area_tillers'] = df_elt.green_area * df_elt.nb_replications.fillna(1.)
        photo_y = df_elt.groupby(['t'], as_index=False).agg({'Photosynthesis_tillers': 'sum', 'PARa_tot_tillers': 'sum', 'green_area_tillers': 'sum'})
        photo_y['Photosynthetic_efficiency_plant'] = photo_y.Photosynthesis_tillers / photo_y.PARa_tot_tillers

        fig, ax = plt.subplots()
        ax.plot(photo_y.t, photo_y.Photosynthetic_efficiency_plant)
        ax.set_ylim(bottom=0.)
        ax.set_xlabel('t')
        ax.set_ylabel(u'Photosynthetic efficiency (µmol C/µmol PARa)')
        ax.set_title('Photosynthetic efficiency of the plant')
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Photosynthetic_efficiency_plant.PNG'), dpi=200, format='PNG', bbox_inches='tight')
        plt.close()

        PARa2 = df_elt.groupby(['day'])['PARa_tot_tillers'].agg('sum')
        PARa2_cum = np.cumsum(PARa2)
        Photosynthesis = df_elt.groupby(['day'])['Photosynthesis_tillers'].agg('sum')
        Photosynthesis_cum = np.cumsum(Photosynthesis)

        avg_photo_y = np.polyfit(PARa2_cum, Photosynthesis_cum, 1)[0]

        fig, ax = plt.subplots()
        ax.plot(PARa2_cum, Photosynthesis_cum)
        ax.set_xlabel(u'Cumulative absorbed PAR (µmol.s$^{-1}$)')
        ax.set_ylabel(u'Cumulative photosynthesis (µmol C.s$^{-1}$)')
        ax.set_title('Average photosynthetic efficiency of the plant: {:.2%}'.format(round(avg_photo_y, 3)))
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'Photosynthetic_efficiency_plant2.PNG'), dpi=200, format='PNG', bbox_inches='tight')
        plt.close()

        # 4) Sum thermal time
        df_SAM = df_SAM[df_SAM['axis'] == 'MS']
        fig, ax = plt.subplots()
        ax.plot(df_SAM['t'], df_SAM['sum_TT'])
        ax.set_xlabel('Hours')
        ax.set_ylabel('Thermal Time')
        ax.set_title('Thermal Time')
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'SumTT.PNG'), dpi=200, format='PNG', bbox_inches='tight')
        plt.close()

        # 5) Residual N : ratio_N_mstruct_max
        df_elt_outputs = pd.read_csv(os.path.join(OUTPUTS_DIRPATH, ELEMENTS_OUTPUTS_FILENAME))
        df_elt_outputs = df_elt_outputs.loc[df_elt_outputs.axis == 'MS']
        df_elt_outputs = df_elt_outputs.loc[df_elt_outputs.mstruct != 0]
        df_elt_outputs['N_content_total'] = df_elt_outputs['N_content_total'] * 100
        x_name = 't'
        x_label = 'Time (Hour)'
        graph_variables_ph_elements = {'N_content_total': u'N content in green + senesced tissues (% mstruct)'}
        for org_ph in (['blade'], ['sheath'], ['internode'], ['peduncle', 'ear']):
            for variable_name, variable_label in graph_variables_ph_elements.items():
                graph_name = variable_name + '_' + '_'.join(org_ph) + '.PNG'
                cnwheat_tools.plot_cnwheat_ouputs(df_elt_outputs,
                                                  x_name=x_name,
                                                  y_name=variable_name,
                                                  x_label=x_label,
                                                  y_label=variable_label,
                                                  colors=[colors[i - 1] for i in df_elt_outputs.metamer.unique().tolist()],
                                                  filters={'organ': org_ph},
                                                  plot_filepath=os.path.join(GRAPHS_DIRPATH, graph_name),
                                                  explicit_label=False)

        # 6) C usages relatif to Net Photosynthesis
        df_org = postprocessing_df_dict[organs_postprocessing_file_basename]
        df_roots = df_org[df_org['organ'] == 'roots'].copy()
        df_roots['day'] = df_roots['t'] // 24 + 1
        df_phloem = df_org[df_org['organ'] == 'phloem'].copy()
        df_phloem['day'] = df_phloem['t'] // 24 + 1

        # --- C usages relatif to Net Photosynthesis
        AMINO_ACIDS_C_RATIO = 4.15  #: Mean number of mol of C in 1 mol of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly)
        AMINO_ACIDS_N_RATIO = 1.25  #: Mean number of mol of N in 1 mol of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly)

        # Photosynthesis
        df_elt['Photosynthesis_tillers'] = df_elt['Photosynthesis'].fillna(0) * df_elt['nb_replications'].fillna(1.)
        Tillers_Photosynthesis_Ag = df_elt.groupby(['t'], as_index=False).agg({'Photosynthesis_tillers': 'sum'})
        C_usages = pd.DataFrame({'t': Tillers_Photosynthesis_Ag['t']})
        C_usages['C_produced'] = np.cumsum(Tillers_Photosynthesis_Ag.Photosynthesis_tillers)

        # Respiration
        C_usages['Respi_roots'] = np.cumsum(df_axe.C_respired_roots)
        C_usages['Respi_shoot'] = np.cumsum(df_axe.C_respired_shoot)

        # Exudation
        C_usages['exudation'] = np.cumsum(df_axe.C_exudated.fillna(0))

        # Structural growth
        C_consumption_mstruct_roots = df_roots.sucrose_consumption_mstruct.fillna(0) + df_roots.AA_consumption_mstruct.fillna(0) * AMINO_ACIDS_C_RATIO / AMINO_ACIDS_N_RATIO
        C_usages['Structure_roots'] = np.cumsum(C_consumption_mstruct_roots.reset_index(drop=True))

        df_hz['C_consumption_mstruct'] = df_hz.sucrose_consumption_mstruct.fillna(0) + df_hz.AA_consumption_mstruct.fillna(0) * AMINO_ACIDS_C_RATIO / AMINO_ACIDS_N_RATIO
        df_hz['C_consumption_mstruct_tillers'] = df_hz['C_consumption_mstruct'] * df_hz['nb_replications']
        C_consumption_mstruct_shoot = df_hz.groupby(['t'])['C_consumption_mstruct_tillers'].sum()
        C_usages['Structure_shoot'] = np.cumsum(C_consumption_mstruct_shoot.reset_index(drop=True))

        # Non structural C
        df_phloem['C_NS'] = df_phloem.sucrose.fillna(0) + df_phloem.amino_acids.fillna(0) * AMINO_ACIDS_C_RATIO / AMINO_ACIDS_N_RATIO
        C_NS_phloem_init = df_phloem.C_NS - df_phloem.C_NS[0]
        C_usages['NS_phloem'] = C_NS_phloem_init.reset_index(drop=True)

        df_elt['C_NS'] = df_elt.sucrose.fillna(0) + df_elt.fructan.fillna(0) + df_elt.starch.fillna(0) + (
                df_elt.amino_acids.fillna(0) + df_elt.proteins.fillna(0)) * AMINO_ACIDS_C_RATIO / AMINO_ACIDS_N_RATIO
        df_elt['C_NS_tillers'] = df_elt['C_NS'] * df_elt['nb_replications'].fillna(1.)
        C_elt = df_elt.groupby(['t']).agg({'C_NS_tillers': 'sum'})

        df_hz['C_NS'] = df_hz.sucrose.fillna(0) + df_hz.fructan.fillna(0) + (df_hz.amino_acids.fillna(0) + df_hz.proteins.fillna(0)) * AMINO_ACIDS_C_RATIO / AMINO_ACIDS_N_RATIO
        df_hz['C_NS_tillers'] = df_hz['C_NS'] * df_hz['nb_replications'].fillna(1.)
        C_hz = df_hz.groupby(['t']).agg({'C_NS_tillers': 'sum'})

        df_roots['C_NS'] = df_roots.sucrose.fillna(0) + df_roots.amino_acids.fillna(0) * AMINO_ACIDS_C_RATIO / AMINO_ACIDS_N_RATIO

        C_NS_autre = df_roots.C_NS.reset_index(drop=True) + C_elt.C_NS_tillers + C_hz.C_NS_tillers
        C_NS_autre_init = C_NS_autre - C_NS_autre[0]
        C_usages['NS_other'] = C_NS_autre_init.reset_index(drop=True)

        # Total
        C_usages['C_budget'] = (C_usages.Respi_roots + C_usages.Respi_shoot + C_usages.exudation + C_usages.Structure_roots + C_usages.Structure_shoot + C_usages.NS_phloem + C_usages.NS_other) / \
                               C_usages.C_produced

        #  Graph
        fig, ax = plt.subplots()
        ax.plot(C_usages.t, C_usages.Structure_shoot / C_usages.C_produced * 100,
                label=u'Structural mass - Shoot', color='g')
        ax.plot(C_usages.t, C_usages.Structure_roots / C_usages.C_produced * 100,
                label=u'Structural mass - Roots', color='r')
        ax.plot(C_usages.t, (C_usages.NS_phloem + C_usages.NS_other) / C_usages.C_produced * 100, label=u'Non-structural C', color='darkorange')
        ax.plot(C_usages.t, (C_usages.Respi_roots + C_usages.Respi_shoot) / C_usages.C_produced * 100, label=u'C loss by respiration', color='b')
        ax.plot(C_usages.t, C_usages.exudation / C_usages.C_produced * 100, label=u'C loss by exudation', color='c')

        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax.set_xlabel('Time (h)')
        ax.set_ylabel(u'Carbon usages : Photosynthesis (%)')
        ax.set_ylim(bottom=0, top=100.)

        fig.suptitle(u'Total cumulated usages are ' + str(round(C_usages.C_budget.tail(1) * 100, 0)) + u' % of Photosynthesis')

        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'C_usages_cumulated.PNG'), format='PNG', bbox_inches='tight')
        plt.close()

        # 7) N shoot in function of above-ground biomass
        fig, ax = plt.subplots()
        ax.plot(df_axe.sum_dry_mass_shoot, df_axe.N_content_shoot)
        ax.set_xlabel(u'Shoot biomass (g)')
        ax.set_ylabel(u'N shoot (% DM)')
        ax.set_ylim(bottom=0)
        ax.set_xlim(left=0)
        plt.savefig(os.path.join(GRAPHS_DIRPATH, 'N_dilution.PNG'), format='PNG', bbox_inches='tight')
        plt.close()


if __name__ == '__main__':
    main(simulation_length=2000, forced_start_time=0,
         run_simu=True, run_postprocessing=True, generate_graphs=True, run_from_outputs=False,
         show_3Dplant=False, heterogeneous_canopy=True, METEO_FILENAME="meteo_PAR500.csv")
