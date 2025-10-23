# -*- coding: latin-1 -*-

import os

import numpy as np
import pandas as pd
import random

from openalea.adel.adel_dynamic import AdelDyn
from openalea.adel.echap_leaf import echap_leaves

from openalea.fspmwheat import caribu_facade
from openalea.fspmwheat import cnwheat_facade
from openalea.fspmwheat import elongwheat_facade
from openalea.fspmwheat import farquharwheat_facade
from openalea.fspmwheat import growthwheat_facade
from openalea.fspmwheat import senescwheat_facade
from openalea.fspmwheat import fspmwheat_facade

from openalea.cnwheat import tools as cnwheat_tools
from openalea.cnwheat import simulation as cnwheat_simulation

"""
    test_cnwheat
    ~~~~~~~~~~~~

    Test the CN-Wheat model.

    You must first install :mod:`cnwheat` (and add it to your PYTHONPATH)
    before running this script with the command `python`.

    To get a coverage report of :mod:`cnwheat`, use the command:
    `nosetests --with-coverage --cover-package=cnwheat test_cnwheat.py`.

    CSV files must contain only ASCII characters and ',' as separator.

"""

random.seed(1234)
np.random.seed(1234)

# Number of seconds in 1 hour
HOUR_TO_SECOND_CONVERSION_FACTOR = 3600

# the precision to use for quantitative comparison test
PRECISION = 4


def test_run(overwrite_desired_data=False):
    # ---------------------------------------------
    # ----- CONFIGURATION OF THE SIMULATION -------
    # ---------------------------------------------

    # -- INPUTS CONFIGURATION --

    # Path of the directory which contains the inputs of the model
    INPUTS_DIRPATH = 'inputs'

    # Name of the CSV files which describes the initial state of the system
    AXES_INITIAL_STATE_FILENAME = 'axes_initial_state.csv'
    ORGANS_INITIAL_STATE_FILENAME = 'organs_initial_state.csv'
    HIDDENZONES_INITIAL_STATE_FILENAME = 'hiddenzones_initial_state.csv'
    ELEMENTS_INITIAL_STATE_FILENAME = 'elements_initial_state.csv'
    SOILS_INITIAL_STATE_FILENAME = 'soils_initial_state.csv'
    METEO_FILENAME = 'meteo_Ljutovac2002.csv'

    # Read the inputs from CSV files and create inputs dataframes
    inputs_dataframes = {}
    for inputs_filename in (AXES_INITIAL_STATE_FILENAME,
                            ORGANS_INITIAL_STATE_FILENAME,
                            HIDDENZONES_INITIAL_STATE_FILENAME,
                            ELEMENTS_INITIAL_STATE_FILENAME,
                            SOILS_INITIAL_STATE_FILENAME):
        inputs_dataframe = pd.read_csv(os.path.join(INPUTS_DIRPATH, inputs_filename))
        inputs_dataframes[inputs_filename] = inputs_dataframe.where(inputs_dataframe.notnull(), None)

    # -- SIMULATION PARAMETERS --
    START_TIME = 0
    SIMULATION_LENGTH = 24
    PLANT_DENSITY = {1: 410}

    # define the time step in hours for each simulator
    CARIBU_TIMESTEP = 4
    SENESCWHEAT_TIMESTEP = 2
    FARQUHARWHEAT_TIMESTEP = 2
    ELONGWHEAT_TIMESTEP = 1
    GROWTHWHEAT_TIMESTEP = 1
    CNWHEAT_TIMESTEP = 1

    # Name of the CSV files which contains the meteo data
    meteo = pd.read_csv(os.path.join(INPUTS_DIRPATH, METEO_FILENAME), index_col='t')

    # -- OUTPUTS CONFIGURATION --

    # Path of the directory which contains the inputs of the model
    OUTPUTS_DIRPATH = 'outputs'

    # Name of the CSV files which will contain the outputs of the model
    DESIRED_AXES_OUTPUTS_FILENAME = 'desired_axes_outputs.csv'
    DESIRED_ORGANS_OUTPUTS_FILENAME = 'desired_organs_outputs.csv'
    DESIRED_HIDDENZONES_OUTPUTS_FILENAME = 'desired_hiddenzones_outputs.csv'
    DESIRED_ELEMENTS_OUTPUTS_FILENAME = 'desired_elements_outputs.csv'
    DESIRED_SOILS_OUTPUTS_FILENAME = 'desired_soils_outputs.csv'
    ACTUAL_AXES_OUTPUTS_FILENAME = 'actual_axes_outputs.csv'
    ACTUAL_ORGANS_OUTPUTS_FILENAME = 'actual_organs_outputs.csv'
    ACTUAL_HIDDENZONES_OUTPUTS_FILENAME = 'actual_hiddenzones_outputs.csv'
    ACTUAL_ELEMENTS_OUTPUTS_FILENAME = 'actual_elements_outputs.csv'
    ACTUAL_SOILS_OUTPUTS_FILENAME = 'actual_soils_outputs.csv'

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

    phytoT = os.path.join(INPUTS_DIRPATH, 'phytoT.csv')

    # Update some parameters
    update_cnwheat_parameters = {'SL_ratio_d': 0.25}

    # Facade initialisation
    elongwheat_facade_ = elongwheat_facade.ElongWheatFacade(g,
                                                            ELONGWHEAT_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                            elongwheat_axes_initial_state,
                                                            elongwheat_hiddenzones_initial_state,
                                                            elongwheat_elements_initial_state,
                                                            shared_axes_inputs_outputs_df,
                                                            shared_hiddenzones_inputs_outputs_df,
                                                            shared_elements_inputs_outputs_df,
                                                            adel_wheat, phytoT, update_cnwheat_parameters)

    # -- CARIBU --
    caribu_facade_ = caribu_facade.CaribuFacade(g,
                                                shared_elements_inputs_outputs_df,
                                                adel_wheat)

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

    # Update some parameters
    update_cnwheat_parameters = {'AGE_EFFECT_SENESCENCE': 10000}

    # Facade initialisation
    senescwheat_facade_ = senescwheat_facade.SenescWheatFacade(g,
                                                               SENESCWHEAT_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                               senescwheat_roots_initial_state,
                                                               senescwheat_axes_initial_state,
                                                               senescwheat_elements_initial_state,
                                                               shared_organs_inputs_outputs_df,
                                                               shared_axes_inputs_outputs_df,
                                                               shared_elements_inputs_outputs_df, update_cnwheat_parameters)

    # -- FARQUHARWHEAT --
    # Initial states
    farquharwheat_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
        farquharwheat_facade.converter.ELEMENT_TOPOLOGY_COLUMNS +
        [i for i in farquharwheat_facade.converter.FARQUHARWHEAT_ELEMENTS_INPUTS if i in inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns]].copy()

    farquharwheat_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
        farquharwheat_facade.converter.AXIS_TOPOLOGY_COLUMNS +
        [i for i in farquharwheat_facade.converter.FARQUHARWHEAT_AXES_INPUTS if i in inputs_dataframes[AXES_INITIAL_STATE_FILENAME].columns]].copy()

    # Facade initialisation
    farquharwheat_facade_ = farquharwheat_facade.FarquharWheatFacade(g,
                                                                     farquharwheat_elements_initial_state,
                                                                     farquharwheat_axes_initial_state,
                                                                     shared_elements_inputs_outputs_df)

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

    # Update some parameters
    update_cnwheat_parameters = {'VMAX_ROOTS_GROWTH_PREFLO': 0.02885625}

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
                                                               shared_axes_inputs_outputs_df, update_cnwheat_parameters)

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

    # Update some parameters
    update_cnwheat_parameters = {'roots': {'K_AMINO_ACIDS_EXPORT': 3E-5,
                                           'K_NITRATE_EXPORT': 1E-6}}

    # Facade initialisation
    cnwheat_facade_ = cnwheat_facade.CNWheatFacade(g,
                                                   CNWHEAT_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                   PLANT_DENSITY,
                                                   update_cnwheat_parameters,
                                                   cnwheat_organs_initial_state,
                                                   cnwheat_hiddenzones_initial_state,
                                                   cnwheat_elements_initial_state,
                                                   cnwheat_soils_initial_state,
                                                   shared_axes_inputs_outputs_df,
                                                   shared_organs_inputs_outputs_df,
                                                   shared_hiddenzones_inputs_outputs_df,
                                                   shared_elements_inputs_outputs_df,
                                                   shared_soils_inputs_outputs_df)
    # -- FSPMWHEAT --
    # Facade initialisation
    fspmwheat_facade_ = fspmwheat_facade.FSPMWheatFacade(g)

    # Update geometry
    adel_wheat.update_geometry(g)

    # ---------------------------------------------
    # -----      RUN OF THE SIMULATION      -------
    # ---------------------------------------------

    for t_caribu in range(START_TIME, SIMULATION_LENGTH, CARIBU_TIMESTEP):
        # run Caribu
        PARi = meteo.loc[t_caribu, ['PARi_MA4']].iloc[0]
        DOY = meteo.loc[t_caribu, ['DOY']].iloc[0]
        hour = meteo.loc[t_caribu, ['hour']].iloc[0]
        PARi_next_hours = meteo.loc[range(t_caribu, t_caribu + CARIBU_TIMESTEP), ['PARi']].sum().values[0]

        if (t_caribu % CARIBU_TIMESTEP == 0) and (PARi_next_hours > 0):
            run_caribu = True
        else:
            run_caribu = False

        caribu_facade_.run(run_caribu, energy=PARi, DOY=DOY, hourTU=hour, latitude=48.85, sun_sky_option='sky', heterogeneous_canopy=True, plant_density=PLANT_DENSITY[1])

        for t_senescwheat in range(t_caribu, t_caribu + CARIBU_TIMESTEP, SENESCWHEAT_TIMESTEP):
            # run SenescWheat
            senescwheat_facade_.run()

            # Run the rest of the model if the plant is alive
            for t_farquharwheat in range(t_senescwheat, t_senescwheat + SENESCWHEAT_TIMESTEP, FARQUHARWHEAT_TIMESTEP):
                # get the meteo of the current step
                Ta, ambient_CO2, RH, Ur = meteo.loc[t_farquharwheat, ['air_temperature_MA2', 'ambient_CO2_MA2', 'humidity_MA2', 'Wind_MA2']]

                # run FarquharWheat
                farquharwheat_facade_.run(Ta, ambient_CO2, RH, Ur)

                for t_elongwheat in range(t_farquharwheat, t_farquharwheat + FARQUHARWHEAT_TIMESTEP, ELONGWHEAT_TIMESTEP):
                    # run ElongWheat
                    Tair, Tsoil = meteo.loc[t_elongwheat, ['air_temperature', 'soil_temperature']]
                    elongwheat_facade_.run(Tair, Tsoil, option_static=False)

                    # Update geometry
                    adel_wheat.update_geometry(g)

                    for t_growthwheat in range(t_elongwheat, t_elongwheat + ELONGWHEAT_TIMESTEP, GROWTHWHEAT_TIMESTEP):
                        # run GrowthWheat
                        growthwheat_facade_.run()

                        for t_cnwheat in range(t_growthwheat, t_growthwheat + GROWTHWHEAT_TIMESTEP, CNWHEAT_TIMESTEP):
                            if t_cnwheat > 0:

                                # run CNWheat
                                Tair = meteo.loc[t_elongwheat, 'air_temperature']
                                Tsoil = meteo.loc[t_elongwheat, 'soil_temperature']
                                cnwheat_facade_.run(Tair, Tsoil)

                            # append outputs at current step to global lists
                            axes_outputs, elements_outputs, hiddenzones_outputs, organs_outputs, soils_outputs = fspmwheat_facade_.build_outputs_df_from_MTG()

                            all_simulation_steps.append(t_cnwheat)
                            axes_all_data_list.append(axes_outputs)
                            organs_all_data_list.append(organs_outputs)
                            hiddenzones_all_data_list.append(hiddenzones_outputs)
                            elements_all_data_list.append(elements_outputs)
                            soils_all_data_list.append(soils_outputs)

    # compare actual to desired outputs at each scale level (an exception is raised if the test failed)
    for (outputs_df_list,
         desired_outputs_filename,
         actual_outputs_filename,
         index_columns,
         state_variables_names) \
            in ((axes_all_data_list, DESIRED_AXES_OUTPUTS_FILENAME,
                 ACTUAL_AXES_OUTPUTS_FILENAME, cnwheat_simulation.Simulation.AXES_T_INDEXES, cnwheat_simulation.Simulation.AXES_STATE),
                (organs_all_data_list, DESIRED_ORGANS_OUTPUTS_FILENAME, ACTUAL_ORGANS_OUTPUTS_FILENAME,
                 cnwheat_simulation.Simulation.ORGANS_T_INDEXES, cnwheat_simulation.Simulation.ORGANS_STATE),
                (hiddenzones_all_data_list, DESIRED_HIDDENZONES_OUTPUTS_FILENAME, ACTUAL_HIDDENZONES_OUTPUTS_FILENAME,
                 cnwheat_simulation.Simulation.HIDDENZONE_T_INDEXES, cnwheat_simulation.Simulation.HIDDENZONE_STATE),
                (elements_all_data_list, DESIRED_ELEMENTS_OUTPUTS_FILENAME, ACTUAL_ELEMENTS_OUTPUTS_FILENAME,
                 cnwheat_simulation.Simulation.ELEMENTS_T_INDEXES, cnwheat_simulation.Simulation.ELEMENTS_STATE),
                (soils_all_data_list, DESIRED_SOILS_OUTPUTS_FILENAME, ACTUAL_SOILS_OUTPUTS_FILENAME,
                 cnwheat_simulation.Simulation.SOILS_T_INDEXES, cnwheat_simulation.Simulation.SOILS_STATE)):
        outputs_df = pd.concat(outputs_df_list, keys=all_simulation_steps, sort=False)
        outputs_df.reset_index(0, inplace=True)
        outputs_df.rename({'level_0': 't'}, axis=1, inplace=True)
        outputs_df = outputs_df.reindex(index_columns + outputs_df.columns.difference(index_columns).tolist(), axis=1, copy=False)
        outputs_df = outputs_df.loc[:, state_variables_names]  # compare only the values of the compartments
        cnwheat_tools.compare_actual_to_desired(OUTPUTS_DIRPATH, outputs_df, desired_outputs_filename,
                                                actual_outputs_filename, precision=PRECISION, overwrite_desired_data=overwrite_desired_data)


if __name__ == '__main__':
    test_run(overwrite_desired_data=False)
