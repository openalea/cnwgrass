# -*- coding: latin-1 -*-

import glob
import os
import logging
import warnings

import pandas as pd

from openalea.cnwheat import simulation as cnwheat_simulation, converter as cnwheat_converter, \
    tools as cnwheat_tools, postprocessing as cnwheat_postprocessing
from openalea.respiwheat import model as respiwheat_model

"""
    test_cnwheat
    ~~~~~~~~~~~~

    Test:

        * the run of a simulation with/without interpolation of the forcings,
        * the logging,
        * the postprocessing,
        * and the graphs generation.

    You must first install model CN-Wheat before running this script with the command `python`. See `README.md` at the
    root directory of the project.

    To get a coverage report of the code of the model, use the command:
    `nosetests --with-coverage --cover-package=cnwheat test_cnwheat.py`.

    CSV files must contain only ASCII characters and ',' as separator.

"""

# Number of seconds in 1 hour
HOUR_TO_SECOND_CONVERSION_FACTOR = 3600

# Precision of floats used to write and format the output CSV files
OUTPUTS_PRECISION = 6

# the precision to use for quantitative comparison test
PRECISION = 4


def force_senescence_and_photosynthesis(t, population, senescence_roots_data_grouped, senescence_elements_data_grouped,
                                        photosynthesis_elements_data_grouped):
    """Force the senescence and photosynthesis data of the population at `t` from input grouped dataframes"""
    for plant in population.plants:
        for axis in plant.axes:
            # Root growth and senescence
            group = senescence_roots_data_grouped.get_group((t, plant.index, axis.label))
            senescence_data_to_use = group.loc[group.first_valid_index(), group.columns.intersection(
                cnwheat_simulation.Simulation.ORGANS_STATE)].dropna().to_dict()
            axis.roots.__dict__.update(senescence_data_to_use)
            for phytomer in axis.phytomers:
                for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):
                    if organ is None:
                        continue
                    for element in (organ.exposed_element, organ.enclosed_element):
                        if element is None:
                            continue
                        # Element senescence
                        group_senesc = senescence_elements_data_grouped.get_group(
                            (t, plant.index, axis.label, phytomer.index, organ.label, element.label))
                        senescence_data_to_use = group_senesc.loc[
                            group_senesc.first_valid_index(), group_senesc.columns.intersection(
                                cnwheat_simulation.Simulation.ELEMENTS_STATE)].dropna().to_dict()
                        element.__dict__.update(senescence_data_to_use)
                        # Element photosynthesis
                        group_photo = photosynthesis_elements_data_grouped.get_group(
                            (t, plant.index, axis.label, phytomer.index, organ.label, element.label))
                        photosynthesis_elements_data_to_use = group_photo.loc[
                            group_photo.first_valid_index(), group_photo.columns.intersection(
                                cnwheat_simulation.Simulation.ELEMENTS_STATE)].dropna().to_dict()
                        element.__dict__.update(photosynthesis_elements_data_to_use)


def test_simulation_run(overwrite_desired_data=False):
    """Test the run of a simulation, without interpolation of the forcings."""

    TEST_DIR_PATH = 'simulation_run'

    # Inputs of the test
    INPUTS_DIRPATH = os.path.join(TEST_DIR_PATH, 'inputs')
    ORGANS_INITIAL_STATE_FILENAME = 'organs_initial_state.csv'
    HIDDENZONES_INITIAL_STATE_FILENAME = 'hiddenzones_initial_state.csv'
    ELEMENTS_INITIAL_STATE_FILENAME = 'elements_initial_state.csv'
    SOILS_INITIAL_STATE_FILENAME = 'soils_initial_state.csv'
    ELEMENTS_PHOTOSYNTHESIS_FORCINGS_FILENAME = 'elements_photosynthesis_forcings.csv'
    ROOTS_SENESCENCE_FORCINGS_FILENAME = 'roots_senescence_forcings.csv'
    ELEMENTS_SENESCENCE_FORCINGS_FILENAME = 'elements_senescence_forcings.csv'

    # Outputs of the test
    OUTPUTS_DIRPATH = os.path.join(TEST_DIR_PATH, 'outputs')
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

    # Simulation parameters
    START_TIME = 0
    SIMULATION_LENGTH = 48
    TIME_STEP = 1
    CULM_DENSITY = {1: 410}

    time_step_seconds = TIME_STEP * HOUR_TO_SECOND_CONVERSION_FACTOR

    # Read the inputs from CSV files and create inputs dataframes
    inputs_dataframes = {}
    for inputs_filename in (
            ORGANS_INITIAL_STATE_FILENAME, HIDDENZONES_INITIAL_STATE_FILENAME, ELEMENTS_INITIAL_STATE_FILENAME,
            SOILS_INITIAL_STATE_FILENAME):
        inputs_dataframes[inputs_filename] = pd.read_csv(os.path.join(INPUTS_DIRPATH, inputs_filename))

    # Convert the inputs dataframes to a population of plants and a dictionary of soils
    population, soils = cnwheat_converter.from_dataframes(inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME],
                                                          inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME],
                                                          inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME],
                                                          inputs_dataframes[SOILS_INITIAL_STATE_FILENAME])

    # Create the simulation
    simulation_ = cnwheat_simulation.Simulation(respiration_model=respiwheat_model, delta_t=time_step_seconds, culm_density=CULM_DENSITY)

    # Initialize the simulation from the population of plants and the dictionary of soils created previously
    simulation_.initialize(population, soils)

    # Read photosynthesis and senescence forcings from CSV files, create dataframes, and group the dataframes by object index
    photosynthesis_elements_data_filepath = os.path.join(INPUTS_DIRPATH, ELEMENTS_PHOTOSYNTHESIS_FORCINGS_FILENAME)
    photosynthesis_elements_data_df = pd.read_csv(photosynthesis_elements_data_filepath)
    photosynthesis_elements_data_grouped = photosynthesis_elements_data_df.groupby(cnwheat_simulation.Simulation.ELEMENTS_T_INDEXES)
    senescence_roots_data_filepath = os.path.join(INPUTS_DIRPATH, ROOTS_SENESCENCE_FORCINGS_FILENAME)
    senescence_roots_data_df = pd.read_csv(senescence_roots_data_filepath)
    senescence_roots_data_grouped = senescence_roots_data_df.groupby(cnwheat_simulation.Simulation.AXES_T_INDEXES)
    senescence_elements_data_filepath = os.path.join(INPUTS_DIRPATH, ELEMENTS_SENESCENCE_FORCINGS_FILENAME)
    senescence_elements_data_df = pd.read_csv(senescence_elements_data_filepath)
    senescence_elements_data_grouped = senescence_elements_data_df.groupby(cnwheat_simulation.Simulation.ELEMENTS_T_INDEXES)

    # Force the senescence and photosynthesis of the population
    force_senescence_and_photosynthesis(0, population, senescence_roots_data_grouped, senescence_elements_data_grouped, photosynthesis_elements_data_grouped)

    # Reinitialize the simulation from forced population and soils
    simulation_.initialize(population, soils)

    # Define the time grid of the simulation
    time_grid = range(START_TIME, SIMULATION_LENGTH + TIME_STEP, TIME_STEP)

    # Create empty lists of dataframes to store the outputs at each step of the simulation
    axes_outputs_df_list = []
    organs_outputs_df_list = []
    hiddenzones_outputs_df_list = []
    elements_outputs_df_list = []
    soils_outputs_df_list = []

    for t in time_grid:

        if t > 0:
            # Run the model of CN exchanges ; the population is internally updated by the model
            simulation_.run()

        # Convert the model outputs to dataframes
        _, axes_outputs_df, _, organs_outputs_df, hiddenzones_outputs_df, elements_outputs_df, soils_outputs_df = cnwheat_converter.to_dataframes(simulation_.population, simulation_.soils)

        # Append the outputs dataframes at current t to the global lists of dataframes
        for df, list_ in ((axes_outputs_df, axes_outputs_df_list), (organs_outputs_df, organs_outputs_df_list),
                          (hiddenzones_outputs_df, hiddenzones_outputs_df_list), (elements_outputs_df, elements_outputs_df_list),
                          (soils_outputs_df, soils_outputs_df_list)):
            df.insert(0, 't', t)
            list_.append(df)

        if 0 < t < SIMULATION_LENGTH:
            # Force the senescence and photosynthesis of the population
            force_senescence_and_photosynthesis(t, population, senescence_roots_data_grouped, senescence_elements_data_grouped, photosynthesis_elements_data_grouped)
            # Reinitialize the simulation from forced population and soils
            simulation_.initialize(population, soils)

    # compare actual to desired outputs at each scale level (an exception is raised if the test failed)
    for (outputs_df_list,
         desired_outputs_filename,
         actual_outputs_filename,
         state_variables_names) \
            in ((axes_outputs_df_list, DESIRED_AXES_OUTPUTS_FILENAME, ACTUAL_AXES_OUTPUTS_FILENAME,
                 cnwheat_simulation.Simulation.AXES_T_INDEXES + cnwheat_simulation.Simulation.AXES_STATE),
                (organs_outputs_df_list, DESIRED_ORGANS_OUTPUTS_FILENAME, ACTUAL_ORGANS_OUTPUTS_FILENAME,
                 cnwheat_simulation.Simulation.ORGANS_T_INDEXES + cnwheat_simulation.Simulation.ORGANS_STATE),
                (hiddenzones_outputs_df_list, DESIRED_HIDDENZONES_OUTPUTS_FILENAME,
                 ACTUAL_HIDDENZONES_OUTPUTS_FILENAME,
                 cnwheat_simulation.Simulation.HIDDENZONE_T_INDEXES + cnwheat_simulation.Simulation.HIDDENZONE_STATE),
                (elements_outputs_df_list, DESIRED_ELEMENTS_OUTPUTS_FILENAME, ACTUAL_ELEMENTS_OUTPUTS_FILENAME,
                 cnwheat_simulation.Simulation.ELEMENTS_T_INDEXES + cnwheat_simulation.Simulation.ELEMENTS_STATE),
                (soils_outputs_df_list, DESIRED_SOILS_OUTPUTS_FILENAME, ACTUAL_SOILS_OUTPUTS_FILENAME,
                 cnwheat_simulation.Simulation.SOILS_T_INDEXES + cnwheat_simulation.Simulation.SOILS_STATE)):
        outputs_df = pd.concat(outputs_df_list, ignore_index=True)
        outputs_df = outputs_df.loc[:, state_variables_names]  # compare only the values of the compartments
        cnwheat_tools.compare_actual_to_desired(OUTPUTS_DIRPATH, outputs_df, desired_outputs_filename,
                                                actual_outputs_filename, precision=PRECISION, overwrite_desired_data=overwrite_desired_data)


def test_simulation_run_with_interpolation(overwrite_desired_data=False):
    """Test the run of a simulation, with interpolation of the forcings."""

    TEST_DIR_PATH = 'simulation_run_with_interpolation'

    # Inputs of the test
    INPUTS_DIRPATH = os.path.join(TEST_DIR_PATH, 'inputs')
    ORGANS_INITIAL_STATE_FILENAME = 'organs_initial_state.csv'
    HIDDENZONES_INITIAL_STATE_FILENAME = 'hiddenzones_initial_state.csv'
    ELEMENTS_INITIAL_STATE_FILENAME = 'elements_initial_state.csv'
    SOILS_INITIAL_STATE_FILENAME = 'soils_initial_state.csv'
    ELEMENTS_PHOTOSYNTHESIS_FORCINGS_FILENAME = 'elements_photosynthesis_forcings.csv'
    ROOTS_SENESCENCE_FORCINGS_FILENAME = 'roots_senescence_forcings.csv'
    ELEMENTS_SENESCENCE_FORCINGS_FILENAME = 'elements_senescence_forcings.csv'

    # Outputs of the test
    OUTPUTS_DIRPATH = os.path.join(TEST_DIR_PATH, 'outputs')
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

    # Simulation parameters
    START_TIME = 0
    SIMULATION_LENGTH = 5
    TIME_STEP = 1
    CULM_DENSITY = {1: 410}

    time_step_seconds = TIME_STEP * HOUR_TO_SECOND_CONVERSION_FACTOR

    # Read the inputs from CSV files and create inputs dataframes
    inputs_dataframes = {}
    for inputs_filename in (
            ORGANS_INITIAL_STATE_FILENAME, HIDDENZONES_INITIAL_STATE_FILENAME, ELEMENTS_INITIAL_STATE_FILENAME,
            SOILS_INITIAL_STATE_FILENAME):
        inputs_dataframes[inputs_filename] = pd.read_csv(os.path.join(INPUTS_DIRPATH, inputs_filename))

    # Convert the inputs dataframes to a population of plants and a dictionary of soils
    population, soils = cnwheat_converter.from_dataframes(inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME],
                                                          inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME],
                                                          inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME],
                                                          inputs_dataframes[SOILS_INITIAL_STATE_FILENAME])

    # Create the simulation
    simulation_ = cnwheat_simulation.Simulation(respiration_model=respiwheat_model, delta_t=time_step_seconds, culm_density=CULM_DENSITY,
                                                interpolate_forcings=True,
                                                senescence_forcings_delta_t=time_step_seconds,
                                                photosynthesis_forcings_delta_t=time_step_seconds)

    # Initialize the simulation from the population of plants and the dictionary of soils created previously
    simulation_.initialize(population, soils)

    # Read photosynthesis and senescence forcings from CSV files, create dataframes, and group the dataframes by object index
    photosynthesis_elements_data_filepath = os.path.join(INPUTS_DIRPATH, ELEMENTS_PHOTOSYNTHESIS_FORCINGS_FILENAME)
    photosynthesis_elements_data_df = pd.read_csv(photosynthesis_elements_data_filepath)
    photosynthesis_elements_data_grouped = photosynthesis_elements_data_df.groupby(cnwheat_simulation.Simulation.ELEMENTS_T_INDEXES)
    senescence_roots_data_filepath = os.path.join(INPUTS_DIRPATH, ROOTS_SENESCENCE_FORCINGS_FILENAME)
    senescence_roots_data_df = pd.read_csv(senescence_roots_data_filepath)
    senescence_roots_data_grouped = senescence_roots_data_df.groupby(cnwheat_simulation.Simulation.AXES_T_INDEXES)
    senescence_elements_data_filepath = os.path.join(INPUTS_DIRPATH, ELEMENTS_SENESCENCE_FORCINGS_FILENAME)
    senescence_elements_data_df = pd.read_csv(senescence_elements_data_filepath)
    senescence_elements_data_grouped = senescence_elements_data_df.groupby(cnwheat_simulation.Simulation.ELEMENTS_T_INDEXES)

    # Force the senescence and photosynthesis of the population
    force_senescence_and_photosynthesis(0, population, senescence_roots_data_grouped, senescence_elements_data_grouped, photosynthesis_elements_data_grouped)

    # Reinitialize the simulation from forced population and soils
    simulation_.initialize(population, soils)

    # Define the time grid of the simulation
    time_grid = range(START_TIME, SIMULATION_LENGTH + TIME_STEP, TIME_STEP)

    # Create empty lists of dataframes to store the outputs at each step of the simulation
    axes_outputs_df_list = []
    organs_outputs_df_list = []
    hiddenzones_outputs_df_list = []
    elements_outputs_df_list = []
    soils_outputs_df_list = []

    for t in time_grid:

        if t > 0:
            # Run the model of CN exchanges ; the population is internally updated by the model
            simulation_.run()

        # Convert the model outputs to dataframes
        _, axes_outputs_df, _, organs_outputs_df, hiddenzones_outputs_df, elements_outputs_df, soils_outputs_df = cnwheat_converter.to_dataframes(simulation_.population, simulation_.soils)

        # Append the outputs dataframes at current t to the global lists of dataframes
        for df, list_ in ((axes_outputs_df, axes_outputs_df_list), (organs_outputs_df, organs_outputs_df_list),
                          (hiddenzones_outputs_df, hiddenzones_outputs_df_list), (elements_outputs_df, elements_outputs_df_list),
                          (soils_outputs_df, soils_outputs_df_list)):
            df.insert(0, 't', t)
            list_.append(df)

        if 0 < t < SIMULATION_LENGTH:
            # Force the senescence and photosynthesis of the population
            force_senescence_and_photosynthesis(t, population, senescence_roots_data_grouped, senescence_elements_data_grouped, photosynthesis_elements_data_grouped)
            # Reinitialize the simulation from forced population and soils
            simulation_.initialize(population, soils)

    # compare actual to desired outputs at each scale level (an exception is raised if the test failed)

    for (outputs_df_list,
         desired_outputs_filename,
         actual_outputs_filename,
         state_variables_names) \
            in ((axes_outputs_df_list, DESIRED_AXES_OUTPUTS_FILENAME, ACTUAL_AXES_OUTPUTS_FILENAME,
                 cnwheat_simulation.Simulation.AXES_T_INDEXES + cnwheat_simulation.Simulation.AXES_STATE),
                (organs_outputs_df_list, DESIRED_ORGANS_OUTPUTS_FILENAME, ACTUAL_ORGANS_OUTPUTS_FILENAME,
                 cnwheat_simulation.Simulation.ORGANS_T_INDEXES + cnwheat_simulation.Simulation.ORGANS_STATE),
                (hiddenzones_outputs_df_list, DESIRED_HIDDENZONES_OUTPUTS_FILENAME,
                 ACTUAL_HIDDENZONES_OUTPUTS_FILENAME,
                 cnwheat_simulation.Simulation.HIDDENZONE_T_INDEXES + cnwheat_simulation.Simulation.HIDDENZONE_STATE),
                (elements_outputs_df_list, DESIRED_ELEMENTS_OUTPUTS_FILENAME, ACTUAL_ELEMENTS_OUTPUTS_FILENAME,
                 cnwheat_simulation.Simulation.ELEMENTS_T_INDEXES + cnwheat_simulation.Simulation.ELEMENTS_STATE),
                (soils_outputs_df_list, DESIRED_SOILS_OUTPUTS_FILENAME, ACTUAL_SOILS_OUTPUTS_FILENAME,
                 cnwheat_simulation.Simulation.SOILS_T_INDEXES + cnwheat_simulation.Simulation.SOILS_STATE)):
        outputs_df = pd.concat(outputs_df_list, ignore_index=True)
        outputs_df = outputs_df.loc[:, state_variables_names]  # compare only the values of the compartments
        cnwheat_tools.compare_actual_to_desired(OUTPUTS_DIRPATH, outputs_df, desired_outputs_filename,
                                                actual_outputs_filename, precision=PRECISION, overwrite_desired_data=overwrite_desired_data)


def test_simulation_logging(overwrite_desired_data=False):
    """Test the logging of a simulation."""

    TEST_DIR_PATH = 'simulation_logging'

    # Inputs of the test
    INPUTS_DIRPATH = os.path.join(TEST_DIR_PATH, 'inputs')
    ORGANS_INITIAL_STATE_FILENAME = 'organs_initial_state.csv'
    HIDDENZONES_INITIAL_STATE_FILENAME = 'hiddenzones_initial_state.csv'
    ELEMENTS_INITIAL_STATE_FILENAME = 'elements_initial_state.csv'
    SOILS_INITIAL_STATE_FILENAME = 'soils_initial_state.csv'
    ELEMENTS_PHOTOSYNTHESIS_FORCINGS_FILENAME = 'elements_photosynthesis_forcings.csv'
    ROOTS_SENESCENCE_FORCINGS_FILENAME = 'roots_senescence_forcings.csv'
    ELEMENTS_SENESCENCE_FORCINGS_FILENAME = 'elements_senescence_forcings.csv'

    # Outputs of the test
    LOGS_DIRPATH = os.path.join(TEST_DIR_PATH, 'logs')
    DESIRED_AXES_COMPARTMENTS_FILENAME = 'desired_axes_compartments.csv'
    DESIRED_ORGANS_COMPARTMENTS_FILENAME = 'desired_organs_compartments.csv'
    DESIRED_HIDDENZONES_COMPARTMENTS_FILENAME = 'desired_hiddenzones_compartments.csv'
    DESIRED_ELEMENTS_COMPARTMENTS_FILENAME = 'desired_elements_compartments.csv'
    DESIRED_SOILS_COMPARTMENTS_FILENAME = 'desired_soils_compartments.csv'
    DESIRED_AXES_DERIVATIVES_FILENAME = 'desired_axes_derivatives.csv'
    DESIRED_ORGANS_DERIVATIVES_FILENAME = 'desired_organs_derivatives.csv'
    DESIRED_HIDDENZONES_DERIVATIVES_FILENAME = 'desired_hiddenzones_derivatives.csv'
    DESIRED_ELEMENTS_DERIVATIVES_FILENAME = 'desired_elements_derivatives.csv'
    DESIRED_SOILS_DERIVATIVES_FILENAME = 'desired_soils_derivatives.csv'
    ACTUAL_AXES_COMPARTMENTS_FILENAME = 'actual_axes_compartments.csv'
    ACTUAL_ORGANS_COMPARTMENTS_FILENAME = 'actual_organs_compartments.csv'
    ACTUAL_HIDDENZONES_COMPARTMENTS_FILENAME = 'actual_hiddenzones_compartments.csv'
    ACTUAL_ELEMENTS_COMPARTMENTS_FILENAME = 'actual_elements_compartments.csv'
    ACTUAL_SOILS_COMPARTMENTS_FILENAME = 'actual_soils_compartments.csv'
    ACTUAL_AXES_DERIVATIVES_FILENAME = 'actual_axes_derivatives.csv'
    ACTUAL_ORGANS_DERIVATIVES_FILENAME = 'actual_organs_derivatives.csv'
    ACTUAL_HIDDENZONES_DERIVATIVES_FILENAME = 'actual_hiddenzones_derivatives.csv'
    ACTUAL_ELEMENTS_DERIVATIVES_FILENAME = 'actual_elements_derivatives.csv'
    ACTUAL_SOILS_DERIVATIVES_FILENAME = 'actual_soils_derivatives.csv'

    # Config file path for logging
    LOGGING_CONFIG_FILEPATH = os.path.join(TEST_DIR_PATH, 'logging.json')

    # Simulation parameters
    START_TIME = 0
    SIMULATION_LENGTH = 5
    TIME_STEP = 1
    CULM_DENSITY = {1: 410}

    time_step_seconds = TIME_STEP * HOUR_TO_SECOND_CONVERSION_FACTOR

    # Remove actual logs files
    for logs_file in glob.glob(os.path.join(LOGS_DIRPATH, "actual*.csv")):
        os.remove(logs_file)

    # Setup the logging (without removing the desired logs since we need them for the comparison test)
    cnwheat_tools.setup_logging(config_filepath=LOGGING_CONFIG_FILEPATH, level=logging.DEBUG,
                                log_model=True, log_compartments=True, log_derivatives=True, remove_old_logs=False)

    # Read the inputs from CSV files and create inputs dataframes
    inputs_dataframes = {}
    for inputs_filename in (
            ORGANS_INITIAL_STATE_FILENAME, HIDDENZONES_INITIAL_STATE_FILENAME, ELEMENTS_INITIAL_STATE_FILENAME,
            SOILS_INITIAL_STATE_FILENAME):
        inputs_dataframes[inputs_filename] = pd.read_csv(os.path.join(INPUTS_DIRPATH, inputs_filename))

    # Convert the inputs dataframes to a population of plants and a dictionary of soils
    population, soils = cnwheat_converter.from_dataframes(inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME],
                                                          inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME],
                                                          inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME],
                                                          inputs_dataframes[SOILS_INITIAL_STATE_FILENAME])

    # Create the simulation
    simulation_ = cnwheat_simulation.Simulation(respiration_model=respiwheat_model, delta_t=time_step_seconds, culm_density=CULM_DENSITY)

    # Initialize the simulation from the population of plants and the dictionary of soils created previously
    simulation_.initialize(population, soils)

    # Read photosynthesis and senescence forcings from CSV files, create dataframes, and group the dataframes by object index
    photosynthesis_elements_data_filepath = os.path.join(INPUTS_DIRPATH, ELEMENTS_PHOTOSYNTHESIS_FORCINGS_FILENAME)
    photosynthesis_elements_data_df = pd.read_csv(photosynthesis_elements_data_filepath)
    photosynthesis_elements_data_grouped = photosynthesis_elements_data_df.groupby(cnwheat_simulation.Simulation.ELEMENTS_T_INDEXES)
    senescence_roots_data_filepath = os.path.join(INPUTS_DIRPATH, ROOTS_SENESCENCE_FORCINGS_FILENAME)
    senescence_roots_data_df = pd.read_csv(senescence_roots_data_filepath)
    senescence_roots_data_grouped = senescence_roots_data_df.groupby(cnwheat_simulation.Simulation.AXES_T_INDEXES)
    senescence_elements_data_filepath = os.path.join(INPUTS_DIRPATH, ELEMENTS_SENESCENCE_FORCINGS_FILENAME)
    senescence_elements_data_df = pd.read_csv(senescence_elements_data_filepath)
    senescence_elements_data_grouped = senescence_elements_data_df.groupby(cnwheat_simulation.Simulation.ELEMENTS_T_INDEXES)

    # Force the senescence and photosynthesis of the population
    force_senescence_and_photosynthesis(0, population, senescence_roots_data_grouped, senescence_elements_data_grouped, photosynthesis_elements_data_grouped)

    # Reinitialize the simulation from forced population and soils
    simulation_.initialize(population, soils)

    # Define the time grid of the simulation
    time_grid = range(START_TIME, SIMULATION_LENGTH + TIME_STEP, TIME_STEP)

    # Create empty lists of dataframes to store the outputs at each step of the simulation
    axes_outputs_df_list = []
    organs_outputs_df_list = []
    hiddenzones_outputs_df_list = []
    elements_outputs_df_list = []
    soils_outputs_df_list = []

    for t in time_grid:

        if t > 0:
            # Run the model of CN exchanges ; the population is internally updated by the model
            simulation_.run()

        # Convert the model outputs to dataframes
        _, axes_outputs_df, _, organs_outputs_df, hiddenzones_outputs_df, elements_outputs_df, soils_outputs_df = cnwheat_converter.to_dataframes(simulation_.population, simulation_.soils)

        # Append the outputs dataframes at current t to the global lists of dataframes
        for df, list_ in ((axes_outputs_df, axes_outputs_df_list), (organs_outputs_df, organs_outputs_df_list),
                          (hiddenzones_outputs_df, hiddenzones_outputs_df_list), (elements_outputs_df, elements_outputs_df_list),
                          (soils_outputs_df, soils_outputs_df_list)):
            df.insert(0, 't', t)
            list_.append(df)

        if 0 < t < SIMULATION_LENGTH:
            # Force the senescence and photosynthesis of the population
            force_senescence_and_photosynthesis(t, population, senescence_roots_data_grouped, senescence_elements_data_grouped, photosynthesis_elements_data_grouped)
            # Reinitialize the simulation from forced population and soils
            simulation_.initialize(population, soils)

    # compare actual to desired logs at each scale level (an exception is raised if the test failed)
    # Compartments logs
    for (desired_compartments_filename,
         actual_compartments_filename) \
            in ((DESIRED_AXES_COMPARTMENTS_FILENAME, ACTUAL_AXES_COMPARTMENTS_FILENAME),
                (DESIRED_ORGANS_COMPARTMENTS_FILENAME, ACTUAL_ORGANS_COMPARTMENTS_FILENAME),
                (DESIRED_HIDDENZONES_COMPARTMENTS_FILENAME, ACTUAL_HIDDENZONES_COMPARTMENTS_FILENAME),
                (DESIRED_ELEMENTS_COMPARTMENTS_FILENAME, ACTUAL_ELEMENTS_COMPARTMENTS_FILENAME),
                (DESIRED_SOILS_COMPARTMENTS_FILENAME, ACTUAL_SOILS_COMPARTMENTS_FILENAME)):
        try:
            actual_compartments_df = pd.read_csv(os.path.join(LOGS_DIRPATH, actual_compartments_filename))
        except pd.errors.EmptyDataError:
            continue  # This file is empty: ignore it.
        cnwheat_tools.compare_actual_to_desired(LOGS_DIRPATH, actual_compartments_df, desired_compartments_filename, precision=PRECISION, overwrite_desired_data=overwrite_desired_data)
    # Derivatives logs
    for (desired_derivatives_filename,
         actual_derivatives_filename) \
            in ((DESIRED_AXES_DERIVATIVES_FILENAME, ACTUAL_AXES_DERIVATIVES_FILENAME),
                (DESIRED_ORGANS_DERIVATIVES_FILENAME, ACTUAL_ORGANS_DERIVATIVES_FILENAME),
                (DESIRED_HIDDENZONES_DERIVATIVES_FILENAME, ACTUAL_HIDDENZONES_DERIVATIVES_FILENAME),
                (DESIRED_ELEMENTS_DERIVATIVES_FILENAME, ACTUAL_ELEMENTS_DERIVATIVES_FILENAME),
                (DESIRED_SOILS_DERIVATIVES_FILENAME, ACTUAL_SOILS_DERIVATIVES_FILENAME)):
        try:
            actual_derivatives_df = pd.read_csv(os.path.join(LOGS_DIRPATH, actual_derivatives_filename))
        except pd.errors.EmptyDataError:
            continue  # This file is empty: ignore it.
        cnwheat_tools.compare_actual_to_desired(LOGS_DIRPATH, actual_derivatives_df, desired_derivatives_filename,
                                                precision=PRECISION, overwrite_desired_data=overwrite_desired_data)


def test_postprocessing(overwrite_desired_data=False):
    """Test the postprocessing."""

    TEST_DIR_PATH = 'postprocessing'

    # Inputs of the test
    OUTPUTS_DIRPATH = os.path.join(TEST_DIR_PATH, 'outputs')
    AXES_OUTPUTS_FILENAME = 'axes_outputs.csv'
    ORGANS_OUTPUTS_FILENAME = 'organs_outputs.csv'
    HIDDENZONES_OUTPUTS_FILENAME = 'hiddenzones_outputs.csv'
    ELEMENTS_OUTPUTS_FILENAME = 'elements_outputs.csv'
    SOILS_OUTPUTS_FILENAME = 'soils_outputs.csv'

    # Outputs of the test
    POSTPROCESSING_DIRPATH = os.path.join(TEST_DIR_PATH, 'postprocessing')
    DESIRED_AXES_POSTPROCESSING_FILENAME = 'desired_axes_postprocessing.csv'
    DESIRED_ORGANS_POSTPROCESSING_FILENAME = 'desired_organs_postprocessing.csv'
    DESIRED_HIDDENZONES_POSTPROCESSING_FILENAME = 'desired_hiddenzones_postprocessing.csv'
    DESIRED_ELEMENTS_POSTPROCESSING_FILENAME = 'desired_elements_postprocessing.csv'
    DESIRED_SOILS_POSTPROCESSING_FILENAME = 'desired_soils_postprocessing.csv'
    ACTUAL_AXES_POSTPROCESSING_FILENAME = 'actual_axes_postprocessing.csv'
    ACTUAL_ORGANS_POSTPROCESSING_FILENAME = 'actual_organs_postprocessing.csv'
    ACTUAL_HIDDENZONES_POSTPROCESSING_FILENAME = 'actual_hiddenzones_postprocessing.csv'
    ACTUAL_ELEMENTS_POSTPROCESSING_FILENAME = 'actual_elements_postprocessing.csv'
    ACTUAL_SOILS_POSTPROCESSING_FILENAME = 'actual_soils_postprocessing.csv'

    # Retrieve outputs dataframes
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

    time_grid = list(outputs_df_dict.values())[0].t
    delta_t = (time_grid.loc[1] - time_grid.loc[0]) * HOUR_TO_SECOND_CONVERSION_FACTOR

    # Compute the post-processing
    axes_postprocessing_file_basename = DESIRED_AXES_POSTPROCESSING_FILENAME.split('.')[0]
    hiddenzones_postprocessing_file_basename = DESIRED_HIDDENZONES_POSTPROCESSING_FILENAME.split('.')[0]
    organs_postprocessing_file_basename = DESIRED_ORGANS_POSTPROCESSING_FILENAME.split('.')[0]
    elements_postprocessing_file_basename = DESIRED_ELEMENTS_POSTPROCESSING_FILENAME.split('.')[0]
    soils_postprocessing_file_basename = DESIRED_SOILS_POSTPROCESSING_FILENAME.split('.')[0]

    postprocessing_df_dict = {}

    try:
        (postprocessing_df_dict[axes_postprocessing_file_basename],
         postprocessing_df_dict[hiddenzones_postprocessing_file_basename],
         postprocessing_df_dict[organs_postprocessing_file_basename],
         postprocessing_df_dict[elements_postprocessing_file_basename],
         postprocessing_df_dict[soils_postprocessing_file_basename]) \
            = cnwheat_postprocessing.postprocessing(axes_df=outputs_df_dict[AXES_OUTPUTS_FILENAME.split('.')[0]],
                                                    hiddenzones_df=outputs_df_dict[HIDDENZONES_OUTPUTS_FILENAME.split('.')[0]],
                                                    organs_df=outputs_df_dict[ORGANS_OUTPUTS_FILENAME.split('.')[0]],
                                                    elements_df=outputs_df_dict[ELEMENTS_OUTPUTS_FILENAME.split('.')[0]],
                                                    soils_df=outputs_df_dict[SOILS_OUTPUTS_FILENAME.split('.')[0]],
                                                    delta_t=delta_t)
    except KeyError as ke:
        warnings.warn(str(ke))
        return

    # Write the postprocessing to CSV files
    for postprocessing_file_basename, postprocessing_filename in ((axes_postprocessing_file_basename, ACTUAL_AXES_POSTPROCESSING_FILENAME),
                                                                  (hiddenzones_postprocessing_file_basename, ACTUAL_HIDDENZONES_POSTPROCESSING_FILENAME),
                                                                  (organs_postprocessing_file_basename, ACTUAL_ORGANS_POSTPROCESSING_FILENAME),
                                                                  (elements_postprocessing_file_basename, ACTUAL_ELEMENTS_POSTPROCESSING_FILENAME),
                                                                  (soils_postprocessing_file_basename, ACTUAL_SOILS_POSTPROCESSING_FILENAME)):
        postprocessing_filepath = os.path.join(POSTPROCESSING_DIRPATH, postprocessing_filename)
        postprocessing_df_dict[postprocessing_file_basename].to_csv(postprocessing_filepath, na_rep='NA', index=False, float_format='%.{}f'.format(OUTPUTS_PRECISION))

    # compare actual to desired postprocessing at each scale level (an exception is raised if the test failed)
    for (postprocessing_file_basename,
         desired_postprocessing_filename,
         actual_postprocessing_filename,
         postprocessing_variables_names) \
            in ((axes_postprocessing_file_basename, DESIRED_AXES_POSTPROCESSING_FILENAME, ACTUAL_AXES_POSTPROCESSING_FILENAME,
                 cnwheat_postprocessing.AXES_T_INDEXES + cnwheat_postprocessing.AXES_POSTPROCESSING_VARIABLES),
                (organs_postprocessing_file_basename, DESIRED_ORGANS_POSTPROCESSING_FILENAME, ACTUAL_ORGANS_POSTPROCESSING_FILENAME,
                 cnwheat_postprocessing.ORGANS_T_INDEXES + cnwheat_postprocessing.ORGANS_POSTPROCESSING_VARIABLES),
                (hiddenzones_postprocessing_file_basename, DESIRED_HIDDENZONES_POSTPROCESSING_FILENAME, ACTUAL_HIDDENZONES_POSTPROCESSING_FILENAME,
                 cnwheat_postprocessing.HIDDENZONE_T_INDEXES + cnwheat_postprocessing.HIDDENZONE_POSTPROCESSING_VARIABLES),
                (elements_postprocessing_file_basename, DESIRED_ELEMENTS_POSTPROCESSING_FILENAME, ACTUAL_ELEMENTS_POSTPROCESSING_FILENAME,
                 cnwheat_postprocessing.ELEMENTS_T_INDEXES + cnwheat_postprocessing.ELEMENTS_POSTPROCESSING_VARIABLES),
                (soils_postprocessing_file_basename, DESIRED_SOILS_POSTPROCESSING_FILENAME, ACTUAL_SOILS_POSTPROCESSING_FILENAME,
                 cnwheat_postprocessing.SOILS_T_INDEXES + cnwheat_postprocessing.SOILS_POSTPROCESSING_VARIABLES)):
        actual_postprocessing_df = postprocessing_df_dict[postprocessing_file_basename]
        actual_postprocessing_df = actual_postprocessing_df.loc[:, actual_postprocessing_df.columns.intersection(postprocessing_variables_names)]  # compare only the postprocessing values
        cnwheat_tools.compare_actual_to_desired(POSTPROCESSING_DIRPATH, actual_postprocessing_df, desired_postprocessing_filename,
                                                actual_postprocessing_filename, precision=PRECISION, overwrite_desired_data=overwrite_desired_data)


def test_graphs_generation():
    """Test the graphs generation."""

    TEST_DIR_PATH = 'graphs_generation'

    # Inputs of the test
    POSTPROCESSING_DIRPATH = os.path.join(TEST_DIR_PATH, 'postprocessing')
    ORGANS_POSTPROCESSING_FILENAME = 'organs_postprocessing.csv'
    HIDDENZONES_POSTPROCESSING_FILENAME = 'hiddenzones_postprocessing.csv'
    ELEMENTS_POSTPROCESSING_FILENAME = 'elements_postprocessing.csv'
    SOILS_POSTPROCESSING_FILENAME = 'soils_postprocessing.csv'

    # Outputs of the test
    GRAPHS_DIRPATH = os.path.join(TEST_DIR_PATH, 'graphs')

    # Retrieve post-processing dataframes
    postprocessing_df_dict = {}

    for postprocessing_filename in (ORGANS_POSTPROCESSING_FILENAME,
                                    HIDDENZONES_POSTPROCESSING_FILENAME,
                                    ELEMENTS_POSTPROCESSING_FILENAME,
                                    SOILS_POSTPROCESSING_FILENAME):
        postprocessing_filepath = os.path.join(POSTPROCESSING_DIRPATH, postprocessing_filename)
        postprocessing_df = pd.read_csv(postprocessing_filepath)
        postprocessing_file_basename = postprocessing_filename.split('.')[0]
        postprocessing_df_dict[postprocessing_file_basename] = postprocessing_df

    try:
        # Generate graphs for validation
        cnwheat_postprocessing.generate_graphs(hiddenzones_df=postprocessing_df_dict[HIDDENZONES_POSTPROCESSING_FILENAME.split('.')[0]],
                                               organs_df=postprocessing_df_dict[ORGANS_POSTPROCESSING_FILENAME.split('.')[0]],
                                               elements_df=postprocessing_df_dict[ELEMENTS_POSTPROCESSING_FILENAME.split('.')[0]],
                                               soils_df=postprocessing_df_dict[SOILS_POSTPROCESSING_FILENAME.split('.')[0]],
                                               graphs_dirpath=GRAPHS_DIRPATH)

    except KeyError as ke:
        warnings.warn(str(ke))
        return


if __name__ == '__main__':
    test_simulation_run(overwrite_desired_data=False)
    print('Simulation Run - OK')

    test_simulation_run_with_interpolation(overwrite_desired_data=False)
    print('Simulation Run with interpolation - OK')

    test_simulation_logging(overwrite_desired_data=False)
    print('Simulation Logging - OK')

    test_postprocessing(overwrite_desired_data=False)
    print('Simulation Postprocessing - OK')

    test_graphs_generation()
    print('Simulation Graphs - OK')
