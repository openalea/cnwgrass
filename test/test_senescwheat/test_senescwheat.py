# -*- coding: latin-1 -*-
import os

import numpy as np
import pandas as pd

from openalea.senescwheat import simulation, converter

"""
    test_senescwheat
    ~~~~~~~~~~~~~~~

    Test the model Senesc-Wheat.

    You must first install :mod:`senescwheat` (and add it to your PYTHONPATH)
    before running this script with the command `python`.

"""

# inputs directory path
INPUTS_DIRPATH = 'inputs'

# the file names of the inputs
ROOTS_INPUTS_FILENAME = 'roots_inputs.csv'
ELEMENTS_INPUTS_FILENAME = 'elements_inputs.csv'
AXES_INPUTS_FILENAME = 'axes_inputs.csv'

# outputs directory path
OUTPUTS_DIRPATH = 'outputs'

# desired outputs filenames
DESIRED_ROOTS_OUTPUTS_FILENAME = 'desired_roots_outputs.csv'
DESIRED_ELEMENTS_OUTPUTS_FILENAME = 'desired_elements_outputs.csv'

# actual outputs filenames
ACTUAL_ROOTS_OUTPUTS_FILENAME = 'actual_roots_outputs.csv'
ACTUAL_ELEMENTS_OUTPUTS_FILENAME = 'actual_elements_outputs.csv'


PRECISION = 6
RELATIVE_TOLERANCE = 10 ** -PRECISION
ABSOLUTE_TOLERANCE = RELATIVE_TOLERANCE


def compare_actual_to_desired(data_dirpath, actual_data_df, desired_data_filename, actual_data_filename=None, overwrite_desired_data=False):
    # read desired data
    desired_data_filepath = os.path.join(data_dirpath, desired_data_filename)
    desired_data_df = pd.read_csv(desired_data_filepath)

    if actual_data_filename is not None:
        actual_data_filepath = os.path.join(data_dirpath, actual_data_filename)
        actual_data_df.to_csv(actual_data_filepath, na_rep='NA', index=False)

    if overwrite_desired_data:
        desired_data_filepath = os.path.join(data_dirpath, desired_data_filename)
        actual_data_df.to_csv(desired_data_filepath, na_rep='NA', index=False)
    else:
        # keep only numerical data
        for column in ('axis', 'organ', 'element', 'is_over', 'is_growing'):
            if column in desired_data_df.columns:
                del desired_data_df[column]
                del actual_data_df[column]

        # compare to the desired data
        np.testing.assert_allclose(actual_data_df.values, desired_data_df.values, RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE)


def test_run(overwrite_desired_data=False):
    # create a simulation
    simulation_ = simulation.Simulation(delta_t=3600)

    # read inputs from Pandas dataframe
    roots_inputs_df = pd.read_csv(os.path.join(INPUTS_DIRPATH, ROOTS_INPUTS_FILENAME))
    elements_inputs_df = pd.read_csv(os.path.join(INPUTS_DIRPATH, ELEMENTS_INPUTS_FILENAME))
    axes_inputs_df = pd.read_csv(os.path.join(INPUTS_DIRPATH, AXES_INPUTS_FILENAME))

    # convert the dataframe to simulation inputs format
    inputs = converter.from_dataframes(roots_inputs_df, axes_inputs_df, elements_inputs_df)

    # initialize the simulation with the inputs
    simulation_.initialize(inputs)

    # create empty lists of dataframes to store the outputs at each step
    roots_outputs_df_list = []
    elements_outputs_df_list = []

    # define the time grid to run the model on
    start_time = 0
    stop_time = 100
    time_step = 1
    time_grid = range(start_time, stop_time + time_step, time_step)

    # run the model on the time grid
    for t in time_grid:
        simulation_.run()

        # convert outputs to dataframes
        roots_outputs_df, _, elements_outputs_df = converter.to_dataframes(simulation_.outputs)

        # append the outputs at current t to the lists of dataframes
        for df, list_ in ((roots_outputs_df, roots_outputs_df_list),
                          (elements_outputs_df, elements_outputs_df_list)):
            df.insert(0, 't', t)
            list_.append(df)

    # compare actual to desired outputs at each scale level (an exception is raised if the test failed)
    for (outputs_df_list,
         desired_outputs_filename,
         actual_outputs_filename) \
            in ((roots_outputs_df_list, DESIRED_ROOTS_OUTPUTS_FILENAME, ACTUAL_ROOTS_OUTPUTS_FILENAME),
                (elements_outputs_df_list, DESIRED_ELEMENTS_OUTPUTS_FILENAME, ACTUAL_ELEMENTS_OUTPUTS_FILENAME)):
        outputs_df = pd.concat(outputs_df_list, ignore_index=True)
        print('Compare {} to {}'.format(actual_outputs_filename, desired_outputs_filename))
        compare_actual_to_desired(OUTPUTS_DIRPATH, outputs_df, desired_outputs_filename, actual_outputs_filename, overwrite_desired_data)
        print('{} OK!'.format(actual_outputs_filename))


if __name__ == '__main__':
    test_run(overwrite_desired_data=False)
