# -*- coding: latin-1 -*-

import os

import numpy as np
import pandas as pd

from openalea.farquharwheat import simulation, converter

"""
    test_farquhar_wheat
    ~~~~~~~~~~~~~~~~~~~

    Test the Farquhar-Wheat model (standalone).

    You must first install :mod:`farquharwheat` and its dependencies
    before running this script with the command `python`.

"""

# the file names of the inputs
INPUTS_ELEMENT_FILENAME = 'elements_inputs.csv'
INPUTS_AXIS_FILENAME = 'axes_inputs.csv'

# desired outputs filenames
DESIRED_OUTPUTS_FILENAME = 'desired_outputs.csv'

# actual outputs filenames
ACTUAL_OUTPUTS_FILENAME = 'actual_outputs.csv'

PRECISION = 6
RELATIVE_TOLERANCE = 10**-PRECISION
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
        for column in ('axis', 'organ', 'element', 'label'):
            if column in desired_data_df.columns:
                del desired_data_df[column]
                del actual_data_df[column]

        # compare to the desired data
        np.testing.assert_allclose(actual_data_df.values, desired_data_df.values, RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE)


def test_run(overwrite_desired_data=False):

    # create a simulation and a converter
    simulation_ = simulation.Simulation()
    # read inputs from Pandas dataframe
    elements_inputs_df = pd.read_csv(INPUTS_ELEMENT_FILENAME)
    axes_inputs_df = pd.read_csv(INPUTS_AXIS_FILENAME)
    # convert the dataframe to simulation inputs format
    inputs = converter.from_dataframe(elements_inputs_df, axes_inputs_df)
    # initialize the simulation with the inputs
    simulation_.initialize(inputs)
    # run the simulation
    simulation_.run(Ta=18.8, ambient_CO2=360, RH=0.530000, Ur=2.200000)
    # convert the outputs to Pandas dataframe
    outputs_df = converter.to_dataframe(simulation_.outputs)
    # compare outputs
    compare_actual_to_desired('.', outputs_df, DESIRED_OUTPUTS_FILENAME, ACTUAL_OUTPUTS_FILENAME, overwrite_desired_data)


if __name__ == '__main__':
    test_run()
