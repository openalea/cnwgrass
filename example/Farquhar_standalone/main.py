# -*- coding: latin-1 -*-

import pandas as pd

from openalea.gasexchange import simulation, converter

"""
    main
    ~~~~

    An example to show how to initialize and run the model Gas-Exchange in a standalone version
    The example runs Gas-Exchange with two different options : either with or without a coupling with a hydraulic model
"""

INPUTS_ELEMENT_FILENAME = 'elements_inputs.csv'
INPUTS_AXIS_FILENAME = 'axes_inputs.csv'
OUTPUTS_FILENAME = 'outputs.csv'
OUTPUTS_FILENAME_HYDRAULICS = 'outputs_hydraulics.csv'

OUTPUTS_PRECISION = 6

if __name__ == '__main__':

    # read inputs from Pandas dataframe
    elements_inputs_df = pd.read_csv(INPUTS_ELEMENT_FILENAME)
    axes_inputs_df = pd.read_csv(INPUTS_AXIS_FILENAME)

    # convert the dataframe to simulation inputs format
    inputs = converter.from_dataframe(elements_inputs_df, axes_inputs_df)

    #-- No hydraulics coupling

    # create a simulation and a converter
    simulation_ = simulation.Simulation()

    # initialize the simulation with the inputs
    simulation_.initialize(inputs)
    # run the simulation
    simulation_.run(Ta=18.8, ambient_CO2=360, RH=0.530000, Ur=2.200000)
    # convert the outputs to Pandas dataframe
    outputs_df = converter.to_dataframe(simulation_.outputs, simulation_.elements_outputs)
    # write the dataframe to CSV
    outputs_df.to_csv(OUTPUTS_FILENAME, index=False, na_rep='NA', float_format='%.{}f'.format(OUTPUTS_PRECISION))


    #-- Hydraulics coupling
    # create a simulation and a converter
    simulation_ = simulation.Simulation(stomatal_model_name='hydraulics', hydraulics=True)

    # initialize the simulation with the inputs
    simulation_.initialize(inputs)
    # run the simulation
    simulation_.run(Ta=18.8, ambient_CO2=360, RH=0.530000, Ur=2.200000)
    # convert the outputs to Pandas dataframe
    outputs_df = converter.to_dataframe(simulation_.outputs, simulation_.elements_outputs)
    # write the dataframe to CSV
    outputs_df.to_csv(OUTPUTS_FILENAME_HYDRAULICS, index=False, na_rep='NA', float_format='%.{}f'.format(OUTPUTS_PRECISION))