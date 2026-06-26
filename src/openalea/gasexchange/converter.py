# -*- coding: latin-1 -*-

from __future__ import division  # use "//" to do integer division
import pandas as pd

"""
    gasexchange.converter
    ~~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`gasexchange.converter` defines functions to convert
    :class:`dataframes <pandas.DataFrame>` to/from Gas-Exchange inputs or outputs format.

"""

#: the columns which define the topology in the input/output elements dataframe
ELEMENT_TOPOLOGY_COLUMNS = ['plant', 'axis', 'metamer', 'organ', 'element']
#: the columns which define the topology in the input/output elements dataframe
AXIS_TOPOLOGY_COLUMNS = ['plant', 'axis']


def from_dataframe(element_inputs, axes_inputs):
    """
    Convert inputs/outputs from Pandas dataframe to Gas-Exchange format.

    :param pandas.DataFrame element_inputs: Emerging and mature element inputs dataframe to convert, with one line by element.
    :param pandas.DataFrame axes_inputs: axes inputs dataframe to convert, with one line per axis  (Shoot Apical Meristem)

    :return: The inputs/outputs in a dictionary.
    :rtype: dict [dict]

    see also:: see :attr:`simulation.Simulation.inputs` and :attr:`simulation.Simulation.outputs`
       for the structure of Gas-Exchange inputs/outputs.
    """
    all_elements_dict = {}
    data_columns = element_inputs.columns.difference(ELEMENT_TOPOLOGY_COLUMNS)
    for current_id, current_group in element_inputs.groupby(ELEMENT_TOPOLOGY_COLUMNS):
        current_series = current_group.loc[current_group.first_valid_index()]
        current_dict = current_series[data_columns].to_dict()
        all_elements_dict[current_id] = current_dict

    all_axes_dict = {}
    data_columns = axes_inputs.columns.difference(AXIS_TOPOLOGY_COLUMNS)
    for current_id, current_group in axes_inputs.groupby(AXIS_TOPOLOGY_COLUMNS):
        current_series = current_group.loc[current_group.first_valid_index()]
        current_dict = current_series[data_columns].to_dict()
        all_axes_dict[current_id] = current_dict

    return {'elements': all_elements_dict, 'axes': all_axes_dict}


def to_dataframe(data_dict, element_outputs):
    """
    Convert inputs/outputs from Gas-Exchange format to Pandas dataframe.

    :param dict data_dict: The inputs/outputs in Gas-Exchange format.
    :param list element_outputs: The list of output names for elements

    :return: one dataframe for element outputs
    :rtype: pandas.DataFrame

    see also:: see :attr:`simulation.Simulation.inputs` and :attr:`simulation.Simulation.outputs`
       for the structure of Gas-Exchange inputs/outputs.
    """
    ids_df = pd.DataFrame(data_dict.keys(), columns=ELEMENT_TOPOLOGY_COLUMNS)
    data_df = pd.DataFrame(data_dict.values())
    df = pd.concat([ids_df, data_df], axis=1)
    df.sort_values(by=ELEMENT_TOPOLOGY_COLUMNS, inplace=True)
    columns_sorted = ELEMENT_TOPOLOGY_COLUMNS + element_outputs
    df = df.reindex(columns_sorted, axis=1, copy=False)
    df.reset_index(drop=True, inplace=True)

    return df
