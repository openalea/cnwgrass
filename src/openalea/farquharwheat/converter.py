# -*- coding: latin-1 -*-

from __future__ import division  # use "//" to do integer division
import pandas as pd

"""
    farquharwheat.converter
    ~~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`farquharwheat.converter` defines functions to convert
    :class:`dataframes <pandas.DataFrame>` to/from FarquharWheat inputs or outputs format.

    :copyright: Copyright 2014-2016 INRA-ECOSYS, see AUTHORS.
    :license: see LICENSE for details.

"""

#: the inputs needed by FarquharWheat at element scale
FARQUHARWHEAT_ELEMENTS_INPUTS = ['width', 'height', 'PARa', 'nitrates', 'amino_acids', 'proteins', 'Nstruct', 'green_area',
                                 'sucrose', 'starch', 'fructan', 'PARa_prim', 'area_prim']
#: the inputs needed by FarquharWheat at axis scale
FARQUHARWHEAT_AXES_INPUTS = ['SAM_temperature', 'height_canopy']

#: the outputs computed by FarquharWheat
FARQUHARWHEAT_ELEMENTS_OUTPUTS = ['Ag', 'An', 'Rd', 'Tr', 'Ts', 'gs', 'width', 'height']

#: the inputs and outputs of FarquharWheat.
FARQUHARWHEAT_ELEMENTS_INPUTS_OUTPUTS = set(FARQUHARWHEAT_ELEMENTS_INPUTS + FARQUHARWHEAT_ELEMENTS_OUTPUTS)

#: the columns which define the topology in the input/output elements dataframe
ELEMENT_TOPOLOGY_COLUMNS = ['plant', 'axis', 'metamer', 'organ', 'element']
#: the columns which define the topology in the input/output elements dataframe
AXIS_TOPOLOGY_COLUMNS = ['plant', 'axis']


def from_dataframe(element_inputs, axes_inputs):
    """
    Convert inputs/outputs from Pandas dataframe to Farquhar-Wheat format.

    :param pandas.DataFrame element_inputs: Emerging and mature element inputs dataframe to convert, with one line by element.
    :param pandas.DataFrame axes_inputs: axes inputs dataframe to convert, with one line per axis  (Shoot Apical Meristem)

    :return: The inputs/outputs in a dictionary.
    :rtype: dict [dict]

    seealso:: see :attr:`simulation.Simulation.inputs` and :attr:`simulation.Simulation.outputs`
       for the structure of Farquhar-Wheat inputs/outputs.
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


def to_dataframe(data_dict):
    """
    Convert inputs/outputs from Farquhar-Wheat format to Pandas dataframe.

    :param dict data_dict: The inputs/outputs in Farquhar-Wheat format.

    :return: one dataframe for element outputs
    :rtype: pandas.DataFrame

    seealso:: see :attr:`simulation.Simulation.inputs` and :attr:`simulation.Simulation.outputs`
       for the structure of Farquhar-Wheat inputs/outputs.
    """
    ids_df = pd.DataFrame(data_dict.keys(), columns=ELEMENT_TOPOLOGY_COLUMNS)
    data_df = pd.DataFrame(data_dict.values())
    df = pd.concat([ids_df, data_df], axis=1)
    df.sort_values(by=ELEMENT_TOPOLOGY_COLUMNS, inplace=True)
    columns_sorted = ELEMENT_TOPOLOGY_COLUMNS + [column_name for column_name in sorted(FARQUHARWHEAT_ELEMENTS_INPUTS_OUTPUTS) if column_name in df.columns]
    df = df.reindex(columns_sorted, axis=1, copy=False)
    df.reset_index(drop=True, inplace=True)
    return df
