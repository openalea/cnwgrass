# -*- coding: latin-1 -*-

from __future__ import division  # use "//" to do integer division
import pandas as pd

"""
    senescwheat.converter
    ~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`senescwheat.converter` defines functions to convert
    :class:`dataframes <pandas.DataFrame>` to/from SenescWheat inputs or outputs format.

    :copyright: Copyright 2014-2016 INRA-ECOSYS, see AUTHORS.
    :license: see LICENSE for details.

"""

#: the inputs needed by SenescWheat at roots scale
SENESCWHEAT_ROOTS_INPUTS = ['sucrose', 'amino_acids', 'mstruct', 'senesced_mstruct', 'Nstruct', 'cytokinins']

#: the inputs needed by SenescWheat at axis scale
SENESCWHEAT_AXES_INPUTS = ['delta_teq', 'delta_teq_roots', 'sum_TT']

#: the inputs needed by SenescWheat at element scale
SENESCWHEAT_ELEMENTS_INPUTS = ['green_area', 'senesced_length_element', 'length', 'proteins', 'mstruct', 'senesced_mstruct', 'max_proteins', 'Nstruct', 'nitrates', 'amino_acids', 'starch', 'fructan',
                               'cytokinins', 'sucrose', 'is_growing', 'max_mstruct', 'Nresidual', 'age']

#: the outputs computed by SenescWheat at roots scale
SENESCWHEAT_ROOTS_OUTPUTS = ['rate_mstruct_death', 'mstruct', 'senesced_mstruct', 'Nstruct', 'cytokinins']

#: the outputs computed by SenescWheat at axis scale
SENESCWHEAT_AXES_OUTPUTS = []

#: the outputs computed by SenescWheat at elements scale
SENESCWHEAT_ELEMENTS_OUTPUTS = ['senesced_length_element', 'green_area', 'mstruct', 'senesced_mstruct', 'Nstruct', 'starch', 'sucrose', 'fructan', 'proteins', 'amino_acids', 'cytokinins', 'nitrates',
                                'N_content_total', 'Nresidual', 'is_over']

#: the inputs and outputs of SenescWheat at roots scale
SENESCWHEAT_ROOTS_INPUTS_OUTPUTS = sorted(list(set(SENESCWHEAT_ROOTS_INPUTS + SENESCWHEAT_ROOTS_OUTPUTS)))

#: the inputs and outputs of SenescWheat at axis scale
SENESCWHEAT_AXES_INPUTS_OUTPUTS = sorted(list(set(SENESCWHEAT_AXES_INPUTS + SENESCWHEAT_AXES_OUTPUTS)))

#: the inputs and outputs of SenescWheat at elements scale
SENESCWHEAT_ELEMENTS_INPUTS_OUTPUTS = sorted(list(set(SENESCWHEAT_ELEMENTS_INPUTS + SENESCWHEAT_ELEMENTS_OUTPUTS)))

#: the columns which define the topology of a roots in the input/output dataframe
ROOTS_TOPOLOGY_COLUMNS = ['plant', 'axis']

#: the columns which define the topology of a roots in the input/output dataframe
AXES_TOPOLOGY_COLUMNS = ['plant', 'axis']

#: the columns which define the topology of an element in the input/output dataframe
ELEMENTS_TOPOLOGY_COLUMNS = ['plant', 'axis', 'metamer', 'organ', 'element']


def from_dataframes(roots_inputs, axes_inputs, elements_inputs):
    """
    Convert inputs/outputs from Pandas dataframes to Senesc-Wheat format.

    :param pandas.DataFrame roots_inputs: Roots inputs dataframe to convert, with one line by roots.
    :param pandas.DataFrame axes_inputs: axes inputs dataframe to convert, with one line by axis.
    :param pandas.DataFrame elements_inputs: Elements inputs dataframe to convert, with one line by element.

    :return: The inputs/outputs in a dictionary.
    :rtype: dict [str, dict]

    .. seealso:: see :attr:`simulation.Simulation.inputs` and :attr:`simulation.Simulation.outputs`
       for the structure of Senesc-Wheat inputs/outputs.

    """
    all_roots_dict = {}
    all_axes_dict = {}
    all_elements_dict = {}
    for (all_current_dict, current_dataframe, current_topology_columns) in ((all_roots_dict, roots_inputs, ROOTS_TOPOLOGY_COLUMNS),
                                                                            (all_axes_dict, axes_inputs, AXES_TOPOLOGY_COLUMNS),
                                                                            (all_elements_dict, elements_inputs, ELEMENTS_TOPOLOGY_COLUMNS)):
        current_columns = current_dataframe.columns.difference(current_topology_columns)
        for current_id, current_group in current_dataframe.groupby(current_topology_columns):
            current_series = current_group.loc[current_group.first_valid_index()]
            current_dict = current_series[current_columns].to_dict()
            all_current_dict[current_id] = current_dict

    return {'roots': all_roots_dict, 'axes': all_axes_dict, 'elements': all_elements_dict}


def to_dataframes(data_dict):
    """
    Convert inputs/outputs from Senesc-Wheat format to Pandas dataframe.

    :param dict data_dict: The inputs/outputs in Senesc-Wheat format.

    :return: One dataframe for roots inputs/outputs, one dataframe for axes inputs/outputs,  one dataframe for elements inputs/outputs.
    :rtype: (pandas.DataFrame, pandas.DataFrame, pandas.DataFrame)

    .. seealso:: see :attr:`simulation.Simulation.inputs` and :attr:`simulation.Simulation.outputs`
       for the structure of Senesc-Wheat inputs/outputs.

    """
    dataframes_dict = {}
    for (current_key, current_topology_columns, current_inputs_outputs_names) in (('roots', ROOTS_TOPOLOGY_COLUMNS, SENESCWHEAT_ROOTS_INPUTS_OUTPUTS),
                                                                                  ('axes', AXES_TOPOLOGY_COLUMNS, SENESCWHEAT_AXES_INPUTS_OUTPUTS),
                                                                                  ('elements', ELEMENTS_TOPOLOGY_COLUMNS, SENESCWHEAT_ELEMENTS_INPUTS_OUTPUTS)):
        current_data_dict = data_dict[current_key]
        current_ids_df = pd.DataFrame(current_data_dict.keys(), columns=current_topology_columns)
        current_data_df = pd.DataFrame(current_data_dict.values())
        current_df = pd.concat([current_ids_df, current_data_df], axis=1)
        current_df.sort_values(by=current_topology_columns, inplace=True)
        current_columns_sorted = current_topology_columns + [input_output for input_output in current_inputs_outputs_names if input_output in current_df.columns]
        current_df = current_df.reindex(current_columns_sorted, axis=1, copy=False)
        current_df.reset_index(drop=True, inplace=True)
        dataframes_dict[current_key] = current_df

    return dataframes_dict['roots'], dataframes_dict['axes'], dataframes_dict['elements']
