# -*- coding: latin-1 -*-

from __future__ import division  # use "//" to do integer division
import pandas as pd

"""
    morphogenesis.converter
    ~~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`morphogenesis.converter` defines functions to convert
    :class:`dataframes <pandas.DataFrame>` to/from Morphogenesis inputs or outputs format.

"""

#: the columns which define the topology in the input/output dataframe
HIDDENZONE_TOPOLOGY_COLUMNS = ['plant', 'axis', 'metamer']
ELEMENT_TOPOLOGY_COLUMNS = ['plant', 'axis', 'metamer', 'organ', 'element']  # Mature + emerging elements
AXIS_TOPOLOGY_COLUMNS = ['plant', 'axis']


def from_dataframes(hiddenzone_inputs, element_inputs, axis_inputs):
    """
    Convert inputs/outputs from Pandas dataframe to Morphogenesis format.

    :param pandas.DataFrame axis_inputs: axis inputs dataframe to convert, with one line by axis
    :param pandas.DataFrame hiddenzone_inputs: Hidden zone inputs dataframe to convert, with one line by Hidden zone.
    :param pandas.DataFrame element_inputs: Emerging and mature element inputs dataframe to convert, with one line by element.

    :return: The inputs in a dictionary.
    :rtype: dict [str, dict]

    see also:: :attr:`simulation.Simulation.inputs` for the structure of Morphogenesis inputs.

    """
    all_hiddenzone_dict = {}
    all_element_dict = {}
    all_axis_dict = {}
    all_length_dict = {}
    cumulated_internode_length = {}

    # -- Convert input dataframe into dictionaries

    hiddenzone_inputs_columns = hiddenzone_inputs.columns.difference(HIDDENZONE_TOPOLOGY_COLUMNS)
    emerging_element_inputs_columns = element_inputs.columns.difference(ELEMENT_TOPOLOGY_COLUMNS)
    axis_inputs_columns = axis_inputs.columns.difference(AXIS_TOPOLOGY_COLUMNS)

    grouped_hiddenzone_inputs = hiddenzone_inputs.groupby(HIDDENZONE_TOPOLOGY_COLUMNS, sort=True)
    grouped_element_inputs = element_inputs.groupby(ELEMENT_TOPOLOGY_COLUMNS, sort=True)

    for axis_inputs_id, axis_inputs_group in axis_inputs.groupby(AXIS_TOPOLOGY_COLUMNS):
        # Axis
        axis_inputs_series = axis_inputs_group.loc[axis_inputs_group.first_valid_index()]
        axis_inputs_dict = axis_inputs_series[axis_inputs_columns].to_dict()
        all_axis_dict[axis_inputs_id] = axis_inputs_dict
        # Complete dict of lengths
        all_length_dict[axis_inputs_id] = {}
        for i in range(axis_inputs_series['nb_leaves']):
            all_length_dict[axis_inputs_id][i+1] = {'sheath': [], 'cumulated_internode': []}
        # For coleoptile
        if (*axis_inputs_id, 0) in grouped_hiddenzone_inputs.indices:  # Not emerged
            all_length_dict[axis_inputs_id][0] = {'sheath': [], 'cumulated_internode': []}

        elif (*axis_inputs_id, 0, 'sheath', 'StemElement') in grouped_element_inputs.indices:  # Emerged
            all_length_dict[axis_inputs_id][0] = {'sheath': [], 'cumulated_internode': []}

        cumulated_internode_length[axis_inputs_id] = []

    for element_inputs_id, element_inputs_group in grouped_element_inputs:
        # Elements
        element_inputs_series = element_inputs_group.loc[element_inputs_group.first_valid_index()]
        element_inputs_dict = element_inputs_series[emerging_element_inputs_columns].to_dict()
        all_element_dict[element_inputs_id] = element_inputs_dict
        # Complete dict of lengths
        plant, axis, phytomer_id, organ, element = element_inputs_id
        axis_id = (plant, axis)
        if organ == 'sheath' and not element_inputs_dict['is_growing']:
            all_length_dict[axis_id][phytomer_id]['sheath'].append(element_inputs_dict['length'])
        elif organ == 'internode' and not element_inputs_dict['is_growing']:  # WARNING: this algo won't copy previous internode length for a phytomer without internode
            cumulated_internode_length[axis_id].append(element_inputs_dict['length'])
            if not all_length_dict[axis_id][phytomer_id]['cumulated_internode']:  # if list is empty for that phytomer, the list of all phytomer lengths is written
                all_length_dict[axis_id][phytomer_id]['cumulated_internode'].extend(cumulated_internode_length[axis_id])
            else:  # only the last internode length is written (case of organs with hidden and visible part)
                all_length_dict[axis_id][phytomer_id]['cumulated_internode'].append(element_inputs_dict['length'])

    for hiddenzone_inputs_id, hiddenzone_inputs_group in grouped_hiddenzone_inputs:
        # hiddenzone
        hiddenzone_inputs_series = hiddenzone_inputs_group.loc[hiddenzone_inputs_group.first_valid_index()]
        hiddenzone_inputs_dict = hiddenzone_inputs_series[hiddenzone_inputs_columns].to_dict()
        all_hiddenzone_dict[hiddenzone_inputs_id] = hiddenzone_inputs_dict
        # Complete dict of length
        plant, axis, phytomer_id = hiddenzone_inputs_id
        axis_id = (plant, axis)
        # Not emerged coleoptile
        if phytomer_id == 0:
            all_length_dict[axis_id][phytomer_id]['sheath'].append(hiddenzone_inputs_dict['leaf_L'])

        if hiddenzone_inputs_dict['leaf_is_emerged'] and hiddenzone_inputs_dict['leaf_is_growing']:
            growing_sheath_length = max(0, hiddenzone_inputs_series['leaf_L'] - hiddenzone_inputs_series['lamina_Lmax'])  # TODO should be moved to model.py
            all_length_dict[axis_id][phytomer_id]['sheath'].append(growing_sheath_length)
        if hiddenzone_inputs_dict['internode_is_growing']:
            cumulated_internode_length[axis_id].append(hiddenzone_inputs_dict['internode_L'])
            all_length_dict[axis_id][phytomer_id]['cumulated_internode'].extend(cumulated_internode_length[axis_id])

        elif not hiddenzone_inputs_dict['internode_is_growing'] and (*hiddenzone_inputs_id, 'internode') not in element_inputs.groupby(ELEMENT_TOPOLOGY_COLUMNS[:-1]).groups.keys():
            all_length_dict[axis_id][phytomer_id]['cumulated_internode'].extend(cumulated_internode_length[axis_id])

    return {'hiddenzone': all_hiddenzone_dict, 'elements': all_element_dict, 'axes': all_axis_dict, 'sheath_internode_lengths': all_length_dict}


def to_dataframes(data_dict, axis_outputs, hiddenzone_outputs, element_outputs):
    """
    Convert outputs from Morphogenesis format to Pandas dataframe.

    :param dict data_dict: The outputs in Morphogenesis format.
    :param list axis_outputs: The list of output names for axes
    :param list hiddenzone_outputs: The list of output names for hiddenzones
    :param list element_outputs: The list of output names for elements

    :return: One dataframe for hiddenzone outputs, one dataframe for element outputs and one dataframe for axis outputs.
    :rtype: (pandas.DataFrame, pandas.DataFrame, pandas.DataFrame)

    """
    dataframes_dict = {}
    for (current_key, current_topology_columns, current_outputs_names) in (('hiddenzone', HIDDENZONE_TOPOLOGY_COLUMNS, hiddenzone_outputs),
                                                                           ('elements', ELEMENT_TOPOLOGY_COLUMNS, element_outputs),
                                                                           ('axes', AXIS_TOPOLOGY_COLUMNS, axis_outputs)):

        current_data_dict = data_dict[current_key]
        current_ids_df = pd.DataFrame(current_data_dict.keys(), columns=current_topology_columns)
        current_data_df = pd.DataFrame(current_data_dict.values())
        current_df = pd.concat([current_ids_df, current_data_df], axis=1)
        current_df.sort_values(by=current_topology_columns, inplace=True)
        current_columns_sorted = current_topology_columns + current_outputs_names
        current_df = current_df.reindex(current_columns_sorted, axis=1, copy=False)
        current_df.reset_index(drop=True, inplace=True)
        dataframes_dict[current_key] = current_df

    return dataframes_dict['hiddenzone'], dataframes_dict['elements'], dataframes_dict['axes']
