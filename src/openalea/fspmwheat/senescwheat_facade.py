# -*- coding: latin-1 -*-

import numpy as np

from senescwheat import converter, simulation

from fspmwheat import tools

"""
    fspmwheat.senescwheat_facade
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    The module :mod:`fspmwheat.senescwheat_facade` is a facade of the model SenescWheat.
    This module permits to initialize and run the model SenescWheat from a :class:`MTG <openalea.mtg.mtg.MTG>`
    in a convenient and transparent way, wrapping all the internal complexity of the model, and dealing
    with all the tedious initialization and conversion processes.
    :copyright: Copyright 2014-2016 INRA-ECOSYS, see AUTHORS.
    :license: see LICENSE for details.
"""

#: the name of the photosynthetic organs modeled by SenescWheat
PHOTOSYNTHETIC_ORGANS_NAMES = {'internode', 'blade', 'sheath', 'peduncle', 'ear'}

#: the columns which define the topology in the organs scale dataframe shared between all models
SHARED_AXES_INPUTS_OUTPUTS_INDEXES = ['plant', 'axis']

#: the columns which define the topology in the organs scale dataframe shared between all models
SHARED_ORGANS_INPUTS_OUTPUTS_INDEXES = ['plant', 'axis', 'organ']

#: the columns which define the topology in the elements scale dataframe shared between all models
SHARED_ELEMENTS_INPUTS_OUTPUTS_INDEXES = ['plant', 'axis', 'metamer', 'organ', 'element']


class SenescWheatFacade(object):
    """
    The SenescWheatFacade class permits to initialize, run the model SenescWheat
    from a :class:`MTG <openalea.mtg.mtg.MTG>`, and update the MTG and the dataframes
    shared between all models.
    Use :meth:`run` to run the model.
    """

    def __init__(self, shared_mtg, delta_t,
                 model_roots_inputs_df,
                 model_axes_inputs_df,
                 model_elements_inputs_df,
                 shared_organs_inputs_outputs_df,
                 shared_axes_inputs_outputs_df,
                 shared_elements_inputs_outputs_df,
                 update_parameters=None,
                 update_shared_df=True,
                 cnwheat_roots=True):

        """
        :param openalea.mtg.mtg.MTG shared_mtg: The MTG shared between all models.
        :param int delta_t: The delta between two runs, in seconds.
        :param pandas.DataFrame model_roots_inputs_df: the inputs of the model at roots scale.
        :param pandas.DataFrame model_axes_inputs_df: the inputs of the model at axes scale.
        :param pandas.DataFrame model_elements_inputs_df: the inputs of the model at elements scale.
        :param pandas.DataFrame shared_organs_inputs_outputs_df: the dataframe of inputs and outputs at organs scale shared between all models.
        :param pandas.DataFrame shared_axes_inputs_outputs_df: the dataframe of inputs and outputs at axis scale shared between all models.
        :param pandas.DataFrame shared_elements_inputs_outputs_df: the dataframe of inputs and outputs at element scale shared between all models.
        :param dict update_parameters: A dictionary with the parameters to update, should have the form {'param1': value1, 'param2': value2, ...}.
        :param bool update_shared_df: If `True`  update the shared dataframes at init and at each run (unless stated otherwise)
        """

        self._shared_mtg = shared_mtg  #: the MTG shared between all models
        self.cnwheat_roots = cnwheat_roots
        self._simulation = simulation.Simulation(delta_t=delta_t, update_parameters=update_parameters, cnwheat_roots=cnwheat_roots)  #: the simulator to use to run the model

        all_senescwheat_inputs_dict = converter.from_dataframes(model_roots_inputs_df, model_axes_inputs_df, model_elements_inputs_df)
        self._update_shared_MTG(all_senescwheat_inputs_dict['roots'], all_senescwheat_inputs_dict['axes'], all_senescwheat_inputs_dict['elements'])

        self._shared_organs_inputs_outputs_df = shared_organs_inputs_outputs_df  #: the dataframe at organs scale shared between all models
        self._shared_axes_inputs_outputs_df = shared_axes_inputs_outputs_df  #: the dataframe at axis scale shared between all models
        self._shared_elements_inputs_outputs_df = shared_elements_inputs_outputs_df  #: the dataframe at elements scale shared between all models
        self._update_shared_df = update_shared_df
        if self._update_shared_df:
            self._update_shared_dataframes(model_roots_inputs_df, model_axes_inputs_df, model_elements_inputs_df)

    def run(self, forced_max_protein_elements=None, postflowering_stages=False, update_shared_df=None):
        """
        Run the model and update the MTG and the dataframes shared between all models.

        :param set forced_max_protein_elements: The elements ids with fixed max proteins.
        :param bool postflowering_stages: True to run a simulation with postflo parameter
        :param bool update_shared_df: if 'True', update the shared dataframes at this time step.
        """

        self._initialize_model()
        self._simulation.run(forced_max_protein_elements=forced_max_protein_elements, postflowering_stages=postflowering_stages)
        self._update_shared_MTG(self._simulation.outputs['roots'], self._simulation.outputs['axes'], self._simulation.outputs['elements'])

        if update_shared_df or (update_shared_df is None and self._update_shared_df):
            senescwheat_roots_outputs_df, senescwheat_axes_outputs_df, senescwheat_elements_outputs_df = converter.to_dataframes(self._simulation.outputs)
            self._update_shared_dataframes(senescwheat_roots_outputs_df, senescwheat_axes_outputs_df, senescwheat_elements_outputs_df)

    def _initialize_model(self):
        """
        Initialize the inputs of the model from the MTG shared between all models.
        """
        all_senescwheat_roots_inputs_dict = {}
        all_senescwheat_axes_inputs_dict = {}
        all_senescwheat_elements_inputs_dict = {}

        # traverse the MTG recursively from the top ...
        for mtg_plant_vid in self._shared_mtg.components_iter(self._shared_mtg.root):
            mtg_plant_index = int(self._shared_mtg.index(mtg_plant_vid))
            for mtg_axis_vid in self._shared_mtg.components_iter(mtg_plant_vid):
                mtg_axis_label = self._shared_mtg.label(mtg_axis_vid)
                if mtg_axis_label != 'MS':
                    continue
                axis_id = (mtg_plant_index, mtg_axis_label)
                mtg_axis_properties = self._shared_mtg.get_vertex_property(mtg_axis_vid)
                if set(mtg_axis_properties).issuperset(converter.SENESCWHEAT_AXES_INPUTS):
                    senescwheat_axis_inputs_dict = {}
                    for senescwheat_axis_input_name in converter.SENESCWHEAT_AXES_INPUTS:
                        senescwheat_axis_inputs_dict[senescwheat_axis_input_name] = mtg_axis_properties[senescwheat_axis_input_name]
                    all_senescwheat_axes_inputs_dict[axis_id] = senescwheat_axis_inputs_dict
                if 'roots' in mtg_axis_properties:
                    mtg_roots_properties = mtg_axis_properties['roots']
                    if set(mtg_roots_properties).issuperset(converter.SENESCWHEAT_ROOTS_INPUTS):
                        senescwheat_roots_inputs_dict = {}
                        for senescwheat_roots_input_name in converter.SENESCWHEAT_ROOTS_INPUTS:
                            senescwheat_roots_inputs_dict[senescwheat_roots_input_name] = mtg_roots_properties[senescwheat_roots_input_name]
                        all_senescwheat_roots_inputs_dict[axis_id] = senescwheat_roots_inputs_dict
                for mtg_metamer_vid in self._shared_mtg.components_iter(mtg_axis_vid):
                    mtg_metamer_index = int(self._shared_mtg.index(mtg_metamer_vid))
                    for mtg_organ_vid in self._shared_mtg.components_iter(mtg_metamer_vid):
                        mtg_organ_label = self._shared_mtg.label(mtg_organ_vid)
                        if mtg_organ_label not in PHOTOSYNTHETIC_ORGANS_NAMES:
                            continue
                        # if np.nan_to_num( self._shared_mtg.property('length').get(mtg_organ_vid,0)) == 0: continue
                        for mtg_element_vid in self._shared_mtg.components_iter(mtg_organ_vid):
                            mtg_element_properties = self._shared_mtg.get_vertex_property(mtg_element_vid)
                            mtg_element_label = self._shared_mtg.label(mtg_element_vid)
                            element_id = (mtg_plant_index, mtg_axis_label, mtg_metamer_index, mtg_organ_label, mtg_element_label)
                            if np.nan_to_num(self._shared_mtg.property('length').get(mtg_element_vid, 0)) == 0:
                                continue
                            if set(mtg_element_properties).issuperset(converter.SENESCWHEAT_ELEMENTS_INPUTS):
                                senescwheat_element_inputs_dict = {}
                                for senescwheat_element_input_name in converter.SENESCWHEAT_ELEMENTS_INPUTS:
                                    senescwheat_element_inputs_dict[senescwheat_element_input_name] = mtg_element_properties[senescwheat_element_input_name]
                                all_senescwheat_elements_inputs_dict[element_id] = senescwheat_element_inputs_dict
                                # TODO: temporary ; replace 'SENESCWHEAT_ELEMENT_PROPERTIES_TEMP' by default values
                                SENESCWHEAT_ELEMENT_PROPERTIES_TEMP = {'starch': 0, 'max_proteins': 0, 'amino_acids': 0,
                                                                       'proteins': 0, 'Nstruct': 0, 'mstruct': 0, 'fructan': 0,
                                                                       'sucrose': 0, 'green_area': 0, 'cytokinins': 0}
                                senescwheat_element_inputs_dict = {}
                                for senescwheat_element_input_name in converter.SENESCWHEAT_ELEMENTS_INPUTS:
                                    if senescwheat_element_input_name in mtg_element_properties:
                                        senescwheat_element_inputs_dict[senescwheat_element_input_name] = mtg_element_properties[senescwheat_element_input_name]
                                    else:
                                        senescwheat_element_inputs_dict[senescwheat_element_input_name] = SENESCWHEAT_ELEMENT_PROPERTIES_TEMP[senescwheat_element_input_name]
                                all_senescwheat_elements_inputs_dict[element_id] = senescwheat_element_inputs_dict

        self._simulation.initialize({'roots': all_senescwheat_roots_inputs_dict, 'axes': all_senescwheat_axes_inputs_dict, 'elements': all_senescwheat_elements_inputs_dict})

    def _update_shared_MTG(self, senescwheat_roots_data_dict, senescwheat_axes_data_dict, senescwheat_elements_data_dict):
        """
        Update the MTG shared between all models from the inputs or the outputs of the model.
        :param dict senescwheat_roots_data_dict: Senesc-Wheat outputs at root scale
        :param dict senescwheat_axes_data_dict: Senesc-Wheat outputs at axis scale
        :param dict senescwheat_elements_data_dict: Senesc-Wheat outputs at element scale
        """

        # add the properties if needed
        mtg_property_names = self._shared_mtg.property_names()
        if 'roots' not in mtg_property_names:
            self._shared_mtg.add_property('roots')
        for senescwheat_elements_data_name in converter.SENESCWHEAT_ELEMENTS_INPUTS_OUTPUTS:
            if senescwheat_elements_data_name not in mtg_property_names:
                self._shared_mtg.add_property(senescwheat_elements_data_name)

        # traverse the MTG recursively from top ...
        for mtg_plant_vid in self._shared_mtg.components_iter(self._shared_mtg.root):
            mtg_plant_index = int(self._shared_mtg.index(mtg_plant_vid))
            for mtg_axis_vid in self._shared_mtg.components_iter(mtg_plant_vid):
                mtg_axis_label = self._shared_mtg.label(mtg_axis_vid)
                if mtg_axis_label != 'MS':
                    continue
                
                if self.cnwheat_roots:
                    # update the axis property in the MTG
                    axis_id = (mtg_plant_index, mtg_axis_label)
                    mtg_axis_properties = self._shared_mtg.get_vertex_property(mtg_axis_vid)
                    mtg_axis_properties.update(senescwheat_axes_data_dict.get(axis_id, []))
                    # update the roots in the MTG
                    if axis_id not in senescwheat_roots_data_dict:
                        continue
                    if 'roots' not in self._shared_mtg.get_vertex_property(mtg_axis_vid):
                        self._shared_mtg.property('roots')[mtg_axis_vid] = {}
                    mtg_roots_properties = self._shared_mtg.get_vertex_property(mtg_axis_vid)['roots']
                    mtg_roots_properties.update(senescwheat_roots_data_dict[axis_id])

                for mtg_metamer_vid in self._shared_mtg.components_iter(mtg_axis_vid):
                    mtg_metamer_index = int(self._shared_mtg.index(mtg_metamer_vid))
                    for mtg_organ_vid in self._shared_mtg.components_iter(mtg_metamer_vid):
                        mtg_organ_label = self._shared_mtg.label(mtg_organ_vid)
                        # senesced_length_organ = 0.  # Temporaire
                        if mtg_organ_label not in PHOTOSYNTHETIC_ORGANS_NAMES:
                            continue
                        for mtg_element_vid in self._shared_mtg.components_iter(mtg_organ_vid):
                            mtg_element_label = self._shared_mtg.label(mtg_element_vid)
                            element_id = (mtg_plant_index, mtg_axis_label, mtg_metamer_index, mtg_organ_label, mtg_element_label)
                            if element_id not in senescwheat_elements_data_dict:
                                continue
                            # update the element in the MTG
                            senescwheat_element_data_dict = senescwheat_elements_data_dict[element_id]
                            for senescwheat_element_data_name, senescwheat_element_data_value in senescwheat_element_data_dict.items():
                                self._shared_mtg.property(senescwheat_element_data_name)[mtg_element_vid] = senescwheat_element_data_value
                                # Temporaire avant de trouver une solution pour :
                                # 1) piloter la senescence des feuilles par green_area plutot que par senesced_length,
                                # 2) updater les organes � partir des �l�ments et non l'inverse.
                                if senescwheat_element_data_name == 'senesced_length_element' and mtg_element_label in ['LeafElement1', 'StemElement']:
                                    self._shared_mtg.property('senesced_length')[mtg_organ_vid] = np.nan_to_num(self._shared_mtg.property(senescwheat_element_data_name).get(mtg_element_vid, 0.))

    def _update_shared_dataframes(self, senescwheat_roots_data_df, senescwheat_axes_data_df, senescwheat_elements_data_df):
        """
        Update the dataframes shared between all models from the inputs dataframes or the outputs dataframes of the model.
        :param pandas.DataFrame senescwheat_roots_data_df: Elong-Wheat shared dataframe at root scale
        :param pandas.DataFrame senescwheat_axes_data_df: Elong-Wheat shared dataframe at axis scale
        :param pandas.DataFrame senescwheat_elements_data_df: Elong-Wheat shared dataframe at element scale
        """

        for senescwheat_data_df, \
            shared_inputs_outputs_indexes, \
            shared_inputs_outputs_df in ((senescwheat_roots_data_df, SHARED_ORGANS_INPUTS_OUTPUTS_INDEXES, self._shared_organs_inputs_outputs_df),
                                         (senescwheat_axes_data_df, SHARED_AXES_INPUTS_OUTPUTS_INDEXES, self._shared_axes_inputs_outputs_df),
                                         (senescwheat_elements_data_df, SHARED_ELEMENTS_INPUTS_OUTPUTS_INDEXES, self._shared_elements_inputs_outputs_df)):

            if senescwheat_data_df is senescwheat_roots_data_df:
                senescwheat_data_df = senescwheat_data_df.copy()
                senescwheat_data_df.loc[:, 'organ'] = 'roots'

            tools.combine_dataframes_inplace(senescwheat_data_df, shared_inputs_outputs_indexes, shared_inputs_outputs_df)
