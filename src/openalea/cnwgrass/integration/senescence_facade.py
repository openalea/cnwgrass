# -*- coding: latin-1 -*-

import numpy as np

from openalea.cnwgrass.senescence import converter, simulation

from openalea.cnwgrass.integration import tools

"""
    integration.senescence_facade
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    The module :mod:`integration.senescence_facade` is a facade of the model Senescence.
    This module permits to initialize and run the model Senescence from a :class:`MTG <openalea.mtg.mtg.MTG>`
    in a convenient and transparent way, wrapping all the internal complexity of the model, and dealing
    with all the tedious initialization and conversion processes.

"""

#: the name of the photosynthetic organs modeled by Senescence
PHOTOSYNTHETIC_ORGANS_NAMES = {'internode', 'blade', 'sheath', 'peduncle', 'ear'}

#: the columns which define the topology in the organs scale dataframe shared between all models
SHARED_AXES_INPUTS_OUTPUTS_INDEXES = ['plant', 'axis']

#: the columns which define the topology in the organs scale dataframe shared between all models
SHARED_ORGANS_INPUTS_OUTPUTS_INDEXES = ['plant', 'axis', 'organ']

#: the columns which define the topology in the elements scale dataframe shared between all models
SHARED_ELEMENTS_INPUTS_OUTPUTS_INDEXES = ['plant', 'axis', 'metamer', 'organ', 'element']


class SENESCENCEFacade(object):
    """
    The SENESCENCEFacade class permits to initialize, run the model Senescence
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
                 update_shared_df=True):

        """
        :param openalea.mtg.mtg.MTG shared_mtg: The MTG shared between all models.
        :param int delta_t: The delta between two runs, in seconds.
        :param pandas.DataFrame model_roots_inputs_df: the inputs of the model at roots scale.
        :param pandas.DataFrame model_axes_inputs_df: the inputs of the model at axes scale.
        :param pandas.DataFrame model_elements_inputs_df: the inputs of the model at elements scale.
        :param pandas.DataFrame shared_organs_inputs_outputs_df: the dataframe of inputs and outputs at organs scale shared between all models.
        :param pandas.DataFrame shared_axes_inputs_outputs_df: the dataframe of inputs and outputs at axis scale shared between all models.
        :param pandas.DataFrame shared_elements_inputs_outputs_df: the dataframe of inputs and outputs at element scale shared between all models.
        :param dict or None update_parameters: A dictionary with the parameters to update, should have the form {'param1': value1, 'param2': value2, ...}.
        :param bool update_shared_df: If `True`  update the shared dataframes at init and at each run (unless stated otherwise)
        """

        self._shared_mtg = shared_mtg  #: the MTG shared between all models

        self._simulation = simulation.Simulation(delta_t=delta_t, update_parameters=update_parameters)  #: the simulator to use to run the model

        all_senescence_inputs_dict = converter.from_dataframes(model_roots_inputs_df, model_axes_inputs_df, model_elements_inputs_df)
        self._update_shared_MTG(all_senescence_inputs_dict['roots'], all_senescence_inputs_dict['axes'], all_senescence_inputs_dict['elements'])

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
            senescence_roots_outputs_df, senescence_axes_outputs_df, senescence_elements_outputs_df = converter.to_dataframes(self._simulation.outputs)
            self._update_shared_dataframes(senescence_roots_outputs_df, senescence_axes_outputs_df, senescence_elements_outputs_df)

    def _initialize_model(self):
        """
        Initialize the inputs of the model from the MTG shared between all models.
        """
        all_senescence_roots_inputs_dict = {}
        all_senescence_axes_inputs_dict = {}
        all_senescence_elements_inputs_dict = {}

        # traverse the MTG recursively from the top ...
        for mtg_plant_vid in self._shared_mtg.components_iter(self._shared_mtg.root):
            mtg_plant_index = int(self._shared_mtg.index(mtg_plant_vid))
            for mtg_axis_vid in self._shared_mtg.components_iter(mtg_plant_vid):
                mtg_axis_label = self._shared_mtg.label(mtg_axis_vid)
                if mtg_axis_label != 'MS':
                    continue
                axis_id = (mtg_plant_index, mtg_axis_label)
                mtg_axis_properties = self._shared_mtg.get_vertex_property(mtg_axis_vid)
                if set(mtg_axis_properties).issuperset(converter.SENESCENCE_AXES_INPUTS):
                    senescence_axis_inputs_dict = {}
                    for senescence_axis_input_name in converter.SENESCENCE_AXES_INPUTS:
                        senescence_axis_inputs_dict[senescence_axis_input_name] = mtg_axis_properties[senescence_axis_input_name]
                    all_senescence_axes_inputs_dict[axis_id] = senescence_axis_inputs_dict
                if 'roots' in mtg_axis_properties:
                    mtg_roots_properties = mtg_axis_properties['roots']
                    if set(mtg_roots_properties).issuperset(converter.SENESCENCE_ROOTS_INPUTS):
                        senescence_roots_inputs_dict = {}
                        for senescence_roots_input_name in converter.SENESCENCE_ROOTS_INPUTS:
                            senescence_roots_inputs_dict[senescence_roots_input_name] = mtg_roots_properties[senescence_roots_input_name]
                        all_senescence_roots_inputs_dict[axis_id] = senescence_roots_inputs_dict
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
                            if set(mtg_element_properties).issuperset(converter.SENESCENCE_ELEMENTS_INPUTS):
                                senescence_element_inputs_dict = {}
                                for senescence_element_input_name in converter.SENESCENCE_ELEMENTS_INPUTS:
                                    senescence_element_inputs_dict[senescence_element_input_name] = mtg_element_properties[senescence_element_input_name]
                                all_senescence_elements_inputs_dict[element_id] = senescence_element_inputs_dict
                                # TODO: temporary ; replace 'SENESCENCE_ELEMENT_PROPERTIES_TEMP' by default values
                                SENESCENCE_ELEMENT_PROPERTIES_TEMP = {'starch': 0, 'max_proteins': 0, 'amino_acids': 0,
                                                                       'proteins': 0, 'Nstruct': 0, 'mstruct': 0, 'fructan': 0,
                                                                       'sucrose': 0, 'green_area': 0, 'cytokinins': 0}
                                senescence_element_inputs_dict = {}
                                for senescence_element_input_name in converter.SENESCENCE_ELEMENTS_INPUTS:
                                    if senescence_element_input_name in mtg_element_properties:
                                        senescence_element_inputs_dict[senescence_element_input_name] = mtg_element_properties[senescence_element_input_name]
                                    else:
                                        senescence_element_inputs_dict[senescence_element_input_name] = SENESCENCE_ELEMENT_PROPERTIES_TEMP[senescence_element_input_name]
                                all_senescence_elements_inputs_dict[element_id] = senescence_element_inputs_dict

        self._simulation.initialize({'roots': all_senescence_roots_inputs_dict, 'axes': all_senescence_axes_inputs_dict, 'elements': all_senescence_elements_inputs_dict})

    def _update_shared_MTG(self, senescence_roots_data_dict, senescence_axes_data_dict, senescence_elements_data_dict):
        """
        Update the MTG shared between all models from the inputs or the outputs of the model.
        :param dict senescence_roots_data_dict: Senescence outputs at root scale
        :param dict senescence_axes_data_dict: Senescence outputs at axis scale
        :param dict senescence_elements_data_dict: Senescence outputs at element scale
        """

        # add the properties if needed
        mtg_property_names = self._shared_mtg.property_names()
        for senescence_axes_data_name in converter.SENESCENCE_AXES_INPUTS_OUTPUTS:
            if senescence_axes_data_name not in mtg_property_names:
                self._shared_mtg.add_property(senescence_axes_data_name)
        if 'roots' not in mtg_property_names:
            self._shared_mtg.add_property('roots')
        for senescence_elements_data_name in converter.SENESCENCE_ELEMENTS_INPUTS_OUTPUTS:
            if senescence_elements_data_name not in mtg_property_names:
                self._shared_mtg.add_property(senescence_elements_data_name)

        # traverse the MTG recursively from top ...
        for mtg_plant_vid in self._shared_mtg.components_iter(self._shared_mtg.root):
            mtg_plant_index = int(self._shared_mtg.index(mtg_plant_vid))
            for mtg_axis_vid in self._shared_mtg.components_iter(mtg_plant_vid):
                mtg_axis_label = self._shared_mtg.label(mtg_axis_vid)
                if mtg_axis_label != 'MS':
                    continue

                # update the axis property in the MTG
                axis_id = (mtg_plant_index, mtg_axis_label)

                if axis_id in senescence_axes_data_dict:
                    senescence_axis_data_dict = senescence_axes_data_dict[axis_id]
                    for axis_data_name, axis_data_value in senescence_axis_data_dict.items():
                        self._shared_mtg.property(axis_data_name)[mtg_axis_vid] = axis_data_value

                # update the roots in the MTG
                if axis_id not in senescence_roots_data_dict:
                    continue
                if 'roots' not in self._shared_mtg.get_vertex_property(mtg_axis_vid):
                    self._shared_mtg.property('roots')[mtg_axis_vid] = {}
                mtg_roots_properties = self._shared_mtg.get_vertex_property(mtg_axis_vid)['roots']
                for roots_data_name, roots_data_value in senescence_roots_data_dict.items():
                    self._shared_mtg.property(roots_data_name)[mtg_axis_vid] = roots_data_value

                for mtg_metamer_vid in self._shared_mtg.components_iter(mtg_axis_vid):
                    mtg_metamer_index = int(self._shared_mtg.index(mtg_metamer_vid))
                    for mtg_organ_vid in self._shared_mtg.components_iter(mtg_metamer_vid):
                        mtg_organ_label = self._shared_mtg.label(mtg_organ_vid)
                        # senesced_length_organ = 0.  # Temporary
                        if mtg_organ_label not in PHOTOSYNTHETIC_ORGANS_NAMES:
                            continue
                        for mtg_element_vid in self._shared_mtg.components_iter(mtg_organ_vid):
                            mtg_element_label = self._shared_mtg.label(mtg_element_vid)
                            element_id = (mtg_plant_index, mtg_axis_label, mtg_metamer_index, mtg_organ_label, mtg_element_label)
                            if element_id not in senescence_elements_data_dict:
                                continue
                            # update the element in the MTG
                            senescence_element_data_dict = senescence_elements_data_dict[element_id]
                            for senescence_element_data_name, senescence_element_data_value in senescence_element_data_dict.items():
                                self._shared_mtg.property(senescence_element_data_name)[mtg_element_vid] = senescence_element_data_value
                                # Temporaire avant de trouver une solution pour :
                                # 1) piloter la senescence des feuilles par green_area plutot que par senesced_length,
                                # 2) updater les organes à partir des éléments et non l'inverse.
                                if senescence_element_data_name == 'senesced_length_element' and mtg_element_label in ['LeafElement1', 'StemElement']:
                                    self._shared_mtg.property('senesced_length')[mtg_organ_vid] = np.nan_to_num(self._shared_mtg.property(senescence_element_data_name).get(mtg_element_vid, 0.))

    def _update_shared_dataframes(self, senescence_roots_data_df, senescence_axes_data_df, senescence_elements_data_df):
        """
        Update the dataframes shared between all models from the inputs dataframes or the outputs dataframes of the model.
        :param pandas.DataFrame senescence_roots_data_df: Morphogenesis shared dataframe at root scale
        :param pandas.DataFrame senescence_axes_data_df: Morphogenesis shared dataframe at axis scale
        :param pandas.DataFrame senescence_elements_data_df: Morphogenesis shared dataframe at element scale
        """

        for senescence_data_df, \
            shared_inputs_outputs_indexes, \
            shared_inputs_outputs_df in ((senescence_roots_data_df, SHARED_ORGANS_INPUTS_OUTPUTS_INDEXES, self._shared_organs_inputs_outputs_df),
                                         (senescence_axes_data_df, SHARED_AXES_INPUTS_OUTPUTS_INDEXES, self._shared_axes_inputs_outputs_df),
                                         (senescence_elements_data_df, SHARED_ELEMENTS_INPUTS_OUTPUTS_INDEXES, self._shared_elements_inputs_outputs_df)):

            if senescence_data_df is senescence_roots_data_df:
                senescence_data_df = senescence_data_df.copy()
                senescence_data_df.loc[:, 'organ'] = 'roots'

            tools.combine_dataframes_inplace(senescence_data_df, shared_inputs_outputs_indexes, shared_inputs_outputs_df)
