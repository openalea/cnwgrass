# -*- coding: latin-1 -*-

from growthwheat import converter, simulation, parameters
from fspmwheat import tools

"""
    fspmwheat.growthwheat_facade
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`fspmwheat.growthwheat_facade` is a facade of the model GrowthWheat.

    This module permits to initialize and run the model GrowthWheat from a :class:`MTG <openalea.mtg.mtg.MTG>`
    in a convenient and transparent way, wrapping all the internal complexity of the model, and dealing
    with all the tedious initialization and conversion processes.

    :copyright: Copyright 2014-2016 INRA-ECOSYS, see AUTHORS.
    :license: LICENSE for details.

    .. seealso:: Barillot et al. 2016.
"""

LEAF_LABELS = ['blade', 'sheath']

EMERGED_GROWING_ORGAN_LABELS = ['StemElement', 'LeafElement1']
ELEMENT_LABELS = ['StemElement', 'LeafElement1', 'HiddenElement']

SHARED_AXES_INPUTS_OUTPUTS_INDEXES = ['plant', 'axis']

SHARED_ORGANS_INPUTS_OUTPUTS_INDEXES = ['plant', 'axis', 'organ']

SHARED_HIDDENZONES_INPUTS_OUTPUTS_INDEXES = ['plant', 'axis', 'metamer']

SHARED_ELEMENTS_INPUTS_OUTPUTS_INDEXES = ['plant', 'axis', 'metamer', 'organ', 'element']


class GrowthWheatFacade(object):
    """
    The GrowthWheatFacade class permits to initialize, run the model GrowthWheat
    from a :class:`MTG <openalea.mtg.mtg.MTG>`, and update the MTG and the dataframes
    shared between all models.

    Use :meth:`run` to run the model.

    :Parameters:

"""

    def __init__(self, shared_mtg, delta_t,
                 model_hiddenzones_inputs_df,
                 model_elements_inputs_df,
                 model_roots_inputs_df,
                 model_axes_inputs_df,
                 shared_organs_inputs_outputs_df,
                 shared_hiddenzones_inputs_outputs_df,
                 shared_elements_inputs_outputs_df,
                 shared_axes_inputs_outputs_df,
                 update_parameters=None,
                 update_shared_df=True,
                 cnwheat_roots=True):

        """
        :param openalea.mtg.mtg.MTG shared_mtg: The MTG shared between all models.
        :param int delta_t: The delta between two runs, in seconds.
        :param pandas.DataFrame model_hiddenzones_inputs_df: the inputs of the model at hiddenzones scale.
        :param pandas.DataFrame model_elements_inputs_df: the inputs of the model at elements scale.
        :param pandas.DataFrame model_roots_inputs_df: the inputs of the model at roots scale.
        :param pandas.DataFrame model_axes_inputs_df: the inputs of the model at axes scale.
        :param pandas.DataFrame shared_organs_inputs_outputs_df: the dataframe of inputs and outputs at organs scale shared between all models.
        :param pandas.DataFrame shared_hiddenzones_inputs_outputs_df: the dataframe of inputs and outputs at hiddenzones scale shared between all models.
        :param pandas.DataFrame shared_elements_inputs_outputs_df: the dataframe of inputs and outputs at elements scale shared between all models.
        :param pandas.DataFrame shared_axes_inputs_outputs_df: the dataframe of inputs and outputs at axis scale shared between all models.
        :param dict update_parameters: A dictionary with the parameters to update, should have the form {'param1': value1, 'param2': value2, ...}.
        :param bool update_shared_df: If `True`  update the shared dataframes at init and at each run (unless stated otherwise)
        """
        if update_parameters is None:
            update_parameters = {}

        self._shared_mtg = shared_mtg  #: the MTG shared between all models

        self._simulation = simulation.Simulation(delta_t=delta_t, update_parameters=update_parameters, cnwheat_roots=cnwheat_roots)  #: the simulator to use to run the model

        all_growthwheat_inputs_dict = converter.from_dataframes(model_hiddenzones_inputs_df, model_elements_inputs_df, model_roots_inputs_df, model_axes_inputs_df)

        self._update_shared_MTG(all_growthwheat_inputs_dict['hiddenzone'], all_growthwheat_inputs_dict['elements'], all_growthwheat_inputs_dict['roots'], all_growthwheat_inputs_dict['axes'])

        self._shared_organs_inputs_outputs_df = shared_organs_inputs_outputs_df  #: the dataframe at organs scale shared between all models
        self._shared_hiddenzones_inputs_outputs_df = shared_hiddenzones_inputs_outputs_df  #: the dataframe at hiddenzones scale shared between all models
        self._shared_elements_inputs_outputs_df = shared_elements_inputs_outputs_df  #: the dataframe at elements scale shared between all models
        self._shared_axes_inputs_outputs_df = shared_axes_inputs_outputs_df  #: the dataframe at axis scale shared between all models
        self._update_shared_df = update_shared_df
        if self._update_shared_df:
            self._update_shared_dataframes(model_hiddenzones_inputs_df, model_elements_inputs_df, model_roots_inputs_df, model_axes_inputs_df)

    def run(self, postflowering_stages=False, update_shared_df=None):
        """
        Run the model and update the MTG and the dataframes shared between all models.
        :param bool postflowering_stages: if True the model will calculate root growth with the parameters calibrated for post flowering stages
        :param bool update_shared_df: if 'True', update the shared dataframes at this time step.
        """
        self._initialize_model()
        self._simulation.run(postflowering_stages)
        self._update_shared_MTG(self._simulation.outputs['hiddenzone'], self._simulation.outputs['elements'], self._simulation.outputs['roots'], self._simulation.outputs['axes'])

        if update_shared_df or (update_shared_df is None and self._update_shared_df):
            growthwheat_hiddenzones_outputs_df, growthwheat_elements_outputs_df, growthwheat_roots_outputs_df, growthwheat_axes_outputs_df = converter.to_dataframes(self._simulation.outputs)
            self._update_shared_dataframes(growthwheat_hiddenzones_outputs_df, growthwheat_elements_outputs_df, growthwheat_roots_outputs_df, growthwheat_axes_outputs_df)

    def _initialize_model(self):
        """
        Initialize the inputs of the model from the MTG shared between all models.
        """

        all_growthwheat_hiddenzones_inputs_dict = {}
        all_growthwheat_elements_inputs_dict = {}
        all_growthwheat_roots_inputs_dict = {}
        all_growthwheat_axes_inputs_dict = {}

        for mtg_plant_vid in self._shared_mtg.components_iter(self._shared_mtg.root):
            mtg_plant_index = int(self._shared_mtg.index(mtg_plant_vid))
            for mtg_axis_vid in self._shared_mtg.components_iter(mtg_plant_vid):
                mtg_axis_label = self._shared_mtg.label(mtg_axis_vid)
                if mtg_axis_label != 'MS':
                    continue

                mtg_axis_properties = self._shared_mtg.get_vertex_property(mtg_axis_vid)
                axis_id = (mtg_plant_index, mtg_axis_label)
                if set(mtg_axis_properties).issuperset(simulation.AXIS_INPUTS):
                    growthwheat_axis_inputs_dict = {}
                    for growthwheat_axis_input_name in simulation.AXIS_INPUTS:
                        growthwheat_axis_inputs_dict[growthwheat_axis_input_name] = mtg_axis_properties[growthwheat_axis_input_name]
                    all_growthwheat_axes_inputs_dict[axis_id] = growthwheat_axis_inputs_dict

                # Roots
                if 'roots' in mtg_axis_properties:
                    roots_id = (mtg_plant_index, mtg_axis_label, 'roots')
                    mtg_roots_properties = mtg_axis_properties['roots']

                    if set(mtg_roots_properties).issuperset(simulation.ROOT_INPUTS):
                        growthwheat_roots_inputs_dict = {}
                        for growthwheat_roots_input_name in simulation.ROOT_INPUTS:
                            growthwheat_roots_inputs_dict[growthwheat_roots_input_name] = mtg_roots_properties[growthwheat_roots_input_name]
                        all_growthwheat_roots_inputs_dict[roots_id] = growthwheat_roots_inputs_dict

                for mtg_metamer_vid in self._shared_mtg.components_iter(mtg_axis_vid):

                    mtg_metamer_index = int(self._shared_mtg.index(mtg_metamer_vid))

                    mtg_metamer_properties = self._shared_mtg.get_vertex_property(mtg_metamer_vid)
                    if 'hiddenzone' in mtg_metamer_properties:
                        hiddenzone_id = (mtg_plant_index, mtg_axis_label, mtg_metamer_index)
                        mtg_hiddenzone_properties = mtg_metamer_properties['hiddenzone']

                        if set(mtg_hiddenzone_properties).issuperset(simulation.HIDDENZONE_INPUTS):  # Initial values are set by elongwheat
                            growthwheat_hiddenzone_inputs_dict = {}
                            for growthwheat_hiddenzone_input_name in simulation.HIDDENZONE_INPUTS:
                                growthwheat_hiddenzone_inputs_dict[growthwheat_hiddenzone_input_name] = mtg_hiddenzone_properties[growthwheat_hiddenzone_input_name]
                            all_growthwheat_hiddenzones_inputs_dict[hiddenzone_id] = growthwheat_hiddenzone_inputs_dict

                        # We take only the elements of growing metamers ie. the ones with hiddenzones
                        for mtg_organ_vid in self._shared_mtg.components_iter(mtg_metamer_vid):
                            mtg_organ_label = self._shared_mtg.label(mtg_organ_vid)

                            for mtg_element_vid in self._shared_mtg.components_iter(mtg_organ_vid):
                                mtg_element_label = self._shared_mtg.label(mtg_element_vid)
                                element_id = (mtg_plant_index, mtg_axis_label, mtg_metamer_index, mtg_organ_label, mtg_element_label)
                                mtg_element_properties = self._shared_mtg.get_vertex_property(mtg_element_vid)

                                if mtg_element_label in ELEMENT_LABELS and \
                                        mtg_element_properties.get('length', 0) > 0:  # Note : ADEL puts length to positive value after updates even for HiddenElement.

                                    growthwheat_element_inputs_dict = {}

                                    # Exclude the HiddenElement appart from remobilization cases
                                    remobilisation = False
                                    if mtg_element_label == 'HiddenElement':
                                        if growthwheat_hiddenzone_inputs_dict['leaf_is_remobilizing'] or growthwheat_hiddenzone_inputs_dict['internode_is_remobilizing']:
                                            remobilisation = True
                                        else:
                                            continue

                                    for growthwheat_element_input_name in simulation.ELEMENT_INPUTS:
                                        mtg_element_input = mtg_element_properties.get(growthwheat_element_input_name)
                                        if mtg_element_input is None:
                                            mtg_element_input = parameters.OrganInit().__dict__[growthwheat_element_input_name]
                                        growthwheat_element_inputs_dict[growthwheat_element_input_name] = mtg_element_input
                                        if remobilisation:
                                            # Needed later on for CN Wheat calculation. TODO: Should it be in elongwheat_facade instead ? (MG)
                                            growthwheat_element_inputs_dict['green_area'] = mtg_element_properties.get('area')
                                    all_growthwheat_elements_inputs_dict[element_id] = growthwheat_element_inputs_dict

        self._simulation.initialize({'hiddenzone': all_growthwheat_hiddenzones_inputs_dict, 'elements': all_growthwheat_elements_inputs_dict,
                                     'roots': all_growthwheat_roots_inputs_dict, 'axes': all_growthwheat_axes_inputs_dict})

    def _update_shared_MTG(self, all_growthwheat_hiddenzones_data_dict, all_growthwheat_elements_data_dict, all_growthwheat_roots_data_dict, all_growthwheat_axes_data_dict):
        """
        Update the MTG shared between all models from the inputs or the outputs of the model.

        :param dict all_growthwheat_hiddenzones_data_dict: Growth-Wheat outputs at hidden zone scale
        :param dict all_growthwheat_elements_data_dict: Growth-Wheat outputs at element scale
        :param dict all_growthwheat_roots_data_dict: Growth-Wheat outputs at root scale
        :param dict all_growthwheat_axes_data_dict: Growth-Wheat outputs at axis scale
        """

        # add the properties if needed
        mtg_property_names = self._shared_mtg.property_names()
        if 'roots' not in mtg_property_names:
            self._shared_mtg.add_property('roots')
        for growthwheat_data_name in set(simulation.HIDDENZONE_INPUTS_OUTPUTS + simulation.ELEMENT_INPUTS_OUTPUTS):
            if growthwheat_data_name not in mtg_property_names:
                self._shared_mtg.add_property(growthwheat_data_name)

        # update the properties of the MTG
        for mtg_plant_vid in self._shared_mtg.components_iter(self._shared_mtg.root):
            mtg_plant_index = int(self._shared_mtg.index(mtg_plant_vid))
            for mtg_axis_vid in self._shared_mtg.components_iter(mtg_plant_vid):
                mtg_axis_label = self._shared_mtg.label(mtg_axis_vid)
                axis_id = (mtg_plant_index, mtg_axis_label)

                if mtg_axis_label != 'MS':
                    continue

                growthwheat_axis_data_dict = all_growthwheat_axes_data_dict[axis_id]
                for axis_data_name, axis_data_value in growthwheat_axis_data_dict.items():
                    self._shared_mtg.property(axis_data_name)[mtg_axis_vid] = axis_data_value

                roots_id = (mtg_plant_index, mtg_axis_label, 'roots')
                if roots_id in all_growthwheat_roots_data_dict:
                    growthwheat_roots_data_dict = all_growthwheat_roots_data_dict[roots_id]
                    mtg_axis_properties = self._shared_mtg.get_vertex_property(mtg_axis_vid)
                    if 'roots' not in mtg_axis_properties:
                        self._shared_mtg.property('roots')[mtg_axis_vid] = {}
                    for roots_data_name, roots_data_value in growthwheat_roots_data_dict.items():
                        self._shared_mtg.property('roots')[mtg_axis_vid][roots_data_name] = roots_data_value

                #: Metamer scale
                for mtg_metamer_vid in self._shared_mtg.components_iter(mtg_axis_vid):
                    mtg_metamer_index = int(self._shared_mtg.index(mtg_metamer_vid))
                    hiddenzone_id = (mtg_plant_index, mtg_axis_label, mtg_metamer_index)
                    if hiddenzone_id in all_growthwheat_hiddenzones_data_dict:
                        growthwheat_hiddenzone_data_dict = all_growthwheat_hiddenzones_data_dict[hiddenzone_id]
                        mtg_metamer_properties = self._shared_mtg.get_vertex_property(mtg_metamer_vid)
                        if 'hiddenzone' not in mtg_metamer_properties:  # MG : when is it used ?
                            self._shared_mtg.property('hiddenzone')[mtg_metamer_vid] = {}
                        for hiddenzone_data_name, hiddenzone_data_value in growthwheat_hiddenzone_data_dict.items():
                            self._shared_mtg.property('hiddenzone')[mtg_metamer_vid][hiddenzone_data_name] = hiddenzone_data_value

                    elif 'hiddenzone' in self._shared_mtg.get_vertex_property(mtg_metamer_vid):
                        # remove the 'hiddenzone' property from this metamer
                        del self._shared_mtg.property('hiddenzone')[mtg_metamer_vid]

                    #: Organ scale
                    for mtg_organ_vid in self._shared_mtg.components_iter(mtg_metamer_vid):
                        mtg_organ_label = self._shared_mtg.label(mtg_organ_vid)

                        #: Element scale
                        for mtg_element_vid in self._shared_mtg.components_iter(mtg_organ_vid):
                            mtg_element_label = self._shared_mtg.label(mtg_element_vid)
                            element_id = (mtg_plant_index, mtg_axis_label, mtg_metamer_index, mtg_organ_label, mtg_element_label)

                            if element_id in all_growthwheat_elements_data_dict:
                                growthwheat_element_data_dict = all_growthwheat_elements_data_dict[element_id]

                                for element_data_name, element_data_value in growthwheat_element_data_dict.items():
                                    self._shared_mtg.property(element_data_name)[mtg_element_vid] = element_data_value

    def _update_shared_dataframes(self, growthwheat_hiddenzones_data_df, growthwheat_elements_data_df, growthwheat_roots_data_df, growthwheat_axes_data_df):
        """
        Update the dataframes shared between all models from the inputs dataframes or the outputs dataframes of the model.

        :param pandas.DataFrame growthwheat_hiddenzones_data_df: Growth-Wheat shared dataframe at hidden zone scale
        :param pandas.DataFrame growthwheat_elements_data_df: Growth-Wheat shared dataframe at element scale
        :param pandas.DataFrame growthwheat_roots_data_df: Growth-Wheat shared dataframe at roots scale
        :param pandas.DataFrame growthwheat_axes_data_df: Growth-Wheat shared dataframe at axis scale
        """

        for growthwheat_data_df, \
            shared_inputs_outputs_indexes, \
            shared_inputs_outputs_df in ((growthwheat_hiddenzones_data_df, SHARED_HIDDENZONES_INPUTS_OUTPUTS_INDEXES, self._shared_hiddenzones_inputs_outputs_df),
                                         (growthwheat_elements_data_df, SHARED_ELEMENTS_INPUTS_OUTPUTS_INDEXES, self._shared_elements_inputs_outputs_df),
                                         (growthwheat_roots_data_df, SHARED_ORGANS_INPUTS_OUTPUTS_INDEXES, self._shared_organs_inputs_outputs_df),
                                         (growthwheat_axes_data_df, SHARED_AXES_INPUTS_OUTPUTS_INDEXES, self._shared_axes_inputs_outputs_df)):

            if growthwheat_data_df is growthwheat_roots_data_df:
                growthwheat_data_df = growthwheat_data_df.copy()
                growthwheat_data_df.loc[:, 'organ'] = 'roots'

            tools.combine_dataframes_inplace(growthwheat_data_df, shared_inputs_outputs_indexes, shared_inputs_outputs_df)
