# -*- coding: latin-1 -*-

from openalea.cnwgrass.respiration import model as respiration_model

from openalea.cnwgrass.cnmetabolism import model as cnmetabolism_model, simulation as cnmetabolism_simulation, \
    converter as cnmetabolism_converter, postprocessing as cnmetabolism_postprocessing
from openalea.cnwgrass.cnmetabolism import parameters as cnmetabolism_parameters

from openalea.cnwgrass.integration import tools

import numpy as np
import math

"""
    integration.cnmetabolism_facade
    ~~~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`integration.cnmetabolism_facade` is a facade of the model CN-Metabolism.

    This module permits to initialize and run the model CN-Metabolism from a :class:`MTG <openalea.mtg.mtg.MTG>`
    in a convenient and transparent way, wrapping all the internal complexity of the model, and dealing
    with all the tedious initialization and conversion processes.

"""

#: the mapping of CN-Metabolism organ classes to the attributes in axis and phytomer which represent an organ
CNMETABOLISM_ATTRIBUTES_MAPPING = {cnmetabolism_model.Internode: 'internode', cnmetabolism_model.Lamina: 'lamina',
                              cnmetabolism_model.Sheath: 'sheath', cnmetabolism_model.Peduncle: 'peduncle', cnmetabolism_model.Chaff: 'chaff',
                              cnmetabolism_model.Roots: 'roots', cnmetabolism_model.Grains: 'grains', cnmetabolism_model.Phloem: 'phloem',
                              cnmetabolism_model.HiddenZone: 'hiddenzone', cnmetabolism_model.Endosperm: 'endosperm'}

#: the mapping of organs (which belong to an axis) labels in MTG to organ classes in CN-Metabolism
MTG_TO_CNMETABOLISM_AXES_ORGANS_MAPPING = {'grains': cnmetabolism_model.Grains, 'phloem': cnmetabolism_model.Phloem, 'roots': cnmetabolism_model.Roots, 'endosperm': cnmetabolism_model.Endosperm}

#: the mapping of organs (which belong to a phytomer) labels in MTG to organ classes in CN-Metabolism
MTG_TO_CNMETABOLISM_PHYTOMERS_ORGANS_MAPPING = {'internode': cnmetabolism_model.Internode, 'blade': cnmetabolism_model.Lamina, 'sheath': cnmetabolism_model.Sheath, 'peduncle': cnmetabolism_model.Peduncle,
                                           'ear': cnmetabolism_model.Chaff, 'hiddenzone': cnmetabolism_model.HiddenZone}

# # the mapping of CN-Metabolism photosynthetic organs to CN-Metabolism photosynthetic organ elements
CNMETABOLISM_ORGANS_TO_ELEMENTS_MAPPING = {cnmetabolism_model.Internode: cnmetabolism_model.InternodeElement, cnmetabolism_model.Lamina: cnmetabolism_model.LaminaElement, cnmetabolism_model.Sheath: cnmetabolism_model.SheathElement,
                                      cnmetabolism_model.Peduncle: cnmetabolism_model.PeduncleElement, cnmetabolism_model.Chaff: cnmetabolism_model.ChaffElement}

#: the parameters and variables which define the state of a CN-Metabolism population
POPULATION_STATE_VARIABLE = set(cnmetabolism_simulation.Simulation.PLANTS_STATE + cnmetabolism_simulation.Simulation.AXES_STATE +
                                cnmetabolism_simulation.Simulation.PHYTOMERS_STATE + cnmetabolism_simulation.Simulation.ORGANS_STATE +
                                cnmetabolism_simulation.Simulation.HIDDENZONE_STATE + cnmetabolism_simulation.Simulation.ELEMENTS_STATE)

#: all the variables of a CN-Metabolism population computed during a run step of the simulation
POPULATION_RUN_VARIABLES = set(cnmetabolism_simulation.Simulation.PLANTS_RUN_VARIABLES + cnmetabolism_simulation.Simulation.AXES_RUN_VARIABLES +
                               cnmetabolism_simulation.Simulation.PHYTOMERS_RUN_VARIABLES + cnmetabolism_simulation.Simulation.ORGANS_RUN_VARIABLES +
                               cnmetabolism_simulation.Simulation.HIDDENZONE_RUN_VARIABLES + cnmetabolism_simulation.Simulation.ELEMENTS_RUN_VARIABLES)

#: all the variables to be stored in the MTG
MTG_RUN_VARIABLES = set(list(POPULATION_RUN_VARIABLES) + cnmetabolism_simulation.Simulation.SOILS_RUN_VARIABLES)

# number of seconds in 1 hour
HOUR_TO_SECOND_CONVERSION_FACTOR = 3600


class CNMetabolismFacade(object):
    """
    The CNMetabolismFacade class permits to initialize, run the model CN-Metabolism
    from a :class:`MTG <openalea.mtg.mtg.MTG>`, and update the MTG and the dataframes
    shared between all models.

    Use :meth:`run` to run the model.

    """

    def __init__(self, shared_mtg, delta_t, culm_density, update_parameters,
                 model_axes_inputs_df,
                 model_organs_inputs_df,
                 model_hiddenzones_inputs_df,
                 model_elements_inputs_df,
                 model_soils_inputs_df,
                 shared_axes_inputs_outputs_df,
                 shared_organs_inputs_outputs_df,
                 shared_hiddenzones_inputs_outputs_df,
                 shared_elements_inputs_outputs_df,
                 shared_soils_inputs_outputs_df,
                 tillers_replications=None,
                 external_soil_model=False,
                 update_shared_df=True):
        """
        :param openalea.mtg.mtg.MTG shared_mtg: The MTG shared between all models.
        :param int delta_t: The delta between two runs, in seconds.
        :param dict culm_density: The density of culm. One key per plant.
        :param dict update_parameters: A dictionary with the parameters to update, should have the form {'Organ_label1': {'param1': value1, 'param2': value2}, ...}.
        :param pandas.DataFrame model_axes_inputs_df: the inputs of the model at axis scale.
        :param pandas.DataFrame model_organs_inputs_df: the inputs of the model at organs scale.
        :param pandas.DataFrame model_hiddenzones_inputs_df: the inputs of the model at hiddenzones scale.
        :param pandas.DataFrame model_elements_inputs_df: the inputs of the model at elements scale.
        :param pandas.DataFrame or None model_soils_inputs_df: the inputs of the model at soils scale.
        :param pandas.DataFrame shared_axes_inputs_outputs_df: the dataframe of inputs and outputs at axes scale shared between all models.
        :param pandas.DataFrame shared_organs_inputs_outputs_df: the dataframe of inputs and outputs at organs scale shared between all models.
        :param pandas.DataFrame shared_hiddenzones_inputs_outputs_df: the dataframe of inputs and outputs at hiddenzones scale shared between all models.
        :param pandas.DataFrame shared_elements_inputs_outputs_df: the dataframe of inputs and outputs at elements scale shared between all models.
        :param pandas.DataFrame shared_soils_inputs_outputs_df: the dataframe of inputs and outputs at soils scale shared between all models.
        :param dict [str, float] tillers_replications: a dictionary with tiller id as key, and weight of replication as value.
        :param bool external_soil_model: whether an external soil model is coupled to cnmetabolism. If True, cnmetabolism will skip calculations made in soil and uptake N by roots
        :param bool update_shared_df: If `True`  update the shared dataframes at init and at each run (unless stated otherwise)

        """

        self._shared_mtg = shared_mtg  #: the MTG shared between all models
        self.tillers_replications = tillers_replications
        self.external_soil_model = external_soil_model

        self._simulation = cnmetabolism_simulation.Simulation(respiration_model=respiration_model, delta_t=delta_t, culm_density=culm_density, external_soil_model=external_soil_model)

        self.population, self.soils = cnmetabolism_converter.from_dataframes(model_axes_inputs_df, model_organs_inputs_df, model_hiddenzones_inputs_df, model_elements_inputs_df, model_soils_inputs_df)

        self._update_parameters = update_parameters

        self._simulation.initialize(self.population, self.soils)

        self._update_shared_MTG()

        self._shared_axes_inputs_outputs_df = shared_axes_inputs_outputs_df  #: the dataframe at axes scale shared between all models
        self._shared_organs_inputs_outputs_df = shared_organs_inputs_outputs_df  #: the dataframe at organs scale shared between all models
        self._shared_hiddenzones_inputs_outputs_df = shared_hiddenzones_inputs_outputs_df  #: the dataframe at hiddenzones scale shared between all models
        self._shared_elements_inputs_outputs_df = shared_elements_inputs_outputs_df  #: the dataframe at elements scale shared between all models
        self._shared_soils_inputs_outputs_df = shared_soils_inputs_outputs_df  #: the dataframe at soils scale shared between all models
        self._update_shared_df = update_shared_df
        if self._update_shared_df:
            self._update_shared_dataframes(cnmetabolism_organs_data_df=model_organs_inputs_df,
                                           cnmetabolism_hiddenzones_data_df=model_hiddenzones_inputs_df,
                                           cnmetabolism_elements_data_df=model_elements_inputs_df,
                                           cnmetabolism_soils_data_df=model_soils_inputs_df)

    def run(self, Tair, Tsoil, update_shared_df=None):
        """
        Run the model and update the MTG and the dataframes shared between all models.

        :param update_shared_df:
        :param float Tair: air temperature (°C)
        :param float Tsoil: soil temperature (°C)
        :param bool update_shared_df: if 'True', update the shared dataframes at this time step.
        """

        self._initialize_model(Tair=Tair, Tsoil=Tsoil)
        self._simulation.run()
        self._update_shared_MTG()

        if update_shared_df or (update_shared_df is None and self._update_shared_df):
            _, cnmetabolism_axes_inputs_outputs_df, _, cnmetabolism_organs_inputs_outputs_df, cnmetabolism_hiddenzones_inputs_outputs_df, cnmetabolism_elements_inputs_outputs_df, cnmetabolism_soils_inputs_outputs_df = \
                cnmetabolism_converter.to_dataframes(self._simulation.population, self._simulation.soils)
            self._update_shared_dataframes(cnmetabolism_axes_data_df=cnmetabolism_axes_inputs_outputs_df,
                                           cnmetabolism_organs_data_df=cnmetabolism_organs_inputs_outputs_df,
                                           cnmetabolism_hiddenzones_data_df=cnmetabolism_hiddenzones_inputs_outputs_df,
                                           cnmetabolism_elements_data_df=cnmetabolism_elements_inputs_outputs_df,
                                           cnmetabolism_soils_data_df=cnmetabolism_soils_inputs_outputs_df)

    @staticmethod
    def postprocessing(axes_outputs_df, organs_outputs_df, hiddenzone_outputs_df, elements_outputs_df, soils_outputs_df, delta_t):
        """
        Run the postprocessing.

        :param pandas.DataFrame axes_outputs_df: the outputs of the model at axis scale.
        :param pandas.DataFrame organs_outputs_df: the outputs of the model at organ scale.
        :param pandas.DataFrame hiddenzone_outputs_df: the outputs of the model at hiddenzone scale.
        :param pandas.DataFrame elements_outputs_df: the outputs of the model at element scale.
        :param pandas.DataFrame soils_outputs_df: the outputs of the model at element scale.
        :param int delta_t: The delta between two runs, in seconds.

    :return: post-processing for each scale:
            * plant (see :attr:`PLANTS_RUN_POSTPROCESSING_VARIABLES`)
            * axis (see :attr:`AXES_RUN_POSTPROCESSING_VARIABLES`)
            * metamer (see :attr:`PHYTOMERS_RUN_POSTPROCESSING_VARIABLES`)
            * organ (see :attr:`ORGANS_RUN_POSTPROCESSING_VARIABLES`)
            * hidden zone (see :attr:`HIDDENZONE_RUN_POSTPROCESSING_VARIABLES`)
            * element (see :attr:`ELEMENTS_RUN_POSTPROCESSING_VARIABLES`)
            * and soil (see :attr:`SOILS_RUN_POSTPROCESSING_VARIABLES`)
        depending on the dataframes given as argument.
        For example, if user passes only dataframes `plants_df`, `axes_df` and `metamers_df`,
        then only post-processing dataframes of plants, axes and metamers are returned.

    :rtype: dict ['scale' : pandas.DataFrame]
        """

        (_, _, organs_postprocessing_df,
         elements_postprocessing_df,
         hiddenzones_postprocessing_df,
         axes_postprocessing_df,
         soils_postprocessing_df) = cnmetabolism_postprocessing.postprocessing(axes_df=axes_outputs_df, hiddenzones_df=hiddenzone_outputs_df,
                                                                          organs_df=organs_outputs_df, elements_df=elements_outputs_df,
                                                                          soils_df=soils_outputs_df, delta_t=delta_t)
        return {
            'axes': axes_postprocessing_df,
            'organs': organs_postprocessing_df,
            'hiddenzones': hiddenzones_postprocessing_df,
            'elements': elements_postprocessing_df,
            'soils': soils_postprocessing_df
        }

    @staticmethod
    def graphs(axes_postprocessing_df, hiddenzones_postprocessing_df, organs_postprocessing_df, elements_postprocessing_df, soils_postprocessing_df, meteo_data=None, graphs_dirpath='.'):
        """
        Generate the graphs and save them into `graphs_dirpath`.

        :param pandas.DataFrame axes_postprocessing_df: CN-Metabolism outputs at axis scale
        :param pandas.DataFrame hiddenzones_postprocessing_df: CN-Metabolism outputs at hidden zone scale
        :param pandas.DataFrame organs_postprocessing_df: CN-Metabolism outputs at organ scale
        :param pandas.DataFrame elements_postprocessing_df: CN-Metabolism outputs at element scale
        :param pandas.DataFrame soils_postprocessing_df: CN-Metabolism outputs at soil scale
        :param pandas.DataFrame meteo_data: the meteo dataframe having the mapping between t (hours) and calendar dates
        :param str graphs_dirpath: the path of the directory to save the generated graphs

        """
        cnmetabolism_postprocessing.generate_graphs(axes_df=axes_postprocessing_df,
                                               hiddenzones_df=hiddenzones_postprocessing_df,
                                               organs_df=organs_postprocessing_df,
                                               elements_df=elements_postprocessing_df,
                                               soils_df=soils_postprocessing_df,
                                               meteo_data=meteo_data,
                                               graphs_dirpath=graphs_dirpath)

    def _initialize_model(self, Tair, Tsoil):
        """
        Initialize the inputs of the model from the MTG shared between all models and the soils.

        :param float Tair: air temperature (°C)
        :param float Tsoil: soil temperature (°C)
        """

        # Convert number of replications per tiller into number of replications per cohort
        cohorts_replications = {}
        if self.tillers_replications is not None:
            for tiller_id, replication_weight in self.tillers_replications.items():
                try:
                    tiller_rank = int(tiller_id[1:])
                except ValueError:
                    continue
                cohorts_replications[tiller_rank + 3] = replication_weight

        self.population = cnmetabolism_model.Population()

        # traverse the MTG recursively from top
        for mtg_plant_vid in self._shared_mtg.components_iter(self._shared_mtg.root):
            mtg_plant_index = int(self._shared_mtg.index(mtg_plant_vid))
            # create a new plant
            cnmetabolism_plant = cnmetabolism_model.Plant(mtg_plant_index)
            is_valid_plant = False

            for mtg_axis_vid in self._shared_mtg.components_iter(mtg_plant_vid):
                mtg_axis_label = self._shared_mtg.label(mtg_axis_vid)

                #: Hack to deal with tillering cases : TEMPORARY
                if mtg_axis_label != 'MS':
                    try:
                        tiller_rank = int(mtg_axis_label[1:])
                        cnmetabolism_plant.cohorts.append(tiller_rank + 3)
                        continue
                    except ValueError:
                        continue

                #: Main Stem
                # create a new axis
                cnmetabolism_axis = cnmetabolism_model.Axis(mtg_axis_label)
                mtg_axis_properties = self._shared_mtg.get_vertex_property(mtg_axis_vid)
                cnmetabolism_axis_data_dict = {}
                for cnmetabolism_axis_data_name in cnmetabolism_simulation.Simulation.AXES_STATE:
                    cnmetabolism_axis_data_dict[cnmetabolism_axis_data_name] = mtg_axis_properties[cnmetabolism_axis_data_name]
                cnmetabolism_axis.__dict__.update(cnmetabolism_axis_data_dict)
                is_valid_axis = True
                for cnmetabolism_organ_class in (cnmetabolism_model.Roots, cnmetabolism_model.Phloem, cnmetabolism_model.Grains, cnmetabolism_model.Endosperm):
                    mtg_organ_label = cnmetabolism_converter.CNMETABOLISM_CLASSES_TO_DATAFRAME_ORGANS_MAPPING[cnmetabolism_organ_class]
                    # create a new organ
                    cnmetabolism_organ = cnmetabolism_organ_class(mtg_organ_label)
                    if mtg_organ_label in mtg_axis_properties:
                        mtg_organ_properties = mtg_axis_properties[mtg_organ_label]
                        cnmetabolism_organ_data_names = set(cnmetabolism_simulation.Simulation.ORGANS_STATE).intersection(cnmetabolism_organ.__dict__)
                        if set(mtg_organ_properties).issuperset(cnmetabolism_organ_data_names):
                            cnmetabolism_organ_data_dict = {}
                            for cnmetabolism_organ_data_name in cnmetabolism_organ_data_names:
                                cnmetabolism_organ_data_dict[cnmetabolism_organ_data_name] = mtg_organ_properties[cnmetabolism_organ_data_name]

                                # Debug: Tell if missing input variable
                                if math.isnan(mtg_organ_properties[cnmetabolism_organ_data_name]) or mtg_organ_properties[cnmetabolism_organ_data_name] is None:
                                    print('Missing variable', cnmetabolism_organ_data_name, 'for vertex id', mtg_axis_vid, 'which is', mtg_organ_label)

                            cnmetabolism_organ.__dict__.update(cnmetabolism_organ_data_dict)
                            if mtg_organ_label == 'roots' and self.external_soil_model:
                                cnmetabolism_organ.Uptake_Nitrates = mtg_organ_properties['Uptake_Nitrates']
                                cnmetabolism_organ.HATS_LATS = mtg_organ_properties['HATS_LATS']

                            # Update parameters if specified
                            if mtg_organ_label in self._update_parameters:
                                cnmetabolism_organ.PARAMETERS.__dict__.update(self._update_parameters[mtg_organ_label])

                            cnmetabolism_organ.initialize()
                            # add the new organ to current axis
                            setattr(cnmetabolism_axis, mtg_organ_label, cnmetabolism_organ)

                        elif cnmetabolism_organ_class is not cnmetabolism_model.Grains:
                            is_valid_axis = False
                            break

                    # For the 1st instantiation of the Grains class during a simulation covering vegetative and reproductive stages
                    elif cnmetabolism_organ_class is cnmetabolism_model.Grains:
                        if mtg_axis_properties['status'] != 'reproductive':
                            continue
                        # grains = cnmetabolism_model.Grains(cnmetabolism_converter.CNMETABOLISM_CLASSES_TO_DATAFRAME_ORGANS_MAPPING[cnmetabolism_model.Grains])
                        # grains.initialize()
                        # setattr(cnmetabolism_axis, cnmetabolism_converter.CNMETABOLISM_CLASSES_TO_DATAFRAME_ORGANS_MAPPING[cnmetabolism_model.Grains], grains)

                    elif cnmetabolism_organ_class is cnmetabolism_model.Endosperm:
                        continue

                    else:
                        is_valid_axis = False
                        print('Invalid axis because of {}'.format(cnmetabolism_organ_class))
                        break

                if not is_valid_axis:
                    continue

                has_valid_phytomer = False
                for mtg_metamer_vid in self._shared_mtg.components_iter(mtg_axis_vid):
                    mtg_metamer_index = int(self._shared_mtg.index(mtg_metamer_vid))

                    # create a new phytomer
                    cnmetabolism_phytomer = cnmetabolism_model.Phytomer(mtg_metamer_index, cohorts=cnmetabolism_plant.cohorts, cohorts_replications=cohorts_replications)  #: Hack to treat tillering cases :TEMPORARY

                    mtg_hiddenzone_label = cnmetabolism_converter.CNMETABOLISM_CLASSES_TO_DATAFRAME_ORGANS_MAPPING[cnmetabolism_model.HiddenZone]
                    mtg_metamer_properties = self._shared_mtg.get_vertex_property(mtg_metamer_vid)

                    if mtg_hiddenzone_label in mtg_metamer_properties:
                        mtg_hiddenzone_properties = mtg_metamer_properties[mtg_hiddenzone_label]

                        if set(mtg_hiddenzone_properties).issuperset(cnmetabolism_simulation.Simulation.HIDDENZONE_STATE) and not mtg_hiddenzone_properties['is_over']:
                            has_valid_hiddenzone = True
                            cnmetabolism_hiddenzone_data_dict = {}
                            for cnmetabolism_hiddenzone_data_name in cnmetabolism_simulation.Simulation.HIDDENZONE_STATE:
                                cnmetabolism_hiddenzone_data_dict[cnmetabolism_hiddenzone_data_name] = mtg_hiddenzone_properties[cnmetabolism_hiddenzone_data_name]

                            # create a new hiddenzone
                            cnmetabolism_hiddenzone = cnmetabolism_model.HiddenZone(mtg_hiddenzone_label, cohorts=cnmetabolism_plant.cohorts, cohorts_replications=cohorts_replications, index=cnmetabolism_phytomer.index,
                                                                          **cnmetabolism_hiddenzone_data_dict)

                            # Update parameters if specified
                            if mtg_hiddenzone_label in self._update_parameters:
                                cnmetabolism_hiddenzone.PARAMETERS.__dict__.update(self._update_parameters[mtg_hiddenzone_label])

                            cnmetabolism_hiddenzone.initialize()
                            # add the new hiddenzone to current phytomer
                            setattr(cnmetabolism_phytomer, mtg_hiddenzone_label, cnmetabolism_hiddenzone)
                        else:
                            has_valid_hiddenzone = False
                    else:
                        has_valid_hiddenzone = False

                    has_valid_organ = False
                    for mtg_organ_vid in self._shared_mtg.components_iter(mtg_metamer_vid):
                        mtg_organ_label = self._shared_mtg.label(mtg_organ_vid)
                        if mtg_organ_label not in MTG_TO_CNMETABOLISM_PHYTOMERS_ORGANS_MAPPING or self._shared_mtg.get_vertex_property(mtg_organ_vid)['length'] == 0:
                            continue

                        # create a new organ
                        cnmetabolism_organ_class = MTG_TO_CNMETABOLISM_PHYTOMERS_ORGANS_MAPPING[mtg_organ_label]
                        cnmetabolism_organ = cnmetabolism_organ_class(mtg_organ_label)

                        # Update parameters if specified
                        if 'PhotosyntheticOrgan' in self._update_parameters:
                            cnmetabolism_organ.PARAMETERS.__dict__.update(self._update_parameters['PhotosyntheticOrgan'])

                        cnmetabolism_organ.initialize()
                        has_valid_element = False

                        # Create a new element
                        for mtg_element_vid in self._shared_mtg.components_iter(mtg_organ_vid):
                            mtg_element_properties = self._shared_mtg.get_vertex_property(mtg_element_vid)
                            mtg_element_label = self._shared_mtg.label(mtg_element_vid)
                            if mtg_element_label not in cnmetabolism_converter.DATAFRAME_TO_CNMETABOLISM_ELEMENTS_NAMES_MAPPING \
                                    or (self._shared_mtg.get_vertex_property(mtg_element_vid)['length'] == 0) \
                                    or (self._shared_mtg.get_vertex_property(mtg_element_vid).get('mstruct', 0) == 0) \
                                    or ((mtg_element_label == 'HiddenElement') and (self._shared_mtg.get_vertex_property(mtg_element_vid).get('is_growing', True))) \
                                    or (self._shared_mtg.get_vertex_property(mtg_element_vid).get('green_area', 0) <= 0.25E-6):
                                continue  # TODO: Check that we are not taking out some relevant cases with the condition on mstruct == 0

                            has_valid_element = True
                            cnmetabolism_element_data_dict = {}
                            for cnmetabolism_element_data_name in cnmetabolism_simulation.Simulation.ELEMENTS_STATE:
                                mtg_element_data_value = mtg_element_properties.get(cnmetabolism_element_data_name)
                                # In case the value is None, or the property is not even defined, we take default value from InitCompartment
                                if mtg_element_data_value is None or np.isnan(mtg_element_data_value):
                                    if cnmetabolism_element_data_name == 'Ts':
                                        mtg_element_data_value = Tair
                                    else:
                                        mtg_element_data_value = cnmetabolism_parameters.PhotosyntheticOrganElementInitCompartments().__dict__[cnmetabolism_element_data_name]
                                cnmetabolism_element_data_dict[cnmetabolism_element_data_name] = mtg_element_data_value
                            cnmetabolism_element = CNMETABOLISM_ORGANS_TO_ELEMENTS_MAPPING[cnmetabolism_organ_class](mtg_element_label, cohorts=cnmetabolism_plant.cohorts, cohorts_replications=cohorts_replications,
                                                                                                      index=cnmetabolism_phytomer.index, **cnmetabolism_element_data_dict)
                            # Add parameters from organ scale
                            cnmetabolism_element.PARAMETERS.__dict__.update(cnmetabolism_organ.PARAMETERS.__dict__)

                            # add the new element to current organ
                            setattr(cnmetabolism_organ, cnmetabolism_converter.DATAFRAME_TO_CNMETABOLISM_ELEMENTS_NAMES_MAPPING[mtg_element_label], cnmetabolism_element)

                        if has_valid_element:
                            has_valid_organ = True
                            setattr(cnmetabolism_phytomer, CNMETABOLISM_ATTRIBUTES_MAPPING[cnmetabolism_organ_class], cnmetabolism_organ)

                    if has_valid_organ or has_valid_hiddenzone:
                        cnmetabolism_axis.phytomers.append(cnmetabolism_phytomer)
                        has_valid_phytomer = True

                if not has_valid_phytomer:
                    is_valid_axis = False

                if is_valid_axis:
                    cnmetabolism_plant.axes.append(cnmetabolism_axis)
                    is_valid_plant = True

            if is_valid_plant:
                self.population.plants.append(cnmetabolism_plant)

        self._simulation.initialize(self.population, self.soils, Tsoil=Tsoil)

    def _update_shared_MTG(self):
        """
        Update the MTG shared between all models from the population of CN-Metabolism.
        """
        # add the missing properties
        mtg_property_names = self._shared_mtg.property_names()
        for cnmetabolism_data_name in MTG_RUN_VARIABLES:
            if cnmetabolism_data_name not in mtg_property_names:
                self._shared_mtg.add_property(cnmetabolism_data_name)
        for cnmetabolism_organ_label in list(MTG_TO_CNMETABOLISM_AXES_ORGANS_MAPPING.keys()) + ['soil'] + [cnmetabolism_converter.CNMETABOLISM_CLASSES_TO_DATAFRAME_ORGANS_MAPPING[cnmetabolism_model.HiddenZone]]:
            if cnmetabolism_organ_label not in mtg_property_names:
                self._shared_mtg.add_property(cnmetabolism_organ_label)

        mtg_plants_iterator = self._shared_mtg.components_iter(self._shared_mtg.root)
        # traverse CN-Metabolism population from top
        for cnmetabolism_plant in self.population.plants:
            cnmetabolism_plant_index = cnmetabolism_plant.index
            while True:
                mtg_plant_vid = next(mtg_plants_iterator)
                if int(self._shared_mtg.index(mtg_plant_vid)) == cnmetabolism_plant_index:
                    break
            mtg_axes_iterator = self._shared_mtg.components_iter(mtg_plant_vid)
            for cnmetabolism_axis in cnmetabolism_plant.axes:
                cnmetabolism_axis_label = cnmetabolism_axis.label
                while True:
                    mtg_axis_vid = next(mtg_axes_iterator)
                    if self._shared_mtg.label(mtg_axis_vid) == cnmetabolism_axis_label:
                        break

                cnmetabolism_axis_property_names = [property_name for property_name in cnmetabolism_simulation.Simulation.AXES_RUN_VARIABLES if hasattr(cnmetabolism_axis, property_name)]
                for cnmetabolism_axis_property_name in cnmetabolism_axis_property_names:
                    cnmetabolism_axis_property_value = getattr(cnmetabolism_axis, cnmetabolism_axis_property_name)
                    self._shared_mtg.property(cnmetabolism_axis_property_name)[mtg_axis_vid] = cnmetabolism_axis_property_value

                for mtg_organ_label in MTG_TO_CNMETABOLISM_AXES_ORGANS_MAPPING.keys():
                    cnmetabolism_organ = getattr(cnmetabolism_axis, mtg_organ_label)
                    if cnmetabolism_organ is None:
                        continue
                    elif mtg_organ_label not in self._shared_mtg.get_vertex_property(mtg_axis_vid) and cnmetabolism_organ is not None:
                        # Add a property describing the organ to the current axis of the MTG
                        self._shared_mtg.property(mtg_organ_label)[mtg_axis_vid] = {}
                    # Update the property describing the organ of the current axis in the MTG
                    mtg_organ_properties = self._shared_mtg.get_vertex_property(mtg_axis_vid)[mtg_organ_label]
                    for cnmetabolism_property_name in cnmetabolism_simulation.Simulation.ORGANS_RUN_VARIABLES:
                        if hasattr(cnmetabolism_organ, cnmetabolism_property_name):
                            mtg_organ_properties[cnmetabolism_property_name] = getattr(cnmetabolism_organ, cnmetabolism_property_name)
                mtg_metamers_iterator = self._shared_mtg.components_iter(mtg_axis_vid)
                for cnmetabolism_phytomer in cnmetabolism_axis.phytomers:
                    cnmetabolism_phytomer_index = cnmetabolism_phytomer.index
                    while True:
                        mtg_metamer_vid = next(mtg_metamers_iterator)
                        if int(self._shared_mtg.index(mtg_metamer_vid)) == cnmetabolism_phytomer_index:
                            break
                    if cnmetabolism_phytomer.hiddenzone is not None:
                        mtg_hiddenzone_label = cnmetabolism_converter.CNMETABOLISM_CLASSES_TO_DATAFRAME_ORGANS_MAPPING[cnmetabolism_model.HiddenZone]
                        if mtg_hiddenzone_label not in self._shared_mtg.get_vertex_property(mtg_metamer_vid):
                            # Add a property describing the hiddenzone to the current metamer of the MTG
                            self._shared_mtg.property(mtg_hiddenzone_label)[mtg_metamer_vid] = {}
                        # Update the property describing the hiddenzone of the current metamer in the MTG
                        mtg_hiddenzone_properties = self._shared_mtg.get_vertex_property(mtg_metamer_vid)[mtg_hiddenzone_label]
                        for cnmetabolism_property_name in cnmetabolism_simulation.Simulation.HIDDENZONE_RUN_VARIABLES:
                            if hasattr(cnmetabolism_phytomer.hiddenzone, cnmetabolism_property_name):
                                mtg_hiddenzone_properties[cnmetabolism_property_name] = getattr(cnmetabolism_phytomer.hiddenzone, cnmetabolism_property_name)
                    for mtg_organ_vid in self._shared_mtg.components_iter(mtg_metamer_vid):
                        mtg_organ_label = self._shared_mtg.label(mtg_organ_vid)
                        if mtg_organ_label not in MTG_TO_CNMETABOLISM_PHYTOMERS_ORGANS_MAPPING:
                            continue
                        cnmetabolism_organ = getattr(cnmetabolism_phytomer, CNMETABOLISM_ATTRIBUTES_MAPPING[MTG_TO_CNMETABOLISM_PHYTOMERS_ORGANS_MAPPING[mtg_organ_label]])
                        if cnmetabolism_organ is None:
                            continue
                        cnmetabolism_organ_property_names = [property_name for property_name in cnmetabolism_simulation.Simulation.ORGANS_RUN_VARIABLES if hasattr(cnmetabolism_organ, property_name)]
                        for cnmetabolism_organ_property_name in cnmetabolism_organ_property_names:
                            attribute_value = getattr(cnmetabolism_organ, cnmetabolism_organ_property_name)
                            # TODO: temporary ; replace by inputs at photosynthetic organs scale
                            if attribute_value is not None:
                                self._shared_mtg.property(cnmetabolism_organ_property_name)[mtg_organ_vid] = attribute_value
                            elif cnmetabolism_organ_property_name not in self._shared_mtg.get_vertex_property(mtg_organ_vid):
                                self._shared_mtg.property(cnmetabolism_organ_property_name)[mtg_organ_vid] = attribute_value

                        for mtg_element_vid in self._shared_mtg.components_iter(mtg_organ_vid):
                            mtg_element_label = self._shared_mtg.label(mtg_element_vid)
                            if mtg_element_label not in cnmetabolism_converter.DATAFRAME_TO_CNMETABOLISM_ELEMENTS_NAMES_MAPPING:
                                continue
                            cnmetabolism_element = getattr(cnmetabolism_organ, cnmetabolism_converter.DATAFRAME_TO_CNMETABOLISM_ELEMENTS_NAMES_MAPPING[mtg_element_label])
                            cnmetabolism_element_property_names = [property_name for property_name in cnmetabolism_simulation.Simulation.ELEMENTS_RUN_VARIABLES if hasattr(cnmetabolism_element, property_name)]
                            for cnmetabolism_element_property_name in cnmetabolism_element_property_names:
                                cnmetabolism_element_property_value = getattr(cnmetabolism_element, cnmetabolism_element_property_name)
                                self._shared_mtg.property(cnmetabolism_element_property_name)[mtg_element_vid] = cnmetabolism_element_property_value
                                self._shared_mtg.property(cnmetabolism_element_property_name)[mtg_organ_vid] = cnmetabolism_element_property_value  # Update organ property too

                #: Temporary: Store Soil variables at axis level
                axis_id = (cnmetabolism_plant_index, cnmetabolism_axis_label)
                if axis_id in self.soils.keys():
                    if 'soil' not in self._shared_mtg.get_vertex_property(mtg_axis_vid):
                        # Add a property describing the organ to the current axis of the MTG
                        self._shared_mtg.property('soil')[mtg_axis_vid] = {}
                    # Update the property describing the organ of the current axis in the MTG
                    mtg_soil_properties = self._shared_mtg.get_vertex_property(mtg_axis_vid)['soil']
                    for cnmetabolism_property_name in cnmetabolism_simulation.Simulation.SOILS_RUN_VARIABLES:
                        if hasattr(self.soils[axis_id], cnmetabolism_property_name):
                            mtg_soil_properties[cnmetabolism_property_name] = getattr(self.soils[axis_id], cnmetabolism_property_name)

    def _update_shared_dataframes(self, cnmetabolism_axes_data_df=None, cnmetabolism_organs_data_df=None,
                                  cnmetabolism_hiddenzones_data_df=None, cnmetabolism_elements_data_df=None,
                                  cnmetabolism_soils_data_df=None):
        """
        Update the dataframes shared between all models from the inputs dataframes or the outputs dataframes of the cnmetabolism model.

        :param pandas.DataFrame cnmetabolism_axes_data_df: CN-Metabolism shared dataframe at axis scale
        :param pandas.DataFrame cnmetabolism_organs_data_df: CN-Metabolism shared dataframe at organ scale
        :param pandas.DataFrame cnmetabolism_hiddenzones_data_df: CN-Metabolism shared dataframe hiddenzone scale
        :param pandas.DataFrame cnmetabolism_elements_data_df: CN-Metabolism shared dataframe at element scale
        :param pandas.DataFrame or None cnmetabolism_soils_data_df: CN-Metabolism shared dataframe at soil scale
        """

        for cnmetabolism_data_df, \
            shared_inputs_outputs_indexes, \
            shared_inputs_outputs_df in ((cnmetabolism_axes_data_df, cnmetabolism_simulation.Simulation.AXES_INDEXES, self._shared_axes_inputs_outputs_df),
                                         (cnmetabolism_organs_data_df, cnmetabolism_simulation.Simulation.ORGANS_INDEXES, self._shared_organs_inputs_outputs_df),
                                         (cnmetabolism_hiddenzones_data_df, cnmetabolism_simulation.Simulation.HIDDENZONE_INDEXES, self._shared_hiddenzones_inputs_outputs_df),
                                         (cnmetabolism_elements_data_df, cnmetabolism_simulation.Simulation.ELEMENTS_INDEXES, self._shared_elements_inputs_outputs_df),
                                         (cnmetabolism_soils_data_df, cnmetabolism_simulation.Simulation.SOILS_INDEXES, self._shared_soils_inputs_outputs_df)):

            if cnmetabolism_data_df is None:
                continue

            tools.combine_dataframes_inplace(cnmetabolism_data_df, shared_inputs_outputs_indexes, shared_inputs_outputs_df)
