# -*- coding: latin-1 -*-
import math

from openalea.hydraulics import model as hydraulics_model, simulation as hydraulics_simulation, \
    converter as hydraulics_converter, postprocessing as hydraulics_postprocessing

from openalea.integration import tools

"""
    integration.hydraulics_facade
    ~~~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`integration.hydraulics_facade` is a facade of the model Hydraulics.

    This module permits to initialize and run the model Hydraulics from a :class:`MTG <openalea.mtg.mtg.MTG>`
    in a convenient and transparent way, wrapping all the internal complexity of the model, and dealing
    with all the tedious initialization and conversion processes.

    :license: TODO, see LICENSE for details.

"""

#: the mapping of Hydraulics organ classes to the attributes in axis and phytomer which represent an organ
hydraulics_ATTRIBUTES_MAPPING = {hydraulics_model.Internode: 'internode', hydraulics_model.Lamina: 'lamina', hydraulics_model.Sheath: 'sheath',
                                   hydraulics_model.Roots: 'roots', hydraulics_model.HiddenZone: 'hiddenzone', hydraulics_model.Xylem: 'xylem'}

#: the mapping of organs (which belong to an axis) labels in MTG to organ classes in Hydraulics
MTG_TO_hydraulics_AXES_ORGANS_MAPPING = {'xylem': hydraulics_model.Xylem, 'roots': hydraulics_model.Roots}

#: the mapping of organs (which belong to a phytomer) labels in MTG to organ classes in Hydraulics
MTG_TO_hydraulics_PHYTOMERS_ORGANS_MAPPING = {'internode': hydraulics_model.Internode, 'blade': hydraulics_model.Lamina, 'sheath': hydraulics_model.Sheath,
                                                'hiddenzone': hydraulics_model.HiddenZone}

#: the mapping of Hydraulics photosynthetic organs to Hydraulics photosynthetic organ elements
hydraulics_ORGANS_TO_ELEMENTS_MAPPING = {hydraulics_model.Internode: hydraulics_model.InternodeElement, hydraulics_model.Lamina: hydraulics_model.LaminaElement,
                                           hydraulics_model.Sheath: hydraulics_model.SheathElement}

#: the parameters and variables which define the state of a Hydraulics population
POPULATION_STATE_VARIABLE = set(hydraulics_simulation.Simulation.PLANTS_STATE + hydraulics_simulation.Simulation.AXES_STATE +
                                hydraulics_simulation.Simulation.ORGANS_STATE + hydraulics_simulation.Simulation.PHYTOMERS_STATE +
                                hydraulics_simulation.Simulation.HIDDENZONE_STATE + hydraulics_simulation.Simulation.ELEMENTS_STATE)

#: all the variables of a Hydraulics population computed during a run step of the simulation
POPULATION_RUN_VARIABLES = set(hydraulics_simulation.Simulation.PLANTS_RUN_VARIABLES + hydraulics_simulation.Simulation.AXES_RUN_VARIABLES +
                               hydraulics_simulation.Simulation.PHYTOMERS_RUN_VARIABLES + hydraulics_simulation.Simulation.ORGANS_RUN_VARIABLES +
                               hydraulics_simulation.Simulation.HIDDENZONE_RUN_VARIABLES + hydraulics_simulation.Simulation.ELEMENTS_RUN_VARIABLES)

#: all the variables to be stored in the MTG
MTG_RUN_VARIABLES = set(list(POPULATION_RUN_VARIABLES) + hydraulics_simulation.Simulation.SOILS_RUN_VARIABLES)

# number of seconds in 1 hour
HOUR_TO_SECOND_CONVERSION_FACTOR = 3600


class hydraulicsFacade(object):
    """
    The hydraulicsFacade class permits to initialize, run the model Hydraulics
    from a :class:`MTG <openalea.mtg.mtg.MTG>`, and update the MTG and the dataframes
    shared between all models.

    Use :meth:`run` to run the model.

    """

    def __init__(self, shared_mtg, delta_t, update_parameters,
                 model_axes_inputs_df,
                 model_hiddenzones_inputs_df,
                 model_elements_inputs_df,
                 model_organs_inputs_df,
                 model_soils_inputs_df,
                 shared_axes_inputs_outputs_df,
                 shared_hiddenzones_inputs_outputs_df,
                 shared_elements_inputs_outputs_df,
                 shared_organs_inputs_outputs_df,
                 shared_soils_inputs_outputs_df,
                 update_shared_df=True):
        """
                :param openalea.mtg.mtg.MTG shared_mtg: The MTG shared between all models.
                :param int delta_t: The delta between two runs, in seconds.
                :param dict update_parameters: A dictionary with the parameters to update, should have the form {'Organ_label1': {'param1': value1, 'param2': value2}, ...}.
                :param pandas.DataFrame model_hiddenzones_inputs_df: the inputs of the model at hiddenzones scale.
                :param pandas.DataFrame model_elements_inputs_df: the inputs of the model at elements scale.
                :param pandas.DataFrame model_organs_inputs_df: the inputs of the model at organ scale.
                :param pandas.DataFrame model_soils_inputs_df: the inputs of the model at soils scale.
                :param pandas.DataFrame shared_hiddenzones_inputs_outputs_df: the dataframe of inputs and outputs at hiddenzones scale shared between all models.
                :param pandas.DataFrame shared_elements_inputs_outputs_df: the dataframe of inputs and outputs at elements scale shared between all models.
                :param pandas.DataFrame shared_organs_inputs_outputs_df: the dataframe of inputs and outputs at organ scale shared between all models.
                :param pandas.DataFrame shared_soils_inputs_outputs_df: the dataframe of inputs and outputs at soils scale shared between all models.
                :param bool update_shared_df: If `True`  update the shared dataframes at init and at each run (unless stated otherwise)
        """
        
        self._shared_mtg = shared_mtg  #: the MTG shared between all models

        self._simulation = hydraulics_simulation.Simulation(delta_t=delta_t)

        self.population, self.soils = hydraulics_converter.from_dataframes(model_axes_inputs_df, model_hiddenzones_inputs_df, model_elements_inputs_df, model_organs_inputs_df, model_soils_inputs_df)

        self._update_parameters = update_parameters

        self._simulation.initialize(self.population, self.soils)

        self._update_shared_MTG()

        self._shared_axes_inputs_outputs_df = shared_axes_inputs_outputs_df         #: the dataframe at axes scale shared between all models
        self._shared_hiddenzones_inputs_outputs_df = shared_hiddenzones_inputs_outputs_df         #: the dataframe at hiddenzones scale shared between all models
        self._shared_elements_inputs_outputs_df = shared_elements_inputs_outputs_df               #: the dataframe at elements scale shared between all models
        self._shared_organs_inputs_outputs_df = shared_organs_inputs_outputs_df  #: the dataframe at organs scale shared between all models
        self._shared_soils_inputs_outputs_df = shared_soils_inputs_outputs_df  #: the dataframe at soils scale shared between all models
        self._update_shared_df = update_shared_df
        if self._update_shared_df:
            self._update_shared_dataframes(hydraulics_axes_data_df=model_axes_inputs_df,
                                           hydraulics_hiddenzones_data_df=model_hiddenzones_inputs_df,
                                           hydraulics_elements_data_df=model_elements_inputs_df,
                                           hydraulics_organs_data_df=model_organs_inputs_df,
                                           hydraulics_soils_data_df=model_soils_inputs_df)

    def run(self, update_shared_df=False):

        """
        Run the model and update the MTG and the dataframes shared between all models.
        """

        self._initialize_model()
        self._simulation.run()
        self._update_shared_MTG()

        if update_shared_df or (update_shared_df is None and self._update_shared_df):
            (_, hydraulics_axes_inputs_outputs_df, _, hydraulics_organs_inputs_outputs_df, hydraulics_hiddenzones_inputs_outputs_df, hydraulics_elements_inputs_outputs_df,
             hydraulics_soils_inputs_outputs_df) = hydraulics_converter.to_dataframes(self._simulation.population, self._simulation.soils)

            self._update_shared_dataframes(hydraulics_axes_data_df=hydraulics_axes_inputs_outputs_df,
                                           hydraulics_hiddenzones_data_df=hydraulics_hiddenzones_inputs_outputs_df,
                                           hydraulics_elements_data_df=hydraulics_elements_inputs_outputs_df,
                                           hydraulics_organs_data_df=hydraulics_organs_inputs_outputs_df,
                                           hydraulics_soils_data_df=hydraulics_soils_inputs_outputs_df)

    @staticmethod
    def postprocessing(axes_outputs_df, hiddenzone_outputs_df, elements_outputs_df, organs_outputs_df, soils_outputs_df, delta_t):
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
        depending of the dataframes given as argument.
        For example, if user passes only dataframes `plants_df`, `axes_df` and `metamers_df`,
        then only post-processing dataframes of plants, axes and metamers are returned.

    :rtype: dict ['scale' : pandas.DataFrame]
        """

        (_, axes_postprocessing_df, _, hiddenzones_postprocessing_df, elements_postprocessing_df, organs_postprocessing_df, soils_postprocessing_df) = (
            hydraulics_postprocessing.postprocessing(axes_df=axes_outputs_df, hiddenzones_df=hiddenzone_outputs_df, elements_df=elements_outputs_df, organs_df=organs_outputs_df, soils_df= soils_outputs_df, delta_t=delta_t))

        return {
            'axes': axes_postprocessing_df,
            'organs': organs_postprocessing_df,
            'hiddenzones': hiddenzones_postprocessing_df,
            'elements': elements_postprocessing_df,
            'soils': soils_postprocessing_df
        }

    @staticmethod
    def graphs(axes_postprocessing_df, hiddenzones_postprocessing_df, elements_postprocessing_df, organs_postprocessing_df, soils_postprocessing_df, meteo_data, graphs_dirpath='.'):
        """
        Generate the graphs and save them into `graphs_dirpath`.

        :param pandas.DataFrame axes_postprocessing_df: CN-Metabolism outputs at axis scale
        :param pandas.DataFrame hiddenzones_postprocessing_df: CN-Metabolism outputs at hidden zone scale
        :param pandas.DataFrame organs_postprocessing_df: CN-Metabolism outputs at organ scale
        :param pandas.DataFrame elements_postprocessing_df: CN-Metabolism outputs at element scale
        :param pandas.DataFrame soils_postprocessing_df: CN-Metabolism outputs at soil scale
        :param pandas.DataFrame meteo_data: the meteo dataframe having the mapping between t (hours) and calendar dates
        :param str graphs_dirpath: the path of the directory to save the generated graphs in
        """
        hydraulics_postprocessing.generate_graphs(axes_df=axes_postprocessing_df, hiddenzones_df=hiddenzones_postprocessing_df, elements_df=elements_postprocessing_df,
                                                    organs_df=organs_postprocessing_df, soils_df=soils_postprocessing_df, meteo_data=meteo_data, graphs_dirpath=graphs_dirpath)

    def _initialize_model(self):
        """
        Initialize the inputs of the model from the MTG shared between all models.
        """

        self.population = hydraulics_model.Population()
        mapping_topology = {'predecessor': {}, 'successor': {}}

        # traverse the MTG recursively from top

        for mtg_plant_vid in self._shared_mtg.components_iter(self._shared_mtg.root):
            mtg_plant_index = int(self._shared_mtg.index(mtg_plant_vid))

            # create a new plant
            hydraulics_plant = hydraulics_model.Plant(mtg_plant_index)
            is_valid_plant = False
            for mtg_axis_vid in self._shared_mtg.components_iter(mtg_plant_vid):
                mtg_axis_label = self._shared_mtg.label(mtg_axis_vid)
                #: keep only MS TODO: temporary
                if mtg_axis_label != 'MS':
                    continue

                #: MS
                # create a new axis
                hydraulics_axis = hydraulics_model.Axis(mtg_axis_label)
                mtg_axis_properties = self._shared_mtg.get_vertex_property(mtg_axis_vid)
                hydraulics_axis_data_dict = {}
                for hydraulics_axis_data_name in hydraulics_simulation.Simulation.AXES_STATE:
                    hydraulics_axis_data_dict[hydraulics_axis_data_name] = mtg_axis_properties[hydraulics_axis_data_name]
                hydraulics_axis.__dict__.update(hydraulics_axis_data_dict)
                is_valid_axis = True
                for hydraulics_organ_class in (hydraulics_model.Roots, hydraulics_model.Xylem):
                    mtg_organ_label = hydraulics_converter.hydraulics_CLASSES_TO_DATAFRAME_ORGANS_MAPPING[hydraulics_organ_class]
                    # create a new organ
                    hydraulics_organ = hydraulics_organ_class(mtg_organ_label)
                    mtg_axis_properties = self._shared_mtg.get_vertex_property(mtg_axis_vid)
                    if mtg_organ_label in mtg_axis_properties:
                        mtg_organ_properties = mtg_axis_properties[mtg_organ_label]
                        hydraulics_organ_data_names = set(hydraulics_simulation.Simulation.ORGANS_STATE).intersection(hydraulics_organ.__dict__)
                        if set(mtg_organ_properties).issuperset(hydraulics_organ_data_names):
                            hydraulics_organ_data_dict = {}
                            for hydraulics_organ_data_name in hydraulics_organ_data_names:
                                hydraulics_organ_data_dict[hydraulics_organ_data_name] = mtg_organ_properties[hydraulics_organ_data_name]

                                # Debug: Tell if missing input variable
                                if math.isnan(mtg_organ_properties[hydraulics_organ_data_name]) or mtg_organ_properties[hydraulics_organ_data_name] is None:
                                    print('Missing variable', hydraulics_organ_data_name, 'for vertex id', mtg_axis_vid, 'which is', mtg_organ_label)

                            hydraulics_organ.__dict__.update(hydraulics_organ_data_dict)

                            # Update parameters if specified
                            if mtg_organ_label in self._update_parameters:
                                hydraulics_organ.PARAMETERS.__dict__.update(self._update_parameters[mtg_organ_label])

                            hydraulics_organ.initialize()
                            # add the new organ to current axis
                            setattr(hydraulics_axis, mtg_organ_label, hydraulics_organ)
                if not is_valid_axis:
                    continue

                has_valid_phytomer = False
                for mtg_metamer_vid in self._shared_mtg.components_iter(mtg_axis_vid):
                    mtg_metamer_index = int(self._shared_mtg.index(mtg_metamer_vid))

                    # create a new phytomer
                    hydraulics_phytomer = hydraulics_model.Phytomer(mtg_metamer_index)
                    mtg_hiddenzone_label = hydraulics_converter.hydraulics_CLASSES_TO_DATAFRAME_ORGANS_MAPPING[hydraulics_model.HiddenZone]
                    mtg_metamer_properties = self._shared_mtg.get_vertex_property(mtg_metamer_vid)

                    if mtg_hiddenzone_label in mtg_metamer_properties and mtg_metamer_properties[mtg_hiddenzone_label]['leaf_is_growing']:
                        has_valid_hiddenzone = True
                        hydraulics_hiddenzone = hydraulics_model.HiddenZone(label=mtg_hiddenzone_label)
                        mtg_hiddenzone_properties = mtg_metamer_properties[mtg_hiddenzone_label]

                        # Adding aggregated variables into inputs
                        hydraulics_hiddenzone_data_names = set(hydraulics_simulation.Simulation.HIDDENZONE_RUN_VARIABLES).intersection(hydraulics_hiddenzone.__dict__)

                        if mtg_hiddenzone_properties.get('leaf_pseudo_age') == 0:  # First time hiddenzone passes into hydraulics sub-model
                            missing_initial_hiddenzone_properties = hydraulics_hiddenzone_data_names - set(mtg_hiddenzone_properties)
                            hydraulics_hiddenzone_data_names -= missing_initial_hiddenzone_properties

                        if set(mtg_hiddenzone_properties).issuperset(hydraulics_hiddenzone_data_names):
                            hydraulics_hiddenzone_data_dict = {}
                            for hydraulics_hiddenzone_data_name in hydraulics_hiddenzone_data_names:
                                mtg_hiddenzone_data_value = mtg_hiddenzone_properties.get(hydraulics_hiddenzone_data_name)
                                hydraulics_hiddenzone_data_dict[hydraulics_hiddenzone_data_name] = mtg_hiddenzone_data_value
                            hydraulics_hiddenzone.__dict__.update(hydraulics_hiddenzone_data_dict)

                        # Update parameters if specified
                        if mtg_hiddenzone_label in self._update_parameters:
                            hydraulics_hiddenzone.PARAMETERS.__dict__.update(self._update_parameters[mtg_hiddenzone_label])

                        hydraulics_hiddenzone.initialize()
                        # add the new hiddenzone to current phytomer
                        setattr(hydraulics_phytomer, mtg_hiddenzone_label, hydraulics_hiddenzone)

                    else:
                        has_valid_hiddenzone = False

                    # create a new organ
                    has_valid_organ = False
                    for mtg_organ_vid in self._shared_mtg.components_iter(mtg_metamer_vid):
                        mtg_organ_label = self._shared_mtg.label(mtg_organ_vid)
                        if mtg_organ_label == "internode":  # No internode in Hydraulics model
                            continue
                        if mtg_organ_label not in MTG_TO_hydraulics_PHYTOMERS_ORGANS_MAPPING or self._shared_mtg.get_vertex_property(mtg_organ_vid)['length'] == 0:
                            continue
                        hydraulics_organ_class = MTG_TO_hydraulics_PHYTOMERS_ORGANS_MAPPING[mtg_organ_label]
                        hydraulics_organ = hydraulics_organ_class(mtg_organ_label)

                        # Update parameters if specified
                        if 'PhotosyntheticOrgan' in self._update_parameters:
                            hydraulics_organ.PARAMETERS.__dict__.update(self._update_parameters['PhotosyntheticOrgan'])

                        hydraulics_organ.initialize()
                        has_valid_element = False

                        # create a new element
                        for mtg_element_vid in self._shared_mtg.components_iter(mtg_organ_vid):
                            mtg_element_label = self._shared_mtg.label(mtg_element_vid)
                            if mtg_element_label not in hydraulics_converter.DATAFRAME_TO_hydraulics_ELEMENTS_NAMES_MAPPING \
                                    or (self._shared_mtg.get_vertex_property(mtg_element_vid)['length'] == 0) \
                                    or ((mtg_element_label == 'HiddenElement') and (self._shared_mtg.get_vertex_property(mtg_element_vid).get('is_growing', True)) \
                                    or (self._shared_mtg.get_vertex_property(mtg_element_vid).get('is_over', True))):
                                continue

                            has_valid_element = True
                            hydraulics_element = hydraulics_ORGANS_TO_ELEMENTS_MAPPING[hydraulics_organ_class](label=mtg_element_label)
                            mtg_element_properties = self._shared_mtg.get_vertex_property(mtg_element_vid)
                            hydraulics_element_data_names = set(hydraulics_simulation.Simulation.ELEMENTS_RUN_VARIABLES).intersection(hydraulics_element.__dict__)    #: Adding aggregated variables into inputs
                            if mtg_element_properties.get('age') == 0:  # First time element passes into hydraulics model
                                missing_initial_element_properties = hydraulics_element_data_names - set(mtg_element_properties)
                                hydraulics_element_data_names -= missing_initial_element_properties

                            if set(mtg_element_properties).issuperset(hydraulics_element_data_names):
                                hydraulics_element_data_dict = {}
                                for hydraulics_element_data_name in hydraulics_element_data_names:
                                    mtg_element_data_value = mtg_element_properties.get(hydraulics_element_data_name)
                                    hydraulics_element_data_dict[hydraulics_element_data_name] = mtg_element_data_value
                                hydraulics_element.__dict__.update(hydraulics_element_data_dict)
                                # add element to organ
                                setattr(hydraulics_organ, hydraulics_converter.DATAFRAME_TO_hydraulics_ELEMENTS_NAMES_MAPPING[mtg_element_label], hydraulics_element)

                            #: Update lamina_Lmax & Wmax
                            if mtg_organ_label == 'blade':
                                if has_valid_hiddenzone is True:
                                    mtg_organ_properties = self._shared_mtg.get_vertex_property(mtg_organ_vid)
                                    if mtg_element_properties['length'] >= mtg_hiddenzone_properties['lamina_Lmax']:
                                        mtg_hiddenzone_properties['lamina_Lmax'] = mtg_element_properties['length']
                                    mtg_element_properties['Wmax'] = mtg_hiddenzone_properties['leaf_Wmax']
                                    mtg_organ_properties['shape_max_width'] = mtg_hiddenzone_properties['leaf_Wmax']
                                    mtg_organ_properties['shape_mature_length'] = mtg_hiddenzone_properties['lamina_Lmax']
                                else:
                                    mtg_organ_properties['shape_max_width'] = mtg_element_properties['Wmax']
                                    mtg_organ_properties['shape_mature_length'] = mtg_element_properties['length']

                        if has_valid_element:
                            has_valid_organ = True
                            setattr(hydraulics_phytomer, hydraulics_ATTRIBUTES_MAPPING[hydraulics_organ_class], hydraulics_organ)

                    if has_valid_organ or has_valid_hiddenzone:
                        hydraulics_axis.phytomers.append(hydraulics_phytomer)
                        has_valid_phytomer = True

                if not has_valid_phytomer:
                    is_valid_axis = False

                if is_valid_axis:
                    hydraulics_plant.axes.append(hydraulics_axis)
                    is_valid_plant = True

            if is_valid_plant:
                self.population.plants.append(hydraulics_plant)

        self._simulation.initialize(self.population, self.soils)

    def _update_shared_MTG(self):
        """
        Update the MTG shared between all models from the population of Hydraulics.
        """
        # add the missing properties
        mtg_property_names = self._shared_mtg.property_names()
        for hydraulics_data_name in MTG_RUN_VARIABLES:
            if hydraulics_data_name not in mtg_property_names:
                self._shared_mtg.add_property(hydraulics_data_name)
        for hydraulics_organ_label in list(MTG_TO_hydraulics_AXES_ORGANS_MAPPING.keys()) + [hydraulics_converter.hydraulics_CLASSES_TO_DATAFRAME_ORGANS_MAPPING[hydraulics_model.HiddenZone]]:
            if hydraulics_organ_label not in mtg_property_names:
                self._shared_mtg.add_property(hydraulics_organ_label)

        mtg_plants_iterator = self._shared_mtg.components_iter(self._shared_mtg.root)
        # traverse Turgor_Growth population from top
        for hydraulics_plant in self.population.plants:
            hydraulics_plant_index = hydraulics_plant.index
            while True:
                mtg_plant_vid = next(mtg_plants_iterator)
                if int(self._shared_mtg.index(mtg_plant_vid)) == hydraulics_plant_index:
                    break
            mtg_axes_iterator = self._shared_mtg.components_iter(mtg_plant_vid)
            for hydraulics_axis in hydraulics_plant.axes:
                hydraulics_axis_label = hydraulics_axis.label
                while True:
                    mtg_axis_vid = next(mtg_axes_iterator)
                    if self._shared_mtg.label(mtg_axis_vid) == hydraulics_axis_label:
                        break

                # Xylem
                hydraulics_axis_property_names = [property_name for property_name in hydraulics_simulation.Simulation.AXES_RUN_VARIABLES if hasattr(hydraulics_axis, property_name)]
                for hydraulics_axis_property_name in hydraulics_axis_property_names:
                    hydraulics_axis_property_value = getattr(hydraulics_axis, hydraulics_axis_property_name)
                    self._shared_mtg.property(hydraulics_axis_property_name)[mtg_axis_vid] = hydraulics_axis_property_value
                for mtg_organ_label in MTG_TO_hydraulics_AXES_ORGANS_MAPPING.keys():
                    if mtg_organ_label not in self._shared_mtg.get_vertex_property(mtg_axis_vid):
                        # Add a property describing the organ to the current axis of the MTG
                        self._shared_mtg.property(mtg_organ_label)[mtg_axis_vid] = {}
                    # Update the property describing the organ of the current axis in the MTG
                    hydraulics_organ = getattr(hydraulics_axis, mtg_organ_label)
                    mtg_organ_properties = self._shared_mtg.get_vertex_property(mtg_axis_vid)[mtg_organ_label]
                    for hydraulics_property_name in hydraulics_simulation.Simulation.ORGANS_RUN_VARIABLES:
                        if hasattr(hydraulics_organ, hydraulics_property_name):
                            mtg_organ_properties[hydraulics_property_name] = getattr(hydraulics_organ, hydraulics_property_name)
                mtg_metamers_iterator = self._shared_mtg.components_iter(mtg_axis_vid)

                for hydraulics_phytomer in hydraulics_axis.phytomers:
                    hydraulics_phytomer_index = hydraulics_phytomer.index
                    while True:
                        mtg_metamer_vid = next(mtg_metamers_iterator)
                        if int(self._shared_mtg.index(mtg_metamer_vid)) == hydraulics_phytomer_index:
                            break
                    if hydraulics_phytomer.hiddenzone is not None:
                        mtg_hiddenzone_label = hydraulics_converter.hydraulics_CLASSES_TO_DATAFRAME_ORGANS_MAPPING[hydraulics_model.HiddenZone]
                        if mtg_hiddenzone_label not in self._shared_mtg.get_vertex_property(mtg_metamer_vid):
                            # Add a property describing the hiddenzone to the current metamer of the MTG
                            self._shared_mtg.property(mtg_hiddenzone_label)[mtg_metamer_vid] = {}
                        # Update the property describing the hiddenzone of the current metamer in the MTG
                        mtg_hiddenzone_properties = self._shared_mtg.get_vertex_property(mtg_metamer_vid)[mtg_hiddenzone_label]

                        mtg_hiddenzone_properties.update(hydraulics_phytomer.hiddenzone.__dict__)
                    for mtg_organ_vid in self._shared_mtg.components_iter(mtg_metamer_vid):
                        mtg_organ_label = self._shared_mtg.label(mtg_organ_vid)
                        if mtg_organ_label == "internode":  # No internode in Hydraulics model
                            continue
                        if mtg_organ_label not in MTG_TO_hydraulics_PHYTOMERS_ORGANS_MAPPING:
                            continue
                        hydraulics_organ = getattr(hydraulics_phytomer, hydraulics_ATTRIBUTES_MAPPING[MTG_TO_hydraulics_PHYTOMERS_ORGANS_MAPPING[mtg_organ_label]])

                        if hydraulics_organ is None:
                            continue
                        # element scale
                        for mtg_element_vid in self._shared_mtg.components_iter(mtg_organ_vid):
                            mtg_element_label = self._shared_mtg.label(mtg_element_vid)

                            if mtg_element_label not in hydraulics_converter.DATAFRAME_TO_hydraulics_ELEMENTS_NAMES_MAPPING: continue

                            hydraulics_element = getattr(hydraulics_organ, hydraulics_converter.DATAFRAME_TO_hydraulics_ELEMENTS_NAMES_MAPPING[mtg_element_label])
                            hydraulics_element_property_names = [property_name for property_name in hydraulics_simulation.Simulation.ELEMENTS_RUN_VARIABLES if hasattr(hydraulics_element, property_name)]
                            for hydraulics_element_property_name in hydraulics_element_property_names:
                                hydraulics_element_property_value = getattr(hydraulics_element, hydraulics_element_property_name)
                                self._shared_mtg.property(hydraulics_element_property_name)[mtg_element_vid] = hydraulics_element_property_value

                        # update of organ scale from elements
                        new_mtg_element_labels = {}

                        for new_element_vid in self._shared_mtg.components_iter(mtg_organ_vid):
                            new_element_label = self._shared_mtg.label(new_element_vid)
                            new_mtg_element_labels[new_element_label] = new_element_vid

                        if mtg_organ_label == 'blade' and 'LeafElement1' in new_mtg_element_labels.keys():
                            leaf_element_mtg_properties = self._shared_mtg.get_vertex_property(new_mtg_element_labels['LeafElement1'])
                            organ_visible_length = leaf_element_mtg_properties['length']
                            self._shared_mtg.property('visible_length')[mtg_organ_vid] = organ_visible_length
                            #: Lamina_Lmax & Wmax in Hydraulics
                            if leaf_element_mtg_properties['is_growing'] is True:
                                if leaf_element_mtg_properties['length'] >= mtg_hiddenzone_properties['lamina_Lmax']:
                                    self._shared_mtg.property('lamina_Lmax')[mtg_metamer_vid] = leaf_element_mtg_properties['length']
                                self._shared_mtg.property('Wmax')[new_mtg_element_labels['LeafElement1']] = mtg_hiddenzone_properties['leaf_Wmax']
                                self._shared_mtg.property('shape_mature_length')[mtg_organ_vid] = mtg_hiddenzone_properties['lamina_Lmax']
                                self._shared_mtg.property('shape_max_width')[mtg_organ_vid] = mtg_hiddenzone_properties['leaf_Wmax']
                            else:
                                self._shared_mtg.property('shape_mature_length')[mtg_organ_vid] = leaf_element_mtg_properties['length']
                                self._shared_mtg.property('shape_max_width')[mtg_organ_vid] = leaf_element_mtg_properties['Wmax']

                        elif mtg_organ_label == 'sheath' and 'StemElement' in new_mtg_element_labels.keys():
                            organ_visible_length = self._shared_mtg.property('length')[new_mtg_element_labels['StemElement']]
                            self._shared_mtg.property('visible_length')[mtg_organ_vid] = organ_visible_length
                        elif mtg_organ_label == 'internode' and 'StemElement' in new_mtg_element_labels.keys():
                            organ_visible_length = self._shared_mtg.property('length')[new_mtg_element_labels['StemElement']]
                            self._shared_mtg.property('visible_length')[mtg_organ_vid] = organ_visible_length
                        else:
                            organ_visible_length = 0

                        #: Internode length
                        if 'HiddenElement' in new_mtg_element_labels.keys():
                            organ_hidden_length = self._shared_mtg.property('length')[new_mtg_element_labels['HiddenElement']]
                        else:
                            organ_hidden_length = 0

                        total_organ_length = organ_visible_length + organ_hidden_length
                        self._shared_mtg.property('length')[mtg_organ_vid] = total_organ_length

                #: Temporary: Store Soil variables at axis level
                axis_id = (hydraulics_plant_index, hydraulics_axis_label)
                if axis_id in self.soils.keys():
                    if 'soil' not in self._shared_mtg.get_vertex_property(mtg_axis_vid):
                        # Add a property describing the organ to the current axis of the MTG
                        self._shared_mtg.property('soil')[mtg_axis_vid] = {}
                    # Update the property describing the organ of the current axis in the MTG
                    mtg_soil_properties = self._shared_mtg.get_vertex_property(mtg_axis_vid)['soil']
                    for hydraulics_property_name in hydraulics_simulation.Simulation.SOILS_RUN_VARIABLES:
                        if hasattr(self.soils[axis_id], hydraulics_property_name):
                            mtg_soil_properties[hydraulics_property_name] = getattr(self.soils[axis_id], hydraulics_property_name)

    def _update_shared_dataframes(self, hydraulics_axes_data_df=None, hydraulics_organs_data_df=None,
                                  hydraulics_hiddenzones_data_df=None, hydraulics_elements_data_df=None,
                                  hydraulics_soils_data_df=None):
        """
        Update the dataframes shared between all models from the inputs dataframes or the outputs dataframes of the cnmetabolism model.

        :param pandas.DataFrame hydraulics_axes_data_df: CN-Metabolism shared dataframe at axis scale
        :param pandas.DataFrame hydraulics_organs_data_df: CN-Metabolism shared dataframe at organ scale
        :param pandas.DataFrame hydraulics_hiddenzones_data_df: CN-Metabolism shared dataframe hiddenzone scale
        :param pandas.DataFrame hydraulics_elements_data_df: CN-Metabolism shared dataframe at element scale
        :param pandas.DataFrame hydraulics_soils_data_df: CN-Metabolism shared dataframe at soil scale
        :param pandas.DataFrame cnmetabolism_soils_data_df: CN-Metabolism shared dataframe at soil scale
        """

        for hydraulics_data_df, \
            shared_inputs_outputs_indexes, \
            shared_inputs_outputs_df in ((hydraulics_axes_data_df, hydraulics_simulation.Simulation.AXES_INDEXES, self._shared_axes_inputs_outputs_df),
                                         (hydraulics_hiddenzones_data_df, hydraulics_simulation.Simulation.HIDDENZONES_INDEXES, self._shared_hiddenzones_inputs_outputs_df),
                                         (hydraulics_elements_data_df, hydraulics_simulation.Simulation.ELEMENTS_INDEXES, self._shared_elements_inputs_outputs_df),
                                         (hydraulics_organs_data_df, hydraulics_simulation.Simulation.ORGANS_INDEXES, self._shared_organs_inputs_outputs_df),
                                         (hydraulics_soils_data_df, hydraulics_simulation.Simulation.SOILS_INDEXES, self._shared_soils_inputs_outputs_df)):

            if hydraulics_data_df is None: continue
            tools.combine_dataframes_inplace(hydraulics_data_df, shared_inputs_outputs_indexes, shared_inputs_outputs_df)
