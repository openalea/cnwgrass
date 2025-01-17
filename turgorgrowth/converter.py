"""
    turgorgrowth.converter
    ~~~~~~~~~~~~~~~~~

    The module :mod:`turgorgrowth.converter` defines functions to convert
    :class:`dataframes <pandas.DataFrame>` to/from Turgor-Growth inputs or outputs format.

    :license: CeCILL-C, see LICENSE for details.
"""

import numpy as np
import pandas as pd

from turgorgrowth import model, simulation

#: the columns of the outputs dataframe at PLANT scale
PLANTS_VARIABLES = simulation.Simulation.PLANTS_INDEXES + simulation.Simulation.PLANTS_RUN_VARIABLES

#: the columns of the outputs dataframe at AXIS scale
AXES_VARIABLES = simulation.Simulation.AXES_INDEXES + simulation.Simulation.AXES_RUN_VARIABLES

#: the columns of the outputs dataframe at PHYTOMER scale
PHYTOMERS_VARIABLES = simulation.Simulation.PHYTOMERS_INDEXES + simulation.Simulation.PHYTOMERS_RUN_VARIABLES

#: the columns of the outputs dataframe at HIDDENZONE scale
HIDDENZONE_VARIABLES = simulation.Simulation.HIDDENZONE_INDEXES + simulation.Simulation.HIDDENZONE_RUN_VARIABLES
HIDDENZONE_OUTPUTS_VARIABLES = HIDDENZONE_VARIABLES
HIDDENZONE_OUTPUTS_RUN_VARIABLES = simulation.Simulation.HIDDENZONE_RUN_VARIABLES

#: the columns of the outputs dataframe at ELEMENTS scale
ELEMENTS_VARIABLES = simulation.Simulation.ELEMENTS_INDEXES + simulation.Simulation.ELEMENTS_RUN_VARIABLES
ELEMENTS_OUTPUTS_VARIABLES = ELEMENTS_VARIABLES
ELEMENTS_OUTPUTS_RUN_VARIABLES = simulation.Simulation.ELEMENTS_RUN_VARIABLES

#: the columns of the outputs dataframe at ORGANS scale
ORGANS_VARIABLES = simulation.Simulation.ORGANS_INDEXES + simulation.Simulation.ORGANS_RUN_VARIABLES
ORGANS_OUTPUTS_VARIABLES = ORGANS_VARIABLES
ORGANS_OUTPUTS_RUN_VARIABLES = simulation.Simulation.ORGANS_RUN_VARIABLES

#: the mapping of the Turgor-Growth organs classes to organs names in MTG
TURGORGROWTH_CLASSES_TO_DATAFRAME_ORGANS_MAPPING = {model.Organ: 'organs', model.Internode: 'internode', model.Lamina: 'blade', model.Sheath: 'sheath',  model.HiddenZone: 'hiddenzone', model.Roots: 'roots', model.Xylem: 'xylem'}

#: the mapping of the name of each element, from Dataframe to Turgor-Growth
DATAFRAME_TO_TURGORGROWTH_ELEMENTS_NAMES_MAPPING = {'HiddenElement': 'enclosed_element', 'StemElement': 'exposed_element', 'LeafElement1': 'exposed_element'}
DATAFRAME_TO_TURGORGROWTH_LEAF_NAMES_MAPPING = {'LeafElement1': 'exposed_element'}


def from_dataframes(axes_inputs = None, hiddenzones_inputs=None, elements_inputs=None, organs_inputs=None):
    """
    If `elements_inputs` and `hiddenzones_inputs` are not `None`, converts `elements_inputs` and `hiddenzones_inputs` to a :class:`population <model.Population>`.

    :param pandas.DataFrame hiddenzones_inputs: Hidden zone inputs, with one line per hidden zone.
    :param pandas.DataFrame elements_inputs: Element inputs, with one line per element.
    :param pandas.DataFrame organs_inputs: Organs (xylem and roots) inputs, with one line per organ.

    :return: If `elements_inputs` and `hiddenzones_inputs` are not `None`, returns a :class:`population <model.Population>`
    :rtype: (model.Population, dict)
    """

    convert_dataframes_to_population = axes_inputs is not None and elements_inputs is not None and hiddenzones_inputs is not None and organs_inputs is not None

    if convert_dataframes_to_population:
        population = model.Population()
        mapping_topology = {'predecessor': {}, 'successor': {}}

        for plant_index in elements_inputs.plant.unique():
            # create a new plant
            plant = model.Plant(plant_index)
            population.plants.append(plant)

            # curr_axes_labels = elements_inputs[elements_inputs['plant'] == plant_index].axis.unique()
            curr_axes_labels = organs_inputs[organs_inputs['plant'] == plant_index].axis.unique()

            for axis_label in curr_axes_labels:
                # create a new axis
                axis = model.Axis(axis_label)
                curr_organs_inputs = organs_inputs[(organs_inputs['plant'] == plant_index) & (organs_inputs['axis'] == axis_label)]
                for axis_attribute_name, axis_attribute_class in (('roots', model.Roots), ('xylem', model.Xylem)):
                    organ_label = TURGORGROWTH_CLASSES_TO_DATAFRAME_ORGANS_MAPPING[axis_attribute_class]
                    organs_inputs = curr_organs_inputs[curr_organs_inputs['organ'] == organ_label]
                    if not organs_inputs.empty:
                        # create a new organs
                        organs = axis_attribute_class(organ_label)
                        organs_attributes_names = [state_var_name for state_var_name in simulation.Simulation.ORGANS_STATE if hasattr(organs, state_var_name)]
                        organs_row = organs_inputs.loc[organs_inputs.first_valid_index()]
                        organs_attributes_values = organs_row[organs_attributes_names].tolist()
                        organs_attributes = dict(zip(organs_attributes_names, organs_attributes_values))
                        organs.__dict__.update(organs_attributes)
                        # Update parameters if specified
                        # if organs_label in update_parameters:
                            # organs.PARAMETERS.__dict__.update(update_parameters[organs_label])

                        organs.initialize()
                        setattr(axis, axis_attribute_name, organs)

                #Topology of organs
                roots = model.Roots(label='roots')
                # last_elongated_internode = roots
                mapping_topology['successor'][roots] = []
                xylem = model.Xylem(label='xylem')
                mapping_topology['successor'][xylem] = []

                curr_metamers_indexes_for_hiddenzones = hiddenzones_inputs[(hiddenzones_inputs['plant'] == plant_index) & (hiddenzones_inputs['axis'] == axis_label)].metamer.unique()
                curr_metamers_indexes_for_elements = elements_inputs[(elements_inputs['plant'] == plant_index) & (elements_inputs['axis'] == axis_label)].metamer.unique()
                curr_metamers_indexes = np.unique(np.concatenate((curr_metamers_indexes_for_hiddenzones,
                                                                  curr_metamers_indexes_for_elements)))

                for metamer_index in sorted(curr_metamers_indexes):
                    # create a new phytomer
                    phytomer = model.Phytomer(metamer_index)
                    axis.phytomers.append(phytomer)

                    for phytomer_attribute_name, phytomer_attribute_class, phytomer_attribute_element_class in \
                        (('lamina', model.Lamina, model.LaminaElement),
                         ('internode', model.Internode, model.InternodeElement),
                         ('sheath', model.Sheath, model.SheathElement)):

                        organs_label = TURGORGROWTH_CLASSES_TO_DATAFRAME_ORGANS_MAPPING[phytomer_attribute_class]

                        if metamer_index in curr_metamers_indexes_for_elements:
                            curr_elements_inputs = elements_inputs[(elements_inputs['plant'] == plant_index) & (elements_inputs['axis'] == axis_label) & (elements_inputs['metamer'] == metamer_index) &
                                                                   (elements_inputs['organ'] == organs_label)]
                            if organs_label not in curr_elements_inputs.organ.values:
                                continue
                            # create a new organs
                            organs = phytomer_attribute_class(label=organs_label)
                            setattr(phytomer, phytomer_attribute_name, organs)

                            for mtg_element_label, turgorgrowth_element_name in DATAFRAME_TO_TURGORGROWTH_ELEMENTS_NAMES_MAPPING.items():
                                element_inputs = curr_elements_inputs[curr_elements_inputs['element'] == mtg_element_label]
                                if len(element_inputs) == 0:
                                    continue
                                element_inputs = element_inputs.loc[:, simulation.Simulation.ELEMENTS_STATE]
                                element_dict = element_inputs.loc[element_inputs.first_valid_index()].dropna().to_dict()
                                # create a new element
                                element = phytomer_attribute_element_class(label=mtg_element_label, **element_dict)

                                # Add parameters from organ scale
                                element.PARAMETERS.__dict__.update(organs.PARAMETERS.__dict__)

                                setattr(organs, turgorgrowth_element_name, element)

                    #: Hidden zones
                    if metamer_index in curr_metamers_indexes_for_hiddenzones:
                        hiddenzone_inputs = hiddenzones_inputs[(hiddenzones_inputs['plant'] == plant_index) & (hiddenzones_inputs['axis'] == axis_label) &
                                                               (hiddenzones_inputs['metamer'] == metamer_index)]
                        if len(hiddenzone_inputs) == 0:
                            continue
                        hiddenzone_inputs = hiddenzone_inputs.loc[:, simulation.Simulation.HIDDENZONE_STATE]
                        hiddenzone_dict = hiddenzone_inputs.loc[hiddenzone_inputs.first_valid_index()].dropna().to_dict()
                        # create a new hidden zone
                        hiddenzone = model.HiddenZone(TURGORGROWTH_CLASSES_TO_DATAFRAME_ORGANS_MAPPING[model.HiddenZone], **hiddenzone_dict)
                        phytomer.hiddenzone = hiddenzone

                plant.axes.append(axis)

    return population, mapping_topology


def to_dataframes(population=None):
    """
    Convert a Turgor-Growth :class:`population <model.Population>` to Pandas dataframes.
    If `population` is not None, convert `population` to Pandas dataframes.

    :param model.Population population: The Turgor-Growth population to convert.

    :return: If `population` is not None, return :class:`dataframes <pandas.DataFrame>` describing the internal state and compartments of the population at each scale:
                 * hidden zones: plant index, axis id, phytomer index, state parameters, state variables, intermediate variables,
                 fluxes and integrative variables of each hidden zone (see :mod:`HIDDENZONE_VARIABLES`)
                 * element scale: plant index, axis id, phytomer index, organs type, element type, state parameters, state variables, intermediate variables,
                 fluxes and integrative variables of each element (see :mod:`ELEMENTS_VARIABLES`)
                 * xylem scale: xylem index, state parameters, state variables, intermediate variables,
                 fluxes and integrative variables of xylem (see :mod:`XYLEM_VARIABLES`)


    :rtype: (pandas.DataFrame, pandas.DataFrame)
    """

    convert_population_to_dataframes = population is not None

    def append_row(model_object, indexes, attributes_names, inputs_df):
        # function to append a row to a dataframe
        attributes_values = []
        for attribute_name in attributes_names:
            attributes_values.append(getattr(model_object, attribute_name, np.nan))
        inputs_df.loc[len(inputs_df), :] = indexes + attributes_values

    if convert_population_to_dataframes:
        # initialize the dataframes
        all_plants_df = pd.DataFrame(columns=PLANTS_VARIABLES)
        all_axes_df = pd.DataFrame(columns=AXES_VARIABLES)
        all_phytomers_df = pd.DataFrame(columns=PHYTOMERS_VARIABLES)
        all_organs_df = pd.DataFrame(columns=ORGANS_VARIABLES)
        all_hiddenzones_df = pd.DataFrame(columns=HIDDENZONE_VARIABLES)
        all_elements_df = pd.DataFrame(columns=ELEMENTS_VARIABLES)

        # run through the population tree and fill the dataframes
        for plant in population.plants:
            append_row(plant, [plant.index], simulation.Simulation.PLANTS_RUN_VARIABLES, all_plants_df)
            for axis in plant.axes:
                append_row(axis, [plant.index, axis.label], simulation.Simulation.AXES_RUN_VARIABLES, all_axes_df)
                for organs in (axis.roots, axis.xylem):
                    if organs is not None:
                        append_row(organs, [plant.index, axis.label, organs.label], simulation.Simulation.ORGANS_RUN_VARIABLES, all_organs_df)
                for phytomer in axis.phytomers:
                    append_row(phytomer, [plant.index, axis.label, phytomer.index], simulation.Simulation.PHYTOMERS_RUN_VARIABLES, all_phytomers_df)
                    if phytomer.hiddenzone is not None:
                        append_row(phytomer.hiddenzone, [plant.index, axis.label, phytomer.index], simulation.Simulation.HIDDENZONE_RUN_VARIABLES, all_hiddenzones_df)
                    for organs in (phytomer.lamina, phytomer.internode, phytomer.sheath):
                        if organs is None:
                            continue
                        for element in (organs.exposed_element, organs.enclosed_element):
                            if element is None:
                                continue
                            append_row(element, [plant.index, axis.label, phytomer.index, organs.label, element.label], simulation.Simulation.ELEMENTS_RUN_VARIABLES, all_elements_df)

        # sort the rows of the dataframes by columns
        all_plants_df.sort_values(by=PLANTS_VARIABLES, inplace=True)
        all_axes_df.sort_values(by=AXES_VARIABLES, inplace=True)
        all_phytomers_df.sort_values(by=PHYTOMERS_VARIABLES, inplace=True)
        all_organs_df.sort_values(by=ORGANS_VARIABLES, inplace=True)
        all_hiddenzones_df.sort_values(by=HIDDENZONE_VARIABLES, inplace=True)
        all_elements_df.sort_values(by=ELEMENTS_VARIABLES, inplace=True)

        # convert the indexes of plants, metamers and elements to integers in the dataframes
        all_plants_df['plant'] = all_plants_df['plant'].astype(int)
        all_axes_df['plant'] = all_axes_df['plant'].astype(int)
        all_phytomers_df[['plant', 'metamer']] = all_phytomers_df[['plant', 'metamer']].astype(int)
        all_organs_df['plant'] = all_organs_df['plant'].astype(int)
        all_hiddenzones_df[['plant', 'metamer']] = all_hiddenzones_df[['plant', 'metamer']].astype(int)
        all_elements_df[['plant', 'metamer']] = all_elements_df[['plant', 'metamer']].astype(int)

        all_plants_df.reset_index(drop=True, inplace=True)
        all_axes_df.reset_index(drop=True, inplace=True)
        all_phytomers_df.reset_index(drop=True, inplace=True)
        all_organs_df.reset_index(drop=True, inplace=True)
        all_hiddenzones_df.reset_index(drop=True, inplace=True)
        all_elements_df.reset_index(drop=True, inplace=True)


    if convert_population_to_dataframes :
        return all_plants_df, all_axes_df, all_phytomers_df, all_organs_df, all_hiddenzones_df, all_elements_df

