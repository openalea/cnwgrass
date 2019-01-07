# -*- coding: latin-1 -*-
import numpy as np
import pandas as pd

from turgorgrowth import model, simulation

"""
    turgorgrowth.converter
    ~~~~~~~~~~~~~~~~~

    The module :mod:`turgorgrowth.converter` defines functions to convert
    :class:`dataframes <pandas.DataFrame>` to/from Turgor-Growth inputs or outputs format.

    :license: CeCILL-C, see LICENSE for details.
"""

#: the columns of the outputs dataframe at PLANT scale
PLANTS_VARIABLES = simulation.Simulation.PLANTS_INDEXES + simulation.Simulation.PLANTS_RUN_VARIABLES

#: the columns of the outputs dataframe at AXIS scale
AXES_VARIABLES = simulation.Simulation.AXES_INDEXES + simulation.Simulation.AXES_RUN_VARIABLES

#: the columns of the outputs dataframe at PHYTOMER scale
PHYTOMERS_VARIABLES = simulation.Simulation.PHYTOMERS_INDEXES + simulation.Simulation.PHYTOMERS_RUN_VARIABLES

#: the columns of the outputs dataframe at HIDDENZONE scale
HIDDENZONE_VARIABLES = simulation.Simulation.HIDDENZONE_INDEXES + simulation.Simulation.HIDDENZONE_RUN_VARIABLES

#: the columns of the outputs dataframe at ORGAN scale
ELEMENTS_VARIABLES = simulation.Simulation.ELEMENTS_INDEXES + simulation.Simulation.ELEMENTS_RUN_VARIABLES

#: the mapping of the Turgor-Growth organ classes to organ names in MTG
TURGORGROWTH_CLASSES_TO_DATAFRAME_ORGANS_MAPPING = {model.Internode: 'internode', model.Lamina: 'blade', model.Sheath: 'sheath',  model.HiddenZone: 'hiddenzone'}

#: the mapping of the name of each element, from Dataframe to Turgor-Growth
DATAFRAME_TO_TURGORGROWTH_ELEMENTS_NAMES_MAPPING = {'HiddenElement': 'enclosed_element', 'StemElement': 'exposed_element', 'LeafElement1': 'exposed_element'}


def from_dataframes(hiddenzones_inputs=None, elements_inputs=None):
    """
    If `elements_inputs` and `hiddenzones_inputs` are not `None`, converts `elements_inputs` and `hiddenzones_inputs` to a :class:`population <model.Population>`.

    :Parameters:
        - `hiddenzones_inputs` (:class:`pandas.DataFrame`) - Hidden zone inputs, with one line per hidden zone.
        - `elements_inputs` (:class:`pandas.DataFrame`) - Element inputs, with one line per element.

    :Returns:
        If `elements_inputs` and `hiddenzones_inputs` are not `None`, returns a :class:`population <model.Population>`,
    :Returns Type:
        :class:`tuple`.

    """

    convert_dataframes_to_population = elements_inputs is not None and hiddenzones_inputs is not None

    if convert_dataframes_to_population:
        population = model.Population()
        mapping_topology = {'predecessor': {}, 'successor': {}}

        for plant_index in elements_inputs.plant.unique():
            # create a new plant
            plant = model.Plant(plant_index)
            population.plants.append(plant)
            curr_axes_labels = elements_inputs[elements_inputs['plant'] == plant_index].axis.unique()
            for axis_label in curr_axes_labels:
                roots = model.Roots(label='roots')  # TODO: temporary. Here, roots have a constant water potential fixed at -0.1 MPa
                last_elongated_internode = roots
                mapping_topology['successor'][roots] = []
                # create a new axis
                axis = model.Axis(axis_label)
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

                        organ_label = TURGORGROWTH_CLASSES_TO_DATAFRAME_ORGANS_MAPPING[phytomer_attribute_class]

                        if metamer_index in curr_metamers_indexes_for_elements:
                            curr_elements_inputs = elements_inputs[(elements_inputs['plant'] == plant_index) & (elements_inputs['axis'] == axis_label) & (elements_inputs['metamer'] == metamer_index) &
                                                                   (elements_inputs['organ'] == organ_label)]
                            if organ_label not in curr_elements_inputs.organ.values:
                                continue
                            # create a new organ
                            organ = phytomer_attribute_class(label=organ_label)
                            setattr(phytomer, phytomer_attribute_name, organ)

                            for mtg_element_label, turgorgrowth_element_name in DATAFRAME_TO_TURGORGROWTH_ELEMENTS_NAMES_MAPPING.items():
                                element_inputs = curr_elements_inputs[curr_elements_inputs['element'] == mtg_element_label]
                                if len(element_inputs) == 0:
                                    continue
                                element_inputs = element_inputs.loc[:, simulation.Simulation.ELEMENTS_STATE]
                                element_dict = element_inputs.loc[element_inputs.first_valid_index()].dropna().to_dict()
                                # create a new element
                                element = phytomer_attribute_element_class(label=mtg_element_label, **element_dict)
                                setattr(organ, turgorgrowth_element_name, element)

                    #: Hidden zones
                    if metamer_index in curr_metamers_indexes_for_hiddenzones:
                        hiddenzone_inputs = hiddenzones_inputs[(hiddenzones_inputs['plant'] == plant_index) & (hiddenzones_inputs['axis'] == axis_label) &
                                                               (hiddenzones_inputs['metamer'] == metamer_index)]
                        if len(hiddenzone_inputs) == 0:
                            continue
                        hiddenzone_inputs = hiddenzone_inputs.loc[:, simulation.Simulation.HIDDENZONE_STATE]
                        hiddenzone_dict = hiddenzone_inputs.loc[hiddenzone_inputs.first_valid_index()].to_dict()
                        # create a new hidden zone
                        hiddenzone = model.HiddenZone(TURGORGROWTH_CLASSES_TO_DATAFRAME_ORGANS_MAPPING[model.HiddenZone], **hiddenzone_dict)
                        phytomer.hiddenzone = hiddenzone
                        # Topology
                        mapping_topology['predecessor'][phytomer.hiddenzone] = last_elongated_internode
                        mapping_topology['successor'][last_elongated_internode].append(phytomer.hiddenzone)

                    # Topoly of elements
                    if phytomer.lamina and phytomer.lamina.exposed_element:
                        if phytomer.sheath:
                            if phytomer.sheath.exposed_element:
                                mapping_topology['predecessor'][phytomer.lamina.exposed_element] = phytomer.sheath.exposed_element
                                mapping_topology['successor'][phytomer.sheath.exposed_element] = phytomer.lamina.exposed_element
                            elif phytomer.sheath.enclosed_element:
                                mapping_topology['predecessor'][phytomer.lamina.exposed_element] = phytomer.sheath.enclosed_element
                                mapping_topology['successor'][phytomer.sheath.enclosed_element] = phytomer.lamina.exposed_element
                        else:
                            mapping_topology['predecessor'][phytomer.lamina.exposed_element] = phytomer.hiddenzone
                            mapping_topology['successor'][phytomer.hiddenzone] = phytomer.lamina.exposed_element

                    if phytomer.internode:
                        if phytomer.internode.enclosed_element:
                            mapping_topology['predecessor'][phytomer.internode.enclosed_element] = last_elongated_internode
                            mapping_topology['successor'][last_elongated_internode].append(phytomer.internode.enclosed_element)
                            last_elongated_internode = phytomer.internode.enclosed_element
                            mapping_topology['successor'][phytomer.internode.enclosed_element] = []

                        if phytomer.internode.exposed_element:
                            if phytomer.internode.enclosed_element:
                                mapping_topology['predecessor'][phytomer.internode.exposed_element] = phytomer.internode.enclosed_element
                                mapping_topology['successor'][phytomer.internode.enclosed_element] = phytomer.internode.exposed_element
                            else:
                                mapping_topology['predecessor'][phytomer.internode.exposed_element] = phytomer.hiddenzone
                                mapping_topology['successor'][phytomer.hiddenzone] = phytomer.internode.exposed_element

                    if phytomer.sheath:
                        if phytomer.sheath.exposed_element:
                            if phytomer.sheath.enclosed_element:
                                mapping_topology['predecessor'][phytomer.sheath.exposed_element] = phytomer.sheath.enclosed_element
                                mapping_topology['successor'][phytomer.sheath.enclosed_element] = phytomer.sheath.exposed_element
                            elif phytomer.hiddenzone:
                                mapping_topology['predecessor'][phytomer.sheath.exposed_element] = phytomer.hiddenzone
                                mapping_topology['successor'][phytomer.hiddenzone] = phytomer.sheath.exposed_element
                            else:
                                mapping_topology['predecessor'][phytomer.sheath.exposed_element] = last_elongated_internode
                                mapping_topology['successor'][last_elongated_internode].append(phytomer.sheath.exposed_element)
                        if phytomer.sheath.enclosed_element:
                            if phytomer.internode:
                                if phytomer.internode.exposed_element:
                                    mapping_topology['predecessor'][phytomer.sheath.enclosed_element] = phytomer.internode.exposed_element
                                    mapping_topology['successor'][phytomer.internode.exposed_element] = phytomer.sheath.enclosed_element
                                else:
                                    mapping_topology['predecessor'][phytomer.sheath.enclosed_element] = phytomer.internode.enclosed_element
                                    mapping_topology['successor'][phytomer.internode.enclosed_element] = phytomer.sheath.enclosed_element

                            elif phytomer.hiddenzone:
                                mapping_topology['predecessor'][phytomer.sheath.enclosed_element] = phytomer.hiddenzone
                                mapping_topology['successor'][phytomer.hiddenzone] = phytomer.sheath.enclosed_element
                            else:
                                mapping_topology['predecessor'][phytomer.sheath.enclosed_element] = last_elongated_internode
                            mapping_topology['successor'][last_elongated_internode].append(phytomer.sheath.enclosed_element)

                plant.axes.append(axis)

        return population, mapping_topology


def to_dataframes(population=None):
    """
    Convert a Turgor-Growth :class:`population <model.Population>` to Pandas dataframes.

    If `population` is not None, convert `population` to Pandas dataframes.

    :Parameters:

        - `population` (:class:`model.Population`) - The Turgor-Growth population to convert.

    :Returns:
        If `population` is not None, return :class:`dataframes <pandas.DataFrame>` describing the internal state and compartments of the population at each scale:

            * hidden zones: plant index, axis id, phytomer index, state parameters, state variables, intermediate variables, fluxes and integrative variables of each hidden zone (see :mod:`HIDDENZONE_VARIABLES`)
            * element scale: plant index, axis id, phytomer index, organ type, element type, state parameters, state variables, intermediate variables, fluxes and integrative variables of each element (see :mod:`ELEMENTS_VARIABLES`)

    :Returns Type:
        :class:`tuple` of :class:`pandas.DataFrame`
        or
        :class:`pandas.DataFrame`

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
        all_hiddenzones_df = pd.DataFrame(columns=HIDDENZONE_VARIABLES)
        all_elements_df = pd.DataFrame(columns=ELEMENTS_VARIABLES)

        # run through the population tree and fill the dataframes
        for plant in population.plants:
            for axis in plant.axes:
                # for organ in (axis.roots, axis.phloem, axis.grains):
                #     if organ is not None:
                #         append_row(organ, [plant.index, axis.label, organ.label], simulation.Simulation.ORGANS_RUN_VARIABLES, all_organs_df)
                for phytomer in axis.phytomers:
                    if phytomer.hiddenzone is not None:
                        append_row(phytomer.hiddenzone, [plant.index, axis.label, phytomer.index], simulation.Simulation.HIDDENZONE_RUN_VARIABLES, all_hiddenzones_df)
                    for organ in (phytomer.lamina, phytomer.internode, phytomer.sheath):
                        if organ is None:
                            continue
                        for element in (organ.exposed_element, organ.enclosed_element):
                            if element is None:
                                continue
                            append_row(element, [plant.index, axis.label, phytomer.index, organ.label, element.label], simulation.Simulation.ELEMENTS_RUN_VARIABLES, all_elements_df)

        # sort the rows of the dataframes by columns
        all_hiddenzones_df.sort_values(by=HIDDENZONE_VARIABLES, inplace=True)
        all_elements_df.sort_values(by=ELEMENTS_VARIABLES, inplace=True)

        # convert the indexes of plants, metamers and elements to integers in the dataframes
        all_hiddenzones_df[['plant', 'metamer']] = all_hiddenzones_df[['plant', 'metamer']].astype(int)
        all_elements_df[['plant', 'metamer']] = all_elements_df[['plant', 'metamer']].astype(int)

        all_hiddenzones_df.reset_index(drop=True, inplace=True)
        all_elements_df.reset_index(drop=True, inplace=True)

        return all_hiddenzones_df, all_elements_df
