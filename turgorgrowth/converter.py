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

#: the columns of the outputs dataframe at SOIL scale
SOIL_VARIABLES = simulation.Simulation.SOIL_INDEXES + simulation.Simulation.SOIL_RUN_VARIABLES
SOIL_OUTPUTS_VARIABLES = SOIL_VARIABLES
SOIL_OUTPUTS_RUN_VARIABLES = simulation.Simulation.SOIL_RUN_VARIABLES

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


def from_dataframes(hiddenzones_inputs=None, elements_inputs=None, organs_inputs=None, soil_inputs=None):
    """
    If `elements_inputs` and `hiddenzones_inputs` are not `None`, converts `elements_inputs` and `hiddenzones_inputs` to a :class:`population <model.Population>`.

    :param pandas.DataFrame hiddenzones_inputs: Hidden zone inputs, with one line per hidden zone.
    :param pandas.DataFrame elements_inputs: Element inputs, with one line per element.
    :param pandas.DataFrame organs_inputs: Organs (xylem and roots) inputs, with one line per organ.
    :param pandas.DataFrame soil_inputs: Soil inputs, with one line per soil.

    :return: If `elements_inputs` and `hiddenzones_inputs` are not `None`, returns a :class:`population <model.Population>`
    :rtype: (model.Population, dict)
    """

    convert_dataframes_to_population = elements_inputs is not None and hiddenzones_inputs is not None and organs_inputs is not None
    ##convert_dataframe_to_soils_dict = soil_inputs is not None

    if convert_dataframes_to_population:
        population = model.Population()
        mapping_topology = {'predecessor': {}, 'successor': {}}

        for plant_index in elements_inputs.plant.unique():
            # create a new plant
            plant = model.Plant(plant_index)
            population.plants.append(plant)
            curr_axes_labels = elements_inputs[elements_inputs['plant'] == plant_index].axis.unique()
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
                roots = model.Roots(label='roots')  # TODO: temporary. Here, roots have a constant water potential fixed at -0.1 MPa
                last_elongated_internode = roots
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
                                    mapping_topology['successor'][phytomer.internode.enclosed_element].append(phytomer.sheath.enclosed_element)

                            elif phytomer.hiddenzone:
                                mapping_topology['predecessor'][phytomer.sheath.enclosed_element] = phytomer.hiddenzone
                                mapping_topology['successor'][phytomer.hiddenzone] = phytomer.sheath.enclosed_element
                            else:
                                mapping_topology['predecessor'][phytomer.sheath.enclosed_element] = last_elongated_internode
                                mapping_topology['successor'][last_elongated_internode].append(phytomer.sheath.enclosed_element)

                plant.axes.append(axis)

    # if convert_dataframe_to_soils_dict:
    #     soils = {}
    #     for soil_id, soil_group in soil_inputs.groupby(simulation.Simulation.SOIL_INDEXES):
    #         # create a new soil
    #         soil_attributes = soil_group.loc[soil_group.first_valid_index(), simulation.Simulation.SOIL_STATE].to_dict()
    #         soil = model.Soil(**soil_attributes) # arguments paquetes en kward qui se comporte comme un dictionnaire
    #         soils[soil_id] = soil
    #
    #
    # if convert_dataframes_to_population and convert_dataframe_to_soils_dict:
    #     return population, soils, mapping_topology
    # elif convert_dataframes_to_population:
    #     return population, mapping_topology
    # else:
    #     return soils

    # return population, mapping_topology, soils
    return population, mapping_topology


def to_dataframes(population=None, soils=None):
    """
    Convert a Turgor-Growth :class:`population <model.Population>` to Pandas dataframes.
    If `population` is not None, convert `population` to Pandas dataframes.
    If `soil` is not None, convert `soil` to Pandas dataframe.


    :param model.Population population: The Turgor-Growth population to convert.
    :param dict soil: The soil to convert.

    :return: If `population` is not None, return :class:`dataframes <pandas.DataFrame>` describing the internal state and compartments of the population at each scale:
                 * hidden zones: plant index, axis id, phytomer index, state parameters, state variables, intermediate variables,
                 fluxes and integrative variables of each hidden zone (see :mod:`HIDDENZONE_VARIABLES`)
                 * element scale: plant index, axis id, phytomer index, organs type, element type, state parameters, state variables, intermediate variables,
                 fluxes and integrative variables of each element (see :mod:`ELEMENTS_VARIABLES`)
                 * xylem scale: xylem index, state parameters, state variables, intermediate variables,
                 fluxes and integrative variables of xylem (see :mod:`XYLEM_VARIABLES`)

            and/or

            if `soil` is not None, return a :class:`dataframe <pandas.DataFrame>` describing internal state and compartments of the soil, with one line per soil:
                * plant index, axis id, state parameters, state variables, intermediate variables, fluxes and integrative variables of each soil (see :mod:`SOIL_RUN_VARIABLES`)


    :rtype: (pandas.DataFrame, pandas.DataFrame)
    """

    convert_population_to_dataframes = population is not None
    convert_soils_to_dataframe = soils is not None

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

    if convert_soils_to_dataframe:
        all_soil_df = pd.DataFrame(columns=SOIL_VARIABLES)
        for soil_id, soil in soils.items():
            append_row(soil, list(soil_id), simulation.Simulation.SOIL_RUN_VARIABLES, all_soil_df)
        all_soil_df.sort_values(by=SOIL_VARIABLES, inplace=True)
        all_soil_df['plant'] = all_soil_df['plant'].astype(int)
        all_soil_df.reset_index(drop=True, inplace=True)
    
    if convert_population_to_dataframes and convert_soils_to_dataframe:
        return all_plants_df, all_axes_df, all_phytomers_df, all_organs_df, all_hiddenzones_df, all_elements_df, all_soil_df
    elif convert_population_to_dataframes:
        return all_plants_df, all_axes_df, all_phytomers_df, all_organs_df, all_hiddenzones_df, all_elements_df
    else:
        return all_soil_df

    ##return all_plants_df, all_axes_df, all_phytomers_df, all_organs_df, all_hiddenzones_df, all_elements_df
