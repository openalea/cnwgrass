# -*- coding: latin-1 -*-
from cnwheat import simulation as cnwheat_simulation
from elongwheat import simulation as elongwheat_simulation
from farquharwheat import converter as farquharwheat_converter
from growthwheat import simulation as growthwheat_simulation
from senescwheat import converter as senescwheat_converter
import numpy as np
import pandas as pd

"""
    fspmwheat.fspmwheat_facade
    ~~~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`fspmwheat.fspmwheat_facade` is a facade of the model FSPMWeat.

    :copyright: Copyright 2014-2016 INRA-ECOSYS, see AUTHORS.
    :license: see LICENSE for details.

"""

"""
    Information about this versioned file:
        $LastChangedBy: mngauthier $
        $LastChangedDate: 2020-03-03 16:00:54 +0100 (mar., 03 mars 2020) $
        $LastChangedRevision: 238 $
        $URL: https://subversion.renater.fr/authscm/mngauthier/svn/fspm-wheat/branches/multirun/fspmwheat/cnwheat_facade.py $
        $Id: cnwheat_facade.py 238 2020-03-03 15:00:54Z mngauthier $
"""

#: columns which define the topology in the input/output dataframe
AXES_TOPOLOGY_COLUMNS = ['plant', 'axis']
ELEMENTS_TOPOLOGY_COLUMNS = ['plant', 'axis', 'metamer', 'organ', 'element']  # Mature + emerging elements
HIDDENZONES_TOPOLOGY_COLUMNS = ['plant', 'axis', 'metamer']
ORGANS_TOPOLOGY_COLUMNS = ['plant', 'axis', 'organ']
SOILS_TOPOLOGY_COLUMNS = ['plant', 'axis']

#: variables in each input/output dataframe
AXES_VARIABLES = set(cnwheat_simulation.Simulation.AXES_RUN_VARIABLES +
                     elongwheat_simulation.AXIS_INPUTS_OUTPUTS +
                     list(growthwheat_simulation.AXIS_INPUTS_OUTPUTS) +
                     senescwheat_converter.SENESCWHEAT_AXES_INPUTS_OUTPUTS)
ELEMENTS_VARIABLES = set(cnwheat_simulation.Simulation.ELEMENTS_RUN_VARIABLES +
                         elongwheat_simulation.ELEMENT_INPUTS_OUTPUTS +
                         list(farquharwheat_converter.FARQUHARWHEAT_ELEMENTS_INPUTS_OUTPUTS) +
                         list(growthwheat_simulation.ELEMENT_INPUTS_OUTPUTS) +
                         senescwheat_converter.SENESCWHEAT_ELEMENTS_INPUTS_OUTPUTS)
HIDDENZONES_VARIABLES = set(cnwheat_simulation.Simulation.HIDDENZONE_RUN_VARIABLES +
                            elongwheat_simulation.HIDDENZONE_INPUTS_OUTPUTS +
                            list(growthwheat_simulation.HIDDENZONE_INPUTS_OUTPUTS))
ORGANS_VARIABLES = set(cnwheat_simulation.Simulation.ORGANS_RUN_VARIABLES +
                       list(growthwheat_simulation.ROOT_INPUTS_OUTPUTS) +
                       senescwheat_converter.SENESCWHEAT_ROOTS_INPUTS_OUTPUTS)
SOILS_VARIABLES = set(cnwheat_simulation.Simulation.SOILS_RUN_VARIABLES)

BOTANICAL_ORGANS_AT_AXIS_SCALE = ['roots', 'phloem', 'grains']
BOTANICAL_COMPARTMENTS_AT_AXIS_SCALE = BOTANICAL_ORGANS_AT_AXIS_SCALE + ['soil']


class FSPMWheatFacade(object):
    """
    The FSPMWheatFacade class permits to ...
    from a :class:`MTG <openalea.mtg.mtg.MTG>`, and update the MTG and the dataframes
    shared between all models.

    """

    def __init__(self, shared_mtg):
        # shared_axes_inputs_outputs_df,
        # shared_organs_inputs_outputs_df,
        # shared_hiddenzones_inputs_outputs_df,
        # shared_elements_inputs_outputs_df,
        # shared_soils_inputs_outputs_df,
        # update_shared_df = True):
        """
        :param openalea.mtg.mtg.MTG shared_mtg: The MTG shared between all models.
        :param pandas.DataFrame shared_axes_inputs_outputs_df: the dataframe of inputs and outputs at axes scale shared between all models.
        :param pandas.DataFrame shared_organs_inputs_outputs_df: the dataframe of inputs and outputs at organs scale shared between all models.
        :param pandas.DataFrame shared_hiddenzones_inputs_outputs_df: the dataframe of inputs and outputs at hiddenzones scale shared between all models.
        :param pandas.DataFrame shared_elements_inputs_outputs_df: the dataframe of inputs and outputs at elements scale shared between all models.
        :param pandas.DataFrame shared_soils_inputs_outputs_df: the dataframe of inputs and outputs at soils scale shared between all models.
        :param bool update_shared_df: If `True`  update the shared dataframes at init and at each run (unless stated otherwise)
        """

        self._shared_mtg = shared_mtg  #: the MTG shared between all models

        # self._shared_axes_inputs_outputs_df = shared_axes_inputs_outputs_df  #: the dataframe at axes scale shared between all models
        # self._shared_organs_inputs_outputs_df = shared_organs_inputs_outputs_df  #: the dataframe at organs scale shared between all models
        # self._shared_hiddenzones_inputs_outputs_df = shared_hiddenzones_inputs_outputs_df  #: the dataframe at hiddenzones scale shared between all models
        # self._shared_elements_inputs_outputs_df = shared_elements_inputs_outputs_df  #: the dataframe at elements scale shared between all models
        # self._shared_soils_inputs_outputs_df = shared_soils_inputs_outputs_df  #: the dataframe at soils scale shared between all models
        # self._update_shared_df = update_shared_df
        # if self._update_shared_df:
        #     self._update_shared_dataframes(cnwheat_organs_data_df=model_organs_inputs_df,
        #                                    cnwheat_hiddenzones_data_df=model_hiddenzones_inputs_df,
        #                                    cnwheat_elements_data_df=model_elements_inputs_df,
        #                                    cnwheat_soils_data_df=model_soils_inputs_df)

    def _read_outputs_on_MTG(self):
        """
        Extract the outputs of all sub-models from the MTG shared between all models.
        """

        axes_dict = {}
        elements_dict = {}
        hiddenzones_dict = {}
        organs_dict = {}
        soils_dict = {}

        for mtg_plant_vid in self._shared_mtg.components_iter(self._shared_mtg.root):
            mtg_plant_index = int(self._shared_mtg.index(mtg_plant_vid))

            # Axis scale
            for mtg_axis_vid in self._shared_mtg.components_iter(mtg_plant_vid):
                mtg_axis_label = self._shared_mtg.label(mtg_axis_vid)
                mtg_axis_properties = self._shared_mtg.get_vertex_property(mtg_axis_vid)
                axis_id = (mtg_plant_index, mtg_axis_label)
                axis_dict = {}

                if mtg_axis_properties.get('nb_leaves') is None:
                    continue
                for axis_run_variable in AXES_VARIABLES:
                    if axis_run_variable in mtg_axis_properties:
                        # use the input from the MTG
                        axis_dict[axis_run_variable] = mtg_axis_properties[axis_run_variable]
                axes_dict[axis_id] = axis_dict

                # Botanical organs at axis scale
                for botanical_organ_name in BOTANICAL_ORGANS_AT_AXIS_SCALE:
                    if botanical_organ_name in mtg_axis_properties:
                        organ_id = (mtg_plant_index, mtg_axis_label, botanical_organ_name)
                        mtg_organ_properties = mtg_axis_properties[botanical_organ_name]
                        if mtg_organ_properties.get('sucrose') is None:
                            continue
                        organ_dict = {}
                        for organ_run_variable in ORGANS_VARIABLES:
                            if organ_run_variable in mtg_organ_properties:
                                organ_dict[organ_run_variable] = mtg_organ_properties[organ_run_variable]
                        organs_dict[organ_id] = organ_dict

                # Soil at axis scale
                if 'soil' in mtg_axis_properties:
                    mtg_soil_properties = mtg_axis_properties['soil']
                    soil_dict = {}
                    for soil_run_variable in SOILS_VARIABLES:
                        if soil_run_variable in mtg_soil_properties:
                            soil_dict[soil_run_variable] = mtg_soil_properties[soil_run_variable]
                    soils_dict[axis_id] = soil_dict

                # Metamer scale
                for mtg_metamer_vid in self._shared_mtg.components_iter(mtg_axis_vid):
                    mtg_metamer_index = int(self._shared_mtg.index(mtg_metamer_vid))
                    mtg_metamer_properties = self._shared_mtg.get_vertex_property(mtg_metamer_vid)
                    if 'hiddenzone' in mtg_metamer_properties:
                        hiddenzone_id = (mtg_plant_index, mtg_axis_label, mtg_metamer_index)
                        mtg_hiddenzone_properties = mtg_metamer_properties['hiddenzone']
                        hiddenzone_dict = {}
                        for hiddenzone_run_variable in HIDDENZONES_VARIABLES:
                            if hiddenzone_run_variable in mtg_hiddenzone_properties:
                                # use the input from the MTG
                                hiddenzone_dict[hiddenzone_run_variable] = mtg_hiddenzone_properties[hiddenzone_run_variable]
                        hiddenzones_dict[hiddenzone_id] = hiddenzone_dict

                    # Photosynthetic organ scale
                    for mtg_organ_vid in self._shared_mtg.components_iter(mtg_metamer_vid):
                        mtg_organ_label = self._shared_mtg.label(mtg_organ_vid)
                        # Element scale
                        for mtg_element_vid in self._shared_mtg.components_iter(mtg_organ_vid):
                            mtg_element_label = self._shared_mtg.label(mtg_element_vid)
                            mtg_element_properties = self._shared_mtg.get_vertex_property(mtg_element_vid)
                            if np.nan_to_num(self._shared_mtg.property('length').get(mtg_element_vid, 0)) == 0:
                                continue
                            element_dict = {}
                            for elongwheat_element_run_variable in ELEMENTS_VARIABLES:
                                element_dict[elongwheat_element_run_variable] = mtg_element_properties.get(elongwheat_element_run_variable, np.nan)
                            element_id = (mtg_plant_index, mtg_axis_label, mtg_metamer_index, mtg_organ_label, mtg_element_label)
                            elements_dict[element_id] = element_dict

        return {'axes': axes_dict, 'elements': elements_dict, 'hiddenzones': hiddenzones_dict, 'organs': organs_dict, 'soils': soils_dict}

    @staticmethod
    def _to_dataframes(data_dict):
        """
        Convert outputs from _read_outputs_on_MTG() which are dictionaries to Pandas dataframes.

        :param dict data_dict: outputs from _read_outputs_on_MTG() which are dictionaries

        :return: Five dataframes: axes, elements, hiddenzones, organs, soils
        :rtype: (pandas.DataFrame, pandas.DataFrame, pandas.DataFrame, pandas.DataFrame, pandas.DataFrame)
        """
        dataframes_dict = {}
        for (current_key, current_topology_columns, current_outputs_names) in (('soils', SOILS_TOPOLOGY_COLUMNS, SOILS_VARIABLES),
                                                                               ('organs', ORGANS_TOPOLOGY_COLUMNS, ORGANS_VARIABLES),
                                                                               ('hiddenzones', HIDDENZONES_TOPOLOGY_COLUMNS, HIDDENZONES_VARIABLES),
                                                                               ('elements', ELEMENTS_TOPOLOGY_COLUMNS, ELEMENTS_VARIABLES),
                                                                               ('axes', AXES_TOPOLOGY_COLUMNS, AXES_VARIABLES)):
            current_data_dict = data_dict[current_key]
            current_ids_df = pd.DataFrame(current_data_dict.keys(), columns=current_topology_columns)
            current_data_df = pd.DataFrame(current_data_dict.values())
            current_df = pd.concat([current_ids_df, current_data_df], axis=1)
            current_df.sort_values(by=current_topology_columns, inplace=True)
            current_columns_sorted = current_topology_columns + list(current_outputs_names)
            current_df = current_df.reindex(current_columns_sorted, axis=1, copy=False)
            # Reset dtypes
            current_df = pd.DataFrame(current_df.where(current_df.notnull(), np.nan).values.tolist(), columns=current_df.columns)
            current_df.reset_index(drop=True, inplace=True)
            dataframes_dict[current_key] = current_df

        return dataframes_dict['axes'], dataframes_dict['elements'], dataframes_dict['hiddenzones'], dataframes_dict['organs'], dataframes_dict['soils']

    def build_outputs_df_from_MTG(self):
        outputs_dict = self._read_outputs_on_MTG()
        return self._to_dataframes(outputs_dict)
