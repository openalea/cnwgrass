# -*- coding: latin-1 -*-

import pandas as pd
import numpy as np
import warnings

from alinea.caribu.CaribuScene import CaribuScene
from alinea.caribu.sky_tools import GenSky, GetLight, Gensun, GetLightsSun, spitters_horaire

from fspmwheat import tools

"""
    fspmwheat.caribu_facade
    ~~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`fspmwheat.caribu_facade` is a facade of the model Caribu.

    This module permits to initialize and run the model Caribu from a :class:`MTG <openalea.mtg.mtg.MTG>`
    in a convenient and transparent way, wrapping all the internal complexity of the model, and dealing
    with all the tedious initialization and conversion processes.

    :copyright: Copyright 2014-2016 INRA-ECOSYS, see AUTHORS.
    :license: TODO, see LICENSE for details.

    .. seealso:: Barillot et al. 2016.
"""

#: the columns which define the topology in the elements scale dataframe shared between all models
SHARED_ELEMENTS_INPUTS_OUTPUTS_INDEXES = ['plant', 'axis', 'metamer', 'organ', 'element']

#: the outputs of Caribu
CARIBU_OUTPUTS = ['PARa', 'Erel', 'PARa_prim', 'Erel_prim']


class CaribuFacade(object):
    """
    The CaribuFacade class permits to initialize, run the model Caribu
    from a :class:`MTG <openalea.mtg.mtg.MTG>`, and update the MTG and the dataframes
    shared between all models.

    Use :meth:`run` to run the model.
    """

    def __init__(self,
                 shared_mtg,
                 shared_elements_inputs_outputs_df,
                 geometrical_model,
                 update_shared_df=True):
        """
        :param openalea.mtg.MTG shared_mtg: The MTG shared between all models.
        :param pandas.DataFrame shared_elements_inputs_outputs_df: The dataframe of inputs and outputs at elements scale shared between all models.
        :param alinea.adel.adel_dynamic.AdelWheatDyn geometrical_model: The model which deals with geometry. This model must have an attribute "domain".
        :param bool update_shared_df: If `True`  update the shared dataframes at init and at each run (unless stated otherwise)
        """
        self._shared_mtg = shared_mtg  #: the MTG shared between all models
        self._shared_elements_inputs_outputs_df = shared_elements_inputs_outputs_df  #: the dataframe at elements scale shared between all models
        self._geometrical_model = geometrical_model  #: the model which deals with geometry
        self._alea_canopy = pd.DataFrame()  #: alea table to generate the heterogeneous canopy
        self._update_shared_df = update_shared_df

    def run(self, run_caribu, sun_sky_option='mix', energy=1, DOY=1, hourTU=12, latitude=48.85, diffuse_model='soc', azimuts=4, zenits=5, heterogeneous_canopy=False,
            plant_density=250., inter_row=0.15, update_shared_df=None, prim_scale=False):
        """
        Run the model and update the MTG and the dataframes shared between all models.

        :param bool run_caribu: If 'True', run the CARIBU model to calculate light distribution inside the 3D canopy.
        :param str sun_sky_option: The irradiance model, should be one of 'mix' or 'sun' or 'sky'
        :param float energy: The incident PAR above the canopy (µmol m-2 s-1)
        :param int DOY: Day Of the Year to be used for solar sources
        :param int hourTU: Hour to be used for solar sources (Universal Time)
        :param float latitude: latitude to be used for solar sources (°)
        :param string diffuse_model: The kind of diffuse model, either 'soc' or 'uoc'.
        :param int azimuts: The number of azimutal positions.
        :param int zenits: The number of zenital positions.
        :param bool heterogeneous_canopy: Whether to create a duplicated heterogeneous canopy from the initial mtg.
        :param float plant_density: Number of plant per m2 in the stand (plant m-2).
        :param float inter_row: Inter-row spacing in the stand (m).
        :param bool update_shared_df: if 'True', update the shared dataframes at this time step.
        :param bool prim_scale: If True, light distribution output at primitive scale, if not at organ scale
        """
        c_scene_sky, c_scene_sun, Erel_input, Erel_input_prim = self._initialize_model(run_caribu,
                                                                                       1,
                                                                                       diffuse_model,
                                                                                       azimuts,
                                                                                       zenits,
                                                                                       DOY,
                                                                                       hourTU,
                                                                                       latitude,
                                                                                       heterogeneous_canopy,
                                                                                       plant_density,
                                                                                       inter_row)
        outputs = {}
        if run_caribu:
            #: Diffuse light
            if sun_sky_option == 'sky':
                raw, aggregated_sky = c_scene_sky.run(direct=True, infinite=True)
                Erel_sky = aggregated_sky['par']['Eabs']  #: Erel is the relative surfacic absorbed energy per organ
                PARa_sky = {k: v * energy for k, v in Erel_sky.items()}
                Erel_output = Erel_sky
                PARa_output = PARa_sky

                # Primitive scale
                if prim_scale:
                    Erel_prim = raw['par']['Eabs']
                    raw_Eabs_abs = {k: [Eabs * energy for Eabs in raw['par']['Eabs'][k]] for k in raw['par']['Eabs']}
                    outputs.update({'Erel_prim': Erel_prim, 'PARa_prim': raw_Eabs_abs, 'area_prim': raw['par']['area']})

            #: Direct light
            elif sun_sky_option == 'sun':
                raw, aggregated_sun = c_scene_sun.run(direct=True, infinite=True)
                Erel_sun = aggregated_sun['par']['Eabs']  #: Erel is the relative surfacic absorbed energy per organ
                PARa_sun = {k: v * energy for k, v in Erel_sun.items()}
                Erel_output = Erel_sun
                PARa_output = PARa_sun

                # Primitive scale
                if prim_scale:
                    Erel_prim = raw['par']['Eabs']
                    raw_Eabs_abs = {k: [Eabs * energy for Eabs in raw['par']['Eabs'][k]] for k in raw['par']['Eabs']}
                    outputs.update({'Erel_prim': Erel_prim, 'PARa_prim': raw_Eabs_abs, 'area_prim': raw['par']['area']})

            #: Mix sky-Sun
            elif sun_sky_option == 'mix':
                #: Diffuse
                raw_sky, aggregated_sky = c_scene_sky.run(direct=True, infinite=True)
                Erel_sky = aggregated_sky['par']['Eabs']
                #: Direct
                raw_sun, aggregated_sun = c_scene_sun.run(direct=True, infinite=True)
                Erel_sun = aggregated_sun['par']['Eabs']

                #: Spitters's model estimating for the diffuse:direct ratio
                Rg = energy / 2.02  #: Global Radiation (W.m-2)
                RdRs = spitters_horaire.RdRsH(Rg=Rg, DOY=DOY, heureTU=hourTU, latitude=latitude)  #: Diffuse fraction of the global irradiance
                Erel = {}
                for element_id, Erel_value in Erel_sky.items():
                    Erel[element_id] = RdRs * Erel_value + (1 - RdRs) * Erel_sun[element_id]

                Erel_output = Erel
                PARa_output = {k: v * energy for k, v in Erel_output.items()}

                # Primitive scale
                if prim_scale:
                    raise ValueError("prim_scale not yet implemented for mix sun_sky_option.")

            else:
                raise ValueError("Unknown sun_sky_option : can be either 'mix', 'sun' or 'sky'.")

            # Ouputs
            outputs.update({'PARa': PARa_output, 'Erel': Erel_output})
        else:
            PARa_output = {k: v * energy for k, v in Erel_input.items()}
            raw_Eabs_abs = {k: [Eabs * energy for Eabs in Erel_input_prim[k]] for k in Erel_input_prim}
            outputs.update({'PARa': PARa_output})

            # Primitive scale
            if prim_scale:
                outputs.update({'PARa_prim': raw_Eabs_abs})

        # Updates
        self.update_shared_MTG(outputs)
        if update_shared_df or (update_shared_df is None and self._update_shared_df):
            self.update_shared_dataframes(outputs)

    def _initialize_model(self, run_caribu, energy, diffuse_model, azimuts, zenits, DOY, hourTU, latitude, heterogeneous_canopy, plant_density, inter_row):
        """
        Initialize the inputs of the model from the MTG shared

        :param bool run_caribu: If 'True', run the CARIBU model to calculate light distribution inside the 3D canopy.
        :param float energy: The incident PAR above the canopy (µmol m-2 s-1)
        :param string diffuse_model: The kind of diffuse model, either 'soc' or 'uoc'.
        :param int azimuts: The number of azimutal positions.
        :param int zenits: The number of zenital positions.
        :param int DOY: Day Of the Year to be used for solar sources
        :param int hourTU: Hour to be used for solar sources (Universal Time)
        :param float latitude: latitude to be used for solar sources (°)
        :param bool heterogeneous_canopy: Whether to create a duplicated heterogeneous canopy from the initial mtg.

        :return: A tuple of Caribu scenes instantiated for sky and sun sources, respectively, and two dictionaries with Erel value per vertex id and per primitive.
        :rtype: (CaribuScene, CaribuScene, dict, dict)
        """
        c_scene_sky, c_scene_sun, Erel, Erel_prim = None, None, None, None

        if run_caribu:

            #: Diffuse light sources : Get the energy and positions of the source for each sector as a string
            sky_string = GetLight.GetLight(GenSky.GenSky()(energy, diffuse_model, azimuts, zenits))  #: (Energy, soc/uoc, azimuts, zenits)

            # Convert string to list in order to be compatible with CaribuScene input format
            sky = []
            for string in sky_string.split('\n'):
                if len(string) != 0:
                    string_split = string.split(' ')
                    t = tuple((float(string_split[0]), tuple((float(string_split[1]), float(string_split[2]), float(string_split[3])))))
                    sky.append(t)

            #: Direct light sources (sun positions)
            sun = Gensun.Gensun()(energy, DOY, hourTU, latitude)
            sun = GetLightsSun.GetLightsSun(sun)
            sun_str_split = sun.split(' ')
            sun = [tuple((float(sun_str_split[0]), tuple((float(sun_str_split[1]), float(sun_str_split[2]), float(sun_str_split[3])))))]

            #: Optical properties
            opt = {'par': {}}
            geom = self._shared_mtg.property('geometry')

            for vid in geom.keys():
                if self._shared_mtg.class_name(vid) in ('LeafElement1', 'LeafElement'):
                    opt['par'][vid] = (0.10, 0.05)  #: (reflectance, transmittance) of the adaxial side of the leaves
                elif self._shared_mtg.class_name(vid) == 'StemElement':
                    opt['par'][vid] = (0.10,)  #: (reflectance,) of the stems
                else:
                    warnings.warn('Warning: unknown element type {}, vid={}'.format(self._shared_mtg.class_name(vid), vid))

            #: Generates CaribuScenes
            if not heterogeneous_canopy:  # TODO: adapt the domain to plant_density
                c_scene_sky = CaribuScene(scene=self._shared_mtg, light=sky, pattern=self._geometrical_model.domain, opt=opt)
                c_scene_sun = CaribuScene(scene=self._shared_mtg, light=sun, pattern=self._geometrical_model.domain, opt=opt)
            else:
                duplicated_scene, domain = self._create_heterogeneous_canopy(plant_density=plant_density, inter_row=inter_row)
                c_scene_sky = CaribuScene(scene=duplicated_scene, light=sky, pattern=domain, opt=opt)
                c_scene_sun = CaribuScene(scene=duplicated_scene, light=sun, pattern=domain, opt=opt)

        else:
            Erel = self._shared_mtg.property('Erel')
            Erel_prim = self._shared_mtg.property('Erel_prim')

        return c_scene_sky, c_scene_sun, Erel, Erel_prim

    def _create_heterogeneous_canopy(self, nplants=50, var_plant_position=0.03, var_leaf_inclination=0.157, var_leaf_azimut=1.57, var_stem_azimut=0.157,
                                     plant_density=250, inter_row=0.15):
        """
        Duplicate a plant in order to obtain a heterogeneous canopy.

        :param int nplants: the desired number of duplicated plants
        :param float var_plant_position: variability for plant position (m)
        :param float var_leaf_inclination: variability for leaf inclination (rad)
        :param float var_leaf_azimut: variability for leaf azimut (rad)
        :param float var_stem_azimut: variability for stem azimut (rad)

        :return: duplicated heterogenous scene and its domain
        :rtype: openalea.plantgl.all.Scene, (float)
        """
        from alinea.adel.Stand import AgronomicStand
        import openalea.plantgl.all as plantgl
        import random

        # Load scene
        initial_scene = self._geometrical_model.scene(self._shared_mtg)

        # Planter
        stand = AgronomicStand(sowing_density=plant_density, plant_density=plant_density, inter_row=inter_row, noise=var_plant_position)
        _, domain, positions, _ = stand.smart_stand(nplants=nplants, at=inter_row, convunit=1)

        random.seed(1234)

        # Built alea table if does not exist yet
        if self._alea_canopy.empty:
            elements_vid_list = []
            for mtg_plant_vid in self._shared_mtg.components_iter(self._shared_mtg.root):
                for mtg_axis_vid in self._shared_mtg.components_iter(mtg_plant_vid):
                    for mtg_metamer_vid in self._shared_mtg.components_iter(mtg_axis_vid):
                        for mtg_organ_vid in self._shared_mtg.components_iter(mtg_metamer_vid):
                            for mtg_element_vid in self._shared_mtg.components_iter(mtg_organ_vid):
                                if self._shared_mtg.label(mtg_element_vid) == 'LeafElement1':
                                    elements_vid_list.append(mtg_element_vid)
            elements_vid_df = pd.DataFrame({'vid': elements_vid_list, 'tmp': 1})
            positions_df = pd.DataFrame({'pos': range(len(positions)),
                                         'tmp': 1,
                                         'azimut_leaf': 0,
                                         'inclination_leaf': 0})
            alea = pd.merge(elements_vid_df, positions_df, on=['tmp'])
            alea = alea.drop('tmp', axis=1)
            for vid in elements_vid_list:
                np.random.seed(vid)
                alea.loc[alea['vid'] == vid, 'azimut_leaf'] = np.random.uniform(-var_leaf_azimut, var_leaf_azimut, size=len(positions))
                alea.loc[alea['vid'] == vid, 'inclination_leaf'] = np.random.uniform(-var_leaf_inclination, var_leaf_inclination, size=len(positions))
            self._alea_canopy = alea

        # Duplication and heterogeneity
        duplicated_scene = plantgl.Scene()
        position_number = 0
        for pos in positions:
            azimut_stem = random.uniform(-var_stem_azimut, var_stem_azimut)
            for shp in initial_scene:
                if self._shared_mtg.label(shp.id) == 'StemElement':
                    rotated_geometry = plantgl.EulerRotated(azimut_stem, 0, 0, shp.geometry)
                    translated_geometry = plantgl.Translated(plantgl.Vector3(pos), rotated_geometry)
                    new_shape = plantgl.Shape(translated_geometry, appearance=shp.appearance, id=shp.id)
                    duplicated_scene += new_shape
                elif self._shared_mtg.label(shp.id) == 'LeafElement1':
                    # Add shp.id in alea_canopy if not in yet:
                    if shp.id not in list(self._alea_canopy['vid']):
                        new_vid_df = pd.DataFrame({'vid': shp.id, 'pos': range(len(positions))})
                        np.random.seed(shp.id)
                        new_vid_df['azimut_leaf'] = np.random.uniform(-var_leaf_azimut, var_leaf_azimut, size=len(positions))
                        new_vid_df['inclination_leaf'] = np.random.uniform(-var_leaf_inclination, var_leaf_inclination, size=len(positions))
                        self._alea_canopy = pd.concat([self._alea_canopy.copy(), new_vid_df], ignore_index=True)
                    # Translation to origin
                    anchor_point = self._shared_mtg.get_vertex_property(shp.id)['anchor_point']
                    trans_to_origin = plantgl.Translated(-anchor_point, shp.geometry)
                    # Rotation variability
                    azimut = self._alea_canopy.loc[(self._alea_canopy.pos == position_number) & (self._alea_canopy.vid == shp.id), 'azimut_leaf'].values[0]
                    inclination = self._alea_canopy.loc[(self._alea_canopy.pos == position_number) & (self._alea_canopy.vid == shp.id), 'inclination_leaf'].values[0]
                    rotated_geometry = plantgl.EulerRotated(azimut, inclination, 0, trans_to_origin)
                    # Restore leaf base at initial anchor point
                    translated_geometry = plantgl.Translated(anchor_point, rotated_geometry)
                    # Translate leaf to new plant position
                    translated_geometry = plantgl.Translated(pos, translated_geometry)
                    new_shape = plantgl.Shape(translated_geometry, appearance=shp.appearance, id=shp.id)
                    duplicated_scene += new_shape
            position_number += 1

        return duplicated_scene, domain

    def update_shared_MTG(self, aggregated_outputs):
        """
        Update the MTG shared between all models from the population of Caribu.
        :param dict aggregated_outputs: {'param1': { vid1: , vid2, ...}, 'param2': { vid1: , vid2, ...}}
        """
        # add the missing property
        for param in aggregated_outputs.keys():
            if param not in self._shared_mtg.properties():
                self._shared_mtg.add_property(param)

            # update the MTG
            self._shared_mtg.property(param).update(aggregated_outputs[param])

    def update_shared_dataframes(self, aggregated_outputs):
        """
        Update the dataframes shared between all models from the inputs dataframes or the outputs dataframes of the model.
        :param dict aggregated_outputs: {'param1': { vid1: , vid2, ...}, 'param2': { vid1: , vid2, ...}}
        """
        ids = []
        ids_lidt_built = False
        aggregated_outputs_list = {}

        for param in aggregated_outputs.keys():
            aggregated_outputs_list[param] = []
            for vid in sorted(aggregated_outputs[param].keys()):
                if not ids_lidt_built:
                    ind = int(self._shared_mtg.index(self._shared_mtg.complex_at_scale(vid, 1))), self._shared_mtg.label(self._shared_mtg.complex_at_scale(vid, 2)), \
                          int(self._shared_mtg.index(self._shared_mtg.complex_at_scale(vid, 3))), self._shared_mtg.label(self._shared_mtg.complex_at_scale(vid, 4)), self._shared_mtg.label(vid)
                    ids.append(ind)
                aggregated_outputs_list[param].append(aggregated_outputs[param][vid])
            ids_lidt_built = True

        ids_df = pd.DataFrame(ids, columns=SHARED_ELEMENTS_INPUTS_OUTPUTS_INDEXES)
        data_df = pd.DataFrame(aggregated_outputs_list)
        df = pd.concat([ids_df, data_df], axis=1)
        df.sort_values(['plant', 'axis', 'metamer', 'organ', 'element'], inplace=True)
        tools.combine_dataframes_inplace(df, SHARED_ELEMENTS_INPUTS_OUTPUTS_INDEXES, self._shared_elements_inputs_outputs_df)
