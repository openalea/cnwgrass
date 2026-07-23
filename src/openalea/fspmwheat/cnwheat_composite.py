# -*- coding: latin-1 -*-

import os
import numpy as np
import pandas as pd
import time
import logging

from dataclasses import dataclass, fields
from openalea.metafspm.component import Model, declare


from openalea.adel.adel_dynamic import AdelDyn
from openalea.adel.echap_leaf import echap_leaves
from openalea.fspmwheat import caribu_facade
from openalea.fspmwheat import cnwheat_facade
from openalea.fspmwheat import elongwheat_facade
from openalea.fspmwheat import farquharwheat_facade
from openalea.fspmwheat import fspmwheat_facade
from openalea.fspmwheat import growthwheat_facade
from openalea.fspmwheat import senescwheat_facade


logger_output = logging.getLogger("Simulation_Logger")

debug = False

@dataclass
class WheatFSPM(Model):
    """
    TODO : Add description
    """

    # Inputs expected from bellowground models
    Export_Nitrates: float = declare(default=0., unit="umol.h-1", unit_comment="of nitrate",
                                        min_value="", max_value="", description="", value_comment="", references="", DOI="",
                                        variable_type="input", by="root_nitrogen", state_variable_type="extensive", edit_by="user")
    Export_Amino_Acids: float = declare(default=0., unit="umol.h-1", unit_comment="of N",
                                        min_value="", max_value="", description="", value_comment="", references="", DOI="",
                                        variable_type="input", by="root_nitrogen", state_variable_type="extensive", edit_by="user")
    sucrose_phloem_outside_solve: float = declare(default=10., unit="µmol of C", unit_comment="amount in equivalent C",
                                        min_value="", max_value="", description="", value_comment="", references="", DOI="",
                                        variable_type="input", by="root_carbon", state_variable_type="extensive", edit_by="user")
    amino_acids_phloem_outside_solve: float = declare(default=1., unit="µmol of N", unit_comment="amount in equivalent N",
                                        min_value="", max_value="", description="", value_comment="", references="", DOI="",
                                        variable_type="input", by="root_nitrogen", state_variable_type="extensive", edit_by="user")                         
    cytokinins: float = declare(default=0., unit="AU", unit_comment="",
                                        min_value="", max_value="", description="", value_comment="", references="", DOI="",
                                        variable_type="input", by="root_nitrogen", state_variable_type="extensive", edit_by="user")
    mstruct: float = declare(default=0., unit="g", unit_comment="", 
                                        min_value="", max_value="", description="", value_comment="", references="", DOI="",
                                        variable_type="input", by="root_carbon", state_variable_type="extensive", edit_by="user")

    # State variables condidered as outputs to bellowground models
    Total_Transpiration: float = declare(default=0., unit="mmol.s-1", unit_comment="of water", 
                                        min_value="", max_value="", description="", value_comment="", references="", DOI="",
                                        variable_type="state_variable", by="model_shoot", state_variable_type="extensive", edit_by="user")
    mstruct_axis: float = declare(default=0.05, unit="g", unit_comment="of axis structural mass", 
                                        min_value="", max_value="", description="", value_comment="", references="", DOI="",
                                        variable_type="state_variable", by="model_shoot", state_variable_type="extensive", edit_by="user")
    sucrose_phloem: float = declare(default=10., unit="µmol", unit_comment="of C", 
                                        min_value="", max_value="", description="", value_comment="", references="", DOI="",
                                        variable_type="state_variable", by="model_shoot", state_variable_type="extensive", edit_by="user")
    amino_acids_phloem: float = declare(default=1, unit="µmol", unit_comment="of N", 
                                        min_value="", max_value="", description="", value_comment="", references="", DOI="",
                                        variable_type="state_variable", by="model_shoot", state_variable_type="extensive", edit_by="user")
    Export_cytokinins: float = declare(default=0., unit="AU.h-1", unit_comment="of cytokinins",
                                        min_value="", max_value="", description="", value_comment="", references="", DOI="",
                                        variable_type="state_variable", by="model_shoot", state_variable_type="extensive", edit_by="user")
    adventitious_to_emerge: list = declare(default=None, unit="", unit_comment="", description="List of adventitous roots delays before emergence starting from current time step, length of list indicates the number to emerge", 
                                                    min_value="", max_value="", value_comment="", references="", DOI="", 
                                                    variable_type="state_variable", by="model_shoot", state_variable_type="descriptor", edit_by="user")
    
    # PARAMETERS
    synchronize_adventitious_emergence: bool = declare(default=True, unit="", unit_comment="", description="boolean to choose option where root nodal emergence depend on shoot leaf emergence dynamic", 
                                                    min_value="", max_value="", value_comment="", references="", DOI="", 
                                                    variable_type="parameter", by="model_shoot", state_variable_type="descriptor", edit_by="user")
    nodal_emergence_delay_since_leaf_emerged: float = declare(default=(180 / 20) * 24 * 3600, unit="s", unit_comment="equivalent at 20°C", description="Emergence delay for nodal primordium to emerge since leaf emerged on this node", 
                                                    min_value="", max_value="", value_comment="", references="Klepper 1984", DOI="", 
                                                    variable_type="parameter", by="model_shoot", state_variable_type="descriptor", edit_by="user")


    def __init__(self, root_mtg, meteo, inputs_dataframes,
                 HIDDENZONES_INITIAL_STATE_FILENAME = 'hiddenzones_initial_state.csv', ELEMENTS_INITIAL_STATE_FILENAME = 'elements_initial_state.csv', 
                 AXES_INITIAL_STATE_FILENAME = 'axes_initial_state.csv', update_parameters_all_models = None, HOUR_TO_SECOND_CONVERSION_FACTOR = 3600, 
                 ORGANS_INITIAL_STATE_FILENAME = 'organs_initial_state.csv', SOILS_INITIAL_STATE_FILENAME = 'soils_initial_state.csv', INPUTS_DIRPATH='inputs', 
                 PLANT_DENSITY=None, tillers_replications=None, stored_times = None, N_fertilizations=None, computing_light_interception=True, heterogeneous_canopy=True, show_3Dplant = False, option_static = False, 
                 isolated_roots = False, cnwheat_roots=True, UPDATE_SHARED_DF=False, START_TIME = 0,
                 CARIBU_TIMESTEP = 4, FARQUHARWHEAT_TIMESTEP = 1, ELONGWHEAT_TIMESTEP = 1, GROWTHWHEAT_TIMESTEP = 1, CNWHEAT_TIMESTEP = 1, SENESCWHEAT_TIMESTEP = 1):
        
        # SELF STORAGE FOR LOOP PARAMETERS
        self.meteo = meteo

        # time steps
        self.time_step_in_hours = START_TIME
        self.CARIBU_TIMESTEP = CARIBU_TIMESTEP
        self.SENESCWHEAT_TIMESTEP = SENESCWHEAT_TIMESTEP
        self.FARQUHARWHEAT_TIMESTEP = FARQUHARWHEAT_TIMESTEP
        self.ELONGWHEAT_TIMESTEP = ELONGWHEAT_TIMESTEP
        self.GROWTHWHEAT_TIMESTEP = GROWTHWHEAT_TIMESTEP
        self.CNWHEAT_TIMESTEP = CNWHEAT_TIMESTEP

        # canopy parameters
        self.PLANT_DENSITY = PLANT_DENSITY
        self.N_fertilizations = N_fertilizations

        # plant parameters
        self.tillers_replications = tillers_replications

        # logging and data structures
        self.stored_times = stored_times
        self.shared_elements_inputs_outputs_df = pd.DataFrame()
        self.shared_elements_inputs_outputs_df = pd.DataFrame()
        self.shared_hiddenzones_inputs_outputs_df = pd.DataFrame()
        self.shared_organs_inputs_outputs_df = pd.DataFrame()
        self.shared_soils_inputs_outputs_df = pd.DataFrame()
        self.shared_axes_inputs_outputs_df = pd.DataFrame()
        self.all_simulation_steps = []
        self.axes_all_data_list = []
        self.organs_all_data_list = []
        self.hiddenzones_all_data_list = []
        self.elements_all_data_list = []
        self.soils_all_data_list = []


        # boolean choices
        self.show_3Dplant = show_3Dplant
        self.heterogeneous_canopy = heterogeneous_canopy
        self.computing_light_interception = computing_light_interception
        self.option_static = option_static

        # -- ADEL and MTG CONFIGURATION --

        # read adelwheat inputs at t0
        self.adel_wheat = AdelDyn(seed=1, scene_unit='m', leaves=echap_leaves(xy_model='Soissons_byleafclass'))
        self.g = self.adel_wheat.load(directory=INPUTS_DIRPATH)
        
        # Section specific to coupling with Root-BRIDGES
        self.shoot_props = self.g.properties()

        self.props = root_mtg.properties()
        self.vertices = root_mtg.vertices(scale=root_mtg.max_scale())

        self.link_self_to_mtg()

        # ---------------------------------------------
        # ----- CONFIGURATION OF THE FACADES -------
        # ---------------------------------------------

        # -- ELONGWHEAT (created first because it is the only facade to add new metamers) --
        # Initial states
        elongwheat_hiddenzones_initial_state = inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME][
            elongwheat_facade.converter.HIDDENZONE_TOPOLOGY_COLUMNS + [i for i in elongwheat_facade.simulation.HIDDENZONE_INPUTS if i in
                                                                    inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME].columns]].copy()
        elongwheat_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
            elongwheat_facade.converter.ELEMENT_TOPOLOGY_COLUMNS + [i for i in elongwheat_facade.simulation.ELEMENT_INPUTS if i in
                                                                    inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns]].copy()
        elongwheat_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
            elongwheat_facade.converter.AXIS_TOPOLOGY_COLUMNS + [i for i in elongwheat_facade.simulation.AXIS_INPUTS if i in inputs_dataframes[AXES_INITIAL_STATE_FILENAME].columns]].copy()

        phytoT_df = pd.read_csv(os.path.join(INPUTS_DIRPATH, 'phytoT.csv'))

        # Update parameters if specified
        if update_parameters_all_models and 'elongwheat' in update_parameters_all_models:
            update_parameters_elongwheat = update_parameters_all_models['elongwheat']
        else:
            update_parameters_elongwheat = None

        # Facade initialisation
        self.elongwheat_facade_ = elongwheat_facade.ElongWheatFacade(self.g,
                                                                ELONGWHEAT_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                                elongwheat_axes_initial_state,
                                                                elongwheat_hiddenzones_initial_state,
                                                                elongwheat_elements_initial_state,
                                                                self.shared_axes_inputs_outputs_df,
                                                                self.shared_hiddenzones_inputs_outputs_df,
                                                                self.shared_elements_inputs_outputs_df,
                                                                self.adel_wheat, phytoT_df,
                                                                update_parameters_elongwheat,
                                                                update_shared_df=UPDATE_SHARED_DF)

        # -- CARIBU --
        if self.computing_light_interception:
            self.caribu_facade_ = caribu_facade.CaribuFacade(self.g,
                                                        self.shared_elements_inputs_outputs_df,
                                                        self.adel_wheat,
                                                        update_shared_df=UPDATE_SHARED_DF)

        # -- SENESCWHEAT --
        # Initial states    
        senescwheat_roots_initial_state = inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].loc[inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME]['organ'] == 'roots'][
            senescwheat_facade.converter.ROOTS_TOPOLOGY_COLUMNS +
            [i for i in senescwheat_facade.converter.SENESCWHEAT_ROOTS_INPUTS if i in inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].columns]].copy()

        senescwheat_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
            senescwheat_facade.converter.ELEMENTS_TOPOLOGY_COLUMNS +
            [i for i in senescwheat_facade.converter.SENESCWHEAT_ELEMENTS_INPUTS if i in inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns]].copy()

        senescwheat_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
            senescwheat_facade.converter.AXES_TOPOLOGY_COLUMNS +
            [i for i in senescwheat_facade.converter.SENESCWHEAT_AXES_INPUTS if i in inputs_dataframes[AXES_INITIAL_STATE_FILENAME].columns]].copy()

        # Update parameters if specified
        if update_parameters_all_models and 'senescwheat' in update_parameters_all_models:
            update_parameters_senescwheat = update_parameters_all_models['senescwheat']
        else:
            update_parameters_senescwheat = None

        # Facade initialisation
        self.senescwheat_facade_ = senescwheat_facade.SenescWheatFacade(self.g,
                                                                SENESCWHEAT_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                                senescwheat_roots_initial_state,
                                                                senescwheat_axes_initial_state,
                                                                senescwheat_elements_initial_state,
                                                                self.shared_organs_inputs_outputs_df,
                                                                self.shared_axes_inputs_outputs_df,
                                                                self.shared_elements_inputs_outputs_df,
                                                                update_parameters_senescwheat,
                                                                update_shared_df=UPDATE_SHARED_DF,
                                                                cnwheat_roots=cnwheat_roots)

        # -- FARQUHARWHEAT --
        # Initial states    
        farquharwheat_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
            farquharwheat_facade.converter.ELEMENT_TOPOLOGY_COLUMNS +
            [i for i in farquharwheat_facade.converter.FARQUHARWHEAT_ELEMENTS_INPUTS if i in inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns]].copy()

        farquharwheat_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
            farquharwheat_facade.converter.AXIS_TOPOLOGY_COLUMNS +
            [i for i in farquharwheat_facade.converter.FARQUHARWHEAT_AXES_INPUTS if i in inputs_dataframes[AXES_INITIAL_STATE_FILENAME].columns]].copy()

        # Use the initial version of the photosynthesis sub-model (as in Barillot et al. 2016, and in Gauthier et al. 2020)
        update_parameters_farquharwheat = {'SurfacicProteins': False, 'NSC_Retroinhibition': False}

        # Facade initialisation
        self.farquharwheat_facade_ = farquharwheat_facade.FarquharWheatFacade(self.g,
                                                                        farquharwheat_elements_initial_state,
                                                                        farquharwheat_axes_initial_state,
                                                                        self.shared_elements_inputs_outputs_df,
                                                                        update_parameters_farquharwheat,
                                                                        update_shared_df=UPDATE_SHARED_DF)

        # -- GROWTHWHEAT --
        # Initial states    
        growthwheat_hiddenzones_initial_state = inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME][
            growthwheat_facade.converter.HIDDENZONE_TOPOLOGY_COLUMNS +
            [i for i in growthwheat_facade.simulation.HIDDENZONE_INPUTS if i in inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME].columns]].copy()

        growthwheat_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
            growthwheat_facade.converter.ELEMENT_TOPOLOGY_COLUMNS +
            [i for i in growthwheat_facade.simulation.ELEMENT_INPUTS if i in inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns]].copy()

        growthwheat_root_initial_state = inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].loc[inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME]['organ'] == 'roots'][
            growthwheat_facade.converter.ROOT_TOPOLOGY_COLUMNS +
            [i for i in growthwheat_facade.simulation.ROOT_INPUTS if i in inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].columns]].copy()

        growthwheat_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
            growthwheat_facade.converter.AXIS_TOPOLOGY_COLUMNS +
            [i for i in growthwheat_facade.simulation.AXIS_INPUTS if i in inputs_dataframes[AXES_INITIAL_STATE_FILENAME].columns]].copy()

        # Update parameters if specified
        if update_parameters_all_models and 'growthwheat' in update_parameters_all_models:
            update_parameters_growthwheat = update_parameters_all_models['growthwheat']
        else:
            update_parameters_growthwheat = None

        # Facade initialisation
        self.growthwheat_facade_ = growthwheat_facade.GrowthWheatFacade(self.g,
                                                                GROWTHWHEAT_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                                growthwheat_hiddenzones_initial_state,
                                                                growthwheat_elements_initial_state,
                                                                growthwheat_root_initial_state,
                                                                growthwheat_axes_initial_state,
                                                                self.shared_organs_inputs_outputs_df,
                                                                self.shared_hiddenzones_inputs_outputs_df,
                                                                self.shared_elements_inputs_outputs_df,
                                                                self.shared_axes_inputs_outputs_df,
                                                                update_parameters_growthwheat,
                                                                update_shared_df=UPDATE_SHARED_DF,
                                                                cnwheat_roots=cnwheat_roots)

        # -- CNWHEAT --
        # Initial states    
        cnwheat_organs_initial_state = inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME][
            [i for i in cnwheat_facade.cnwheat_converter.ORGANS_VARIABLES if i in inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].columns]].copy()

        cnwheat_hiddenzones_initial_state = inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME][
            [i for i in cnwheat_facade.cnwheat_converter.HIDDENZONE_VARIABLES if i in inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME].columns]].copy()

        cnwheat_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
            [i for i in cnwheat_facade.cnwheat_converter.ELEMENTS_VARIABLES if i in inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns]].copy()

        cnwheat_soils_initial_state = inputs_dataframes[SOILS_INITIAL_STATE_FILENAME][
            [i for i in cnwheat_facade.cnwheat_converter.SOILS_VARIABLES if i in inputs_dataframes[SOILS_INITIAL_STATE_FILENAME].columns]].copy()

        # Update parameters if specified
        if update_parameters_all_models and 'cnwheat' in update_parameters_all_models:
            update_parameters_cnwheat = update_parameters_all_models['cnwheat']
        else:
            update_parameters_cnwheat = {}

        # Force solver separation if a different root model has been chosen
        if not cnwheat_roots:
            isolated_roots = True

        # Facade initialisation
        self.cnwheat_facade_ = cnwheat_facade.CNWheatFacade(self.g,
                                                    CNWHEAT_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                    PLANT_DENSITY,
                                                    update_parameters_cnwheat,
                                                    cnwheat_organs_initial_state,
                                                    cnwheat_hiddenzones_initial_state,
                                                    cnwheat_elements_initial_state,
                                                    cnwheat_soils_initial_state,
                                                    self.shared_axes_inputs_outputs_df,
                                                    self.shared_organs_inputs_outputs_df,
                                                    self.shared_hiddenzones_inputs_outputs_df,
                                                    self.shared_elements_inputs_outputs_df,
                                                    self.shared_soils_inputs_outputs_df,
                                                    update_shared_df=UPDATE_SHARED_DF,
                                                    isolated_roots=isolated_roots,
                                                    cnwheat_roots=cnwheat_roots)

        # Run cnwheat with constant nitrates concentration in the soil if specified
        if N_fertilizations is not None and 'constant_Conc_Nitrates' in N_fertilizations.keys():
            self.cnwheat_facade_.soils[(1, 'MS')].constant_Conc_Nitrates = True  # TODO: make (1, 'MS') more general
            self.cnwheat_facade_.soils[(1, 'MS')].nitrates = N_fertilizations['constant_Conc_Nitrates'] * self.cnwheat_facade_.soils[(1, 'MS')].volume

        # -- FSPMWHEAT --
        # Facade initialisation
        self.fspmwheat_facade_ = fspmwheat_facade.FSPMWheatFacade(self.g)

        # Update geometry
        self.adel_wheat.update_geometry(self.g)
        if show_3Dplant:
            self.adel_wheat.plot(self.g)
        
        self.cn_wheat_root_props = self.g.get_vertex_property(2)["roots"]

        # TODO : Temporary
        self.g.properties()["Total_Transpiration"][2] = self.props["Total_Transpiration"][1]

        self.g.get_vertex_property(2)['phloem']['sucrose'] = self.props["sucrose_phloem"][1]
        self.g.get_vertex_property(2)['phloem']['amino_acids'] = self.props["amino_acids_phloem"][1]

        
        if self.synchronize_adventitious_emergence:
            # Specific initialization for already emerged leaves at the begining of the simulation
            self.main_stem_vid = [v for v in self.g.vertices(scale=2) if "MS" in str(self.g.node(v).label)][0]
            self.already_emerged_leaves = []
            hiddenzones = self.g.vertices(scale=3)
            for v in hiddenzones:
                n = self.g.node(v)
                if hasattr(n, "hiddenzone"):
                    if n.hiddenzone is not None:
                        hz = n.hiddenzone
                        if hz["leaf_is_emerged"]:
                            self.already_emerged_leaves.append(v)
            nodal_emergence_delays = []
            for vid in self.g.components_at_scale(self.main_stem_vid, scale=5):
                n = self.g.node(vid)
                if "Leaf" in n.label:
                    remaining_time_to_emergence = max(0, self.nodal_emergence_delay_since_leaf_emerged - (n.age / 20) * 24 * 3600)
                    nodal_emergence_delays += [remaining_time_to_emergence,
                                               remaining_time_to_emergence]
                    # Third likely nodal on this node
                    if np.random.random() < 0.8:
                        nodal_emergence_delays.append(remaining_time_to_emergence + (100 / 20) * 24 * 3600 )
                    # Fourth less likely nodal on this node
                    elif np.random.random() < 0.2:
                        nodal_emergence_delays.append(remaining_time_to_emergence + (100 / 20) * 24 * 3600 )
            
            self.props["adventitious_to_emerge"].update({1: nodal_emergence_delays})

        self.sync_shoot_outputs_with_root_mtg()
        

    def sync_shoot_inputs_with_shoot_mtg(self):
        for name in self.inputs:
            if name == 'sucrose_phloem_outside_solve':
                self.g.get_vertex_property(2)['phloem']['sucrose'] == self.props[name][1]
            elif name == 'amino_acids_phloem_outside_solve':
                self.g.get_vertex_property(2)['phloem']['amino_acids'] == self.props[name][1]
            else:
                self.cn_wheat_root_props[name] = self.props[name][1]

    def sync_shoot_outputs_with_root_mtg(self):
        # Link this specific data structure to self for variables exchange, only for outputs that will be read by other models  here.
        # Note : here eval is necessary to ensure intended lambda function definition
        for name in self.state_variables:
            if name == "Total_Transpiration":
                self.props[name].update({1: self.g.get_vertex_property(2)[name]})

            elif name == "mstruct_axis":
                self.props[name].update({1: self.g.get_vertex_property(2)['mstruct']})

            elif name == "amino_acids_phloem":
                self.props[name].update({1: self.g.get_vertex_property(2)['phloem']['amino_acids']})

            elif name == "sucrose_phloem":
                self.props[name].update({1: self.g.get_vertex_property(2)['phloem']['sucrose']})

            elif name == "adventitious_to_emerge" and self.synchronize_adventitious_emergence:
                hiddenzones = self.g.vertices(scale=3)
                main_stem_hz = self.g.components_at_scale(self.main_stem_vid, scale=3)
                nodal_emergence_delays = [] # Several per time step are possible!
                for v in hiddenzones:
                    n = self.g.node(v)
                    
                    if hasattr(n, "hiddenzone"):
                        if n.hiddenzone is not None:
                            hz = n.hiddenzone
                            if hz["leaf_is_emerged"] and v not in self.already_emerged_leaves:
                                nodal_emergence_delays += [self.nodal_emergence_delay_since_leaf_emerged,
                                                        self.nodal_emergence_delay_since_leaf_emerged]
                                if v in main_stem_hz:
                                    # Third likely nodal on this node
                                    if np.random.random() < 0.8:
                                        nodal_emergence_delays.append(self.nodal_emergence_delay_since_leaf_emerged + (100 / 20) * 24 * 3600 )
                                    # Fourth less likely nodal on this node
                                    elif np.random.random() < 0.2:
                                        nodal_emergence_delays.append(self.nodal_emergence_delay_since_leaf_emerged + (100 / 20) * 24 * 3600 )

                                self.props[name].update({1: nodal_emergence_delays})
                                self.already_emerged_leaves.append(v)
                            
            else:
                self.props[name].update({1: self.cn_wheat_root_props[name]})

        


    def __call__(self):
        # SPECIFIC TO COUPLING, syncs the shoot mtg variables with bellowground models
        self.pull_available_inputs()
        self.sync_shoot_inputs_with_shoot_mtg()

        # run Caribu
        PARi = self.meteo.loc[self.time_step_in_hours, ['PARi']].iloc[0]
        DOY = self.meteo.loc[self.time_step_in_hours, ['DOY']].iloc[0]
        hour = self.meteo.loc[self.time_step_in_hours, ['hour']].iloc[0]

        if self.computing_light_interception:
            PARi_next_hours = self.meteo.loc[range(self.time_step_in_hours, self.time_step_in_hours + self.CARIBU_TIMESTEP), ['PARi']].sum().values[0]

            if (self.time_step_in_hours % self.CARIBU_TIMESTEP == 0) and (PARi_next_hours > 0):
                run_caribu = True
            else:
                run_caribu = False

            self.caribu_facade_.run(run_caribu, energy=PARi, DOY=DOY, hourTU=hour, latitude=48.85, sun_sky_option='sky', heterogeneous_canopy=self.heterogeneous_canopy, plant_density=self.PLANT_DENSITY[1])

        for t_senescwheat in range(self.time_step_in_hours, self.time_step_in_hours + self.SENESCWHEAT_TIMESTEP, self.SENESCWHEAT_TIMESTEP):
            # run SenescWheat
            t1 = time.time()
            self.senescwheat_facade_.run()
            # if (1, 'MS', 5, 'internode', 'HiddenElement') not in self.senescwheat_facade_._simulation.outputs["elements"].keys():
            #     logger_output.info("Post SW, Internode not found")
            hidelt = self.g.node(80)
            if (hidelt.length is None) or (np.isnan(hidelt.length)) or (np.isinf(hidelt.length)) or (hidelt.length <= 0.):
                logger_output.info(f"Post SW, Internode element has no proper length in MTG : {hidelt.length}")
            hidorg = self.g.node(70)
            if (hidorg.length is None) or (np.isnan(hidorg.length)) or (np.isinf(hidorg.length)) or (hidorg.length <= 0.):
                logger_output.info(f"Post SW, Internode organ has no proper length in MTG : {hidorg.length}")
            if debug: print("senescwheat took :", time.time() - t1)

            # Test for dead plant # TODO: adapt in case of multiple plants
            if not self.shared_elements_inputs_outputs_df.empty and \
                    np.nansum(self.shared_elements_inputs_outputs_df.loc[self.shared_elements_inputs_outputs_df['element'].isin(['StemElement', 'LeafElement1']), 'green_area']) == 0:
                # append the inputs and outputs at current step to global lists
                self.all_simulation_steps.append(t_senescwheat)
                self.axes_all_data_list.append(self.shared_axes_inputs_outputs_df.copy())
                self.organs_all_data_list.append(self.shared_organs_inputs_outputs_df.copy())
                self.hiddenzones_all_data_list.append(self.shared_hiddenzones_inputs_outputs_df.copy())
                self.elements_all_data_list.append(self.shared_elements_inputs_outputs_df.copy())
                self.soils_all_data_list.append(self.shared_soils_inputs_outputs_df.copy())
                break

            # Run the rest of the model if the plant is alive
            for t_farquharwheat in range(t_senescwheat, t_senescwheat + self.SENESCWHEAT_TIMESTEP, self.FARQUHARWHEAT_TIMESTEP):
                # get the meteo of the current step
                Ta, ambient_CO2, RH, Ur = self.meteo.loc[t_farquharwheat, ['air_temperature', 'ambient_CO2', 'humidity', 'Wind']]

                # run FarquharWheat
                t1 = time.time()
                self.farquharwheat_facade_.run(Ta, ambient_CO2, RH, Ur)
                # if (1, 'MS', 5, 'internode', 'HiddenElement') not in self.farquharwheat_facade_._simulation.outputs.keys():
                #     logger_output.info("Post FW, Internode not found")
                hidelt = self.g.node(80)
                if (hidelt.length is None) or (np.isnan(hidelt.length)) or (np.isinf(hidelt.length)) or (hidelt.length <= 0.):
                    logger_output.info(f"Post FW, Internode element has no proper length in MTG : {hidelt.length}")
                hidorg = self.g.node(70)
                if (hidorg.length is None) or (np.isnan(hidorg.length)) or (np.isinf(hidorg.length)) or (hidorg.length <= 0.):
                    logger_output.info(f"Post FW, Internode organ has no proper length in MTG : {hidorg.length}")
                if debug: print("farquharwheat took :", time.time() - t1)

                for t_elongwheat in range(t_farquharwheat, t_farquharwheat + self.FARQUHARWHEAT_TIMESTEP, self.ELONGWHEAT_TIMESTEP):
                    # run ElongWheat
                    Tair, Tsoil = self.meteo.loc[t_elongwheat, ['air_temperature', 'soil_temperature']]
                    t1 = time.time()
                    self.elongwheat_facade_.run(Tair, Tsoil, option_static=self.option_static)
                    # if (1, 'MS', 5, 'internode', 'HiddenElement') not in self.elongwheat_facade_._simulation.outputs["elements"].keys():
                    #     logger_output.info("Post EW, Internode not found")
                    hidelt = self.g.node(80)
                    if (hidelt.length is None) or (np.isnan(hidelt.length)) or (np.isinf(hidelt.length)) or (hidelt.length <= 0.):
                        logger_output.info(f"Post EW, Internode element has no proper length in MTG : {hidelt.length}")
                    hidorg = self.g.node(70)
                    if (hidorg.length is None) or (np.isnan(hidorg.length)) or (np.isinf(hidorg.length)) or (hidorg.length <= 0.):
                        logger_output.info(f"Post EW, Internode organ has no proper length in MTG : {hidorg.length}")
                    if debug: print("elongwheat took :", time.time() - t1)

                    # Update geometry
                    self.adel_wheat.update_geometry(self.g)
                    if self.show_3Dplant:
                        self.adel_wheat.plot(self.g)
                    hidelt = self.g.node(80)
                    if (hidelt.length is None) or (np.isnan(hidelt.length)) or (np.isinf(hidelt.length)) or (hidelt.length <= 0.):
                        logger_output.info(f"Post AW, Internode element has no proper length in MTG : {hidelt.length}")
                    hidorg = self.g.node(70)
                    if (hidorg.length is None) or (np.isnan(hidorg.length)) or (np.isinf(hidorg.length)) or (hidorg.length <= 0.):
                        logger_output.info(f"Post AW, Internode organ has no proper length in MTG : {hidorg.length}")

                    for t_growthwheat in range(t_elongwheat, t_elongwheat + self.ELONGWHEAT_TIMESTEP, self.GROWTHWHEAT_TIMESTEP):
                        # run GrowthWheat
                        t1 = time.time()
                        self.growthwheat_facade_.run()
                        # if (1, 'MS', 5, 'internode', 'HiddenElement') not in self.growthwheat_facade_._simulation.outputs["elements"].keys():
                        #     logger_output.info("Post GW, Internode not found")
                        hidelt = self.g.node(80)
                        if (hidelt.length is None) or (np.isnan(hidelt.length)) or (np.isinf(hidelt.length)) or (hidelt.length <= 0.):
                            logger_output.info(f"Post GW, Internode element has no proper length in MTG : {hidelt.length}")
                        hidorg = self.g.node(70)
                        if (hidorg.length is None) or (np.isnan(hidorg.length)) or (np.isinf(hidorg.length)) or (hidorg.length <= 0.):
                            logger_output.info(f"Post GW, Internode organ has no proper length in MTG : {hidorg.length}")

                        for t_cnwheat in range(t_growthwheat, t_growthwheat + self.GROWTHWHEAT_TIMESTEP, self.CNWHEAT_TIMESTEP):
                            #print('t cnwheat is {}'.format(t_cnwheat))

                            # N fertilization if any
                            if self.N_fertilizations is not None and len(self.N_fertilizations) > 0:
                                if t_cnwheat in self.N_fertilizations.keys():
                                    self.cnwheat_facade_.soils[(1, 'MS')].nitrates += self.N_fertilizations[t_cnwheat]

                            if t_cnwheat > 0:
                                # run CNWheat
                                Tair = self.meteo.loc[t_elongwheat, 'air_temperature']
                                Tsoil = self.meteo.loc[t_elongwheat, 'soil_temperature']
                                t1 = time.time()
                                self.cnwheat_facade_.run(Tair, Tsoil, self.tillers_replications)
                                if debug: print("cnwheat took :", time.time() - t1)
                                hidelt = self.g.node(80)
                                if (hidelt.length is None) or (np.isnan(hidelt.length)) or (np.isinf(hidelt.length)) or (hidelt.length <= 0.):
                                    logger_output.info(f"Post CNW, Internode element has no proper length in MTG : {hidelt.length}")
                                hidorg = self.g.node(70)
                                if (hidorg.length is None) or (np.isnan(hidorg.length)) or (np.isinf(hidorg.length)) or (hidorg.length <= 0.):
                                    logger_output.info(f"Post CNW, Internode organ has no proper length in MTG : {hidorg.length}")

                            # append outputs at current step to global lists
                            if (self.stored_times == 'all') or (t_cnwheat in self.stored_times):
                                axes_outputs, elements_outputs, hiddenzones_outputs, organs_outputs, soils_outputs = self.fspmwheat_facade_.build_outputs_df_from_MTG()

                                self.all_simulation_steps.append(t_cnwheat)
                                self.axes_all_data_list.append(axes_outputs)
                                self.organs_all_data_list.append(organs_outputs)
                                self.hiddenzones_all_data_list.append(hiddenzones_outputs)
                                self.elements_all_data_list.append(elements_outputs)
                                self.soils_all_data_list.append(soils_outputs)

        self.time_step_in_hours += self.SENESCWHEAT_TIMESTEP

        self.sync_shoot_outputs_with_root_mtg()



def scenario_utility(time_step_in_seconds: int = 3600, INPUTS_DIRPATH = "inputs", OUTPUTS_DIRPATH = "outputs", METEO_FILENAME = "meteo_Ljutovac2002.csv", PLANT_DENSITY = {1:250},
                     forced_start_time = 0, tillers_replications={'T1': 0.5, 'T2': 0.5, 'T3': 0.5, 'T4': 0.5}, N_fertilizations = {2016: 357143, 2520: 1000000},
                     stored_times = None, option_static = False, show_3Dplant = False, run_from_outputs = False, heterogeneous_canopy = True, update_parameters_all_models = None,
                     isolated_roots = False, cnwheat_roots = True):
    scenario = {}

    ### DIRS ###

    scenario["INPUTS_DIRPATH"] = INPUTS_DIRPATH
    
    # Save the outputs with a full scan of the MTG at each time step (or at selected time steps)
    UPDATE_SHARED_DF = False
    if stored_times is None:
        stored_times = 'all'
    if not (stored_times == 'all' or type(stored_times) == list):
        print('stored_times should be either \'all\', a list or an empty list.')
        raise

    scenario["stored_times"] = stored_times

    ### METEO PARAMETER ###
    scenario["meteo"] = pd.read_csv(os.path.join(INPUTS_DIRPATH, METEO_FILENAME), index_col='t')

    AXES_INDEX_COLUMNS = ['t', 'plant', 'axis']
    ELEMENTS_INDEX_COLUMNS = ['t', 'plant', 'axis', 'metamer', 'organ', 'element']
    HIDDENZONES_INDEX_COLUMNS = ['t', 'plant', 'axis', 'metamer']
    ORGANS_INDEX_COLUMNS = ['t', 'plant', 'axis', 'organ']
    SOILS_INDEX_COLUMNS = ['t', 'plant', 'axis']

    # Name of the CSV files which describes the initial state of the system
    AXES_INITIAL_STATE_FILENAME = 'axes_initial_state.csv'
    ORGANS_INITIAL_STATE_FILENAME = 'organs_initial_state.csv'
    HIDDENZONES_INITIAL_STATE_FILENAME = 'hiddenzones_initial_state.csv'
    ELEMENTS_INITIAL_STATE_FILENAME = 'elements_initial_state.csv'
    SOILS_INITIAL_STATE_FILENAME = 'soils_initial_state.csv'

    # Name of the CSV files which will contain the outputs of the model
    AXES_OUTPUTS_FILENAME = 'axes_outputs.csv'
    ORGANS_OUTPUTS_FILENAME = 'organs_outputs.csv'
    HIDDENZONES_OUTPUTS_FILENAME = 'hiddenzones_outputs.csv'
    ELEMENTS_OUTPUTS_FILENAME = 'elements_outputs.csv'
    SOILS_OUTPUTS_FILENAME = 'soils_outputs.csv'

    ### INPUT DATAFRAMES PARAMETER ###
    # Read the inputs from CSV files and create inputs dataframes
    inputs_dataframes = {}
    if run_from_outputs:

        previous_outputs_dataframes = {}

        for initial_state_filename, outputs_filename, index_columns in ((AXES_INITIAL_STATE_FILENAME, AXES_OUTPUTS_FILENAME, AXES_INDEX_COLUMNS),
                                                                        (ORGANS_INITIAL_STATE_FILENAME, ORGANS_OUTPUTS_FILENAME, ORGANS_INDEX_COLUMNS),
                                                                        (HIDDENZONES_INITIAL_STATE_FILENAME, HIDDENZONES_OUTPUTS_FILENAME, HIDDENZONES_INDEX_COLUMNS),
                                                                        (ELEMENTS_INITIAL_STATE_FILENAME, ELEMENTS_OUTPUTS_FILENAME, ELEMENTS_INDEX_COLUMNS),
                                                                        (SOILS_INITIAL_STATE_FILENAME, SOILS_OUTPUTS_FILENAME, SOILS_INDEX_COLUMNS)):

            previous_outputs_dataframe = pd.read_csv(os.path.join(OUTPUTS_DIRPATH, outputs_filename))
            # Convert NaN to None
            previous_outputs_dataframes[outputs_filename] = previous_outputs_dataframe.where(previous_outputs_dataframe.notnull(), None)

            assert 't' in previous_outputs_dataframes[outputs_filename].columns
            if forced_start_time > 0:
                new_start_time = forced_start_time + 1
                previous_outputs_dataframes[outputs_filename] = previous_outputs_dataframes[outputs_filename][previous_outputs_dataframes[outputs_filename]['t'] <= forced_start_time]
            else:
                last_t_step = max(previous_outputs_dataframes[outputs_filename]['t'])
                new_start_time = last_t_step + 1

            if initial_state_filename == ELEMENTS_INITIAL_STATE_FILENAME:
                elements_previous_outputs = previous_outputs_dataframes[outputs_filename]
                new_initial_state = elements_previous_outputs[~elements_previous_outputs.is_over.isnull()]
            else:
                new_initial_state = previous_outputs_dataframes[outputs_filename]
            idx = new_initial_state.groupby([col for col in index_columns if col != 't'])['t'].transform(max) == new_initial_state['t']
            inputs_dataframes[initial_state_filename] = new_initial_state[idx].drop(['t'], axis=1)

        # Make sure boolean columns have either type bool or float
        bool_columns = ['is_over', 'is_growing', 'leaf_is_emerged', 'internode_is_visible', 'leaf_is_growing', 'internode_is_growing', 'leaf_is_remobilizing', 'internode_is_remobilizing']
        for df in [inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME], inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME]]:
            for cln in bool_columns:
                if cln in df.keys():
                    df[cln].replace(to_replace='False', value=0.0, inplace=True)
                    df[cln].replace(to_replace='True', value=1.0, inplace=True)
                    df[cln] = pd.to_numeric(df[cln])
    else:
        new_start_time = -1
        for inputs_filename in (AXES_INITIAL_STATE_FILENAME,
                                ORGANS_INITIAL_STATE_FILENAME,
                                HIDDENZONES_INITIAL_STATE_FILENAME,
                                ELEMENTS_INITIAL_STATE_FILENAME,
                                SOILS_INITIAL_STATE_FILENAME):
            inputs_dataframe = pd.read_csv(os.path.join(INPUTS_DIRPATH, inputs_filename))
            inputs_dataframes[inputs_filename] = inputs_dataframe.where(inputs_dataframe.notnull(), None)

    # Start time of the simulation
    START_TIME = max(0, new_start_time)
    scenario["START_TIME"] = START_TIME

    # Pass a unified time-step
    
    # TODO constrained to interger for now, and using float, so bellow hour time stepping breaks both for loops and input dataframe accessions
    time_step_in_hours = int(time_step_in_seconds / 3600)
    scenario["CARIBU_TIMESTEP"] = time_step_in_hours
    scenario["FARQUHARWHEAT_TIMESTEP"] = time_step_in_hours
    scenario["ELONGWHEAT_TIMESTEP"] = time_step_in_hours
    scenario["GROWTHWHEAT_TIMESTEP"] = time_step_in_hours
    scenario["CNWHEAT_TIMESTEP"] = time_step_in_hours
    scenario["SENESCWHEAT_TIMESTEP"] = time_step_in_hours

    scenario["inputs_dataframes"] = inputs_dataframes

    ### OPTIONS ###
    scenario["PLANT_DENSITY"] = PLANT_DENSITY
    scenario["option_static"] = option_static
    scenario["show_3Dplant"] = show_3Dplant
    scenario["tillers_replications"] = tillers_replications
    scenario["heterogeneous_canopy"] = heterogeneous_canopy
    scenario["N_fertilizations"] = N_fertilizations
    scenario["update_parameters_all_models"] = update_parameters_all_models
    scenario["isolated_roots"] = isolated_roots
    scenario["cnwheat_roots"] = cnwheat_roots


    return scenario
