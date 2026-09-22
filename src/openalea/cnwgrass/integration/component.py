# -*- coding: latin-1 -*-

import os
import numpy as np
import pandas as pd

from dataclasses import dataclass, fields
from openalea.metafspm.component import Model, declare

from openalea.adel.adel_dynamic import AdelDyn
from openalea.adel.echap_leaf import echap_leaves
from openalea.cnwgrass.integration import caribu_facade
from openalea.cnwgrass.integration import cnmetabolism_facade
from openalea.cnwgrass.integration import morphogenesis_facade
from openalea.cnwgrass.integration import gasexchange_facade
from openalea.cnwgrass.integration import build_outputs
from openalea.cnwgrass.integration import growth_facade
from openalea.cnwgrass.integration import senescence_facade
from openalea.cnwgrass.integration import hydraulics_facade


@dataclass
class CNW_Grass(Model):
    """
    CNW-Grass is a Functional Structural Plant Model (FSPM) of grasses which fully integrates shoot morphogenesis and the metabolism of carbon (C) and nitrogen (N) at organ scale within a 3D representation of plant architecture. Plants are described as a collection of tillers, each consisting in individual shoot organs (lamina, sheath, internode, peduncle, chaff), a single root compartment, the grains, and a phloem. CNW-Grass also includes a hydraulic model allowing to compute water flow in the plant and the co-regulation of leaf growth by metabolic and hydraulic processes. In this case, the plants also include a xylem compartment.

    CNW-Grass simulates:

        Organ photosynthesis, temperature and transpiration from light distribution within the 3D canopy.
        Leaf and internode elongation.
        Leaf, internode and root growth in mass.
        N acquisition, synthesis and allocation of C and N metabolites at organ level and among tiller organs.
        Senescence of shoot organs and roots.
        Water fluxes and water potentials.

    Model inputs are the pedoclimatic conditions (temperature, light, humidity, CO2, wind, soil NO3-, Soil Relative Water Content) and initial dimensions, mass and metabolic composition of individual organs.
    """

    # Inputs expected from bellowground models
    Export_Nitrates: float = declare(default=0., unit="umol.h-1", unit_comment="of nitrate",
                                        min_value="", max_value="", description="", value_comment="", references="", DOI="",
                                        variable_type="input", by="root_nitrogen", state_variable_type="extensive", edit_by="user")
    Export_Amino_Acids: float = declare(default=0., unit="umol.h-1", unit_comment="of N",
                                        min_value="", max_value="", description="", value_comment="", references="", DOI="",
                                        variable_type="input", by="root_nitrogen", state_variable_type="extensive", edit_by="user")
    Unloading_Sucrose_phloem: float = declare(default=0.1, unit="umol.h-1", unit_comment="of C",
                                        min_value="", max_value="", description="", value_comment="", references="", DOI="",
                                        variable_type="input", by="root_nitrogen", state_variable_type="extensive", edit_by="user")
    Unloading_Amino_Acids_phloem: float = declare(default=0.1, unit="umol.h-1", unit_comment="of N",
                                        min_value="", max_value="", description="", value_comment="", references="", DOI="",
                                        variable_type="input", by="root_nitrogen", state_variable_type="extensive", edit_by="user")                      
    cytokinins: float = declare(default=0., unit="AU", unit_comment="",
                                        min_value="", max_value="", description="", value_comment="", references="", DOI="",
                                        variable_type="input", by="root_nitrogen", state_variable_type="extensive", edit_by="user")
    mstruct: float = declare(default=0., unit="g", unit_comment="", 
                                        min_value="", max_value="", description="", value_comment="", references="", DOI="",
                                        variable_type="input", by="root_carbon", state_variable_type="extensive", edit_by="user")
    xylem_water_potential: float = declare(default=-0.6, unit="MPa", unit_comment="input shoot xylem uniform total water potential (hydrostatic + osmotic)",
                                        min_value="", max_value="", description="", value_comment="", references="", DOI="",
                                        variable_type="input", by="root_water", state_variable_type="extensive", edit_by="user")

    # State variables condidered as outputs to bellowground models
    root_to_shoot_xylem_water_flow: float = declare(default=0., unit="g H2O . time_step-1", unit_comment="",
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
    Unloading_Sucrose: float = declare(default=30., unit="umol.g-1.h-1", unit_comment="of equivalent C mol? sucrose",
                                        min_value="", max_value="", description="", value_comment="", references="", DOI="",
                                        variable_type="state_variable", by="model_shoot", state_variable_type="extensive", edit_by="user")
    Unloading_Amino_Acids: float = declare(default=1., unit="umol.g-1.h-1", unit_comment="of equivalent amino acids N", 
                                        min_value="", max_value="", description="", value_comment="", references="", DOI="",
                                        variable_type="state_variable", by="model_shoot", state_variable_type="extensive", edit_by="user")
    Unloading_Sucrose_shoot_organs: float = declare(default=30., unit="umol.h-1", unit_comment="of equivalent C mol? sucrose",
                                        min_value="", max_value="", description="", value_comment="", references="", DOI="",
                                        variable_type="state_variable", by="model_shoot", state_variable_type="extensive", edit_by="user")
    Unloading_Amino_Acids_shoot_organs: float = declare(default=1., unit="umol.h-1", unit_comment="of equivalent amino acids N", 
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
                 AXES_INITIAL_STATE_FILENAME = 'axes_initial_state.csv', update_parameters_all_models = None, step_callback=None, HOUR_TO_SECOND_CONVERSION_FACTOR = 3600, 
                 ORGANS_INITIAL_STATE_FILENAME = 'organs_initial_state.csv', SOILS_INITIAL_STATE_FILENAME = 'soils_initial_state.csv', INPUTS_DIRPATH='inputs', 
                 # Canopy parameters
                 single_plant = False, plant_density = {1: 250}, inter_row = 0.15, sowing_depth = 0.025, N_fertilizations = None,
                 # Options
                 tillers_replications=None, stored_times = None, computing_light_interception=True, external_soil_model=False, heterogeneous_canopy=True, show_3Dplant = False, 
                 hydraulics = False, stomatal_model_name='BWB', drought_trigger=None, rehydration_scenario=None, optimal_growth_option=False, option_static = False, 
                 isolated_roots = False, cnwgrass_roots = True, UPDATE_SHARED_DF=False, START_TIME = 0,
                 CARIBU_TIMESTEP = 4, MORPHOGENESIS_TIMESTEP = 1, GROWTH_TIMESTEP = 1, CNMETABOLISM_TIMESTEP = 1, SENESCENCE_TIMESTEP = 1, HYDRAULICS_TIMESTEP = 1):
        """
        :param openalea.mtg.mtg.MTG root_mtg: the MTG of the belowground root model this shoot component is coupled to (e.g. Root-BRIDGES). Its vertex-1 properties are the channel through which the two models
                                            exchange state at each time step: this component reads them as its declared ``input`` variables (``Export_Nitrates``, ``Export_Amino_Acids``, ``Unloading_Sucrose_phloem``,
                                            ``Unloading_Amino_Acids_phloem``, ``cytokinins``, ``mstruct``, ``xylem_water_potential``) and writes back its declared ``state_variable`` outputs
                                            (``water_outflux``, ``mstruct_axis``, ``sucrose_phloem``, ``amino_acids_phloem``, ``Unloading_Sucrose``, ``Unloading_Amino_Acids``,
                                            ``Unloading_Sucrose_shoot_organs``, ``Unloading_Amino_Acids_shoot_organs``, ``Export_cytokinins``, ``adventitious_to_emerge``).
        :param pandas.DataFrame meteo: hourly meteorological forcing, indexed by time step ``t`` (hours), with at least the columns ``PARi``, ``DOY``, ``hour``, ``air_temperature``, ``ambient_CO2``, ``humidity``,
                                            ``Wind`` and ``soil_temperature``.
        :param dict [str, pandas.DataFrame] inputs_dataframes: the initial-state dataframes of the shoot model, keyed by the corresponding ``*_INITIAL_STATE_FILENAME`` argument
                                            (one entry each for hiddenzones, elements, axes, organs and soils). Typically built by :func:`scenario_utility`.
        :param str HIDDENZONES_INITIAL_STATE_FILENAME: key of `inputs_dataframes` holding the initial state of the hiddenzones (growing leaf/internode bases).
        :param str ELEMENTS_INITIAL_STATE_FILENAME: key of `inputs_dataframes` holding the initial state of the photosynthetic and non-photosynthetic organ elements.
        :param str AXES_INITIAL_STATE_FILENAME: key of `inputs_dataframes` holding the initial state at axis scale.
        :param dict or None update_parameters_all_models: a dict to override default model parameters, keyed by sub-model name
                                            {'cnmetabolism': {'organ1': {'param1': 'val1', 'param2': 'val2'}, 'organ2': {...}}, 'morphogenesis': {'param1': 'val1', ...}, ...}.
        :param dict or None step_callback: a dict of functions used to force some external inputs that are natively computed by the model, keyed by function name. Recognised keys are
                                            ``'ADEL_mtg'`` (called as ``step_callback['ADEL_mtg'](adel_wheat, INPUTS_DIRPATH, nff)`` to build the initial MTG instead of loading a serialised one)
                                            and ``'nitrate_uptake'`` (called at each time step as ``step_callback['nitrate_uptake'](t, population, g)`` when `external_soil_model` is True).
        :param int HOUR_TO_SECOND_CONVERSION_FACTOR: number of seconds in one hour, used to convert the `*_TIMESTEP` arguments (in hours) into a `delta_t` in seconds for each facade.
        :param str ORGANS_INITIAL_STATE_FILENAME: key of `inputs_dataframes` holding the initial state of the non-photosynthetic organs (roots, grains, ...).
        :param str SOILS_INITIAL_STATE_FILENAME: key of `inputs_dataframes` holding the initial state at soil scale.
        :param str INPUTS_DIRPATH: path of the directory containing the model inputs (serialised ADEL-Wheat scene/MTG, ``phytoT.csv``, ...).
        :param bool single_plant: if True, build a single isolated plant with :class:`AdelDyn` instead of a stand generated from `plant_density` and `inter_row`.
        :param dict [int, float] plant_density: the density of plants (plants.m-2), one key per plant id.
        :param float inter_row: the inter-row distance (m), used to build the canopy stand when `single_plant` is False.
        :param float sowing_depth: the sowing depth (m, stored as `self.Zsowing`), used in the perceived-temperature computation of the morphogenesis model.
        :param dict or None N_fertilizations: nitrogen fertilization events applied to the soil nitrate pool. Either {time_step_in_hours: amount of nitrates added (umol) at that hour, ...},
                                            or {'constant_Conc_Nitrates': value} to force a constant soil nitrate concentration instead of punctual additions. Only used if `cnwgrass_roots` is True.
        :param dict [str, float] or None tillers_replications: a dictionary with tiller id as key, and weight of replication as value, used to scale up CN-Metabolism fluxes for unreplicated tillers.
        :param str or list or None stored_times: time steps at which the model outputs are stored. Either 'all', a list of time steps (hours), or an empty list.
        :param bool computing_light_interception: whether to run the Caribu facade to compute light interception; if False, organs keep whatever `PARa`/`Erel` they were initialised with.
        :param bool external_soil_model: whether an external soil model is coupled to cnmetabolism. If True, cnmetabolism will skip its own soil calculations and expect root N uptake to be forced
                                            through the `step_callback['nitrate_uptake']` function.
        :param bool heterogeneous_canopy: whether to create a duplicated heterogeneous canopy from the initial MTG (passed to the Caribu facade at each run).
        :param bool show_3Dplant: whether to plot the 3D scene with PlantGL after each geometry update.
        :param bool hydraulics: if True, couple the model to the turgor-driven hydraulics model (xylem water potentials and turgor-driven growth); also enables the 'Tuzet' and 'hydraulics' stomatal
                                            conductance models.
        :param str stomatal_model_name: the model of stomatal conductance. Should be one of 'BWB', 'Leuning', 'Tuzet' or 'hydraulics'. 'Tuzet' and 'hydraulics' require `hydraulics` to be True.
        :param dict or None drought_trigger: a dict to define an external drought-control scenario, {'trigger_variable': value}. For now only implemented for the 'green_area' variable
                                            (value = total green area below which the drought is triggered). Only used if `hydraulics` and `cnwgrass_roots` are True.
        :param dict or None rehydration_scenario: a dict to specify the rehydration scenario following a drought event,
                                            {'stop_drought_SRWC': SRWC at which the drought event stops (%), 'SRWC_target': target SRWC for rehydration (%), 'rehydration_duration': duration of the
                                            rehydration period (days)}. Only used if `hydraulics` and `cnwgrass_roots` are True.
        :param bool optimal_growth_option: if True, the morphogenesis model assumes optimal growth conditions instead of metabolism-limited growth.
        :param bool option_static: whether the model should be run on a static plant architecture (no new organ initiation from morphogenesis).
        :param bool isolated_roots: if True, the root/xylem compartments are integrated by the ODE solvers separately from the shoot compartments, so that roots can be driven by another,
                                            independently-stepped model. Forced to True whenever `cnwgrass_roots` is False.
        :param bool cnwgrass_roots: whether the CNW-Grass root model is used to represent the roots and soil. Set to False when roots and soil are entirely delegated to a coupled belowground
                                            model (e.g. Root-BRIDGES via `root_mtg`), in which case `isolated_roots` is forced to True.
        :param bool UPDATE_SHARED_DF: if True, update the shared inputs/outputs dataframes of every facade at init and after each run; otherwise they are only filled when outputs are stored.
        :param int START_TIME: the first time step of the simulation (hours), used to initialise `self.time_step_in_hours` and, when resuming from previous outputs, as the first index read from `meteo`.
        :param int CARIBU_TIMESTEP: number of hours between two calls to the Caribu light-interception facade.
        :param int MORPHOGENESIS_TIMESTEP: number of hours between two calls to the morphogenesis facade.
        :param int GROWTH_TIMESTEP: number of hours between two calls to the growth facade.
        :param int CNMETABOLISM_TIMESTEP: number of hours between two calls to the CN-metabolism facade.
        :param int SENESCENCE_TIMESTEP: number of hours between two calls to the senescence facade; also used as the increment of `self.time_step_in_hours` at the end of each `__call__`.
        :param int HYDRAULICS_TIMESTEP: number of hours between two calls to the hydraulics facade (used only if `hydraulics` is True).
        """
        # SELF STORAGE FOR LOOP PARAMETERS
        self.meteo = meteo

        # time steps
        self.time_step_in_hours = START_TIME
        self.CARIBU_TIMESTEP = CARIBU_TIMESTEP
        self.SENESCENCE_TIMESTEP = SENESCENCE_TIMESTEP
        self.MORPHOGENESIS_TIMESTEP = MORPHOGENESIS_TIMESTEP
        self.GROWTH_TIMESTEP = GROWTH_TIMESTEP
        self.CNMETABOLISM_TIMESTEP = CNMETABOLISM_TIMESTEP
        self.HYDRAULICS_TIMESTEP = HYDRAULICS_TIMESTEP

        # canopy parameters
        # Management data
        self.plant_density = plant_density
        self.inter_row = inter_row
        self.Zsowing = sowing_depth
        self.N_fertilizations = N_fertilizations

        # plant parameters
        self.tillers_replications = tillers_replications

        # logging and data structures
        self.stored_times = stored_times
        self.shared_axes_inputs_outputs_df = pd.DataFrame()
        self.shared_organs_inputs_outputs_df = pd.DataFrame()
        self.shared_hiddenzones_inputs_outputs_df = pd.DataFrame()
        self.shared_elements_inputs_outputs_df = pd.DataFrame()
        self.shared_soils_inputs_outputs_df = pd.DataFrame()
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
        self.external_soil_model = external_soil_model
        self.hydraulics = hydraulics
        self.option_static = option_static
        self.cnwgrass_roots = cnwgrass_roots

        # -- ADEL and MTG CONFIGURATION --

        # read adelwheat inputs at t0
        if single_plant:
            self.adel_wheat = AdelDyn(seed=1, scene_unit='m', leaves=echap_leaves(xy_model='Soissons_byleafclass'))
        else:
            stand = AgronomicStand(sowing_density=self.plant_density[1], plant_density=self.plant_density[1], inter_row=inter_row, noise=0.) #todo to be adapted if multiple cultivars
            self.adel_wheat = AdelDyn(seed=1, scene_unit='m', leaves=echap_leaves(xy_model='Soissons_byleafclass'), stand=stand)
        
        # MTG generation
        if step_callback is not None and 'ADEL_mtg' in step_callback.keys():
            nff = update_parameters_all_models['morphogenesis']['max_nb_leaves']
            self.g = step_callback['ADEL_mtg'](self.adel_wheat, INPUTS_DIRPATH, nff)  # Create a new MTG
        else:
            self.g = self.adel_wheat.load(directory=INPUTS_DIRPATH)  # read adelwheat inputs at t0 from a serialised MTG
        
        # Section specific to coupling with Root-BRIDGES
        self.shoot_props = self.g.properties()

        self.props = root_mtg.properties()
        self.vertices = root_mtg.vertices(scale=root_mtg.max_scale())

        self.link_self_to_mtg()

        # ---------------------------------------------
        # ----- CONFIGURATION OF THE FACADES -------
        # ---------------------------------------------

        # -- MORPHOGENESIS (created first because it is the only facade to add new metamers) --
        # Initial states
        morphogenesis_hiddenzones_initial_state = inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME]
        morphogenesis_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME]
        morphogenesis_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME]

        phytoT_df = pd.read_csv(os.path.join(INPUTS_DIRPATH, 'phytoT.csv'))

        # Update parameters if specified
        if update_parameters_all_models and 'morphogenesis' in update_parameters_all_models:
            update_parameters_morphogenesis = update_parameters_all_models['morphogenesis']
        else:
            update_parameters_morphogenesis = None

        # Facade initialisation
        self.morphogenesis_facade_ = morphogenesis_facade.MorphogenesisFacade(self.g,
                                                                MORPHOGENESIS_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                                morphogenesis_axes_initial_state,
                                                                morphogenesis_hiddenzones_initial_state,
                                                                morphogenesis_elements_initial_state,
                                                                self.shared_axes_inputs_outputs_df,
                                                                self.shared_hiddenzones_inputs_outputs_df,
                                                                self.shared_elements_inputs_outputs_df,
                                                                self.adel_wheat, phytoT_df,
                                                                hydraulics=hydraulics,
                                                                optimal_growth_option=optimal_growth_option,
                                                                option_static=option_static,
                                                                update_parameters=update_parameters_morphogenesis,
                                                                update_shared_df=UPDATE_SHARED_DF)

        # -- CARIBU --
        if self.computing_light_interception:
            self.caribu_facade_ = caribu_facade.CaribuFacade(self.g,
                                                        self.shared_elements_inputs_outputs_df,
                                                        self.adel_wheat,
                                                        update_shared_df=UPDATE_SHARED_DF)

        # -- SENESCENCE --
        # Initial states    
        senescence_roots_initial_state = inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].loc[inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME]['organ'] == 'roots'][
            senescence_facade.converter.ROOTS_TOPOLOGY_COLUMNS +
            [i for i in senescence_facade.converter.SENESCENCE_ROOTS_INPUTS if i in inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].columns]].copy()

        senescence_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
            senescence_facade.converter.ELEMENTS_TOPOLOGY_COLUMNS +
            [i for i in senescence_facade.converter.SENESCENCE_ELEMENTS_INPUTS if i in inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns]].copy()

        senescence_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
            senescence_facade.converter.AXES_TOPOLOGY_COLUMNS +
            [i for i in senescence_facade.converter.SENESCENCE_AXES_INPUTS if i in inputs_dataframes[AXES_INITIAL_STATE_FILENAME].columns]].copy()

        # Update parameters if specified
        if update_parameters_all_models and 'senescence' in update_parameters_all_models:
            update_parameters_senescence = update_parameters_all_models['senescence']
        else:
            update_parameters_senescence = None

        # Facade initialisation
        self.senescence_facade_ = senescence_facade.SENESCENCEFacade(self.g,
                                                                SENESCENCE_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                                senescence_roots_initial_state,
                                                                senescence_axes_initial_state,
                                                                senescence_elements_initial_state,
                                                                self.shared_organs_inputs_outputs_df,
                                                                self.shared_axes_inputs_outputs_df,
                                                                self.shared_elements_inputs_outputs_df,
                                                                update_parameters=update_parameters_senescence,
                                                                update_shared_df=UPDATE_SHARED_DF,
                                                                cnwgrass_roots=cnwgrass_roots)

        # -- GAS-EXCHANGE --
        # Initial states
        gasexchange_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME]
        gasexchange_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME]

        # Update parameters if specified
        if update_parameters_all_models and 'gasexchange' in update_parameters_all_models:
            update_parameters_gasexchange = update_parameters_all_models['gasexchange']
        else:
            update_parameters_gasexchange = None

        # Facade initialisation
        self.gasexchange_facade_ = gasexchange_facade.GasExchangeFacade(self.g,
                                                                        gasexchange_elements_initial_state,
                                                                        gasexchange_axes_initial_state,
                                                                        self.shared_elements_inputs_outputs_df,
                                                                        stomatal_model_name=stomatal_model_name,
                                                                        hydraulics=hydraulics,
                                                                        update_parameters=update_parameters_gasexchange,
                                                                        update_shared_df=UPDATE_SHARED_DF)

        # -- GROWTH --
        # Initial states
        growth_hiddenzones_initial_state = inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME]
        growth_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME]
        growth_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME]
        growth_root_initial_state = inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME]

        # Update parameters if specified
        if update_parameters_all_models and 'growth' in update_parameters_all_models:
            update_parameters_growth = update_parameters_all_models['growth']
        else:
            update_parameters_growth = None

        # Facade initialisation
        self.growth_facade_ = growth_facade.GrowthFacade(self.g,
                                                                GROWTH_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                                growth_hiddenzones_initial_state,
                                                                growth_elements_initial_state,
                                                                growth_root_initial_state,
                                                                growth_axes_initial_state,
                                                                self.shared_organs_inputs_outputs_df,
                                                                self.shared_hiddenzones_inputs_outputs_df,
                                                                self.shared_elements_inputs_outputs_df,
                                                                self.shared_axes_inputs_outputs_df,
                                                                hydraulics=hydraulics,
                                                                update_parameters=update_parameters_growth,
                                                                update_shared_df=UPDATE_SHARED_DF,
                                                                cnwgrass_roots=cnwgrass_roots)

        # -- CNMETABOLISM --
        # Initial states
        cnmetabolism_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
            [i for i in cnmetabolism_facade.cnmetabolism_converter.AXES_VARIABLES if i in inputs_dataframes[AXES_INITIAL_STATE_FILENAME].columns]].copy()

        cnmetabolism_organs_initial_state = inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME][
            [i for i in cnmetabolism_facade.cnmetabolism_converter.ORGANS_VARIABLES if i in inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].columns]].copy()

        cnmetabolism_hiddenzones_initial_state = inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME][
            [i for i in cnmetabolism_facade.cnmetabolism_converter.HIDDENZONE_VARIABLES if i in inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME].columns]].copy()

        cnmetabolism_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
            [i for i in cnmetabolism_facade.cnmetabolism_converter.ELEMENTS_VARIABLES if i in inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns]].copy()

        cnmetabolism_soils_initial_state = inputs_dataframes[SOILS_INITIAL_STATE_FILENAME][
            [i for i in cnmetabolism_facade.cnmetabolism_converter.SOILS_VARIABLES if i in inputs_dataframes[SOILS_INITIAL_STATE_FILENAME].columns]].copy()

        if cnwgrass_roots: # Not needed if there is a substitute root and soil model
            if not hydraulics and 'SRWC' not in inputs_dataframes[SOILS_INITIAL_STATE_FILENAME].columns:
                cnmetabolism_soils_initial_state['SRWC'] = 100
            elif hydraulics and 'SRWC' not in inputs_dataframes[SOILS_INITIAL_STATE_FILENAME].columns:
                raise(ValueError('Hydraulics option is True but SRWC not found in {}.'.format(SOILS_INITIAL_STATE_FILENAME)))

        # Update parameters if specified
        if update_parameters_all_models and 'cnmetabolism' in update_parameters_all_models:
            update_parameters_cnmetabolism = update_parameters_all_models['cnmetabolism']
        else:
            update_parameters_cnmetabolism = {}

        # Force solver separation if a different root model has been chosen
        if not cnwgrass_roots:
            isolated_roots = True

        # Facade initialisation
        self.cnmetabolism_facade_ = cnmetabolism_facade.CNMetabolismFacade(self.g,
                                                    CNMETABOLISM_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                    plant_density, # TODO GB check why it is needed
                                                    update_parameters_cnmetabolism,
                                                    cnmetabolism_axes_initial_state,
                                                    cnmetabolism_organs_initial_state,
                                                    cnmetabolism_hiddenzones_initial_state,
                                                    cnmetabolism_elements_initial_state,
                                                    cnmetabolism_soils_initial_state,
                                                    self.shared_axes_inputs_outputs_df,
                                                    self.shared_organs_inputs_outputs_df,
                                                    self.shared_hiddenzones_inputs_outputs_df,
                                                    self.shared_elements_inputs_outputs_df,
                                                    self.shared_soils_inputs_outputs_df,
                                                    tillers_replications=tillers_replications,
                                                    external_soil_model=external_soil_model,
                                                    update_shared_df=UPDATE_SHARED_DF,
                                                    isolated_roots=isolated_roots,
                                                    cnwgrass_roots=cnwgrass_roots)

        if cnwgrass_roots:
            # Run cn-metabolism with constant nitrates concentration in the soil if specified
            if N_fertilizations is not None: 
                if 'constant_Conc_Nitrates' in N_fertilizations.keys():
                    self.cnmetabolism_facade_.soils[(1, 'MS')].constant_Conc_Nitrates = True  # TODO: make (1, 'MS') more general
                    self.cnmetabolism_facade_.soils[(1, 'MS')].nitrates = N_fertilizations['constant_Conc_Nitrates'] * self.cnmetabolism_facade_.soils[(1, 'MS')].volume

            # Force root nitrate uptake if specified
            if external_soil_model and step_callback is not None:
                try:
                    step_callback['nitrate_uptake'](0, self.cnmetabolism_facade_.population, self.g)
                except KeyError:
                    print('Function name error in step_callback keys. It should be nitrate_uptake')


        # -- HYDRAULICS --
        drought_ongoing = False  # Is a drought event ongoing (bool)
        drought_passed = False   # Has a drought event occurred (bool)
        rehydration = False      # Is a rehydration period ongoing (bool)
        self.hydraulics_facade_ = None

        if hydraulics:
            # Initial states
            hydraulics_axes_initial_state = inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
                [i for i in hydraulics_facade.hydraulics_converter.AXES_VARIABLES if
                i in inputs_dataframes[AXES_INITIAL_STATE_FILENAME].columns]].copy()

            hydraulics_organs_initial_state = inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME][
                [i for i in hydraulics_facade.hydraulics_converter.ORGANS_VARIABLES if
                i in inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].columns]].copy()

            hydraulics_hiddenzones_initial_state = inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME][
                [i for i in hydraulics_facade.hydraulics_converter.HIDDENZONE_VARIABLES if
                i in inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME].columns]].copy()

            hydraulics_elements_initial_state = inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
                [i for i in hydraulics_facade.hydraulics_converter.ELEMENTS_VARIABLES if
                i in inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns]].copy()

            # Update parameters if specified
            if update_parameters_all_models and 'hydraulics' in update_parameters_all_models:
                update_parameters_hydraulics = update_parameters_all_models['hydraulics']
            else:
                update_parameters_hydraulics = {}

            hydraulics_soils_initial_state = inputs_dataframes[SOILS_INITIAL_STATE_FILENAME][
                [i for i in hydraulics_facade.hydraulics_converter.SOILS_VARIABLES if
                i in inputs_dataframes[SOILS_INITIAL_STATE_FILENAME].columns]].copy()

            # Facade initialisation
            self.hydraulics_facade_ = hydraulics_facade.hydraulicsFacade(self.g,
                                                                    HYDRAULICS_TIMESTEP * HOUR_TO_SECOND_CONVERSION_FACTOR,
                                                                    update_parameters_hydraulics,
                                                                    hydraulics_axes_initial_state,
                                                                    hydraulics_hiddenzones_initial_state,
                                                                    hydraulics_elements_initial_state,
                                                                    hydraulics_organs_initial_state,
                                                                    hydraulics_soils_initial_state,
                                                                    self.shared_axes_inputs_outputs_df,
                                                                    self.shared_hiddenzones_inputs_outputs_df,
                                                                    self.shared_elements_inputs_outputs_df,
                                                                    self.shared_organs_inputs_outputs_df,
                                                                    self.shared_soils_inputs_outputs_df,
                                                                    update_shared_df=UPDATE_SHARED_DF,
                                                                    isolated_roots=isolated_roots,
                                                                    cnwgrass_roots=cnwgrass_roots)

        # -- MODEL INTEGRATION --
        # Facade initialisation
        self.build_outputs_ = build_outputs.BuildOutputs(self.g, self.morphogenesis_facade_, self.growth_facade_, self.gasexchange_facade_, self.hydraulics_facade_)

        # Update geometry
        self.adel_wheat.update_geometry(self.g)
        if show_3Dplant:
            self.adel_wheat.plot(self.g)
        
        self.cnw_grass_root_props = self.g.get_vertex_property(2)["roots"]

        # TODO GB : Temporary inputs initialization at Component's prescribed value
        self.cnw_grass_root_props["Unloading_Sucrose"] = self.props["Unloading_Sucrose"][1]
        self.cnw_grass_root_props["Unloading_Amino_Acids"] = self.props["Unloading_Amino_Acids"][1]
        self.g.get_vertex_property(2)['xylem']['water_potential'] = self.props["xylem_water_potential"][1]

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

        # TODO GB Temporary initialization of model's outputs with component's prescribed values
        self.g.get_vertex_property(2)['phloem']['sucrose'] = self.props["sucrose_phloem"][1]
        self.g.get_vertex_property(2)['phloem']['amino_acids'] = self.props["amino_acids_phloem"][1]
        self.g.get_vertex_property(2)['phloem']['Unloading_Sucrose_shoot_organs'] = 30.
        self.g.get_vertex_property(2)['phloem']['Unloading_Amino_Acids_shoot_organs'] = 1.
        self.g.get_vertex_property(2)['xylem']['root_to_shoot_xylem_water_flow'] = self.props["root_to_shoot_xylem_water_flow"][1]

        self.sync_shoot_outputs_with_root_mtg()
        

    def sync_shoot_inputs_with_shoot_mtg(self):
        for name in self.inputs:
            if name == "Unloading_Sucrose_phloem":
                self.cnw_grass_root_props["Unloading_Sucrose"] = self.props[name][1] / self.props['mstruct'][1]
                # print(self.cnw_grass_root_props["Unloading_Sucrose"])

            elif name == "Unloading_Amino_Acids_phloem":
                self.cnw_grass_root_props["Unloading_Amino_Acids"] = self.props[name][1] / self.props['mstruct'][1]

            elif name == "xylem_water_potential":
                self.g.get_vertex_property(2)['xylem']["water_potential"] = self.props[name][1]

            else:
                self.cnw_grass_root_props[name] = self.props[name][1]

    def sync_shoot_outputs_with_root_mtg(self):
        # Link this specific data structure to self for variables exchange, only for outputs that will be read by other models  here.
        # Note : here eval is necessary to ensure intended lambda function definition
        for name in self.state_variables:
            if name == "root_to_shoot_xylem_water_flow":
                self.props[name].update({1: self.g.get_vertex_property(2)['xylem']['root_to_shoot_xylem_water_flow']})

            elif name == "mstruct_axis":
                self.props[name].update({1: self.g.get_vertex_property(2)['mstruct']})

            elif name == "amino_acids_phloem":
                self.props[name].update({1: self.g.get_vertex_property(2)['phloem']['amino_acids']})

            elif name == "sucrose_phloem":
                self.props[name].update({1: self.g.get_vertex_property(2)['phloem']['sucrose']})
            
            # Was added to Simulation.ORGANS_INTEGRATIVE_VARIABLES to be recognized
            elif name == "Unloading_Sucrose_shoot_organs":
                self.props[name].update({1: self.g.get_vertex_property(2)['phloem']['Unloading_Sucrose_shoot_organs']})

            elif name == "Unloading_Amino_Acids_shoot_organs":
                self.props[name].update({1: self.g.get_vertex_property(2)['phloem']['Unloading_Amino_Acids_shoot_organs']})

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
                self.props[name].update({1: self.cnw_grass_root_props[name]})

        


    def __call__(self):
        # SPECIFIC TO COUPLING, syncs the shoot mtg variables with bellowground models
        self.pull_available_inputs()
        self.sync_shoot_inputs_with_shoot_mtg()

        # run Caribu
        PARi = self.meteo.loc[self.time_step_in_hours, ['PARi']].iloc[0]
        DOY = self.meteo.loc[self.time_step_in_hours, ['DOY']].iloc[0]
        hour = self.meteo.loc[self.time_step_in_hours, ['hour']].iloc[0]

        # Run Caribu, if activated
        if self.computing_light_interception:
            PARi_next_hours = self.meteo.loc[range(self.time_step_in_hours, self.time_step_in_hours + self.CARIBU_TIMESTEP), ['PARi']].sum().values[0]

            if (self.time_step_in_hours % self.CARIBU_TIMESTEP == 0) and (PARi_next_hours > 0):
                run_caribu = True
            else:
                run_caribu = False

            self.caribu_facade_.run(run_caribu, energy=PARi, DOY=DOY, hourTU=hour, latitude=48.85, sun_sky_option='sky', 
                                    heterogeneous_canopy=self.heterogeneous_canopy, plant_density=self.plant_density[1], inter_row=inter_row)

        # Run Senescence
        self.senescence_facade_.run()

        # Test for dead plant # TODO: adapt in case of multiple plants
        if not self.shared_elements_inputs_outputs_df.empty and \
                np.nansum(self.shared_elements_inputs_outputs_df.loc[self.shared_elements_inputs_outputs_df['element'].isin(['StemElement', 'LeafElement1']), 'green_area']) == 0:
            # append the inputs and outputs at current step to global lists
            self.all_simulation_steps.append(self.time_step_in_hours)
            self.axes_all_data_list.append(self.shared_axes_inputs_outputs_df.copy())
            self.organs_all_data_list.append(self.shared_organs_inputs_outputs_df.copy())
            self.hiddenzones_all_data_list.append(self.shared_hiddenzones_inputs_outputs_df.copy())
            self.elements_all_data_list.append(self.shared_elements_inputs_outputs_df.copy())
            self.soils_all_data_list.append(self.shared_soils_inputs_outputs_df.copy())
            raise RuntimeError('Dead plant')

        # Run the rest of the model if the plant is alive
        # get the meteo of the current step
        Ta, ambient_CO2, RH, Ur = self.meteo.loc[self.time_step_in_hours, ['air_temperature', 'ambient_CO2', 'humidity', 'Wind']]

        # run GasExchange
        self.gasexchange_facade_.run(Ta, ambient_CO2, RH, Ur)

        # run Morphogeneis
        Tair, Tsoil = self.meteo.loc[self.time_step_in_hours, ['air_temperature', 'soil_temperature']]
        self.morphogenesis_facade_.run(Tair, Tsoil, self.Zsowing) # NOTE Zsowing needed for perceived temperature computation

        # Update geometry
        self.adel_wheat.update_geometry(self.g)
        if self.show_3Dplant:
            self.adel_wheat.plot(self.g)

        # run hydraulics
        if self.hydraulics:
            if self.cnwgrass_roots:
                turgor_soil = hydraulics_facade_.soils[(1, 'MS')]
                # Trigger drought
                if drought_trigger is not None and 'green_area' in drought_trigger.keys():
                    if  (sum(g.property('green_area').values()) >= drought_trigger['green_area'] or drought_ongoing) and not drought_passed:
                        drought_ongoing = True
                        turgor_soil.constant_water_content = False
                    # Rehydration scenario. Only implemented for a hourly and linear rehydration scenario.
                    if rehydration_scenario is not None:
                        # Maximum of drought, start of rehydration
                        if turgor_soil.SRWC <= rehydration_scenario['stop_drought_SRWC'] and not rehydration:
                            rehydration = True
                            total_irrigation = (rehydration_scenario['SRWC_target'] * turgor_soil.PARAMETERS.AWC) / 100 - turgor_soil.water_content  # Total amount of water to add to the soil in order to reach the target SRWC
                            turgor_soil.hourly_irrigation = total_irrigation / (rehydration_scenario['rehydration_duration'] * 24)  # Amount of water to add each hour to reach the target SRWC at the end of the rehydration period

                        # Ongoing rehydration
                        elif rehydration:
                            # Target SRWC reached after rehydration, end of drought event
                            if turgor_soil.SRWC >= rehydration_scenario['SRWC_target']:
                                rehydration = False
                                drought_ongoing = False
                                drought_passed = True
                                turgor_soil.water_content = (rehydration_scenario['SRWC_target'] * turgor_soil.PARAMETERS.AWC) / 100
                                turgor_soil.SRWC = rehydration_scenario['SRWC_target']
                                turgor_soil.constant_water_content = True
                                turgor_soil.hourly_rehydration = 0

            self.hydraulics_facade_.run()

            # Update geometry
            self.adel_wheat.update_geometry(self.g)
            if self.show_3Dplant:
                self.adel_wheat.plot(self.g)

        # run Growth
        self.growth_facade_.run()

        # run CN-Metabolism
        if self.cnwgrass_roots:
            # N fertilization if any
            if self.N_fertilizations is not None:
                if self.time_step_in_hours in self.N_fertilizations.keys():
                    self.cnmetabolism_facade_.soils[(1, 'MS')].nitrates += self.N_fertilizations[self.time_step_in_hours]

            # Force root nitrate uptake if specified
            if self.external_soil_model and step_callback is not None:
                try:
                    step_callback['nitrate_uptake'](self.time_step_in_hours, cnmetabolism_facade_.population, g)
                except KeyError:
                    print(
                        'Function name error in step_callback keys. It should be nitrate_uptake')

        self.cnmetabolism_facade_.run(Tair, Tsoil, self.tillers_replications)

        # append outputs at current step to global lists
        if (self.stored_times == 'all') or (self.time_step_in_hours in self.stored_times):
            axes_outputs, elements_outputs, hiddenzones_outputs, organs_outputs, soils_outputs = self.build_outputs_.build_outputs_df_from_MTG()

            self.all_simulation_steps.append(self.time_step_in_hours)
            self.axes_all_data_list.append(axes_outputs)
            self.organs_all_data_list.append(organs_outputs)
            self.hiddenzones_all_data_list.append(hiddenzones_outputs)
            self.elements_all_data_list.append(elements_outputs)
            self.soils_all_data_list.append(soils_outputs)

        self.time_step_in_hours += self.SENESCENCE_TIMESTEP

        self.sync_shoot_outputs_with_root_mtg()



def scenario_utility(time_step_in_seconds: int = 3600, INPUTS_DIRPATH = "inputs", OUTPUTS_DIRPATH = "outputs", METEO_FILENAME = "meteo_Ljutovac2002.csv", plant_density = {1:250},
                     forced_start_time = 0, tillers_replications={'T1': 0.5, 'T2': 0.5, 'T3': 0.5, 'T4': 0.5}, N_fertilizations = {2016: 357143, 2520: 1000000},
                     stored_times = None, option_static = False, single_plant = False, show_3Dplant = False, run_from_outputs = False, heterogeneous_canopy = True, update_parameters_all_models = None,
                     isolated_roots = False, cnwgrass_roots = True, hydraulics = False, stomatal_model_name = 'BWB'):
    scenario = {}

    ### DIRS ###

    scenario["INPUTS_DIRPATH"] = INPUTS_DIRPATH
    
    # Save the outputs with a full scan of the MTG at each time step (or at selected time steps)
    UPDATE_SHARED_DF = False
    if stored_times is None:
        stored_times = 'all'
    if not (stored_times == 'all' or isinstance(stored_times, list)):
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
        previous_axes_outputs_dataframe = pd.read_csv(os.path.join(OUTPUTS_DIRPATH, AXES_OUTPUTS_FILENAME))
        assert 't' in previous_axes_outputs_dataframe.columns
        if forced_start_time > 0:
            new_start_time = forced_start_time + 1
        else:
            last_t_step = int(previous_axes_outputs_dataframe['t'].max())
            new_start_time = last_t_step + 1

        previous_outputs_dataframes = {}

        for initial_state_filename, outputs_filename, index_columns in ((AXES_INITIAL_STATE_FILENAME, AXES_OUTPUTS_FILENAME, AXES_INDEX_COLUMNS),
                                                                        (ORGANS_INITIAL_STATE_FILENAME, ORGANS_OUTPUTS_FILENAME, ORGANS_INDEX_COLUMNS),
                                                                        (HIDDENZONES_INITIAL_STATE_FILENAME, HIDDENZONES_OUTPUTS_FILENAME, HIDDENZONES_INDEX_COLUMNS),
                                                                        (ELEMENTS_INITIAL_STATE_FILENAME, ELEMENTS_OUTPUTS_FILENAME, ELEMENTS_INDEX_COLUMNS),
                                                                        (SOILS_INITIAL_STATE_FILENAME, SOILS_OUTPUTS_FILENAME, SOILS_INDEX_COLUMNS)):

            previous_outputs_dataframe = pd.read_csv(os.path.join(OUTPUTS_DIRPATH, outputs_filename))
            # Convert NaN to None
            previous_outputs_dataframes[outputs_filename] = previous_outputs_dataframe.where(previous_outputs_dataframe.notnull(), None)

            # assert 't' in previous_outputs_dataframes[outputs_filename].columns
            if forced_start_time > 0:
                previous_outputs_dataframes[outputs_filename] = previous_outputs_dataframes[outputs_filename][previous_outputs_dataframes[outputs_filename]['t'] <= forced_start_time]

            if initial_state_filename == ELEMENTS_INITIAL_STATE_FILENAME:
                elements_previous_outputs = previous_outputs_dataframes[outputs_filename]
                new_initial_state = elements_previous_outputs[~elements_previous_outputs.is_over.isnull()]
            else:
                new_initial_state = previous_outputs_dataframes[outputs_filename]
            idx = new_initial_state.groupby([col for col in index_columns if col != 't'])['t'].transform(max) == new_initial_state['t']
            inputs_dataframes[initial_state_filename] = new_initial_state[idx].drop(['t'], axis=1)

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
    scenario["MORPHOGENESIS_TIMESTEP"] = time_step_in_hours
    scenario["GROWTH_TIMESTEP"] = time_step_in_hours
    scenario["CNMETABOLISM_TIMESTEP"] = time_step_in_hours
    scenario["SENESCENCE_TIMESTEP"] = time_step_in_hours

    scenario["inputs_dataframes"] = inputs_dataframes

    ### OPTIONS ###
    scenario["plant_density"] = plant_density
    scenario["option_static"] = option_static
    scenario["show_3Dplant"] = show_3Dplant
    scenario["tillers_replications"] = tillers_replications
    scenario["heterogeneous_canopy"] = heterogeneous_canopy
    scenario["N_fertilizations"] = N_fertilizations
    scenario["update_parameters_all_models"] = update_parameters_all_models
    scenario["single_plant"] = single_plant
    scenario["isolated_roots"] = isolated_roots
    scenario["cnwgrass_roots"] = cnwgrass_roots
    scenario["hydraulics"] = hydraulics
    scenario["stomatal_model_name"] = stomatal_model_name


    return scenario

