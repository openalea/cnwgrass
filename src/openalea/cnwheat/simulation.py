# -*- coding: latin-1 -*-

from __future__ import division  # use "//" to do integer division
import logging

import numpy as np
from scipy.integrate import solve_ivp
from scipy import interpolate

from openalea.cnwheat import model
from openalea.cnwheat import tools

"""
    cnwheat.simulation
    ~~~~~~~~~~~~~~~~~~

    The module :mod:`cnwheat.simulation` is the front-end to run the model CN-Wheat.
    The public API consists of methods :meth:`initialize` and :meth:`run`.

"""


class SimulationError(Exception):
    """
    Abstract class for the management of simulation errors. Do not instance it directly.
    """
    pass


class SimulationConstructionError(SimulationError):
    """
    Exception raised when a problem occurs in the constructor, in particular
    when the arguments are not consistent with each other.
    """
    pass


class SimulationInitializationError(SimulationError):
    """
    Exception raised when a problem occurs at initialization time, in particular
    when checking the consistency of inputs `population` and `soils` (see :meth:`initialize`).
    """
    pass


class SimulationRunError(SimulationError):
    """
    Exception raised when running a simulation, for example when a problem occurs
    during the integration of the system of differential equations.
    """
    pass


class Simulation(object):
    """
    The Simulation class permits to initialize and run the model.

    User should use method :meth:`initialize` to initialize the model, and method
    :meth:`run` to run the model.

    :param class respiration_model: the model of respiration to use.
          This model must define a class implementing these functions:
            * R_Nnit_upt(U_Nnit, sucrose): Nitrate uptake respiration.
                * Parameters:
                    - `U_Nnit` (:class:`float`) - uptake of N nitrates (µmol` N)
                    - `sucrose` (:class:`float`) -  amount of C sucrose in organ (µmol` C)
                * Returns: _R_Nnit_upt (µmol` C respired)
                * Returns Type: :class:`float`

            * R_phloem(sucrose_loading, sucrose, mstruct): Phloem loading respiration
                * Parameters:
                    - `sucrose_loading` (:class:`float`) -  Loading flux from the C substrate pool to phloem (µmol` C g-1 mstruct)
                    - `sucrose` (:class:`float`) -  amount of C sucrose in organ (µmol` C)
                    - `mstruct` (:class:`float`) -  structural dry mass of organ (g)
                * Returns: _R_phloem (µmol` C respired)
                * Returns Type: :class:`float`

            * R_Nnit_red(s_amino_acids, sucrose, mstruct, root=False): Nitrate reduction-linked respiration
              Distinction is made between nitrate realised in roots or in shoots where a part of the energy required is derived from ATP
              and reducing power obtained directly from photosynthesis (rather than C substrate)

                * Parameters:
                    - `s_amino_acids` (:class:`float`) - consumption of N for the synthesis of amino acids (µmol` N g-1 mstruct)
                      (in the present version, this is used to approximate nitrate reduction needed in the original model of Thornley and Cannell, 2000)
                    - `sucrose` (:class:`float`) -  amount of C sucrose in organ (µmol` C)
                    - `mstruct` (:class:`float`) -  structural dry mass of organ (g)
                    - `root` (:class:`bool`) - specifies if the nitrate reduction-linked respiration is computed for shoot (False) or root (True) tissues.
                * Returns: _R_Nnit_upt (µmol` C respired)
                * Returns Type: :class:`float`

            * R_residual(sucrose, mstruct, Ntot, delta_t, Ts): Residual maintenance respiration (cost from protein turn-over, cell ion gradients, futile cycles...)
                * Parameters:
                    - `sucrose` (:class:`float`) - amount of C sucrose (µmol` C)
                    - `mstruct` (:class:`float`) - structural dry mass of organ (g)
                    - `Ntot` (:class:`float`) - total N in organ (µmol` N)
                    - `delta_t` (:class:`float`) - timestep (s)
                    - `Ts` (:class:`float`) - organ temperature (°C)
                * Returns: _R_residual (µmol` C respired)
                * Returns Type: :class:`float`

            * R_grain_growth(mstruct_growth, starch_filling, mstruct): Grain growth respiration
                * Parameters:
                    - `mstruct_growth` (:class:`float`) - gross growth of grain structure (µmol` C added in grain structure)
                    - `starch_filling` (:class:`float`) - gross growth of grain starch (µmol` C added in grain starch g-1 mstruct)
                    - `mstruct` (:class:`float`) -  structural dry mass of organ (g)
                * Returns: R_grain_growth (µmol` C respired)
                * Returns Type: :class:`float`

    :param int delta_t: the delta t of the simulation (in seconds) ; default is `1`.
    :param dict [int, int] culm_density: culm density (culm m-2).
    :param bool interpolate_forcings: if True: interpolate senescence and photosynthesis forcings from values of `senescence_forcings_delta_t`and `senescence_forcings_delta_t`.
           Default is `False` (do not interpolate the forcings).
    :param int senescence_forcings_delta_t: the delta t of the senescence forcings (in seconds) ; default is `None`.
           If the user sets `interpolate_forcings` to `True`, then he/she must also set `senescence_forcings_delta_t` to an integer value greater or equal to `delta_t`.
           For example, if `interpolate_forcings` is `True` and `delta_t==3600`, then `senescence_forcings_delta_t` must be greater or equal to `3600`, that is for example `7200`.
     :param int photosynthesis_forcings_delta_t: the delta t of the photosynthesis forcings (in seconds) ; default is `None`.
           If the user sets `interpolate_forcings` to `True`, then he/she must also set `photosynthesis_forcings_delta_t` to an integer value greater or equal to `delta_t`.
           For example, if `interpolate_forcings` is `True` and `delta_t==3600`, then `photosynthesis_forcings_delta_t` must be greater or equal to `3600`, that is for example `7200`.

        - interpolate_forcings (:class:`bool`) - if True: interpolate senescence and photosynthesis forcings from values of `senescence_forcings_delta_t`
          and `senescence_forcings_delta_t`. Default is `False` (do not interpolate the forcings).

        - senescence_forcings_delta_t (:class:`int`) - the delta t of the senescence forcings (in seconds) ; default is `None`.
          If the user sets `interpolate_forcings` to `True`, then he/she must also set `senescence_forcings_delta_t` to an integer value greater or equal to `delta_t`.
          For example, if `interpolate_forcings` is `True` and `delta_t==3600`, then `senescence_forcings_delta_t` must be greater or equal to `3600`, that is for example `7200`.

        - photosynthesis_forcings_delta_t (:class:`int`) - the delta t of the photosynthesis forcings (in seconds) ; default is `None`.
          If the user sets `interpolate_forcings` to `True`, then he/she must also set `photosynthesis_forcings_delta_t` to an integer value greater or equal to `delta_t`.
          For example, if `interpolate_forcings` is `True` and `delta_t==3600`, then `photosynthesis_forcings_delta_t` must be greater or equal to `3600`, that is for example `7200`.

    """

    #: the name of the compartments attributes in the model, for objects of types
    #: :class:`model.Plant`, :class:`model.Axis`, :class:`model.Phytomer`,
    #: :class:`model.Organ`, :class:`model.HiddenZone`, :class:`model.PhotosyntheticOrganElement`,
    #: and :class:`model.Soil`.
    MODEL_COMPARTMENTS_NAMES = {model.Plant: [],
                                model.Axis: ['C_exudated', 'sum_respi_shoot', 'sum_respi_roots'],
                                model.Phytomer: [],
                                model.Organ: ['age_from_flowering', 'amino_acids', 'cytokinins',
                                              'nitrates', 'proteins', 'starch', 'structure', 'sucrose'],
                                model.HiddenZone: ['amino_acids', 'fructan', 'proteins', 'sucrose'],
                                model.PhotosyntheticOrganElement: ['amino_acids', 'cytokinins', 'fructan',
                                                                   'nitrates', 'proteins', 'starch', 'sucrose', 'triosesP'],
                                model.Soil: ['nitrates']}

    #: the time index
    T_INDEX = ['t']

    # ---------------------------------------------------------------------------
    # ---------- DEFINITION OF THE PARAMETERS AND COMPUTED VARIABLES ------------
    # ---------------------------------------------------------------------------

    # ---------- PLANT scale ----------

    #: the index to locate the plants in the modeled system
    PLANTS_INDEXES = ['plant']
    #: concatenation of :attr:`T_INDEX` and :attr:`PLANTS_INDEXES`
    PLANTS_T_INDEXES = T_INDEX + PLANTS_INDEXES
    #: the parameters which define the state of the modeled system at plant scale
    PLANTS_STATE_PARAMETERS = ['Tair']
    #: the variables which define the state of the modeled system at plant scale,
    #: formed be the concatenation of :attr:`PLANTS_STATE_PARAMETERS` and the names
    #: of the compartments associated to each plant (see :attr:`MODEL_COMPARTMENTS_NAMES`)
    PLANTS_STATE = PLANTS_STATE_PARAMETERS + MODEL_COMPARTMENTS_NAMES.get(model.Plant, [])
    #: the variables that we need to compute in order to compute fluxes and/or compartments values at plant scale
    PLANTS_INTERMEDIATE_VARIABLES = []
    #: the fluxes exchanged between the compartments at plant scale
    PLANTS_FLUXES = []
    #: the variables computed by integrating values of plant components parameters/variables recursively
    PLANTS_INTEGRATIVE_VARIABLES = []
    #: all the variables computed during a run step of the simulation at plant scale
    PLANTS_RUN_VARIABLES = PLANTS_STATE + PLANTS_INTERMEDIATE_VARIABLES + PLANTS_FLUXES + PLANTS_INTEGRATIVE_VARIABLES

    # ---------- AXIS scale ----------

    #: the indexes to locate the axes in the modeled system
    AXES_INDEXES = ['plant', 'axis']
    #: concatenation of :attr:`T_INDEX` and :attr:`AXES_INDEXES`
    AXES_T_INDEXES = T_INDEX + AXES_INDEXES
    #: the parameters which define the state of the modeled system at axis scale
    AXES_STATE_PARAMETERS = ['mstruct', 'senesced_mstruct']
    #: the variables which define the state of the modeled system at axis scale,
    #: formed be the concatenation of :attr:`AXES_STATE_PARAMETERS` and the names
    #: of the compartments associated to each axis (see :attr:`MODEL_COMPARTMENTS_NAMES`)
    AXES_STATE = AXES_STATE_PARAMETERS + MODEL_COMPARTMENTS_NAMES.get(model.Axis, [])
    #: the variables that we need to compute in order to compute fluxes and/or compartments values at axis scale
    AXES_INTERMEDIATE_VARIABLES = []
    #: the fluxes exchanged between the compartments at axis scale
    AXES_FLUXES = []
    #: the variables computed by integrating values of axis components parameters/variables recursively
    AXES_INTEGRATIVE_VARIABLES = ['Total_Transpiration']
    #: all the variables computed during a run step of the simulation at axis scale
    AXES_RUN_VARIABLES = AXES_STATE + AXES_INTERMEDIATE_VARIABLES + AXES_FLUXES + AXES_INTEGRATIVE_VARIABLES

    # ---------- PHYTOMER scale ----------

    #: the indexes to locate the phytomers in the modeled system
    PHYTOMERS_INDEXES = ['plant', 'axis', 'metamer']
    #: concatenation of :attr:`T_INDEX` and :attr:`PHYTOMERS_INDEXES`
    PHYTOMERS_T_INDEXES = T_INDEX + PHYTOMERS_INDEXES
    #: the parameters which define the state of the modeled system at phytomer scale
    PHYTOMERS_STATE_PARAMETERS = ['mstruct']
    #: the variables which define the state of the modeled system at phytomer scale,
    #: formed be the concatenation of :attr:`PHYTOMERS_STATE_PARAMETERS` and the names
    #: of the compartments associated to each phytomer (see :attr:`MODEL_COMPARTMENTS_NAMES`)
    PHYTOMERS_STATE = PHYTOMERS_STATE_PARAMETERS + MODEL_COMPARTMENTS_NAMES.get(model.Phytomer, [])
    #: the variables that we need to compute in order to compute fluxes and/or compartments values at phytomer scale
    PHYTOMERS_INTERMEDIATE_VARIABLES = []
    #: the fluxes exchanged between the compartments at phytomer scale
    PHYTOMERS_FLUXES = []
    #: the variables computed by integrating values of phytomer components parameters/variables recursively
    PHYTOMERS_INTEGRATIVE_VARIABLES = []
    #: all the variables computed during a run step of the simulation at phytomer scale
    PHYTOMERS_RUN_VARIABLES = PHYTOMERS_STATE + PHYTOMERS_INTERMEDIATE_VARIABLES + PHYTOMERS_FLUXES + PHYTOMERS_INTEGRATIVE_VARIABLES

    # ---------- ORGAN scale ----------

    #: the indexes to locate the organs in the modeled system
    ORGANS_INDEXES = ['plant', 'axis', 'organ']
    #: concatenation of :attr:`T_INDEX` and :attr:`ORGANS_INDEXES`
    ORGANS_T_INDEXES = T_INDEX + ORGANS_INDEXES
    #: the parameters which define the state of the modeled system at organ scale
    ORGANS_STATE_PARAMETERS = ['mstruct', 'Nstruct', 'senesced_mstruct']
    #: the variables which define the state of the modeled system at organ scale,
    #: formed be the concatenation of :attr:`ORGANS_STATE_PARAMETERS` and the names
    #: of the compartments associated to each organ (see :attr:`MODEL_COMPARTMENTS_NAMES`)
    ORGANS_STATE = ORGANS_STATE_PARAMETERS + MODEL_COMPARTMENTS_NAMES.get(model.Organ, [])
    #: the variables that we need to compute in order to compute fluxes and/or compartments values at organ scale
    ORGANS_INTERMEDIATE_VARIABLES = ['C_exudation', 'HATS_LATS', 'N_exudation', 'RGR_Structure', 'R_Nnit_red', 'R_Nnit_upt', 'Respi_growth',
                                     'R_grain_growth_starch', 'R_grain_growth_struct', 'R_residual', 'regul_transpiration', 'sum_respi']
    #: the fluxes exchanged between the compartments at organ scale
    ORGANS_FLUXES = ['Export_Amino_Acids', 'Export_Nitrates', 'Export_cytokinins', 'S_Amino_Acids', 'S_cytokinins', 'S_grain_starch',
                     'S_grain_structure', 'S_Proteins', 'Unloading_Amino_Acids', 'Unloading_Sucrose', 'Uptake_Nitrates']
    #: the variables computed by integrating values of organ components parameters/variables recursively
    ORGANS_INTEGRATIVE_VARIABLES = ['Total_Organic_Nitrogen']
    #: all the variables computed during a run step of the simulation at organ scale
    ORGANS_RUN_VARIABLES = ORGANS_STATE + ORGANS_INTERMEDIATE_VARIABLES + ORGANS_FLUXES + ORGANS_INTEGRATIVE_VARIABLES

    # ---------- HIDDENZONE scale ----------

    #: the indexes to locate the hidden zones in the modeled system
    HIDDENZONE_INDEXES = ['plant', 'axis', 'metamer']
    #: concatenation of :attr:`T_INDEX` and :attr:`HIDDENZONE_INDEXES`
    HIDDENZONE_T_INDEXES = T_INDEX + HIDDENZONE_INDEXES
    #: the parameters which define the state of the modeled system at hidden zone scale
    HIDDENZONE_STATE_PARAMETERS = ['Nstruct', 'mstruct', 'ratio_DZ']
    #: the variables which define the state of the modeled system at hidden zone scale,
    #: formed be the concatenation of :attr:`HIDDENZONE_STATE_PARAMETERS` and the names
    #: of the compartments associated to each hidden zone (see :attr:`MODEL_COMPARTMENTS_NAMES`)
    HIDDENZONE_STATE = HIDDENZONE_STATE_PARAMETERS + MODEL_COMPARTMENTS_NAMES.get(model.HiddenZone, [])
    #: the variables that we need to compute in order to compute fluxes and/or compartments values at hidden zone scale
    HIDDENZONE_INTERMEDIATE_VARIABLES = ['nb_replications']
    #: the fluxes exchanged between the compartments at hidden zone scale
    HIDDENZONE_FLUXES = ['D_Fructan', 'D_Proteins', 'S_Fructan', 'S_Proteins', 'Unloading_Amino_Acids', 'Unloading_Sucrose']
    #: the variables computed by integrating values of hidden zone components parameters/variables recursively
    HIDDENZONE_INTEGRATIVE_VARIABLES = []
    #: all the variables computed during a run step of the simulation at plnat scale
    HIDDENZONE_RUN_VARIABLES = HIDDENZONE_STATE + HIDDENZONE_INTERMEDIATE_VARIABLES + HIDDENZONE_FLUXES + HIDDENZONE_INTEGRATIVE_VARIABLES

    # ---------- ELEMENT scale ----------

    #: the indexes to locate the elements in the modeled system
    ELEMENTS_INDEXES = ['plant', 'axis', 'metamer', 'organ', 'element']
    #: concatenation of :attr:`T_INDEX` and :attr:`ELEMENTS_INDEXES`
    ELEMENTS_T_INDEXES = T_INDEX + ELEMENTS_INDEXES
    #: the parameters which define the state of the modeled system at element scale
    ELEMENTS_STATE_PARAMETERS = ['Ag', 'Nstruct', 'Tr', 'Ts', 'green_area', 'is_growing', 'mstruct', 'senesced_mstruct']
    #: the variables which define the state of the modeled system at element scale,
    #: formed be the concatenation of :attr:`ELEMENTS_STATE_PARAMETERS` and the names
    #: of the compartments associated to each element (see :attr:`MODEL_COMPARTMENTS_NAMES`)
    ELEMENTS_STATE = ELEMENTS_STATE_PARAMETERS + MODEL_COMPARTMENTS_NAMES.get(model.PhotosyntheticOrganElement, [])
    #: the variables that we need to compute in order to compute fluxes and/or compartments values at element scale
    ELEMENTS_INTERMEDIATE_VARIABLES = ['Photosynthesis', 'R_Nnit_red', 'R_phloem_loading', 'R_residual', 'Transpiration', 'sum_respi', 'nb_replications']
    #: the fluxes exchanged between the compartments at element scale
    ELEMENTS_FLUXES = ['Amino_Acids_import', 'D_Fructan', 'D_Proteins', 'D_Starch', 'D_cytokinins', 'Loading_Amino_Acids', 'Loading_Sucrose',
                       'Nitrates_import', 'Regul_S_Fructan', 'S_Fructan', 'S_Starch', 'S_Sucrose', 'S_Amino_Acids', 'S_Proteins',
                       'cytokinins_import']
    #: the variables computed by integrating values of element components parameters/variables recursively
    ELEMENTS_INTEGRATIVE_VARIABLES = ['Total_Organic_Nitrogen']
    #: all the variables computed during a run step of the simulation at element scale
    ELEMENTS_RUN_VARIABLES = ELEMENTS_STATE + ELEMENTS_INTERMEDIATE_VARIABLES + ELEMENTS_FLUXES + ELEMENTS_INTEGRATIVE_VARIABLES

    # ---------- SOIL scale ----------

    #: the indexes to locate the soils in the modeled system
    SOILS_INDEXES = ['plant', 'axis']
    #: concatenation of :attr:`T_INDEX` and :attr:`SOILS_INDEXES`
    SOILS_T_INDEXES = T_INDEX + SOILS_INDEXES
    #: the parameters which define the state of the modeled system at soil scale
    SOILS_STATE_PARAMETERS = ['Tsoil', 'volume']
    #: the variables which define the state of the modeled system at soil scale,
    #: formed be the concatenation of :attr:`SOILS_STATE_PARAMETERS` and the names
    #: of the compartments associated to each soil (see :attr:`MODEL_COMPARTMENTS_NAMES`)
    SOILS_STATE = SOILS_STATE_PARAMETERS + MODEL_COMPARTMENTS_NAMES.get(model.Soil, [])
    #: the variables that we need to compute in order to compute fluxes and/or compartments values at soil scale
    SOILS_INTERMEDIATE_VARIABLES = ['Conc_Nitrates_Soil', 'mineralisation']
    #: the fluxes exchanged between the compartments at soil scale
    SOILS_FLUXES = []
    #: the variables computed by integrating values of soil components parameters/variables recursively
    SOILS_INTEGRATIVE_VARIABLES = []
    #: all the variables computed during a run step of the simulation at soil scale
    SOILS_RUN_VARIABLES = SOILS_STATE + SOILS_INTERMEDIATE_VARIABLES + SOILS_FLUXES + SOILS_INTEGRATIVE_VARIABLES

    #: a dictionary of all the variables which define the state of the modeled system, for each scale
    ALL_STATE_PARAMETERS = {model.Plant: PLANTS_STATE_PARAMETERS,
                            model.Axis: AXES_STATE_PARAMETERS,
                            model.Phytomer: PHYTOMERS_STATE_PARAMETERS,
                            model.Organ: ORGANS_STATE_PARAMETERS,
                            model.HiddenZone: HIDDENZONE_STATE_PARAMETERS,
                            model.PhotosyntheticOrganElement: ELEMENTS_STATE_PARAMETERS,
                            model.Soil: SOILS_STATE_PARAMETERS}

    #: the names of the roots (scenescence) forcings
    ROOTS_FORCINGS = ('Nstruct', 'mstruct')
    #: the names of the elements photosynthesis forcings
    ELEMENTS_PHOTOSYNTHESIS_FORCINGS = ('Ag', 'Tr', 'Ts')
    #: the names of the elements scenescence forcings
    ELEMENTS_SENESCENCE_FORCINGS = ('Nstruct', 'green_area', 'mstruct')
    #: the names of the elements photosynthesis and scenescence forcings
    ELEMENTS_FORCINGS = ELEMENTS_PHOTOSYNTHESIS_FORCINGS + ELEMENTS_SENESCENCE_FORCINGS

    #: the name of the loggers for compartments and derivatives
    LOGGERS_NAMES = {'compartments': {model.Plant: 'cnwheat.compartments.plants',
                                      model.Axis: 'cnwheat.compartments.axes',
                                      model.Phytomer: 'cnwheat.compartments.phytomers',
                                      model.Organ: 'cnwheat.compartments.organs',
                                      model.HiddenZone: 'cnwheat.compartments.hiddenzones',
                                      model.PhotosyntheticOrganElement: 'cnwheat.compartments.elements',
                                      model.Soil: 'cnwheat.compartments.soils'},
                     'derivatives': {model.Plant: 'cnwheat.derivatives.plants',
                                     model.Axis: 'cnwheat.derivatives.axes',
                                     model.Phytomer: 'cnwheat.derivatives.phytomers',
                                     model.Organ: 'cnwheat.derivatives.organs',
                                     model.HiddenZone: 'cnwheat.derivatives.hiddenzones',
                                     model.PhotosyntheticOrganElement: 'cnwheat.derivatives.elements',
                                     model.Soil: 'cnwheat.derivatives.soils'}}

    def __init__(self, respiration_model, delta_t=1, culm_density=None, interpolate_forcings=False, senescence_forcings_delta_t=None, photosynthesis_forcings_delta_t=None):

        self.respiration_model = respiration_model  #: the model of respiration to use

        self.population = model.Population()  #: the population to simulate on

        #: The inputs of the soils.
        #:
        #: `soils` is a dictionary of objects of type :class:`model.Soil`:
        #:     {(plant_index, axis_label): soil_object, ...}
        self.soils = {}

        self.initial_conditions = []  #: the initial conditions of the compartments in the population and soils
        self.initial_conditions_mapping = {}  #: dictionary to map the compartments to their indexes in :attr:`initial_conditions`

        self.progressbar = tools.ProgressBar(title='Solver progress')  #: progress bar to show the progress of the solver
        self.show_progressbar = False  #: True: show the progress bar ; False: DO NOT show the progress bar

        self.delta_t = delta_t  #: the delta t of the simulation (in seconds)

        self.time_step = self.delta_t / 3600.0  #: time step of the simulation (in hours)

        self.time_grid = np.array([0.0, self.time_step])  #: the time grid of the simulation (in hours)

        self.culm_density = culm_density  #: culm density (culm m-2)

        self.interpolate_forcings = interpolate_forcings  #: a boolean flag which indicates if we want to interpolate or not the forcings (True: interpolate, False: do not interpolate)

        # set the loggers for compartments and derivatives
        compartments_logger = logging.getLogger('cnwheat.compartments')
        derivatives_logger = logging.getLogger('cnwheat.derivatives')
        if compartments_logger.isEnabledFor(logging.DEBUG) or derivatives_logger.isEnabledFor(logging.DEBUG):
            sep = ','
            if compartments_logger.isEnabledFor(logging.DEBUG):
                plants_compartments_logger = logging.getLogger('cnwheat.compartments.plants')
                plants_compartments_logger.debug(sep.join(Simulation.PLANTS_T_INDEXES + Simulation.PLANTS_STATE))
                axes_compartments_logger = logging.getLogger('cnwheat.compartments.axes')
                axes_compartments_logger.debug(sep.join(Simulation.AXES_T_INDEXES + Simulation.AXES_STATE))
                phytomers_compartments_logger = logging.getLogger('cnwheat.compartments.phytomers')
                phytomers_compartments_logger.debug(sep.join(Simulation.PHYTOMERS_T_INDEXES + Simulation.PHYTOMERS_STATE))
                organs_compartments_logger = logging.getLogger('cnwheat.compartments.organs')
                organs_compartments_logger.debug(sep.join(Simulation.ORGANS_T_INDEXES + Simulation.ORGANS_STATE))
                hiddenzones_compartments_logger = logging.getLogger('cnwheat.compartments.hiddenzones')
                hiddenzones_compartments_logger.debug(sep.join(Simulation.HIDDENZONE_T_INDEXES + Simulation.HIDDENZONE_STATE))
                elements_compartments_logger = logging.getLogger('cnwheat.compartments.elements')
                elements_compartments_logger.debug(sep.join(Simulation.ELEMENTS_T_INDEXES + Simulation.ELEMENTS_STATE))
                soils_compartments_logger = logging.getLogger('cnwheat.compartments.soils')
                soils_compartments_logger.debug(sep.join(Simulation.SOILS_T_INDEXES + Simulation.SOILS_STATE))
            if derivatives_logger.isEnabledFor(logging.DEBUG):
                plants_derivatives_logger = logging.getLogger('cnwheat.derivatives.plants')
                plants_derivatives_logger.debug(sep.join(Simulation.PLANTS_T_INDEXES + Simulation.PLANTS_STATE))
                axes_derivatives_logger = logging.getLogger('cnwheat.derivatives.axes')
                axes_derivatives_logger.debug(sep.join(Simulation.AXES_T_INDEXES + Simulation.AXES_STATE))
                phytomers_derivatives_logger = logging.getLogger('cnwheat.derivatives.phytomers')
                phytomers_derivatives_logger.debug(sep.join(Simulation.PHYTOMERS_T_INDEXES + Simulation.PHYTOMERS_STATE))
                organs_derivatives_logger = logging.getLogger('cnwheat.derivatives.organs')
                organs_derivatives_logger.debug(sep.join(Simulation.ORGANS_T_INDEXES + Simulation.ORGANS_STATE))
                hiddenzones_derivatives_logger = logging.getLogger('cnwheat.derivatives.hiddenzones')
                hiddenzones_derivatives_logger.debug(sep.join(Simulation.HIDDENZONE_T_INDEXES + Simulation.HIDDENZONE_STATE))
                elements_derivatives_logger = logging.getLogger('cnwheat.derivatives.elements')
                elements_derivatives_logger.debug(sep.join(Simulation.ELEMENTS_T_INDEXES + Simulation.ELEMENTS_STATE))
                soils_derivatives_logger = logging.getLogger('cnwheat.derivatives.soils')
                soils_derivatives_logger.debug(sep.join(Simulation.SOILS_T_INDEXES + Simulation.SOILS_STATE))

        logger = logging.getLogger(__name__)
        if logger.isEnabledFor(logging.DEBUG):
            self.t_offset = 0.0  #: the absolute time offset elapsed from the beginning of the simulation

        if interpolate_forcings:
            if senescence_forcings_delta_t is not None and photosynthesis_forcings_delta_t is not None and \
                    senescence_forcings_delta_t >= delta_t and photosynthesis_forcings_delta_t >= delta_t:
                self.senescence_forcings_delta_t_ratio = senescence_forcings_delta_t / delta_t  #: the ratio between the delta t of the senescence forcings and the delta t of the simulation
                self.photosynthesis_forcings_delta_t_ratio = photosynthesis_forcings_delta_t / delta_t  #: the ratio between the delta t of the photosynthesis forcings and the delta t of the simulation
            elif senescence_forcings_delta_t is None:
                message = """The value of `interpolate_forcings` passed to the Simulation constructor is `True`, but `senescence_forcings_delta_t` is `None`. 
        Please set `senescence_forcings_delta_t` (through the Simulation constructor) to a not `None` value."""
                logger.exception(message)
                raise SimulationConstructionError(message)
            elif photosynthesis_forcings_delta_t is None:
                message = """The value of `interpolate_forcings` passed to the Simulation constructor is `True`, but `photosynthesis_forcings_delta_t` is `None`. 
        Please set `photosynthesis_forcings_delta_t` (through the Simulation constructor) to a not `None` value."""
                logger.exception(message)
                raise SimulationConstructionError(message)
            elif senescence_forcings_delta_t < delta_t:
                message = """The value of `senescence_forcings_delta_t` passed to the Simulation constructor is lesser than the one of `delta_t`. 
        Please set a `senescence_forcings_delta_t` that is at least equal to `delta_t`."""
                logger.exception(message)
                raise SimulationConstructionError(message)
            elif photosynthesis_forcings_delta_t < delta_t:
                message = """The value of `photosynthesis_forcings_delta_t` passed to the Simulation constructor is lesser than the one of `delta_t`. 
        Please set a `photosynthesis_forcings_delta_t` that is at least equal to `delta_t`."""
                logger.exception(message)
                raise SimulationConstructionError(message)

            self.previous_forcings_values = {}  #: previous values of the forcings
            self.new_forcings_values = {}  #: new values of the forcings
            self.interpolation_functions = {}  #: functions to interpolate the forcings

        self.nfev_total = 0  #: cumulative number of RHS function evaluations

    def initialize(self, population, soils, Tair=12, Tsoil=12):
        """
        Initialize:
            * :attr:`population`,
            * :attr:`soils`,
            * :attr:`initial_conditions_mapping`,
            * and :attr:`initial_conditions`

        from `population` and `soils`.

        :param model.Population population: a population of plants.
        :param dict soils: the soil associated to each axis. `soils` must be a dictionary with the same structure as :attr:`soils`
        :param float Tair: air temperature (°C)
        :param float Tsoil: soil temperature (°C)
        """

        logger = logging.getLogger(__name__)

        logger.info('Initialization of the simulation...')

        # clean the attributes of the simulation
        del self.population.plants[:]
        self.soils.clear()
        del self.initial_conditions[:]
        self.initial_conditions_mapping.clear()

        # create new population and soils
        self.population.plants.extend(population.plants)
        self.soils.update(soils)

        # check the consistency of population and soils
        if len(self.population.plants) != 0:  # population must contain at least 1 plant
            for plant in self.population.plants:
                if len(plant.axes) != 0:  # each plant must contain at least 1 axis
                    for axis in plant.axes:
                        if axis.roots is None:  # each axis must have a "roots"
                            message = 'No roots found in (plant={},axis={})'.format(plant.index, axis.label)
                            logger.exception(message)
                            raise SimulationInitializationError(message)
                        if axis.phloem is None:  # each axis must have a phloem
                            message = 'No phloem found in (plant={},axis={})'.format(plant.index, axis.label)
                            logger.exception(message)
                            raise SimulationInitializationError(message)
                        if len(axis.phytomers) != 0:  # each axis must contain at least 1 phytomer
                            for phytomer in axis.phytomers:
                                phytomer_organs = (phytomer.lamina, phytomer.internode, phytomer.sheath, phytomer.chaff, phytomer.peduncle)
                                # each phytomer must contain at least 1 photosynthetic organ or an hidden growing zone
                                if phytomer_organs.count(None) != len(phytomer_organs) or phytomer.hiddenzone is not None:
                                    for organ in phytomer_organs:
                                        if organ is not None:
                                            organ_elements = (organ.exposed_element, organ.enclosed_element)
                                            # each photosynthetic organ must contain at least 1 element
                                            if organ_elements.count(None) != len(organ_elements):
                                                for element in organ_elements:
                                                    if element is not None:
                                                        # an element must belong to an organ of the same type (e.g. a LaminaElement must belong to a Lamina)
                                                        if organ.__class__.__name__ not in element.__class__.__name__:
                                                            message = 'In (plant={},axis={},phytomer={}), a {} belongs to a {}'.format(plant.index,
                                                                                                                                       axis.label,
                                                                                                                                       phytomer.index,
                                                                                                                                       element.__class__.__name__,
                                                                                                                                       organ.__class__.__name__)
                                                            logger.exception(message)
                                                            raise SimulationInitializationError(message)
                                            else:
                                                message = 'No element found in (plant={},axis={},phytomer={},organ={})'.format(plant.index,
                                                                                                                               axis.label,
                                                                                                                               phytomer.index,
                                                                                                                               organ.label)
                                                logger.exception(message)
                                                raise SimulationInitializationError(message)
                                else:
                                    message = 'Neither photosynthetic organ nor hidden growing zone found in (plant={},axis={},phytomer={})'.format(plant.index,
                                                                                                                                                    axis.label,
                                                                                                                                                    phytomer.index)
                                    logger.exception(message)
                                    raise SimulationInitializationError(message)
                        else:
                            message = 'No phytomer found in (plant={},axis={})'.format(plant.index,
                                                                                       axis.label)
                            logger.exception(message)
                            raise SimulationInitializationError(message)
                        if (plant.index, axis.label) not in self.soils:  # each axis must be associated to a soil
                            message = 'No soil found in (plant={},axis={})'.format(plant.index,
                                                                                   axis.label)
                            logger.exception(message)
                            raise SimulationInitializationError(message)
                else:
                    message = 'No axis found in (plant={})'.format(plant.index)
                    logger.exception(message)
                    raise SimulationInitializationError(message)
        else:
            message = 'No plant found in the population.'
            logger.exception(message)
            raise SimulationInitializationError(message)

        if self.interpolate_forcings:
            # Save the new value of each forcing and set the state parameters to the previous forcing values.
            self.new_forcings_values.clear()
            for plant in self.population.plants:
                for axis in plant.axes:
                    if axis.roots is not None:
                        roots_id = (plant.index, axis.label)
                        self.new_forcings_values[roots_id] = {}
                        for forcing_label in Simulation.ROOTS_FORCINGS:
                            self.new_forcings_values[roots_id][forcing_label] = getattr(axis.roots, forcing_label)
                            if roots_id in self.previous_forcings_values:
                                setattr(axis.roots, forcing_label, self.previous_forcings_values[roots_id][forcing_label])
                    for phytomer in axis.phytomers:
                        for organ in (phytomer.lamina, phytomer.sheath):
                            if organ is None:
                                continue
                            for element in (organ.exposed_element, organ.enclosed_element):
                                if element is not None:
                                    element_id = (plant.index, axis.label, phytomer.index, organ.label, element.label)
                                    self.new_forcings_values[element_id] = {}
                                    for forcing_label in Simulation.ELEMENTS_FORCINGS:
                                        self.new_forcings_values[element_id][forcing_label] = getattr(element, forcing_label)
                                        if element_id in self.previous_forcings_values:
                                            setattr(element, forcing_label, self.previous_forcings_values[element_id][forcing_label])

        # Update soil and air temperature using weather data
        for soil_id, soil_inputs in self.soils.items():
            self.soils[soil_id].Tsoil = Tsoil
        for plant in self.population.plants:
            plant.Tair = Tair

        # initialize initial conditions
        def _init_initial_conditions(model_object, index):
            class_ = model_object.__class__
            if issubclass(class_, model.HiddenZone):
                class_ = model.HiddenZone
            elif issubclass(class_, model.Organ):
                class_ = model.Organ
            elif issubclass(class_, model.PhotosyntheticOrganElement):
                class_ = model.PhotosyntheticOrganElement
            compartments_names = Simulation.MODEL_COMPARTMENTS_NAMES[class_]
            self.initial_conditions_mapping[model_object] = {}
            for compartment_name in compartments_names:
                if hasattr(model_object, compartment_name):
                    self.initial_conditions_mapping[model_object][compartment_name] = index
                    self.initial_conditions.append(0)
                    index += 1
            return index

        i = 0

        for soil in soils.values():
            i = _init_initial_conditions(soil, i)

        for plant in self.population.plants:
            i = _init_initial_conditions(plant, i)
            for axis in plant.axes:
                i = _init_initial_conditions(axis, i)
                for organ in (axis.roots, axis.phloem, axis.grains):
                    if organ is None:
                        continue
                    i = _init_initial_conditions(organ, i)
                for phytomer in axis.phytomers:
                    i = _init_initial_conditions(phytomer, i)
                    for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath, phytomer.hiddenzone):
                        if organ is None:
                            continue
                        i = _init_initial_conditions(organ, i)
                        if organ is phytomer.hiddenzone:
                            continue
                        for element in (organ.exposed_element, organ.enclosed_element):
                            if element is None:
                                continue
                            i = _init_initial_conditions(element, i)

        self.population.calculate_aggregated_variables()

        logger.info('Initialization of the simulation DONE')

    def run(self, show_progressbar=False):
        """
        Compute CN exchanges which occurred in :attr:`population` and :attr:`soils` over :attr:`delta_t`.

        :param bool show_progressbar: True: show the progress bar of the solver ; False: do not show the progress bar (default).
        """
        logger = logging.getLogger(__name__)
        logger.info('Run of CN-Wheat...')

        if self.interpolate_forcings:
            # interpolate the forcings
            self._interpolate_forcings()

        # set the progress-bar
        self.show_progressbar = show_progressbar
        if self.show_progressbar:
            self.progressbar.set_t_max(self.time_step)

        self._update_initial_conditions()

        if logger.isEnabledFor(logging.DEBUG):
            logger.debug("Run the solver with delta_t = %s", self.time_step)

        # call :func:`scipy.integrate.solve_ivp` to integrate the system during 1 time step ;
        # :func:`scipy.integrate.solve_ivp` computes the derivatives of each function by calling :meth:`_calculate_all_derivatives`
        sol = solve_ivp(fun=self._calculate_all_derivatives, t_span=self.time_grid, y0=self.initial_conditions,
                        method='BDF', t_eval=np.array([self.time_step]), dense_output=False)

        self.nfev_total += sol.nfev

        if logger.isEnabledFor(logging.DEBUG):
            logger.debug("Run of the solver DONE")

        # check the integration ; raise an exception if the integration failed
        if not sol.success:
            message = "Integration failed: {}".format(sol.message)
            logger.exception(message)
            raise SimulationRunError(message)

        # Re-compute integrative variables
        self.population.calculate_aggregated_variables()

        if logger.isEnabledFor(logging.DEBUG):
            self.t_offset += self.time_step

        logger.info('Run of CN-Wheat DONE')

    def _update_initial_conditions(self):
        """Update the compartments values in :attr:`initial_conditions` from the compartments values of :attr:`population` and :attr:`soils`.
        """
        # Update the compartments values
        for model_object, compartments in self.initial_conditions_mapping.items():
            for compartment_name, compartment_index in compartments.items():
                self.initial_conditions[compartment_index] = getattr(model_object, compartment_name)

    def _interpolate_forcings(self):
        """Create functions to interpolate the forcings of the model to any time inside the time grid (see `self.time_grid`).

        If this is the first run of the model, then we consider that the forcings are constant.
        The interpolation functions are stored in :attr:`interpolation_functions`, and will be used later on and as needed by the SciPy solver.
        """
        self.interpolation_functions.clear()
        next_forcings_values = {}
        for plant in self.population.plants:
            for axis in plant.axes:
                if axis.roots is not None:
                    roots_id = (plant.index, axis.label)
                    self.interpolation_functions[roots_id] = {}
                    next_forcings_values[roots_id] = {}
                    for forcing_label in Simulation.ROOTS_FORCINGS:
                        if roots_id in self.previous_forcings_values and \
                                self.previous_forcings_values[roots_id][forcing_label] != self.new_forcings_values[roots_id][forcing_label]:
                            prev_forcing_value = self.previous_forcings_values[roots_id][forcing_label]
                            next_forcing_value = prev_forcing_value + (self.new_forcings_values[roots_id][forcing_label] - prev_forcing_value) / self.senescence_forcings_delta_t_ratio
                        else:
                            next_forcing_value = self.new_forcings_values[roots_id][forcing_label]
                            prev_forcing_value = next_forcing_value
                        self.interpolation_functions[roots_id][forcing_label] = interpolate.interp1d(self.time_grid, [prev_forcing_value, next_forcing_value], assume_sorted=True)
                        next_forcings_values[roots_id][forcing_label] = next_forcing_value
                for phytomer in axis.phytomers:
                    for organ in (phytomer.lamina, phytomer.sheath):
                        if organ is None:
                            continue
                        for element in (organ.exposed_element, organ.enclosed_element):
                            if element is not None:
                                element_id = (plant.index, axis.label, phytomer.index, organ.label, element.label)
                                self.interpolation_functions[element_id] = {}
                                next_forcings_values[element_id] = {}
                                for (forcing_labels, forcings_delta_t_ratio) in ((Simulation.ELEMENTS_PHOTOSYNTHESIS_FORCINGS, self.photosynthesis_forcings_delta_t_ratio),
                                                                                 (Simulation.ELEMENTS_SENESCENCE_FORCINGS, self.senescence_forcings_delta_t_ratio)):
                                    for forcing_label in forcing_labels:
                                        if element_id in self.previous_forcings_values and \
                                                self.previous_forcings_values[element_id][forcing_label] != self.new_forcings_values[element_id][forcing_label]:
                                            prev_forcing_value = self.previous_forcings_values[element_id][forcing_label]
                                            next_forcing_value = prev_forcing_value + (self.new_forcings_values[element_id][forcing_label] - prev_forcing_value) / forcings_delta_t_ratio
                                        else:
                                            next_forcing_value = self.new_forcings_values[element_id][forcing_label]
                                            prev_forcing_value = next_forcing_value
                                        self.interpolation_functions[element_id][forcing_label] = interpolate.interp1d(self.time_grid, [prev_forcing_value, next_forcing_value], assume_sorted=True)
                                        next_forcings_values[element_id][forcing_label] = next_forcing_value

        self.previous_forcings_values.clear()
        self.previous_forcings_values.update(next_forcings_values)

    def _log_compartments(self, t, y, loggers_names):
        """Log the values in `y` to the loggers in `loggers_names`.
        """

        def update_rows(model_object, indexes, rows, i):
            """Update list `rows` appending a new row corresponding to the compartment
            values associated to object `model_object` located at indexes `indexes`.
            `i` is used to reach the values associated to object `model_object`
            from array `y`.
            """
            row = []
            class_ = model_object.__class__
            if issubclass(class_, model.HiddenZone):
                class_ = model.HiddenZone
            elif issubclass(class_, model.Organ):
                class_ = model.Organ
            elif issubclass(class_, model.PhotosyntheticOrganElement):
                class_ = model.PhotosyntheticOrganElement
            parameters_names = Simulation.ALL_STATE_PARAMETERS[class_]
            for parameter_name in parameters_names:
                if hasattr(model_object, parameter_name):
                    row.append(str(getattr(model_object, parameter_name)))
                else:
                    row.append('NA')
            compartments_names = Simulation.MODEL_COMPARTMENTS_NAMES[class_]
            for compartment_name in compartments_names:
                if hasattr(model_object, compartment_name):
                    row.append(str(y[i]))
                    i += 1
                else:
                    row.append('NA')
            rows.append([str(index) for index in indexes] + row)
            return i

        i = 0
        all_rows = dict([(class_, []) for class_ in loggers_names])

        for soil_id, soil in self.soils.items():
            i = update_rows(soil, (t,) + soil_id, all_rows[model.Soil], i)

        for plant in self.population.plants:
            i = update_rows(plant, [t, plant.index], all_rows[model.Plant], i)
            for axis in plant.axes:
                i = update_rows(axis, [t, plant.index, axis.label], all_rows[model.Axis], i)
                for organ in (axis.roots, axis.phloem, axis.grains):
                    if organ is None:
                        continue
                    i = update_rows(organ, [t, plant.index, axis.label, organ.label], all_rows[model.Organ], i)
                for phytomer in axis.phytomers:
                    i = update_rows(phytomer, [t, plant.index, axis.label, phytomer.index], all_rows[model.Phytomer], i)
                    for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath, phytomer.hiddenzone):
                        if organ is None:
                            continue
                        if organ is phytomer.hiddenzone:
                            i = update_rows(organ, [t, plant.index, axis.label, phytomer.index], all_rows[model.HiddenZone], i)
                            continue
                        for element in (organ.exposed_element, organ.enclosed_element):
                            if element is None:
                                continue
                            i = update_rows(element, [t, plant.index, axis.label, phytomer.index, organ.label, element.label], all_rows[model.PhotosyntheticOrganElement], i)

        row_sep = '\n'
        column_sep = ','
        for class_, logger_name in loggers_names.items():
            compartments_logger = logging.getLogger(logger_name)
            formatted_initial_conditions = row_sep.join([column_sep.join(row) for row in all_rows[class_]])
            compartments_logger.debug(formatted_initial_conditions)

    def _calculate_all_derivatives(self, t, y):
        """Compute the derivative of `y` at `t`.

        :meth:`_calculate_all_derivatives` is passed as **func** argument to
        :func:`solve_ivp(fun, t_span, y0,...) <scipy.integrate.solve_ivp>`.
        :meth:`_calculate_all_derivatives` is called automatically by
        :func:`scipy.integrate.solve_ivp <scipy.integrate.solve_ivp>`.

        First call to :meth:`_calculate_all_derivatives` uses `y` = **y0** and
        `t` = **t_span** [0], where **y0** and **t_span** are arguments passed to :func:`solve_ivp(fun, t_span, y0,...) <scipy.integrate.solve_ivp>`.

        Following calls to :meth:`_calculate_all_derivatives` use `t` in [**t_span** [0], **t_span** [1]]. 

        :param float t: The current t at which we want to compute the derivatives.
              Values of `t` are chosen automatically by :func:`scipy.integrate.solve_ivp`.
              At first call to :meth:`_calculate_all_derivatives` by :func:`scipy.integrate.solve_ivp`,
              `t` = **t_span** [0], where **t_span** is one of the arguments passed to :func:`solve_ivp(fun, t_span, y0,...) <scipy.integrate.solve_ivp>`.
              For each following call to :meth:`_calculate_all_derivatives`, `t` belongs
              to the interval [**t_span** [0], **t_span** [1]].
        :param list [float] y: The current values of y.
              At first call to :meth:`_calculate_all_derivatives` by :func:`scipy.integrate.solve_ivp`, `y` = **y0**
              where **y0** is one of the arguments passed to :func:`solve_ivp(fun, t_span, y0,...) <scipy.integrate.solve_ivp>`.
              Then, values of `y` are chosen automatically by :func:`scipy.integrate.solve_ivp`.

        :return: The derivatives of `y` at `t`.
        :rtype: list [float]
        """
        logger = logging.getLogger(__name__)

        if logger.isEnabledFor(logging.DEBUG):
            t_abs = t + self.t_offset
            logger.debug('t = {}'.format(t_abs))

        if self.interpolate_forcings:
            # Update state parameters using interpolation functions
            for plant in self.population.plants:
                for axis in plant.axes:
                    if axis.roots is not None:
                        roots_id = (plant.index, axis.label)
                        for forcing_label in Simulation.ROOTS_FORCINGS:
                            setattr(axis.roots, forcing_label, float(self.interpolation_functions[roots_id][forcing_label](t)))
                    for phytomer in axis.phytomers:
                        for organ in (phytomer.lamina, phytomer.sheath):
                            if organ is None:
                                continue
                            for element in (organ.exposed_element, organ.enclosed_element):
                                if element is not None:
                                    element_id = (plant.index, axis.label, phytomer.index, organ.label, element.label)
                                    for forcing_label in Simulation.ELEMENTS_FORCINGS:
                                        setattr(element, forcing_label, float(self.interpolation_functions[element_id][forcing_label](t)))

            # Compute integrative variables
            self.population.calculate_aggregated_variables()

        compartments_logger = logging.getLogger('cnwheat.compartments')
        if logger.isEnabledFor(logging.DEBUG) and compartments_logger.isEnabledFor(logging.DEBUG):
            self._log_compartments(t_abs, y, Simulation.LOGGERS_NAMES['compartments'])

        # check that the solver is not crashed
        y_isnan = np.isnan(y)
        if y_isnan.any():
            message = 'The solver did not manage to compute a compartment. See the logs. NaN found in y'
            logger.exception(message)
            raise SimulationRunError(message)

        y_derivatives = np.zeros_like(y)

        # TODO: TEMP !!!!
        soil_contributors = []
        soil = self.soils[(1, 'MS')]
        soil.nitrates = y[self.initial_conditions_mapping[soil]['nitrates']]
        soil.Conc_Nitrates_Soil = soil.calculate_Conc_Nitrates(soil.nitrates)

        soil.T_effect_Vmax = soil.calculate_temperature_effect_on_Vmax(soil.Tsoil)
        soil.T_effect_conductivity = soil.calculate_temperature_effect_on_conductivity(soil.Tsoil)

        for plant in self.population.plants:

            plant.T_effect_conductivity = plant.calculate_temperature_effect_on_conductivity(plant.Tair)
            plant.T_effect_Vmax = plant.calculate_temperature_effect_on_Vmax(plant.Tair)

            for axis in plant.axes:
                sum_respi_shoot = 0.0
                axis.C_exudated = y[self.initial_conditions_mapping[axis]['C_exudated']]
                axis.sum_respi_shoot = y[self.initial_conditions_mapping[axis]['sum_respi_shoot']]
                axis.sum_respi_roots = y[self.initial_conditions_mapping[axis]['sum_respi_roots']]

                # Phloem
                phloem_contributors = []
                axis.phloem.sucrose = y[self.initial_conditions_mapping[axis.phloem]['sucrose']]
                axis.phloem.amino_acids = y[self.initial_conditions_mapping[axis.phloem]['amino_acids']]
                # Roots
                axis.roots.nitrates = y[self.initial_conditions_mapping[axis.roots]['nitrates']]
                axis.roots.amino_acids = y[self.initial_conditions_mapping[axis.roots]['amino_acids']]
                axis.roots.sucrose = y[self.initial_conditions_mapping[axis.roots]['sucrose']]
                axis.roots.cytokinins = y[self.initial_conditions_mapping[axis.roots]['cytokinins']]
                phloem_contributors.append(axis.roots)

                # compute total transpiration at t_inf
                axis.Total_Transpiration = 0.0  # mmol s-1
                total_green_area = 0.0  # m2
                for phytomer in axis.phytomers:
                    for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):
                        if organ is not None:
                            for element in (organ.exposed_element, organ.enclosed_element):
                                if element is not None and element.green_area > 0:
                                    element.Transpiration = element.calculate_Total_Transpiration(element.Tr, element.green_area)
                                    axis.Total_Transpiration += (element.Transpiration * element.nb_replications)
                                    total_green_area += (element.green_area * element.nb_replications)

                # Compute the regulating factor of root exports by shoot transpiration
                axis.roots.regul_transpiration = axis.roots.calculate_regul_transpiration(axis.Total_Transpiration)

                # compute the flows from/to the roots to/from photosynthetic organs
                axis.roots.Uptake_Nitrates, axis.roots.HATS_LATS = axis.roots.calculate_Uptake_Nitrates(soil.Conc_Nitrates_Soil, axis.roots.nitrates, axis.roots.sucrose,
                                                                                                        soil.T_effect_Vmax)
                soil_contributors.append((axis.roots.Uptake_Nitrates, plant.index))  #: TODO TEMP!!!
                axis.roots.R_Nnit_upt = self.respiration_model.RespirationModel.R_Nnit_upt(axis.roots.Uptake_Nitrates, axis.roots.sucrose)
                axis.roots.Export_Nitrates = axis.roots.calculate_Export_Nitrates(axis.roots.nitrates, axis.roots.regul_transpiration)
                axis.roots.Export_Amino_Acids = axis.roots.calculate_Export_Amino_Acids(axis.roots.amino_acids, axis.roots.regul_transpiration)
                axis.roots.Export_cytokinins = axis.roots.calculate_Export_cytokinins(axis.roots.cytokinins, axis.roots.regul_transpiration)

                # compute the derivative of each photosynthetic organ element compartment
                for phytomer in axis.phytomers:
                    # Hidden zone
                    hiddenzone = phytomer.hiddenzone
                    if phytomer.hiddenzone is not None:
                        if hiddenzone.mstruct == 0:
                            continue
                        hiddenzone.sucrose = y[self.initial_conditions_mapping[hiddenzone]['sucrose']]
                        hiddenzone.fructan = y[self.initial_conditions_mapping[hiddenzone]['fructan']]
                        hiddenzone.amino_acids = y[self.initial_conditions_mapping[hiddenzone]['amino_acids']]
                        hiddenzone.proteins = y[self.initial_conditions_mapping[hiddenzone]['proteins']]
                        phloem_contributors.append(hiddenzone)

                    hiddenzone_Loading_Sucrose_contribution = 0
                    hiddenzone_Loading_Amino_Acids_contribution = 0
                    for organ in (phytomer.chaff, phytomer.peduncle, phytomer.lamina, phytomer.internode, phytomer.sheath):
                        if organ is None:
                            continue

                        for element in (organ.exposed_element, organ.enclosed_element):
                            if element is None or element.green_area <= 0.25E-6 or element.mstruct <= 0.0:
                                continue

                            element.starch = y[self.initial_conditions_mapping[element]['starch']]
                            element.sucrose = y[self.initial_conditions_mapping[element]['sucrose']]
                            element.triosesP = y[self.initial_conditions_mapping[element]['triosesP']]
                            element.fructan = y[self.initial_conditions_mapping[element]['fructan']]
                            element.nitrates = y[self.initial_conditions_mapping[element]['nitrates']]
                            element.amino_acids = y[self.initial_conditions_mapping[element]['amino_acids']]
                            element.proteins = y[self.initial_conditions_mapping[element]['proteins']]
                            element.cytokinins = y[self.initial_conditions_mapping[element]['cytokinins']]

                            # intermediate variables
                            element.Photosynthesis = element.calculate_total_Photosynthesis(element.Ag, element.green_area)

                            # flows
                            if element.is_growing:  #: Export of sucrose and amino acids towards the HZ. Several growing elements might export toward the HZ at the same time (leaf and internode)
                                element.Loading_Sucrose = element.calculate_export_sucrose(element.sucrose, hiddenzone.sucrose, hiddenzone.mstruct, plant.T_effect_conductivity)
                                hiddenzone_Loading_Sucrose_contribution += element.Loading_Sucrose
                                element.Loading_Amino_Acids = element.calculate_Export_Amino_Acids(element.amino_acids, hiddenzone.amino_acids, hiddenzone.mstruct, plant.T_effect_conductivity)
                                hiddenzone_Loading_Amino_Acids_contribution += element.Loading_Amino_Acids

                            else:  #: Loading of sucrose and amino acids towards the phloem
                                phloem_contributors.append(element)
                                element.Loading_Sucrose = element.calculate_Loading_Sucrose(element.sucrose, axis.phloem.sucrose, axis.mstruct, plant.T_effect_conductivity)
                                element.Loading_Amino_Acids = element.calculate_Loading_Amino_Acids(element.amino_acids, axis.phloem.amino_acids, axis.mstruct, plant.T_effect_conductivity)

                            element.Regul_S_Fructan = element.calculate_Regul_S_Fructan(element.Loading_Sucrose)
                            element.S_Fructan = element.calculate_S_Fructan(element.sucrose, element.Regul_S_Fructan, plant.T_effect_Vmax)
                            element.D_Fructan = element.calculate_D_Fructan(element.sucrose, element.fructan, plant.T_effect_Vmax)
                            element.S_Starch = element.calculate_S_Starch(element.triosesP, plant.T_effect_Vmax)
                            element.D_Starch = element.calculate_D_Starch(element.starch, plant.T_effect_Vmax)
                            element.S_Sucrose = element.calculate_S_Sucrose(element.triosesP, plant.T_effect_Vmax)
                            element.R_phloem_loading, element.Loading_Sucrose = self.respiration_model.RespirationModel.R_phloem(element.Loading_Sucrose,
                                                                                                                                 element.mstruct * element.__class__.PARAMETERS.ALPHA)
                            element.Nitrates_import = element.calculate_Nitrates_import(axis.roots.Export_Nitrates, element.Transpiration, axis.Total_Transpiration)
                            element.Amino_Acids_import = element.calculate_Amino_Acids_import(axis.roots.Export_Amino_Acids, element.Transpiration, axis.Total_Transpiration)
                            element.S_Amino_Acids = element.calculate_S_amino_acids(element.nitrates, element.triosesP, plant.T_effect_Vmax)
                            element.R_Nnit_red, element.S_Amino_Acids = self.respiration_model.RespirationModel.R_Nnit_red(element.S_Amino_Acids, element.sucrose,
                                                                                                                           element.mstruct * element.__class__.PARAMETERS.ALPHA)
                            element.S_Proteins = element.calculate_S_proteins(element.amino_acids, plant.T_effect_Vmax)
                            element.D_Proteins = element.calculate_D_Proteins(element.proteins, element.cytokinins, plant.T_effect_Vmax)
                            element.cytokinins_import = element.calculate_cytokinins_import(axis.roots.Export_cytokinins, element.Transpiration, axis.Total_Transpiration)
                            element.D_cytokinins = element.calculate_D_cytokinins(element.cytokinins, plant.T_effect_Vmax)

                            # compartments derivatives
                            starch_derivative = element.calculate_starch_derivative(element.S_Starch, element.D_Starch)
                            element.R_residual = self.respiration_model.RespirationModel.R_residual(element.sucrose, element.mstruct * element.__class__.PARAMETERS.ALPHA,
                                                                                                                           element.Total_Organic_Nitrogen, element.Ts)
                            element.sum_respi = element.R_phloem_loading + element.R_Nnit_red + element.R_residual
                            sum_respi_shoot += element.sum_respi * element.nb_replications
                            sucrose_derivative = element.calculate_sucrose_derivative(element.S_Sucrose, element.D_Starch, element.Loading_Sucrose, element.S_Fructan,
                                                                                      element.D_Fructan, element.sum_respi)
                            triosesP_derivative = element.calculate_triosesP_derivative(element.Photosynthesis, element.S_Sucrose, element.S_Starch, element.S_Amino_Acids)
                            fructan_derivative = element.calculate_fructan_derivative(element.S_Fructan, element.D_Fructan)
                            nitrates_derivative = element.calculate_nitrates_derivative(element.Nitrates_import, element.S_Amino_Acids)
                            amino_acids_derivative = element.calculate_amino_acids_derivative(element.Amino_Acids_import, element.S_Amino_Acids, element.S_Proteins, element.D_Proteins,
                                                                                              element.Loading_Amino_Acids)
                            proteins_derivative = element.calculate_proteins_derivative(element.S_Proteins, element.D_Proteins)
                            cytokinins_derivative = element.calculate_cytokinins_derivative(element.cytokinins_import, element.D_cytokinins)

                            y_derivatives[self.initial_conditions_mapping[element]['starch']] = starch_derivative
                            y_derivatives[self.initial_conditions_mapping[element]['sucrose']] = sucrose_derivative
                            y_derivatives[self.initial_conditions_mapping[element]['triosesP']] = triosesP_derivative
                            y_derivatives[self.initial_conditions_mapping[element]['fructan']] = fructan_derivative
                            y_derivatives[self.initial_conditions_mapping[element]['nitrates']] = nitrates_derivative
                            y_derivatives[self.initial_conditions_mapping[element]['amino_acids']] = amino_acids_derivative
                            y_derivatives[self.initial_conditions_mapping[element]['proteins']] = proteins_derivative
                            y_derivatives[self.initial_conditions_mapping[element]['cytokinins']] = cytokinins_derivative

                    if phytomer.hiddenzone is not None:
                        # Unloading of sucrose from phloem
                        hiddenzone.Unloading_Sucrose = hiddenzone.calculate_Unloading_Sucrose(hiddenzone.sucrose, axis.phloem.sucrose, axis.mstruct, plant.T_effect_conductivity)

                        # Unloading of AA from phloem
                        hiddenzone.Unloading_Amino_Acids = hiddenzone.calculate_Unloading_Amino_Acids(hiddenzone.amino_acids, axis.phloem.amino_acids, axis.mstruct, plant.T_effect_conductivity)

                        # Fructan synthesis
                        Regul_Sfructanes = hiddenzone.calculate_Regul_S_Fructan(hiddenzone.Unloading_Sucrose)
                        hiddenzone.S_Fructan = hiddenzone.calculate_S_Fructan(hiddenzone.sucrose, Regul_Sfructanes, plant.T_effect_Vmax)

                        # Fructan degradation
                        hiddenzone.D_Fructan = hiddenzone.calculate_D_Fructan(hiddenzone.sucrose, hiddenzone.fructan, plant.T_effect_Vmax)

                        # Synthesis proteins
                        hiddenzone.S_Proteins = hiddenzone.calculate_S_proteins(hiddenzone.amino_acids, plant.T_effect_Vmax)

                        # Degradation proteins
                        hiddenzone.D_Proteins = hiddenzone.calculate_D_Proteins(hiddenzone.proteins, plant.T_effect_Vmax)

                        # Residual respiration
                        hiddenzone.R_residual = self.respiration_model.RespirationModel.R_residual(hiddenzone.sucrose,
                                                                                                                             hiddenzone.mstruct * hiddenzone.__class__.PARAMETERS.ALPHA,
                                                                                                                             hiddenzone.Total_Organic_Nitrogen,
                                                                                                                             plant.Tair)
                        sum_respi_shoot += hiddenzone.R_residual * hiddenzone.nb_replications

                        # compute the derivatives of the hidden zone
                        y_derivatives[self.initial_conditions_mapping[hiddenzone]['sucrose']] = hiddenzone.calculate_sucrose_derivative(hiddenzone.Unloading_Sucrose, hiddenzone.S_Fructan,
                                                                                                                                        hiddenzone.D_Fructan,
                                                                                                                                        hiddenzone_Loading_Sucrose_contribution, hiddenzone.R_residual)
                        y_derivatives[self.initial_conditions_mapping[hiddenzone]['amino_acids']] = hiddenzone.calculate_amino_acids_derivative(hiddenzone.Unloading_Amino_Acids, hiddenzone.S_Proteins,
                                                                                                                                                hiddenzone.D_Proteins,
                                                                                                                                                hiddenzone_Loading_Amino_Acids_contribution)
                        y_derivatives[self.initial_conditions_mapping[hiddenzone]['fructan']] = hiddenzone.calculate_fructan_derivative(hiddenzone.S_Fructan, hiddenzone.D_Fructan)
                        y_derivatives[self.initial_conditions_mapping[hiddenzone]['proteins']] = hiddenzone.calculate_proteins_derivative(hiddenzone.S_Proteins, hiddenzone.D_Proteins)

                if axis.grains is not None:
                    phloem_contributors.append(axis.grains)
                    # compute the derivative of each compartment of grains
                    axis.grains.structure = y[self.initial_conditions_mapping[axis.grains]['structure']]
                    axis.grains.starch = y[self.initial_conditions_mapping[axis.grains]['starch']]
                    axis.grains.proteins = y[self.initial_conditions_mapping[axis.grains]['proteins']]
                    axis.grains.age_from_flowering = y[self.initial_conditions_mapping[axis.grains]['age_from_flowering']]

                    # intermediate variables
                    T_effect_growth = axis.grains.calculate_temperature_effect_on_growth(plant.Tair)
                    axis.grains.RGR_Structure = axis.grains.calculate_RGR_Structure(axis.phloem.sucrose, axis.mstruct, T_effect_growth)
                    axis.grains.structural_dry_mass = axis.grains.calculate_structural_dry_mass(axis.grains.structure)

                    # flows
                    axis.grains.S_grain_structure = axis.grains.calculate_S_grain_structure(axis.grains.structure, axis.grains.RGR_Structure)
                    axis.grains.S_grain_starch = axis.grains.calculate_S_grain_starch(axis.phloem.sucrose, axis.mstruct, plant.T_effect_Vmax)
                    axis.grains.S_Proteins = axis.grains.calculate_S_proteins(axis.grains.S_grain_structure, axis.grains.S_grain_starch, axis.phloem.amino_acids, axis.phloem.sucrose,
                                                                              axis.grains.structural_dry_mass)
                    # compartments derivatives
                    axis.grains.R_grain_growth_struct, axis.grains.R_grain_growth_starch = self.respiration_model.RespirationModel.R_grain_growth(axis.grains.S_grain_structure,
                                                                                                                                                  axis.grains.S_grain_starch,
                                                                                                                                                  axis.grains.structural_dry_mass)
                    sum_respi_shoot += axis.grains.R_grain_growth_struct + axis.grains.R_grain_growth_starch
                    structure_derivative = axis.grains.calculate_structure_derivative(axis.grains.S_grain_structure, axis.grains.R_grain_growth_struct)
                    starch_derivative = axis.grains.calculate_starch_derivative(axis.grains.S_grain_starch, axis.grains.structural_dry_mass, axis.grains.R_grain_growth_starch)
                    proteins_derivative = axis.grains.calculate_proteins_derivative(axis.grains.S_Proteins)
                    y_derivatives[self.initial_conditions_mapping[axis.grains]['structure']] = structure_derivative
                    y_derivatives[self.initial_conditions_mapping[axis.grains]['starch']] = starch_derivative
                    y_derivatives[self.initial_conditions_mapping[axis.grains]['proteins']] = proteins_derivative
                    y_derivatives[self.initial_conditions_mapping[axis.grains]['age_from_flowering']] += (self.delta_t * T_effect_growth) #TODO: create a function

                # compute the derivative of each compartment of roots
                # flows
                axis.roots.Unloading_Sucrose = axis.roots.calculate_Unloading_Sucrose(axis.roots.sucrose, axis.phloem.sucrose, axis.mstruct, plant.T_effect_conductivity)
                axis.roots.Unloading_Amino_Acids = axis.roots.calculate_Unloading_Amino_Acids(axis.roots.Unloading_Sucrose, axis.phloem.sucrose, axis.phloem.amino_acids)
                axis.roots.S_Amino_Acids = axis.roots.calculate_S_amino_acids(axis.roots.nitrates, axis.roots.sucrose, soil.T_effect_Vmax)
                axis.roots.R_Nnit_red, axis.roots.S_Amino_Acids = self.respiration_model.RespirationModel.R_Nnit_red(axis.roots.S_Amino_Acids, axis.roots.sucrose,
                                                                                                                     axis.roots.mstruct * model.Roots.PARAMETERS.ALPHA, root=True)
                axis.roots.C_exudation, axis.roots.N_exudation = axis.roots.calculate_exudation(axis.roots.Unloading_Sucrose, axis.roots.sucrose, axis.roots.amino_acids, axis.phloem.amino_acids)
                axis.roots.S_cytokinins = axis.roots.calculate_S_cytokinins(axis.roots.sucrose, axis.roots.nitrates, soil.T_effect_Vmax)

                # compartments derivatives
                axis.roots.R_residual = self.respiration_model.RespirationModel.R_residual(axis.roots.sucrose, axis.roots.mstruct * model.Roots.PARAMETERS.ALPHA, axis.roots.Total_Organic_Nitrogen,
                                                                                              soil.Tsoil)
                axis.roots.sum_respi = axis.roots.R_Nnit_upt + axis.roots.R_Nnit_red + axis.roots.R_residual
                sucrose_derivative = axis.roots.calculate_sucrose_derivative(axis.roots.Unloading_Sucrose, axis.roots.S_Amino_Acids, axis.roots.C_exudation, axis.roots.sum_respi)
                nitrates_derivative = axis.roots.calculate_nitrates_derivative(axis.roots.Uptake_Nitrates, axis.roots.Export_Nitrates, axis.roots.S_Amino_Acids)
                amino_acids_derivative = axis.roots.calculate_amino_acids_derivative(axis.roots.Unloading_Amino_Acids, axis.roots.S_Amino_Acids, axis.roots.Export_Amino_Acids, axis.roots.N_exudation)
                cytokinins_derivative = axis.roots.calculate_cytokinins_derivative(axis.roots.S_cytokinins, axis.roots.Export_cytokinins)

                y_derivatives[self.initial_conditions_mapping[axis.roots]['sucrose']] = sucrose_derivative
                y_derivatives[self.initial_conditions_mapping[axis.roots]['nitrates']] = nitrates_derivative
                y_derivatives[self.initial_conditions_mapping[axis.roots]['amino_acids']] = amino_acids_derivative
                y_derivatives[self.initial_conditions_mapping[axis.roots]['cytokinins']] = cytokinins_derivative

                # compute the derivative of each compartment of phloem
                sucrose_phloem_derivative = axis.phloem.calculate_sucrose_derivative(phloem_contributors)
                amino_acids_phloem_derivative = axis.phloem.calculate_amino_acids_derivative(phloem_contributors)
                y_derivatives[self.initial_conditions_mapping[axis.phloem]['sucrose']] = sucrose_phloem_derivative
                y_derivatives[self.initial_conditions_mapping[axis.phloem]['amino_acids']] = amino_acids_phloem_derivative

                # compute the derivative of each compartment of axis
                C_exudated = axis.calculate_C_exudated(axis.roots.C_exudation, axis.roots.N_exudation, axis.roots.mstruct)
                y_derivatives[self.initial_conditions_mapping[axis]['C_exudated']] += C_exudated
                y_derivatives[self.initial_conditions_mapping[axis]['sum_respi_roots']] += axis.roots.sum_respi
                y_derivatives[self.initial_conditions_mapping[axis]['sum_respi_shoot']] += sum_respi_shoot

        # compute the derivative of each compartment of soil
        soil.mineralisation = soil.calculate_mineralisation(soil.T_effect_Vmax)
        y_derivatives[self.initial_conditions_mapping[soil]['nitrates']] = soil.calculate_nitrates_derivative(soil.mineralisation, soil_contributors, self.culm_density, soil.constant_Conc_Nitrates)

        if self.show_progressbar:
            self.progressbar.update(t)

        derivatives_logger = logging.getLogger('cnwheat.derivatives')
        if logger.isEnabledFor(logging.DEBUG) and derivatives_logger.isEnabledFor(logging.DEBUG):
            self._log_compartments(t_abs, y_derivatives, Simulation.LOGGERS_NAMES['derivatives'])

        return y_derivatives
