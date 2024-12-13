"""
    turgorgrowth.simulation
    ~~~~~~~~~~~~~~~~~~

    The module :mod:`turgorgrowth.simulation` is the front-end to run the model Turgor-Growth.
    The public API consists of methods :meth:`initialize` and :meth:`run`.

    :license: CeCILL-C, see LICENSE for details.

"""

from __future__ import division  # use "//" to do integer division
import logging

import numpy as np
from scipy.integrate import solve_ivp
from scipy import interpolate


from turgorgrowth import model
from turgorgrowth import parameters


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
    when checking the consistency of inputs `population` (see :meth:`initialize`).
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
    The Simulation class allows to initialize and run the model.

    User should use method :meth:`initialize` to initialize the model, and method
    :meth:`run` to run the model.

    """

    #: the name of the compartments attributes in the model
    MODEL_COMPARTMENTS_NAMES = {model.Plant: [],
                                model.Axis: [],
                                model.Roots: [],
                                model.Phytomer: [],
                                model.Organ: [],
                                model.HiddenZone: ['leaf_L', 'turgor_water_potential', 'water_content', 'width', 'thickness'],
                                model.PhotosyntheticOrganElement: ['length', 'turgor_water_potential', 'water_content', 'width', 'thickness']
                                }

    #: the time index
    T_INDEX = ['t']

    # --------------------------------------------------------------------------
    #           DEFINITION OF THE PARAMETERS AND COMPUTED VARIABLES
    # --------------------------------------------------------------------------

    # ---------- PLANT scale ----------

    #: the index to locate the plants in the modeled system
    PLANTS_INDEXES = ['plant']
    #: concatenation of :attr:`T_INDEX` and :attr:`PLANTS_INDEXES`
    PLANTS_T_INDEXES = T_INDEX + PLANTS_INDEXES
    #: the parameters which define the state of the modeled system at plant scale
    PLANTS_STATE_PARAMETERS = []
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
    AXES_STATE_PARAMETERS = ['Tr', 'green_area', 'mstruct']
    #: the variables which define the state of the modeled system at axis scale,
    #: formed be the concatenation of :attr:`AXES_STATE_PARAMETERS` and the names
    #: of the compartments associated to each axis (see :attr:`MODEL_COMPARTMENTS_NAMES`)
    AXES_STATE = AXES_STATE_PARAMETERS + MODEL_COMPARTMENTS_NAMES.get(model.Axis, [])
    #: the variables that we need to compute in order to compute fluxes and/or compartments values at axis scale
    AXES_INTERMEDIATE_VARIABLES = []
    #: the fluxes exchanged between the compartments at axis scale
    AXES_FLUXES = []
    #: the variables computed by integrating values of axis components parameters/variables recursively
    AXES_INTEGRATIVE_VARIABLES = ['Total_Transpiration_turgor', 'Growth', 'total_water_influx', 'xylem_water_potential', 'delta_plant_water_content', 'plant_water_content', 'plant_WC_DM', 'LAI_turgor']
    #: all the variables computed during a run step of the simulation at axis scale
    AXES_RUN_VARIABLES = AXES_STATE + AXES_INTERMEDIATE_VARIABLES + AXES_FLUXES + AXES_INTEGRATIVE_VARIABLES

    # ---------- PHYTOMER scale ----------

    #: the indexes to locate the phytomers in the modeled system
    PHYTOMERS_INDEXES = ['plant', 'axis', 'metamer']
    #: concatenation of :attr:`T_INDEX` and :attr:`PHYTOMERS_INDEXES`
    PHYTOMERS_T_INDEXES = T_INDEX + PHYTOMERS_INDEXES
    #: the parameters which define the state of the modeled system at phytomer scale
    PHYTOMERS_STATE_PARAMETERS = ['Tr', 'green_area']
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

    #: the indexes to locate the organ in the modeled system
    ORGANS_INDEXES = ['plant', 'axis', 'organ']
    #: concatenation of :attr:`T_INDEX` and :attr:`ORGAN_INDEXES`
    ORGANS_T_INDEXES = T_INDEX + ORGANS_INDEXES
    #: the parameters which define the state of the modeled system at organ scale
    ORGANS_STATE_PARAMETERS = ['age', 'amino_acids', 'proteins', 'sucrose', 'mstruct', 'Tr', 'green_area', 'SRWC', 'Tsoil']
    #: the variables which define the state of the modeled system at organ scale,
    #: formed be the concatenation of :attr:`ORGAN_STATE_PARAMETERS` and the names
    #: of the compartments associated to each organ (see :attr:`MODEL_COMPARTMENTS_NAMES`)
    ORGANS_STATE = ORGANS_STATE_PARAMETERS + MODEL_COMPARTMENTS_NAMES.get(model.Organ, [])
    #: the variables that we need to compute in order to compute fluxes and/or compartments values at organ scale
    ORGANS_INTERMEDIATE_VARIABLES = ['soil_water_potential', 'total_water_potential', 'Ksoil']
    #: the fluxes exchanged between the compartments at organ scale
    ORGANS_FLUXES = []
    #: the variables computed by integrating values of xylem components parameters/variables recursively
    ORGANS_INTEGRATIVE_VARIABLES = ['xylem_water_potential']
    #: all the variables computed during a run step of the simulation at organ scale
    ORGANS_RUN_VARIABLES = ORGANS_STATE + ORGANS_INTERMEDIATE_VARIABLES + ORGANS_FLUXES + ORGANS_INTEGRATIVE_VARIABLES

    # ---------- HIDDENZONE scale ----------

    #: the indexes to locate the hidden zones in the modeled system
    HIDDENZONE_INDEXES = ['plant', 'axis', 'metamer']
    #: concatenation of :attr:`T_INDEX` and :attr:`HIDDENZONE_INDEXES`
    HIDDENZONE_T_INDEXES = T_INDEX + HIDDENZONE_INDEXES
    #: the parameters which define the state of the modeled system at hidden zone scale
    HIDDENZONE_STATE_PARAMETERS = ['delta_teq', 'leaf_pseudo_age', 'leaf_pseudostem_length', 'fructan', 'amino_acids', 'proteins', 'sucrose', 'temperature', 'mstruct', 'hiddenzone_age', 'leaf_enclosed_mstruct', 'width_prev', 'thickness_prev', 'leaf_Wmax', 'init_leaf_L']
    #: the variables which define the state of the modeled system at hidden zone scale,
    #: formed be the concatenation of :attr:`HIDDENZONE_STATE_PARAMETERS` and the names
    #: of the compartments associated to each hidden zone (see :attr:`MODEL_COMPARTMENTS_NAMES`)
    HIDDENZONE_STATE = HIDDENZONE_STATE_PARAMETERS + MODEL_COMPARTMENTS_NAMES.get(model.HiddenZone, [])
    #: the variables that we need to compute in order to compute fluxes and/or compartments values at hidden zone scale
    HIDDENZONE_INTERMEDIATE_VARIABLES = ['osmotic_water_potential', 'resistance', 'total_water_potential', 'volume', 'length', 'phi_width', 'phi_thickness', 'phi_length', 'phi_volume', 'epsilon_volume', 'organ_volume', 'WC_mstruct', 'omega', 'leaf_Lmax', 'delta_hiddenzone_dimensions_plastic']
    #: the fluxes exchanged between the compartments at hidden zone scale
    HIDDENZONE_FLUXES = ['water_influx', 'water_outflow', 'Growth']
    #: the variables computed by integrating values of hidden zone components parameters/variables recursively
    HIDDENZONE_INTEGRATIVE_VARIABLES = []
    #: all the variables computed during a run step of the simulation at plant scale
    HIDDENZONE_RUN_VARIABLES = HIDDENZONE_STATE + HIDDENZONE_INTERMEDIATE_VARIABLES + HIDDENZONE_FLUXES + HIDDENZONE_INTEGRATIVE_VARIABLES


    # ---------- ELEMENT scale ----------

    #: the indexes to locate the ELEMENTS in the modeled system
    ELEMENTS_INDEXES = ['plant', 'axis', 'metamer', 'organ', 'element']
    #: concatenation of :attr:`T_INDEX` and :attr:`ELEMENTS_INDEXES`
    ELEMENTS_T_INDEXES = T_INDEX + ELEMENTS_INDEXES
    #: the parameters which define the state of the modeled system at organ scale
    ELEMENTS_STATE_PARAMETERS = ['amino_acids', 'green_area', 'mstruct', 'proteins', 'sucrose', 'fructan', 'Ts', 'Tr', 'age', 'is_growing', 'Wmax']
    #: the variables which define the state of the modeled system at organ scale,
    #: formed be the concatenation of :attr:`ELEMENTS_STATE_PARAMETERS` and the names
    #: of the compartments associated to each organ (see :attr:`MODEL_COMPARTMENTS_NAMES`)
    ELEMENTS_STATE = ELEMENTS_STATE_PARAMETERS + MODEL_COMPARTMENTS_NAMES.get(model.PhotosyntheticOrganElement, [])
    #: the variables that we need to compute in order to compute fluxes and/or compartments values at organ scale
    ELEMENTS_INTERMEDIATE_VARIABLES = ['osmotic_water_potential', 'resistance', 'total_water_potential', 'volume', 'epsilon_volume', 'organ_volume', 'WC_mstruct']
    #: the fluxes exchanged between the compartments at organ scale
    ELEMENTS_FLUXES = ['water_influx']
    #: the variables computed by integrating values of organ components parameters/variables recursively
    ELEMENTS_INTEGRATIVE_VARIABLES = ['Total_Transpiration_turgor', 'total_water_influx']
    #: all the variables computed during a run step of the simulation at organ scale
    ELEMENTS_RUN_VARIABLES = ELEMENTS_STATE + ELEMENTS_INTERMEDIATE_VARIABLES + ELEMENTS_FLUXES + ELEMENTS_INTEGRATIVE_VARIABLES

    #: a dictionary of all the variables which define the state of the modeled system, for each scale
    ALL_STATE_PARAMETERS = {model.Plant: PLANTS_STATE_PARAMETERS,
                            model.Axis: AXES_STATE_PARAMETERS,
                            model.Organ: ORGANS_STATE_PARAMETERS,
                            model.Phytomer: PHYTOMERS_STATE_PARAMETERS,
                            model.HiddenZone: HIDDENZONE_STATE_PARAMETERS,
                            model.PhotosyntheticOrganElement: ELEMENTS_STATE_PARAMETERS,
                            }

    #: the names of the elements forcings
    ELEMENTS_FORCINGS = ('green_area', 'Tr')

    #: the name of the loggers for compartments and derivatives
    LOGGERS_NAMES = {'compartments': {model.Plant: 'turgorgrowth.compartments.plants',
                                      model.Axis: 'turgorgrowth.compartments.axes',
                                      model.Phytomer: 'turgorgrowth.compartments.phytomers',
                                      model.Organ: 'turgorgrowth.compartments.organs',
                                      model.HiddenZone: 'turgorgrowth.compartments.hiddenzones',
                                      model.PhotosyntheticOrganElement: 'turgorgrowth.compartments.elements'},
                     'derivatives': {model.Plant: 'turgorgrowth.derivatives.plants',
                                     model.Axis: 'turgorgrowth.derivatives.axes',
                                     model.Phytomer: 'turgorgrowth.derivatives.phytomers',
                                     model.Organ: 'turgorgrowth.derivatives.organs',
                                     model.HiddenZone: 'turgorgrowth.derivatives.hiddenzones',
                                     model.PhotosyntheticOrganElement: 'turgorgrowth.derivatives.elements'}}

    def __init__(self, delta_t=1, interpolate_forcings=False, elements_forcings_delta_t=None, hiddenzone_forcings_delta_t=None):

        self.population = model.Population()  #: the population to simulate on
        self.mapping_topology = {}  #: a dict describing plant topology at element scale

        self.initial_conditions = []  #: the initial conditions of the compartments in the population and soil
        self.initial_conditions_mapping = {}  #: dictionary to map the compartments to their indexes in :attr:`initial_conditions`

        self.delta_t = delta_t  #: the delta t of the simulation (in seconds)

        self.time_step = self.delta_t / 3600.0  #: time step of the simulation (in hours)

        self.time_grid = np.array([0.0, self.time_step])  #: the time grid of the simulation (in hours)

        self.interpolate_forcings = interpolate_forcings  #: a boolean flag which indicates if we want to interpolate or not the forcings (True: interpolate, False: do not interpolate)

        # set the loggers for compartments and derivatives
        compartments_logger = logging.getLogger('turgorgrowth.compartments')
        derivatives_logger = logging.getLogger('turgorgrowth.derivatives')
        if compartments_logger.isEnabledFor(logging.DEBUG) or derivatives_logger.isEnabledFor(logging.DEBUG):
            sep = ','
            if compartments_logger.isEnabledFor(logging.DEBUG):
                plants_compartments_logger = logging.getLogger('turgorgrowth.compartments.plants')
                plants_compartments_logger.debug(sep.join(Simulation.PLANTS_T_INDEXES + Simulation.PLANTS_STATE))
                axes_compartments_logger = logging.getLogger('turgorgrowth.compartments.axes')
                axes_compartments_logger.debug(sep.join(Simulation.AXES_T_INDEXES + Simulation.AXES_STATE))
                phytomers_compartments_logger = logging.getLogger('turgorgrowth.compartments.phytomers')
                phytomers_compartments_logger.debug(sep.join(Simulation.PHYTOMERS_T_INDEXES + Simulation.PHYTOMERS_STATE))
                organs_compartments_logger = logging.getLogger('turgorgrowth.compartments.organs')
                organs_compartments_logger.debug(sep.join(Simulation.ORGANS_T_INDEXES + Simulation.ORGANS_STATE))
                hiddenzones_compartments_logger = logging.getLogger('turgorgrowth.compartments.hiddenzones')
                hiddenzones_compartments_logger.debug(sep.join(Simulation.HIDDENZONE_T_INDEXES + Simulation.HIDDENZONE_STATE + Simulation.HIDDENZONE_INTERMEDIATE_VARIABLES))
                # hiddenzones_compartments_logger.debug(sep.join(Simulation.HIDDENZONE_T_INDEXES + Simulation.HIDDENZONE_RUN_VARIABLES))
                elements_compartments_logger = logging.getLogger('turgorgrowth.compartments.elements')
                elements_compartments_logger.debug(sep.join(Simulation.ELEMENTS_T_INDEXES + Simulation.ELEMENTS_STATE + Simulation.ELEMENTS_INTERMEDIATE_VARIABLES))
                # elements_compartments_logger.debug(sep.join(Simulation.ELEMENTS_T_INDEXES + Simulation.ELEMENTS_RUN_VARIABLES))
                organs_compartments_logger = logging.getLogger('turgorgrowth.compartments.organs')
                organs_compartments_logger.debug(sep.join(Simulation.ORGANS_T_INDEXES + Simulation.ORGANS_STATE))
            if derivatives_logger.isEnabledFor(logging.DEBUG):
                plants_derivatives_logger = logging.getLogger('turgorgrowth.derivatives.plants')
                plants_derivatives_logger.debug(sep.join(Simulation.PLANTS_T_INDEXES + Simulation.PLANTS_STATE))
                axes_derivatives_logger = logging.getLogger('turgorgrowth.derivatives.axes')
                axes_derivatives_logger.debug(sep.join(Simulation.AXES_T_INDEXES + Simulation.AXES_STATE))
                phytomers_derivatives_logger = logging.getLogger('turgorgrowth.derivatives.phytomers')
                phytomers_derivatives_logger.debug(sep.join(Simulation.PHYTOMERS_T_INDEXES + Simulation.PHYTOMERS_STATE))
                organs_derivatives_logger = logging.getLogger('turgorgrowth.derivatives.organs')
                organs_derivatives_logger.debug(sep.join(Simulation.ORGANS_T_INDEXES + Simulation.ORGANS_STATE))
                hiddenzones_derivatives_logger = logging.getLogger('turgorgrowth.derivatives.hiddenzones')
                hiddenzones_derivatives_logger.debug(sep.join(Simulation.HIDDENZONE_T_INDEXES + Simulation.HIDDENZONE_STATE + Simulation.HIDDENZONE_INTERMEDIATE_VARIABLES))
                # hiddenzones_derivatives_logger.debug(sep.join(Simulation.HIDDENZONE_T_INDEXES + Simulation.HIDDENZONE_RUN_VARIABLES))
                elements_derivatives_logger = logging.getLogger('turgorgrowth.derivatives.elements')
                elements_derivatives_logger.debug(sep.join(Simulation.ELEMENTS_T_INDEXES + Simulation.ELEMENTS_STATE + Simulation.ELEMENTS_INTERMEDIATE_VARIABLES))
                # elements_derivatives_logger.debug(sep.join(Simulation.ELEMENTS_T_INDEXES + Simulation.HIDDENZONE_RUN_VARIABLES))
                organs_derivatives_logger = logging.getLogger('turgorgrowth.derivatives.organs')
                organs_derivatives_logger.debug(sep.join(Simulation.ORGANS_T_INDEXES + Simulation.ORGANS_STATE))

        logger = logging.getLogger(__name__)
        if logger.isEnabledFor(logging.DEBUG):
            self.t_offset = 0.0  #: the absolute time offset elapsed from the beginning of the simulation

        if interpolate_forcings:
            if (elements_forcings_delta_t is not None and hiddenzone_forcings_delta_t is not None):
                self.elements_forcings_delta_t_ratio = elements_forcings_delta_t / delta_t  #: the ratio between the delta t of the elements forcings and the delta t of the simulation
                self.hiddenzone_forcings_delta_t_ratio = hiddenzone_forcings_delta_t / delta_t  #: the ratio between the delta t of the hiddenzone forcings and the delta t of the simulation
            elif elements_forcings_delta_t is None:
                message = """The value of `interpolate_forcings` passed to the Simulation constructor is `True`, but `elements_forcings_delta_t` is `None`. 
        Please set `elements_forcings_delta_t` (through the Simulation constructor) to a not `None` value."""
                logger.exception(message)
                raise SimulationConstructionError(message)
            elif hiddenzone_forcings_delta_t is None:
                message = """The value of `interpolate_forcings` passed to the Simulation constructor is `True`, but `hiddenzone_forcings_delta_t` is `None`. 
                Please set `hiddenzone_forcings_delta_t` (through the Simulation constructor) to a not `None` value."""
                logger.exception(message)
                raise SimulationConstructionError(message)
            elif elements_forcings_delta_t < delta_t:
                message = """The value of `elements_forcings_delta_t` passed to the Simulation constructor is lesser than the one of `delta_t`. Please set a `elements_forcings_delta_t` that is at least equal to `delta_t`."""
                logger.exception(message)
                raise SimulationConstructionError(message)
            elif hiddenzone_forcings_delta_t < delta_t:
                message = """The value of `hiddenzone_forcings_delta_t` passed to the Simulation constructor is lesser than the one of `delta_t`. Please set a `hiddenzone_forcings_delta_t` that is at least equal to `delta_t`."""
                logger.exception(message)
                raise SimulationConstructionError(message)

            self.previous_forcings_values = {}  #: previous values of the forcings
            self.new_forcings_values = {}  #: new values of the forcings
            self.interpolation_functions = {}  #: functions to interpolate the forcings

        self.nfev_total = 0  #: cumulative number of RHS function evaluations

    def initialize(self, population, mapping_topology, SRWC=80):
        """
        Initialize:
            * :attr:`population`,
            * :attr:`initial_conditions_mapping`,
            * and :attr:`initial_conditions`

        from `population`

        :param model.Population population: a population of plants.
        :param dict mapping_topology: a dictionary with the topological mapping of the predecessor and successor of each organ
        """

        logger = logging.getLogger(__name__)

        logger.info('Initialization of the simulation...')

        # clean the attributes of the simulation
        del self.population.plants[:]

        del self.initial_conditions[:]
        self.initial_conditions_mapping.clear()

        self.mapping_topology = mapping_topology

        # create new population
        self.population.plants.extend(population.plants)

        # check the consistency of population
        if len(self.population.plants) != 0:  # population must contain at least 1 plant
            for plant in self.population.plants:
                if len(plant.axes) != 0:  # each plant must contain at least 1 axis
                    for axis in plant.axes:
                        if axis.roots is None:  # each axis must have a "roots"
                            message = 'No roots found in (plant={},axis={})'.format(plant.index, axis.label)
                            logger.exception(message)
                            raise SimulationInitializationError(message)
                        if axis.xylem is None:  # each axis must have a xylem
                            message = 'No xylem found in (plant={},axis={})'.format(plant.index, axis.label)
                            logger.exception(message)
                            raise SimulationInitializationError(message)
                        if len(axis.phytomers) != 0:  # each axis must contain at least 1 phytomer
                            for phytomer in axis.phytomers:
                                # phytomer_organs = (phytomer.lamina, phytomer.internode, phytomer.sheath)
                                phytomer_organs = (phytomer.lamina, phytomer.sheath)
                                # each phytomer must contain at least 1 photosynthetic organ or an hidden zone
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
                    for organ in (axis.xylem, axis.roots):
                        xylem_id = (plant.index, axis.label, organ.label)
                        self.new_forcings_values[xylem_id] = {}
                        forcing_labels = {}
                        self.new_forcings_values[xylem_id][forcing_labels] = getattr(axis.xylem, forcing_labels)
                        if axis.xylem in self.previous_forcings_values:
                            setattr(axis.xylem, forcing_labels, self.previous_forcings_values[xylem_id][forcing_labels])
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

        # Update SRWC using weather data
        for plant in self.population.plants:
            plant.SRWC = SRWC

        # # Update SRWC using weather data
        # for plant in self.population.plants:
        #     for axis in plant.axes:
        #         for organ in (axis.xylem, axis.roots):
        #             if organ == axis.xylem:
        #                 axis.xylem.SRWC = SRWC

        # initialize initial conditions
        def _init_initial_conditions(model_object, i):
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
                    self.initial_conditions_mapping[model_object][compartment_name] = i
                    self.initial_conditions.append(0)
                    i += 1
            return i

        i = 0

        for plant in self.population.plants:
            i = _init_initial_conditions(plant, i)
            for axis in plant.axes:
                i = _init_initial_conditions(axis, i)
                for organ in (axis.roots, axis.xylem):
                    if organ is None:
                        continue
                    i = _init_initial_conditions(organ, i)
                for phytomer in axis.phytomers:
                    i = _init_initial_conditions(phytomer, i)
                    if phytomer.hiddenzone:
                        i = _init_initial_conditions(phytomer.hiddenzone, i)
                    # for organ in (phytomer.lamina, phytomer.internode, phytomer.sheath):
                    for organ in (phytomer.lamina, phytomer.sheath):
                        if organ:
                            for element in (organ.exposed_element, organ.enclosed_element):
                                if element:
                                    i = _init_initial_conditions(element, i)

        self.population.calculate_aggregated_variables()

        logger.info('Initialization of the simulation DONE')

    def run(self):
        """
        Compute turgor pressure driven growth in :attr:`population` over :attr:`delta_t`.
        """
        logger = logging.getLogger(__name__)
        logger.info('Run of Turgor-Growth...')

        if self.interpolate_forcings:
            # interpolate the forcings
            self._interpolate_forcings()

        self._update_initial_conditions()

        if logger.isEnabledFor(logging.DEBUG):
            logger.debug("Run the solver with delta_t = %s", self.time_step)

        #: Call :func:`scipy.integrate.solve_ivp` to integrate the system over self.time_grid.
        sol = solve_ivp(fun=self._calculate_all_derivatives, t_span=self.time_grid, y0=self.initial_conditions,
                        method='LSODA', t_eval=np.array([self.time_step]), dense_output=False)

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

        logger.info('Run of Turgor-Growth DONE')

    def _update_initial_conditions(self):
        """Update the compartments values in :attr:`initial_conditions` from the compartments values of :attr:`population`.
        """
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
                for phytomer in axis.phytomers:
                    for organ in (phytomer.lamina, phytomer.sheath):
                        if organ is None:
                            continue
                        for element in (organ.exposed_element, organ.enclosed_element):
                            if element is not None:
                                element_id = (plant.index, axis.label, phytomer.index, organ.label, element.label)
                                self.interpolation_functions[element_id] = {}
                                next_forcings_values[element_id] = {}
                                forcing_labels = Simulation.ELEMENTS_FORCINGS
                                forcings_delta_t_ratio = self.elements_forcings_delta_t_ratio
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
            intermediate_variables = None
            if issubclass(class_, model.HiddenZone):
                class_ = model.HiddenZone
                intermediate_variables = Simulation.HIDDENZONE_INTERMEDIATE_VARIABLES
            elif issubclass(class_, model.PhotosyntheticOrganElement):
                class_ = model.PhotosyntheticOrganElement
                intermediate_variables = Simulation.ELEMENTS_INTERMEDIATE_VARIABLES
            elif issubclass(class_, model.Organ):
                class_ = model.Organ
                intermediate_variables = Simulation.ORGANS_INTERMEDIATE_VARIABLES
            parameters_names = Simulation.ALL_STATE_PARAMETERS[class_] + intermediate_variables
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

        for plant in self.population.plants:
            for axis in plant.axes:
                for organ in (axis.roots, axis.xylem):
                    if organ is None:
                        continue
                    i = update_rows(organ, [t, plant.index, axis.label, organ.label], all_rows[model.Organ], i)
                for phytomer in axis.phytomers:
                    # for organ in (phytomer.lamina, phytomer.internode, phytomer.sheath, phytomer.hiddenzone):
                    for organ in (phytomer.lamina, phytomer.sheath, phytomer.hiddenzone):
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
              At first call to :meth:`_calculate_all_derivatives` by :func:`scipy.integrate.solve_ivp`
              `t` = **t_span** [0], where **t_span** is one of the arguments passed to :func:`solve_ivp(fun, t_span, y0,...) <scipy.integrate.solve_ivp>`.
              For each following call to :meth:`_calculate_all_derivatives`, `t` belongs
              to the interval [**t_span** [0], **t_span** [1]].
        :param list y: The current values of y.
              At first call to :meth:`_calculate_all_derivatives` by :func:`scipy.integrate.solve_ivp`, `y` = **y0**
              where **y0** is one of the arguments passed to :func:`solve_ivp(fun, t_span, y0,...) <scipy.integrate.solve_ivp>`.
              Then, values of `y` are chosen automatically by :func:`scipy.integrate.solve_ivp`.

        :return: The derivatives of `y` at `t`.
        :rtype: list
        """
        logger = logging.getLogger(__name__)

        if logger.isEnabledFor(logging.DEBUG):
            t_abs = t + self.t_offset
            logger.debug('t = {}'.format(t_abs))

        if self.interpolate_forcings:
            # Update state parameters using interpolation functions
            for plant in self.population.plants:
                for axis in plant.axes:
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

        compartments_logger = logging.getLogger('turgorgrowth.compartments')
        if logger.isEnabledFor(logging.DEBUG) and compartments_logger.isEnabledFor(logging.DEBUG):
            self._log_compartments(t_abs, y, Simulation.LOGGERS_NAMES['compartments'])

        # Check that the solver is not crashed
        y_isnan = np.isnan(y)
        if y_isnan.any():
            message = 'The solver did not manage to compute a compartment. See the logs. NaN found in y'
            logger.exception(message)
            raise SimulationRunError(message)

        y_derivatives = np.zeros_like(y)

        #: Water flux with xylem and organs
        for plant in self.population.plants:
            for axis in plant.axes:
                # Xylem
                #: Soil water potential
                axis.xylem.Ksoil = axis.xylem.calculate_soil_conductivity(axis.LAI_turgor)
                axis.xylem.soil_water_potential = axis.xylem.calculate_soil_water_potential(plant.SRWC)
                #: SRWC
                axis.xylem.SRWC = plant.SRWC
                #: Total water potential
                axis.xylem.total_water_potential = axis.xylem.calculate_xylem_water_potential(axis.xylem.soil_water_potential, axis.total_water_influx, axis.Growth, axis.xylem.Ksoil, self.delta_t)

                for phytomer in axis.phytomers:
                    # Hidden zone
                    hiddenzone = phytomer.hiddenzone
                    if hiddenzone is not None:
                        hiddenzone.leaf_L = y[self.initial_conditions_mapping[hiddenzone]['leaf_L']]
                        hiddenzone.thickness = y[self.initial_conditions_mapping[hiddenzone]['thickness']]
                        hiddenzone.width = y[self.initial_conditions_mapping[hiddenzone]['width']]
                        # Update of leaf_Lmax
                        hiddenzone.leaf_Lmax = hiddenzone.leaf_L

                        if hiddenzone.leaf_pseudo_age == 0:  #: First time after previous leaf emergence
                            #: Width and thickness of the previous leaf
                            # hiddenzone.width = hiddenzone.width_prev
                            # hiddenzone.thickness = hiddenzone.thickness_prev

                            # TEST 11.2024 VICTORIA
                            # v4
                            # width_ratio = 0.6
                            # thickness_ratio = 0.125
                            # v5
                            # width_ratio = 0.675
                            # thickness_ratio = 0.13

                            # v6
                            # # width_ratio = 0.60
                            # # width_ratio = 0.865
                            # width_ratio = 0.85
                            # # thickness_ratio = 0.125
                            # thickness_ratio = 0.225

                            # v7
                            # width_ratio = 0.75  # (1)
                            # thickness_ratio = 0.16  # (1)
                            # width_ratio = 0.6   # (2)
                            # thickness_ratio = 0.10  # (2)
                            # width_ratio = 0.625   # (3)
                            # thickness_ratio = 0.1  # (3)

                            # width_ratio = 0.8   # (3)
                            width_ratio = 0.775   # (3)
                            thickness_ratio = 0.1  # (3)

                            hiddenzone.width = hiddenzone.leaf_L * width_ratio
                            hiddenzone.thickness = hiddenzone.leaf_L * thickness_ratio

                            #: Volume & water content as function of dimensions
                            hiddenzone.volume = hiddenzone.length * hiddenzone.width * hiddenzone.thickness
                            hiddenzone.water_content = hiddenzone.volume * parameters.RHO_WATER
                            #: Osmotic water potential
                            hiddenzone.omega = hiddenzone.calculate_solutes_contribution(hiddenzone.fructan, hiddenzone.sucrose, hiddenzone.amino_acids, hiddenzone.volume)
                            hiddenzone.osmotic_water_potential = hiddenzone.calculate_osmotic_water_potential(hiddenzone.fructan, hiddenzone.sucrose, hiddenzone.amino_acids, hiddenzone.volume, hiddenzone.temperature, hiddenzone.omega)
                            #: Total water potential
                            hiddenzone.total_water_potential = axis.xylem.total_water_potential
                            #: Turgor water potential
                            hiddenzone.turgor_water_potential = hiddenzone.total_water_potential - hiddenzone.osmotic_water_potential
                            #: Length
                            hiddenzone.length = hiddenzone.calculate_hiddenzone_length(hiddenzone.leaf_L, hiddenzone.leaf_pseudostem_length)
                        elif hiddenzone.leaf_pseudo_age > 0:    #: After previous leaf emergence
                            #: Turgor water potential
                            hiddenzone.turgor_water_potential = y[self.initial_conditions_mapping[hiddenzone]['turgor_water_potential']]
                            #: Volume and water content
                            hiddenzone.water_content = y[self.initial_conditions_mapping[hiddenzone]['water_content']]
                            hiddenzone.volume = hiddenzone.calculate_volume(hiddenzone.water_content)
                            #: Osmotic water potential
                            hiddenzone.omega = hiddenzone.calculate_solutes_contribution(hiddenzone.fructan, hiddenzone.sucrose, hiddenzone.amino_acids, hiddenzone.volume)
                            hiddenzone.osmotic_water_potential = hiddenzone.calculate_osmotic_water_potential(hiddenzone.fructan, hiddenzone.sucrose, hiddenzone.amino_acids, hiddenzone.volume, hiddenzone.temperature, hiddenzone.omega)
                            #: Total water potential
                            hiddenzone.total_water_potential = hiddenzone.calculate_water_potential(hiddenzone.turgor_water_potential, hiddenzone.osmotic_water_potential)
                            #: Length
                            hiddenzone.length = hiddenzone.calculate_hiddenzone_length(hiddenzone.leaf_L, hiddenzone.leaf_pseudostem_length)
                        else:   #: Before previous leaf emergence (calculation in elong-wheat)
                            continue

                        #: Resistance to water flow
                        hiddenzone_dimensions = {'length': hiddenzone.length, 'thickness': hiddenzone.thickness, 'width': hiddenzone.width}
                        hiddenzone.resistance = hiddenzone.calculate_resistance(hiddenzone_dimensions)
                        #: Flows with xylem
                        hiddenzone.water_influx = hiddenzone.calculate_water_flux(hiddenzone.total_water_potential, axis.xylem.total_water_potential, hiddenzone.resistance, self.delta_t)
                        hiddenzone.water_outflow = 0    #: No water flow between hiddenzone and element


                    # Photosynthetic Organ Elements
                    # for organ in (phytomer.lamina, phytomer.internode, phytomer.sheath):
                    for organ in (phytomer.lamina, phytomer.sheath):
                        if organ is None:
                            continue
                        for element in (organ.exposed_element, organ.enclosed_element):
                            if element is None:
                                continue
                            element.length = y[self.initial_conditions_mapping[element]['length']]
                            element.width = y[self.initial_conditions_mapping[element]['width']]
                            element.thickness = y[self.initial_conditions_mapping[element]['thickness']]
                            element.organ_dimensions = {'length': element.length, 'width': element.width, 'thickness': element.thickness}

                            if element.age == 0:    #: First time after element emergence
                                #: Width and thickness of the hiddenzone
                                if organ.label == "blade":
                                    element.thickness = hiddenzone.thickness
                                    element.width = hiddenzone.width
                                #: Volume and water content as function of dimensions
                                element.volume = element.calculate_organ_volume(element.organ_dimensions)
                                element.water_content = element.volume * parameters.RHO_WATER
                                # : Osmotic water potential
                                element.osmotic_water_potential = element.calculate_osmotic_water_potential(element.sucrose, element.amino_acids, element.volume, element.temperature, element.fructan)
                                #: Total water potential
                                element.total_water_potential = axis.xylem.total_water_potential
                                # : Turgor water potential
                                element.turgor_water_potential = element.total_water_potential - element.osmotic_water_potential

                            elif element.age > 0:   #: Emerged element
                                #: Water content
                                element.water_content = y[self.initial_conditions_mapping[element]['water_content']]
                                # : Volume
                                element.volume = element.calculate_volume(element.water_content)
                                # : Turgor water potential
                                element.turgor_water_potential = y[self.initial_conditions_mapping[element]['turgor_water_potential']]
                                #: Osmotic water potential
                                element.osmotic_water_potential = element.calculate_osmotic_water_potential(element.sucrose, element.amino_acids, element.volume, element.temperature, element.fructan)
                                #: Total water potential
                                element.total_water_potential = element.calculate_water_potential(element.turgor_water_potential, element.osmotic_water_potential)

                            #: Resistance to water flow
                            element.resistance = element.calculate_resistance(element.organ_dimensions)
                            #: Water fluxes with xylem
                            element.water_influx = element.calculate_water_flux(element.total_water_potential, axis.xylem.total_water_potential, element.resistance, self.delta_t)


        #: compute the derivative of each compartment of element
        for plant in self.population.plants:
            for axis in plant.axes:
                for phytomer in axis.phytomers:
                    # Hidden zone
                    hiddenzone = phytomer.hiddenzone
                    if hiddenzone is not None:
                        #: Delta water content
                        delta_water_content_hz = hiddenzone.calculate_delta_water_content(hiddenzone.water_influx, hiddenzone.water_outflow)
                        #: Extensibility
                        # phi = hiddenzone.calculate_extensibility_init(hiddenzone.leaf_pseudo_age, self.delta_t)
                        # phi = hiddenzone.calculate_extensibility_dimensions(hiddenzone.leaf_pseudo_age, self.delta_t)
                        phi = hiddenzone.calculate_extensibility_temperature(hiddenzone.leaf_pseudo_age, hiddenzone.delta_teq, self.delta_t)
                        hiddenzone.phi_length = phi['z']  # extensibility for length
                        hiddenzone.phi_width = phi['x']  # extensibility for length
                        hiddenzone.phi_thickness = phi['y']  # extensibility for length
                        hiddenzone.phi_volume = phi['x'] + phi['y'] + phi['z']  # volumetric extensibility
                        #: Elasticity
                        epsilon_x, epsilon_y, epsilon_z = hiddenzone.PARAMETERS.epsilon['x'], hiddenzone.PARAMETERS.epsilon['y'], hiddenzone.PARAMETERS.epsilon['z']
                        hiddenzone.epsilon_volume = (epsilon_x * epsilon_y * epsilon_z) / (epsilon_z * epsilon_x + epsilon_z * epsilon_y + epsilon_x * epsilon_y)  #: Elastic reversible growth (MPa)

                        if hiddenzone.leaf_pseudo_age > 0:  #: After previous leaf emergence
                            #: Derivatives
                            #: Delta turgor pressure
                            delta_turgor_water_potential = hiddenzone.calculate_delta_turgor_water_potential(phi, hiddenzone.turgor_water_potential, hiddenzone.volume, delta_water_content_hz)
                            #: Dimensions with plastic and elastic deformation
                            hiddenzone_dimensions = {'length': hiddenzone.length, 'width': hiddenzone.width, 'thickness': hiddenzone.thickness}
                            delta_hiddenzone_dimensions_plastic = hiddenzone.calculate_delta_organ_dimensions_plastic(hiddenzone.turgor_water_potential, phi, hiddenzone_dimensions)
                            delta_hiddenzone_dimensions_elastic = hiddenzone.calculate_delta_organ_dimensions_elastic(delta_turgor_water_potential, hiddenzone_dimensions)

                            # Saving into outputs
                            hiddenzone.delta_hiddenzone_dimensions_plastic = delta_hiddenzone_dimensions_plastic['length'] * 1000

                            if hiddenzone.leaf_L >= hiddenzone.leaf_pseudostem_length:  # Emerged blade
                                # Growing leaf with emerged blade
                                y_derivatives[self.initial_conditions_mapping[hiddenzone]['water_content']] = 0     # Transfer of water content to emerged element
                                y_derivatives[self.initial_conditions_mapping[hiddenzone]['turgor_water_potential']] = delta_turgor_water_potential
                                # Elastic deformation
                                y_derivatives[self.initial_conditions_mapping[hiddenzone]['width']] = delta_hiddenzone_dimensions_elastic['width']
                                y_derivatives[self.initial_conditions_mapping[hiddenzone]['thickness']] = delta_hiddenzone_dimensions_elastic['thickness']
                                y_derivatives[self.initial_conditions_mapping[hiddenzone]['leaf_L']] = delta_hiddenzone_dimensions_elastic['length']

                            else:  # Enclosed blade
                                y_derivatives[self.initial_conditions_mapping[hiddenzone]['water_content']] = delta_water_content_hz
                                y_derivatives[self.initial_conditions_mapping[hiddenzone]['turgor_water_potential']] = delta_turgor_water_potential
                                #  Plastic deformation
                                y_derivatives[self.initial_conditions_mapping[hiddenzone]['thickness']] = delta_hiddenzone_dimensions_elastic['thickness'] + delta_hiddenzone_dimensions_plastic['thickness']
                                y_derivatives[self.initial_conditions_mapping[hiddenzone]['width']] = delta_hiddenzone_dimensions_elastic['width'] + delta_hiddenzone_dimensions_plastic['width']
                                y_derivatives[self.initial_conditions_mapping[hiddenzone]['leaf_L']] = delta_hiddenzone_dimensions_elastic['length'] + delta_hiddenzone_dimensions_plastic['length']

                            hiddenzone.organ_volume = hiddenzone.calculate_organ_volume(hiddenzone_dimensions)
                            hiddenzone.WC_mstruct = hiddenzone.water_content / (hiddenzone.water_content + hiddenzone.mstruct) * 100

                        else:  #: Before previous leaf emergence
                            continue


                    # Photosynthetic Organ Elements
                    # for organ in (phytomer.lamina, phytomer.internode, phytomer.sheath):
                    for organ in (phytomer.lamina, phytomer.sheath):
                        if organ is None:
                            continue
                        for element in (organ.exposed_element, organ.enclosed_element):
                            if element is None:
                                continue
                            epsilon_x, epsilon_y, epsilon_z = element.PARAMETERS.epsilon['x'], element.PARAMETERS.epsilon['y'], element.PARAMETERS.epsilon['z']
                            element.epsilon_volume = (epsilon_z * epsilon_x * epsilon_y) / (epsilon_z * epsilon_x + epsilon_z * epsilon_y + epsilon_x * epsilon_y)  #: Elastic reversible growth (MPa)
                            element.organ_dimensions = {'length': element.length, 'width': element.width, 'thickness': element.thickness}
                            #: Delta water content
                            delta_water_content_ele = element.calculate_delta_water_content(element.water_influx, element.Total_Transpiration_turgor)
                            #: Delta turgor pressure
                            delta_turgor_water_potential = element.calculate_delta_turgor_water_potential(element.volume, delta_water_content_ele)
                            #: Dimensions
                            delta_element_dimensions = element.calculate_delta_organ_dimensions(delta_turgor_water_potential, element.organ_dimensions)
                            #: Derivatives
                            y_derivatives[self.initial_conditions_mapping[element]['turgor_water_potential']] = delta_turgor_water_potential

                            if hiddenzone is not None:  # Growing leaf
                                if element.is_growing == True:  # Growing element
                                    y_derivatives[self.initial_conditions_mapping[element]['water_content']] = delta_water_content_ele + delta_water_content_hz  # Transfer of water content from hiddenzone
                                    # Plastic deformation
                                    y_derivatives[self.initial_conditions_mapping[hiddenzone]['leaf_L']] = delta_hiddenzone_dimensions_plastic['length'] + delta_element_dimensions['length']
                                    y_derivatives[self.initial_conditions_mapping[element]['length']] = delta_hiddenzone_dimensions_plastic['length'] + delta_element_dimensions['length']
                                    y_derivatives[self.initial_conditions_mapping[element]['width']] = delta_hiddenzone_dimensions_plastic['width'] + delta_element_dimensions['width']
                                    y_derivatives[self.initial_conditions_mapping[element]['thickness']] = delta_hiddenzone_dimensions_plastic['thickness'] + delta_element_dimensions['thickness']
                                    y_derivatives[self.initial_conditions_mapping[element]['length']] = delta_hiddenzone_dimensions_plastic['length'] + delta_element_dimensions['length']

                                    # # Test meteo Gauthier et al. (2020) - 11.2024
                                    # if phytomer.index == 3:
                                    #     if hiddenzone.leaf_L >= 0.14:
                                    #         hiddenzone.leaf_L = 0.14
                                    #         element_L = hiddenzone.leaf_L - hiddenzone.leaf_pseudostem_length
                                    #         element.length = element_L
                                    # if phytomer.index == 4:
                                    #     if hiddenzone.leaf_L >= 0.16:
                                    #         hiddenzone.leaf_L = 0.16
                                    #         element_L = hiddenzone.leaf_L - hiddenzone.leaf_pseudostem_length
                                    #         element.length = element_L
                                    # if phytomer.index == 5:
                                    #     if hiddenzone.leaf_L >= 0.18:
                                    #         hiddenzone.leaf_L = 0.18
                                    #         element_L = hiddenzone.leaf_L - hiddenzone.leaf_pseudostem_length
                                    #         element.length = element_L
                                    # if phytomer.index == 6:
                                    #     if hiddenzone.leaf_L >= 0.20:
                                    #         hiddenzone.leaf_L = 0.20
                                    #         element_L = hiddenzone.leaf_L - hiddenzone.leaf_pseudostem_length
                                    #         element.length = element_L

                                elif element.is_growing == False:   # End of leaf elongation
                                    #: Derivatives
                                    y_derivatives[self.initial_conditions_mapping[element]['water_content']] = delta_water_content_ele
                                    # Elastic deformation
                                    y_derivatives[self.initial_conditions_mapping[element]['length']] = delta_element_dimensions['length']
                                    y_derivatives[self.initial_conditions_mapping[element]['width']] = delta_element_dimensions['width']
                                    y_derivatives[self.initial_conditions_mapping[element]['thickness']] = delta_element_dimensions['thickness']
                                    y_derivatives[self.initial_conditions_mapping[hiddenzone]['leaf_L']] = delta_element_dimensions['length']

                            else:   # Mature leaf
                                #: Derivatives
                                y_derivatives[self.initial_conditions_mapping[element]['water_content']] = delta_water_content_ele
                                # Elastic deformation
                                y_derivatives[self.initial_conditions_mapping[element]['length']] = delta_element_dimensions['length']
                                y_derivatives[self.initial_conditions_mapping[element]['width']] = delta_element_dimensions['width']
                                y_derivatives[self.initial_conditions_mapping[element]['thickness']] = delta_element_dimensions['thickness']

                            #: Dimensions volume of element
                            element.organ_volume = element.calculate_organ_volume(element.organ_dimensions)
                            element.WC_mstruct = element.water_content / (element.water_content + element.mstruct) * 100


        derivatives_logger = logging.getLogger('turgorgrowth.derivatives')
        if logger.isEnabledFor(logging.DEBUG) and derivatives_logger.isEnabledFor(logging.DEBUG):
            self._log_compartments(t_abs, y_derivatives, Simulation.LOGGERS_NAMES['derivatives'])

        return y_derivatives