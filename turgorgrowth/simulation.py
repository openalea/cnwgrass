# -*- coding: latin-1 -*-

from __future__ import division  # use "//" to do integer division
import logging

import numpy as np
from scipy.integrate import solve_ivp

import model

"""
    turgorgrowth.simulation
    ~~~~~~~~~~~~~~~~~~

    The module :mod:`turgorgrowth.simulation` is the front-end to run the model Turgor-Growth.
    The public API consists of methods :meth:`initialize` and :meth:`run`.

    :license: CeCILL-C, see LICENSE for details.

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

    """

    #: the name of the compartments attributes in the model
    MODEL_COMPARTMENTS_NAMES = {model.Plant: [],
                                model.Axis: [],
                                model.Phytomer: [],
                                model.HiddenZone: ['age', 'length', 'radius', 'turgor_water_potential', 'water_content'],
                                model.PhotosyntheticOrganElement: ['age', 'length', 'radius', 'thickness', 'turgor_water_potential', 'water_content', 'width']}

    #: the time index
    T_INDEX = ['t']

    ############################################################################
    ########### DEFINITION OF THE PARAMETERS AND COMPUTED VARIABLES ############
    ############################################################################

    ########### PLANT scale ############

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

    ########### AXIS scale ############

    #: the indexes to locate the axes in the modeled system
    AXES_INDEXES = ['plant', 'axis']
    #: concatenation of :attr:`T_INDEX` and :attr:`AXES_INDEXES`
    AXES_T_INDEXES = T_INDEX + AXES_INDEXES
    #: the parameters which define the state of the modeled system at axis scale
    AXES_STATE_PARAMETERS = []
    #: the variables which define the state of the modeled system at axis scale,
    #: formed be the concatenation of :attr:`AXES_STATE_PARAMETERS` and the names
    #: of the compartments associated to each axis (see :attr:`MODEL_COMPARTMENTS_NAMES`)
    AXES_STATE = AXES_STATE_PARAMETERS + MODEL_COMPARTMENTS_NAMES.get(model.Axis, [])
    #: the variables that we need to compute in order to compute fluxes and/or compartments values at axis scale
    AXES_INTERMEDIATE_VARIABLES = []
    #: the fluxes exchanged between the compartments at axis scale
    AXES_FLUXES = []
    #: the variables computed by integrating values of axis components parameters/variables recursively
    AXES_INTEGRATIVE_VARIABLES = []
    #: all the variables computed during a run step of the simulation at axis scale
    AXES_RUN_VARIABLES = AXES_STATE + AXES_INTERMEDIATE_VARIABLES + AXES_FLUXES + AXES_INTEGRATIVE_VARIABLES

    ########### PHYTOMER scale ############

    #: the indexes to locate the phytomers in the modeled system
    PHYTOMERS_INDEXES = ['plant', 'axis', 'metamer']
    #: concatenation of :attr:`T_INDEX` and :attr:`PHYTOMERS_INDEXES`
    PHYTOMERS_T_INDEXES = T_INDEX + PHYTOMERS_INDEXES
    #: the parameters which define the state of the modeled system at phytomer scale
    PHYTOMERS_STATE_PARAMETERS = []
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

    ########### HIDDENZONE scale ############

    #: the indexes to locate the hidden zones in the modeled system
    HIDDENZONE_INDEXES = ['plant', 'axis', 'metamer']
    #: concatenation of :attr:`T_INDEX` and :attr:`HIDDENZONE_INDEXES`
    HIDDENZONE_T_INDEXES = T_INDEX + HIDDENZONE_INDEXES
    #: the parameters which define the state of the modeled system at hidden zone scale
    HIDDENZONE_STATE_PARAMETERS = ['amino_acids', 'proteins', 'sucrose', 'temperature']
    #: the variables which define the state of the modeled system at hidden zone scale,
    #: formed be the concatenation of :attr:`HIDDENZONE_STATE_PARAMETERS` and the names
    #: of the compartments associated to each hidden zone (see :attr:`MODEL_COMPARTMENTS_NAMES`)
    HIDDENZONE_STATE = HIDDENZONE_STATE_PARAMETERS + MODEL_COMPARTMENTS_NAMES.get(model.HiddenZone, [])
    #: the variables that we need to compute in order to compute fluxes and/or compartments values at hidden zone scale
    HIDDENZONE_INTERMEDIATE_VARIABLES = ['osmotic_water_potential', 'resistance', 'total_water_potential']
    #: the fluxes exchanged between the compartments at hidden zone scale
    HIDDENZONE_FLUXES = ['water_influx']
    #: the variables computed by integrating values of hidden zone components parameters/variables recursively
    HIDDENZONE_INTEGRATIVE_VARIABLES = []
    #: all the variables computed during a run step of the simulation at plnat scale
    HIDDENZONE_RUN_VARIABLES = HIDDENZONE_STATE + HIDDENZONE_INTERMEDIATE_VARIABLES + HIDDENZONE_FLUXES + HIDDENZONE_INTEGRATIVE_VARIABLES

    ########### ELEMENT scale ############

    #: the indexes to locate the ELEMENTS in the modeled system
    ELEMENTS_INDEXES = ['plant', 'axis', 'metamer', 'organ', 'element']
    #: concatenation of :attr:`T_INDEX` and :attr:`ELEMENTS_INDEXES`
    ELEMENTS_T_INDEXES = T_INDEX + ELEMENTS_INDEXES
    #: the parameters which define the state of the modeled system at organ scale
    ELEMENTS_STATE_PARAMETERS = ['amino_acids', 'proteins', 'sucrose', 'temperature', 'transpiration']
    #: the variables which define the state of the modeled system at organ scale,
    #: formed be the concatenation of :attr:`ELEMENTS_STATE_PARAMETERS` and the names
    #: of the compartments associated to each organ (see :attr:`MODEL_COMPARTMENTS_NAMES`)
    ELEMENTS_STATE = ELEMENTS_STATE_PARAMETERS + MODEL_COMPARTMENTS_NAMES.get(model.PhotosyntheticOrganElement, [])
    #: the variables that we need to compute in order to compute fluxes and/or compartments values at organ scale
    ELEMENTS_INTERMEDIATE_VARIABLES = ['osmotic_water_potential', 'resistance', 'total_water_potential']
    #: the fluxes exchanged between the compartments at organ scale
    ELEMENTS_FLUXES = ['water_influx']
    #: the variables computed by integrating values of organ components parameters/variables recursively
    ELEMENTS_INTEGRATIVE_VARIABLES = []
    #: all the variables computed during a run step of the simulation at organ scale
    ELEMENTS_RUN_VARIABLES = ELEMENTS_STATE + ELEMENTS_INTERMEDIATE_VARIABLES + ELEMENTS_FLUXES + ELEMENTS_INTEGRATIVE_VARIABLES

    #: a dictionary of all the variables which define the state of the modeled system, for each scale
    ALL_STATE_PARAMETERS = {model.Plant: PLANTS_STATE_PARAMETERS,
                            model.Axis: AXES_STATE_PARAMETERS,
                            model.Phytomer: PHYTOMERS_STATE_PARAMETERS,
                            model.HiddenZone: HIDDENZONE_STATE_PARAMETERS,
                            model.PhotosyntheticOrganElement: ELEMENTS_STATE_PARAMETERS}

    #: the name of the loggers for compartments and derivatives
    LOGGERS_NAMES = {'compartments': {model.HiddenZone: 'turgorgrowth.compartments.hiddenzones',
                                      model.PhotosyntheticOrganElement: 'turgorgrowth.compartments.elements'},
                     'derivatives': {model.HiddenZone: 'turgorgrowth.derivatives.hiddenzones',
                                     model.PhotosyntheticOrganElement: 'turgorgrowth.derivatives.elements'}}

    def __init__(self, delta_t=1):
        
        self.population = model.Population()  #: the population to simulate on
        self.mapping_topology = {}  #: a dict describing plant topology at element scale

        self.initial_conditions = []  #: the initial conditions of the compartments in the population and soils
        self.initial_conditions_mapping = {}  #: dictionary to map the compartments to their indexes in :attr:`initial_conditions`

        self.delta_t = delta_t  #: the delta t of the simulation (in seconds)

        self.time_step = self.delta_t / 3600.0  #: time step of the simulation (in hours)
        
        self.time_grid = np.array([0.0, self.time_step])  #: the time grid of the simulation (in hours)

        # set the loggers for compartments and derivatives
        compartments_logger = logging.getLogger('turgorgrowth.compartments')
        derivatives_logger = logging.getLogger('turgorgrowth.derivatives')
        if compartments_logger.isEnabledFor(logging.DEBUG) or derivatives_logger.isEnabledFor(logging.DEBUG):
            sep = ','
            if compartments_logger.isEnabledFor(logging.DEBUG):
                hiddenzones_compartments_logger = logging.getLogger('turgorgrowth.compartments.hiddenzones')
                hiddenzones_compartments_logger.debug(sep.join(Simulation.HIDDENZONE_T_INDEXES + Simulation.HIDDENZONE_STATE))
                elements_compartments_logger = logging.getLogger('turgorgrowth.compartments.elements')
                elements_compartments_logger.debug(sep.join(Simulation.ELEMENTS_T_INDEXES + Simulation.ELEMENTS_STATE))
            if derivatives_logger.isEnabledFor(logging.DEBUG):
                hiddenzones_derivatives_logger = logging.getLogger('turgorgrowth.derivatives.hiddenzones')
                hiddenzones_derivatives_logger.debug(sep.join(Simulation.HIDDENZONE_T_INDEXES + Simulation.HIDDENZONE_STATE))
                elements_derivatives_logger = logging.getLogger('turgorgrowth.derivatives.elements')
                elements_derivatives_logger.debug(sep.join(Simulation.ELEMENTS_T_INDEXES + Simulation.ELEMENTS_STATE))

        logger = logging.getLogger(__name__)
        if logger.isEnabledFor(logging.DEBUG):
            self.t_offset = 0.0  #: the absolute time offset elapsed from the beginning of the simulation
        
        self.nfev_total = 0  #: cumulative number of RHS function evaluations

    def initialize(self, population, mapping_topology):
        """
        Initialize:

            * :attr:`population`,
            * :attr:`initial_conditions_mapping`,
            * and :attr:`initial_conditions`

        from `population`.

        :Parameters:

            - `population` (:class:`model.Population`) - a population of plants.
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

        # check the consistency of population and soils
        if len(self.population.plants) != 0:  # population must contain at least 1 plant
            for plant in self.population.plants:
                if len(plant.axes) != 0:  # each plant must contain at least 1 axis
                    for axis in plant.axes:
                        # if axis.roots is None:  # each axis must have a "roots"
                        #     message = 'No roots found in (plant={},axis={})'.format(plant.index, axis.label)
                        #     logger.exception(message)
                        #     raise SimulationInitializationError(message)
                        if len(axis.phytomers) != 0:  # each axis must contain at least 1 phytomer
                            for phytomer in axis.phytomers:
                                phytomer_organs = (phytomer.lamina, phytomer.internode, phytomer.sheath)
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
        
        # initialize initial conditions
        def _init_initial_conditions(model_object, i):
            class_ = model_object.__class__
            if issubclass(class_, model.HiddenZone):
                class_ = model.HiddenZone
            # elif issubclass(class_, model.Organ):
            #     class_ = model.Organ
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
                # for organ in (axis.roots, axis.phloem, axis.grains):
                #     if organ is None:
                #         continue
                #     i = _init_initial_conditions(organ, i)
                for phytomer in axis.phytomers:
                    i = _init_initial_conditions(phytomer, i)
                    if phytomer.hiddenzone:
                        i = _init_initial_conditions(phytomer.hiddenzone, i)
                    for organ in (phytomer.lamina, phytomer.internode, phytomer.sheath):
                        if organ:
                            for element in (organ.exposed_element, organ.enclosed_element):
                                if element:
                                    i = _init_initial_conditions(element, i)

        logger.info('Initialization of the simulation DONE')

    def run(self):
        """
        Compute turgor pressure driven growth in :attr:`population` over :attr:`delta_t`.

        """
        logger = logging.getLogger(__name__)
        logger.info('Run of Turgor-Growth...')
        
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
        
        if logger.isEnabledFor(logging.DEBUG):
            self.t_offset += self.time_step
            
        logger.info('Run of Turgor-Growth DONE')

    def _update_initial_conditions(self):
        """Update the compartments values in :attr:`initial_conditions` from the compartments values of :attr:`population`.
        """
        for model_object, compartments in self.initial_conditions_mapping.items():
            for compartment_name, compartment_index in compartments.items():
                self.initial_conditions[compartment_index] = getattr(model_object, compartment_name)
                
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

        for plant in self.population.plants:
            for axis in plant.axes:
                for phytomer in axis.phytomers:
                    for organ in (phytomer.lamina, phytomer.internode, phytomer.sheath, phytomer.hiddenzone):
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

        :Parameters:

            - `t` (:class:`float`) - The current t at which we want to compute the derivatives.
              Values of `t` are chosen automatically by :func:`scipy.integrate.solve_ivp`.
              At first call to :meth:`_calculate_all_derivatives` by :func:`scipy.integrate.solve_ivp`,
              `t` = **t_span** [0], where **t_span** is one of the arguments passed to :func:`solve_ivp(fun, t_span, y0,...) <scipy.integrate.solve_ivp>`.
              For each following call to :meth:`_calculate_all_derivatives`, `t` belongs
              to the interval [**t_span** [0], **t_span** [1]].
              
            - `y` (:class:`list`) - The current values of y. 
              At first call to :meth:`_calculate_all_derivatives` by :func:`scipy.integrate.solve_ivp`, `y` = **y0**
              where **y0** is one of the arguments passed to :func:`solve_ivp(fun, t_span, y0,...) <scipy.integrate.solve_ivp>`.
              Then, values of `y` are chosen automatically by :func:`scipy.integrate.solve_ivp`.

        :Returns:
            The derivatives of `y` at `t`.

        :Returns Type:
            :class:`list`


        """
        logger = logging.getLogger(__name__)
        
        if logger.isEnabledFor(logging.DEBUG):
            t_abs = t + self.t_offset
            logger.debug('t = {}'.format(t_abs))

        compartments_logger = logging.getLogger('turgorgrowth.compartments')
        if logger.isEnabledFor(logging.DEBUG) and compartments_logger.isEnabledFor(logging.DEBUG):
            self._log_compartments(t_abs, y, Simulation.LOGGERS_NAMES['compartments'])

        # check that the solver is not crashed
        y_isnan = np.isnan(y)
        if y_isnan.any():
            message = 'The solver did not manage to compute a compartment. See the logs. NaN found in y'
            logger.exception(message)
            raise SimulationRunError(message)

        y_derivatives = np.zeros_like(y)

        for plant in self.population.plants:
            for axis in plant.axes:
                # compute the derivative of each organ compartment
                for phytomer in axis.phytomers:
                    # Hidden zone
                    hiddenzone = phytomer.hiddenzone
                    if hiddenzone:
                        hiddenzone.age = y[self.initial_conditions_mapping[hiddenzone]['age']]
                        hiddenzone.length = y[self.initial_conditions_mapping[hiddenzone]['length']]
                        hiddenzone.length = y[self.initial_conditions_mapping[hiddenzone]['length']]
                        hiddenzone.radius = y[self.initial_conditions_mapping[hiddenzone]['radius']]
                        hiddenzone.turgor_water_potential = y[self.initial_conditions_mapping[hiddenzone]['turgor_water_potential']]
                        hiddenzone.water_content = y[self.initial_conditions_mapping[hiddenzone]['water_content']]

                        #: Hidden zone volume
                        hiddenzone_volume = hiddenzone.calculate_volume(hiddenzone.length, hiddenzone.radius)

                        #: Osmotic water potential
                        hiddenzone.osmotic_water_potential = hiddenzone.calculate_osmotic_water_potential(hiddenzone.sucrose, hiddenzone.amino_acids, hiddenzone.proteins,
                                                                                                          hiddenzone_volume, hiddenzone.temperature)

                        #: Total water potential
                        hiddenzone.total_water_potential = hiddenzone.calculate_water_potential(hiddenzone.turgor_water_potential, hiddenzone.osmotic_water_potential)

                        #: Resistance to water flow
                        hiddenzone_dimensions = {'length': hiddenzone.length, 'radius': hiddenzone.radius}
                        previous_organ = self.mapping_topology['predecessor'][hiddenzone]
                        hiddenzone.resistance = hiddenzone.calculate_resistance(hiddenzone_dimensions, previous_organ)

                    for organ in (phytomer.lamina, phytomer.internode, phytomer.sheath):
                        if organ is None:
                            continue
                        for element in (organ.exposed_element, organ.enclosed_element):
                            if element is None:
                                continue
                            element.age = y[self.initial_conditions_mapping[element]['age']]
                            element.length = y[self.initial_conditions_mapping[element]['length']]
                            element.turgor_water_potential = y[self.initial_conditions_mapping[element]['turgor_water_potential']]
                            element.water_content = y[self.initial_conditions_mapping[element]['water_content']]

                            #: Lamina
                            if organ.label == 'blade':
                                element.length = y[self.initial_conditions_mapping[element]['length']]
                                element.width = y[self.initial_conditions_mapping[element]['width']]
                                element.thickness = y[self.initial_conditions_mapping[element]['thickness']]
        
                                #: Lamina volume
                                element_volume = element.calculate_volume(element.length, element.width, element.thickness)
        
                                #: Osmotic water potential
                                element.osmotic_water_potential = element.calculate_osmotic_water_potential(element.sucrose, element.amino_acids, element.proteins,
                                                                                                            element_volume, element.temperature)
                                
                                #: Total water potential
                                element.total_water_potential = element.calculate_water_potential(element.turgor_water_potential, element.osmotic_water_potential)
        
                                #: Resistance to water flow
                                element.dimensions = {'length': element.length, 'width': element.width, 'thickness': element.thickness}
                                previous_organ = self.mapping_topology['predecessor'][element]
                                element.resistance = element.calculate_resistance(element.dimensions, previous_organ)
        
                            # Internode or Sheath
                            else:
                                element.radius = y[self.initial_conditions_mapping[element]['radius']]
        
                                #: Hidden zone volume
                                element_volume = element.calculate_volume(element.length, element.radius)
        
                                #: Osmotic water potential
                                element.osmotic_water_potential = element.calculate_osmotic_water_potential(element.sucrose, element.amino_acids, element.proteins,
                                                                                                            element_volume, element.temperature)
        
                                #: Total water potential
                                element.total_water_potential = element.calculate_water_potential(element.turgor_water_potential, element.osmotic_water_potential)
        
                                #: Resistance to water flow
                                element.dimensions = {'length': element.length, 'radius': element.radius}
                                previous_organ = self.mapping_topology['predecessor'][element]
                                element.resistance = element.calculate_resistance(element.dimensions, previous_organ)
        #: Water influx
        for plant in self.population.plants:
            for axis in plant.axes:
                # compute the derivative of each organ compartment
                for phytomer in axis.phytomers:
                    # Hidden zone
                    hiddenzone = phytomer.hiddenzone
                    if phytomer.hiddenzone:
                        previous_organ = self.mapping_topology['predecessor'][hiddenzone]
                        hiddenzone.water_influx = hiddenzone.calculate_water_influx(hiddenzone.total_water_potential, previous_organ.total_water_potential, hiddenzone.resistance, self.delta_t)

                    for organ in (phytomer.lamina, phytomer.internode, phytomer.sheath):
                        if organ is None:
                            continue
                        for element in (organ.exposed_element, organ.enclosed_element):
                            if element is None:
                                continue
                            previous_organ = self.mapping_topology['predecessor'][element]
                            element.water_influx = element.calculate_water_influx(element.total_water_potential, previous_organ.total_water_potential, element.resistance, self.delta_t)

        for plant in self.population.plants:
            for axis in plant.axes:
                # compute the derivative of each organ compartment
                for phytomer in axis.phytomers:
                    # Hidden zone
                    hiddenzone = phytomer.hiddenzone
                    if phytomer.hiddenzone:
                        #: Delta water content
                        if hiddenzone in self.mapping_topology['successor'].keys():
                            next_organ = self.mapping_topology['successor'][hiddenzone]
                            delta_water_content = hiddenzone.calculate_delta_water_content(hiddenzone.water_influx, next_organ.water_influx)
                        else:
                            delta_water_content = hiddenzone.calculate_delta_water_content(hiddenzone.water_influx)

                        #: Delta turgor pressure
                        phi = hiddenzone.calculate_extensibility(hiddenzone.age, self.delta_t)
                        hiddenzone_dimensions = {'length': hiddenzone.length, 'radius': hiddenzone.radius}
                        delta_turgor_water_potential = hiddenzone.calculate_delta_turgor_water_potential(phi, hiddenzone.turgor_water_potential, hiddenzone_dimensions, delta_water_content)

                        #: Store compartments' derivative
                        delta_organ_dimensions = hiddenzone.calculate_delta_organ_dimensions(delta_turgor_water_potential, hiddenzone.turgor_water_potential, phi, hiddenzone_dimensions)
                        delta_age = hiddenzone.calculate_growing_degree_days(self.delta_t)
                        y_derivatives[self.initial_conditions_mapping[hiddenzone]['age']] = delta_age
                        y_derivatives[self.initial_conditions_mapping[hiddenzone]['length']] = delta_organ_dimensions['length']
                        y_derivatives[self.initial_conditions_mapping[hiddenzone]['radius']] = delta_organ_dimensions['radius']
                        y_derivatives[self.initial_conditions_mapping[hiddenzone]['turgor_water_potential']] = delta_turgor_water_potential
                        y_derivatives[self.initial_conditions_mapping[hiddenzone]['water_content']] = delta_water_content

                    for organ in (phytomer.lamina, phytomer.internode, phytomer.sheath):
                        if organ is None:
                            continue
                        for element in (organ.exposed_element, organ.enclosed_element):
                            if element is None:
                                continue

                            # Lamina
                            if organ.label == 'blade':
                                delta_water_content = element.calculate_delta_water_content(element.water_influx, element.transpiration, self.delta_t)
                                #: Delta turgor pressure
                                phi = element.calculate_extensibility(element.age, self.delta_t)
                                delta_turgor_water_potential = element.calculate_delta_turgor_water_potential(phi, element.turgor_water_potential, element.dimensions, delta_water_content)

                                #: Store compartments' derivative
                                delta_organ_dimensions = element.calculate_delta_organ_dimensions(delta_turgor_water_potential, element.turgor_water_potential, phi, element.dimensions)
                                delta_age = element.calculate_growing_degree_days(self.delta_t)
                                y_derivatives[self.initial_conditions_mapping[element]['age']] = delta_age
                                y_derivatives[self.initial_conditions_mapping[element]['length']] = delta_organ_dimensions['length']
                                y_derivatives[self.initial_conditions_mapping[element]['width']] = delta_organ_dimensions['width']
                                y_derivatives[self.initial_conditions_mapping[element]['thickness']] = delta_organ_dimensions['thickness']
                                y_derivatives[self.initial_conditions_mapping[element]['turgor_water_potential']] = delta_turgor_water_potential
                                y_derivatives[self.initial_conditions_mapping[element]['water_content']] = delta_water_content

                            else:
                                next_organ = self.mapping_topology['successor'][element]
                                #: Delta water content
                                delta_water_content = element.calculate_delta_water_content(element.water_influx, next_organ.water_influx, element.transpiration, self.delta_t)

                                #: Delta turgor pressure
                                phi = element.calculate_extensibility(element.age, self.delta_t)
                                element_dimensions = {'length': element.length, 'radius': element.radius}
                                delta_turgor_water_potential = element.calculate_delta_turgor_water_potential(phi, element.turgor_water_potential, element_dimensions, delta_water_content)

                                #: Store compartments' derivative
                                delta_organ_dimensions = element.calculate_delta_organ_dimensions(delta_turgor_water_potential, element.turgor_water_potential, phi, element_dimensions)
                                delta_age = element.calculate_growing_degree_days(self.delta_t)
                                y_derivatives[self.initial_conditions_mapping[element]['age']] = delta_age
                                y_derivatives[self.initial_conditions_mapping[element]['length']] = delta_organ_dimensions['length']
                                y_derivatives[self.initial_conditions_mapping[element]['radius']] = delta_organ_dimensions['radius']
                                y_derivatives[self.initial_conditions_mapping[element]['turgor_water_potential']] = delta_turgor_water_potential
                                y_derivatives[self.initial_conditions_mapping[element]['water_content']] = delta_water_content

        derivatives_logger = logging.getLogger('turgorgrowth.derivatives')
        if logger.isEnabledFor(logging.DEBUG) and derivatives_logger.isEnabledFor(logging.DEBUG):
            self._log_compartments(t_abs, y_derivatives, Simulation.LOGGERS_NAMES['derivatives'])

        return y_derivatives
