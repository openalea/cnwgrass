# -*- coding: latin-1 -*-
"""
    turgorgrowth.model
    ~~~~~~~~~~~~~~~~~~

    The module :mod:`turgorgrowth.model` defines the equations of wtaer flow, turgor pressure and growth.

    :license: CeCILL-C, see LICENSE for details.

"""

from __future__ import division  # use "//" to do integer division
import numpy as np
from math import exp
from math import ceil

from turgorgrowth import parameters


class Population(object):
    """
    The class :class:`Population`.

    A :class:`population <Population>` must have at least one :class:`plant <Plant>`.
    """

    PARAMETERS = parameters.POPULATION_PARAMETERS  #: the internal parameters of the population

    def __init__(self, plants=None):
        """
        :param list [Plant] plants: the list of Plant objects
        """
        if plants is None:
            plants = []
        self.plants = plants

    def calculate_aggregated_variables(self):
        """Calculate the integrative variables of the population recursively.
        """
        for plant in self.plants:
            plant.calculate_aggregated_variables()


class Soil(object):
    """
    The class :class:`Soil` defines the soil water potential as function of the soil relative water content.
    """

    PARAMETERS = parameters.SOIL_PARAMETERS  #: the internal parameters of the soil
    INIT_COMPARTMENTS = parameters.SOIL_INIT_COMPARTMENTS  #: the initial values of compartments and state parameters


class Plant(object):
    """
    The class :class:`Plant` defines the water flow at plant scale.

    A :class:`plant <Plant>` must have at least one :class:`axis <Axis>`.
    """

    PARAMETERS = parameters.PLANT_PARAMETERS  #: the internal parameters of the plants

    def __init__(self, index=None, axes=None):
        """
        :param int index: plant index
        :param list [Axis] axes: the list of Axis objects
        """
        self.index = index
        if axes is None:
            axes = []
        self.axes = axes  #: the list of axes
        self.cohorts = []  #: list of cohort values - Hack to treat tillering cases : TEMPORARY

    def calculate_aggregated_variables(self):
        """Calculate the integrative variables of the plant recursively.
        """
        for axis in self.axes:
            axis.calculate_aggregated_variables()


class Axis(object):
    """
    The class :class:`Axis`.

    An :class:`axis <Axis>` must have:
        * one :class:`root compartment <Roots>`,
        * one :class:`xylem <Xylem>`,
        * at least one :class:`phytomer<Phytomer>`.
    """

    PARAMETERS = parameters.AXIS_PARAMETERS  #: the internal parameters of the axes

    def __init__(self, label=None, roots=None, xylem=None, phytomers=None):
        """
        :param str label: the label of the axis
        :param Root roots: Root object
        :param Xylem xylem: Xylem object
        :param list [Phytomer] phytomers: list of Phytomer objects
        """
        self.label = label
        self.roots = roots
        self.xylem = xylem
        if phytomers is None:
            phytomers = []
        self.phytomers = phytomers  #: the list of phytomers

        # integrative variables
        self.Total_Transpiration = None  #: the total transpiration (mmol s-1)
        self.Growth = None  #: the hiddenzone growth (g H2O)
        self.total_water_influx = None  #: water influx between xylem and photosynthetic organ element (g H2O)
        self.xylem_water_potential = None  #: xylem water potential (Mpa)
        self.plant_water_content = None

    def calculate_aggregated_variables(self):
        """Calculate the integrative variables of the axis recursively.
        """
        self.Total_Transpiration = 0
        self.Growth = 0
        self.total_water_influx = 0
        self.xylem_water_potential = 0
        self.plant_water_content = 0

        for phytomer in self.phytomers:
            phytomer.calculate_aggregated_variables()
            self.Total_Transpiration += phytomer.Total_Transpiration
            if phytomer != phytomer.hiddenzone:
                self.total_water_influx += phytomer.total_water_influx
            if phytomer.hiddenzone is not None:
                self.Growth += phytomer.hiddenzone.water_influx * phytomer.hiddenzone.nb_replications

        self.xylem_water_potential = self.xylem.total_water_potential
        self.plant_water_content = self.total_water_influx + self.Growth - self.Total_Transpiration     # TEST


class Phytomer(object):
    """
    The class :class:`Phytomer`.

    A :class:`phytomer <Phytomer>` must have at least:
        * 1 photosynthetic organ: :class:`lamina <Lamina>`, :class:`internode <Internode>`,
                                  or :class:`sheath <Sheath>`.
        * or 1 :class:`hiddenzone <HiddenZone>`.
    """

    PARAMETERS = parameters.PHYTOMER_PARAMETERS  #: the internal parameters of the phytomers
    INIT_COMPARTMENTS = parameters.PHYTOMER_INIT_COMPARTMENTS  #: the initial values of compartments and state parameters

    def __init__(self, index=None, lamina=None, internode=None, sheath=None, hiddenzone=None, cohorts=None, cohorts_replications=None):
        """
        :param int index: index of the phytomer
        :param Lamina lamina: Lamina object
        :param Internode internode: Internode object
        :param Sheath sheath: Sheath object
        :param HiddenZone hiddenzone: HiddenZone object
        """
        self.index = index
        self.lamina = lamina
        self.internode = internode
        self.sheath = sheath
        self.hiddenzone = hiddenzone
        if cohorts is None:
            cohorts = []
        self.cohorts = cohorts  #: list of cohort values - Hack to treat tillering cases : TEMPORARY. Devrait être porté à l'échelle de la plante uniquement mais je ne vois pas comment faire mieux
        self.cohorts_replications = cohorts_replications  #: dictionary of number of replications per cohort rank

        # integrative variables
        self.Tr = None  #: mmol H2O m-2 s-1
        self.green_area = None  #: m2
        self.Total_Transpiration = None #: g H20
        self.water_influx = None    #: g H20
        self.Growth = None  #: g H20
        self.total_water_influx = None  #: g H20

    def calculate_aggregated_variables(self):
        """Calculate the integrative variables of the phytomer recursively.
        """
        self.Total_Transpiration = 0
        self.water_influx = 0
        self.total_water_influx = 0
        for organ_ in (self.lamina, self.internode, self.sheath, self.hiddenzone):
            if organ_ is not None:
                organ_.calculate_aggregated_variables()
                if hasattr(organ_, 'Total_Transpiration'):
                    self.Total_Transpiration += organ_.Total_Transpiration
                    self.total_water_influx += organ_.total_water_influx
                if organ_ == 'hiddenzone':
                    self.water_influx += organ_.water_influx

    @property
    def nb_replications(self):
        return sum(int(v <= self.index) * self.cohorts_replications.get(v, 0) for v in self.cohorts) + 1

class Organ(object):
    """
    The class :class:`Organ`.

    :class:`Organ` is the base class of all organs. DO NOT INSTANTIATE IT.
    """
    PARAMETERS = parameters.ORGAN_PARAMETERS  #: the internal parameters of the organ

    def __init__(self, label):
        """
        :param str label: the label of the organ
        """
        self.label = label

    def initialize(self):
        """Initialize the derived attributes of the organ.
        """
        pass

    def calculate_aggregated_variables(self):
        """Calculate the integrative variables of the organ recursively.
        """
        pass


class Roots(Organ):
    """
    The class :class:`Roots`.
    """

    PARAMETERS = parameters.ROOTS_PARAMETERS  #: the internal parameters of the roots
    INIT_COMPARTMENTS = parameters.ROOTS_INIT_COMPARTMENTS  #: the initial values of compartments and state parameters

    def __init__(self, label='roots'):
        """
        :param str label: root label
        """
        super(Roots, self).__init__(label)


class Xylem(Organ):
    """
       The class :class:`Xylem` defines the water exchanges in a xylem.
       """

    PARAMETERS = parameters.XYLEM_PARAMETERS  #: the internal parameters of the xylem
    INIT_COMPARTMENTS = parameters.XYLEM_INIT_COMPARTMENTS  #: the initial values of compartments and state parameters

    def __init__(self, label='xylem', total_water_potential=INIT_COMPARTMENTS.total_water_potential, soil_water_potential=INIT_COMPARTMENTS.soil_water_potential,
                 total_water_influx=INIT_COMPARTMENTS.total_water_influx, Total_Transpiration=INIT_COMPARTMENTS.Total_Transpiration, Growth=INIT_COMPARTMENTS.Growth):

        super(Xylem, self).__init__(label)

        # state parameters
        self.soil_water_potential = soil_water_potential  #: MPa

        # # fluxes from xylem
        self.total_water_influx = total_water_influx  #: g H2O
        self.Total_Transpiration = Total_Transpiration  #: g H2O
        self.Growth = Growth  #: g H2O

        # other fluxes
        self.total_water_potential = total_water_potential    #: MPa

        # integrative variables
        self.delta_t = 3600     #: the delta t of the simulation (in seconds)


    def calculate_aggregated_variables(self):
        """Calculate the integrative variables of the phytomer recursively.
        """
        self.xylem_water_potential = self.calculate_xylem_water_potential(self.soil_water_potential, self.total_water_influx, self.Growth, self.Total_Transpiration, self.delta_t)

    #:  Model equations for water flux
    @staticmethod
    def calculate_soil_water_potential(SRWC):
        """Total water potential of the xylem
        Equation from Chen et al., 1998 (Genotypic variation in drought tolerance of poplar in relation to abscisic acid)

        :param float SRWC: %

        :return: Total water potential (MPa)
        :rtype: float
        """
        #: TO DO declaration of parameters into parameters.py
        soil_water_potential = - exp((-SRWC + 39.765) / 18.902)
        return soil_water_potential

    @staticmethod
    def calculate_xylem_water_potential(soil_water_potential, total_water_influx, Growth, Total_Transpiration, delta_t):
        """Total water potential of the xylem

        :param float soil_water_potential: MPa
        :param float total_water_influx: g H2O
        :param float Growth: g H2O
        :param float delta_t: time step of the simulation (s)

        :return: Total water potential (MPa)
        :rtype: float
        """
        #: Axial resistance between soil and xylem is a fixed parameter : R_soil

        # xylem_water_potential = soil_water_potential - ((Total_Transpiration * 1E-3 * 18) + Growth) * Xylem.PARAMETERS.R_soil * delta_t
        total_water_potential = soil_water_potential - ((Growth + total_water_influx) * Xylem.PARAMETERS.R_soil * delta_t)

        return total_water_potential


class HiddenZone(Organ):
    """
    The class :class:`HiddenZone`.
    """

    PARAMETERS = parameters.HIDDEN_ZONE_PARAMETERS                #: the internal parameters of the hidden zone
    INIT_COMPARTMENTS = parameters.HIDDEN_ZONE_INIT_COMPARTMENTS  #: the initial values of compartments and state parameters

    def __init__(self, label='hiddenzone', leaf_pseudo_age=INIT_COMPARTMENTS.age, age=INIT_COMPARTMENTS.age, amino_acids=INIT_COMPARTMENTS.amino_acids, proteins=INIT_COMPARTMENTS.proteins, sucrose=INIT_COMPARTMENTS.sucrose,
                 temperature=INIT_COMPARTMENTS.temperature, mstruct=INIT_COMPARTMENTS.mstruct, osmotic_water_potential=INIT_COMPARTMENTS.osmotic_water_potential,
                 total_water_potential=INIT_COMPARTMENTS.total_water_potential, leaf_pseudostem_length=INIT_COMPARTMENTS.leaf_pseudostem_length,
                 leaf_L=INIT_COMPARTMENTS.leaf_L, thickness=INIT_COMPARTMENTS.thickness, width=INIT_COMPARTMENTS.width,
                 turgor_water_potential=INIT_COMPARTMENTS.turgor_water_potential, water_content=INIT_COMPARTMENTS.water_content,
                 Tr=INIT_COMPARTMENTS.Tr, green_area=INIT_COMPARTMENTS.green_area, water_influx=INIT_COMPARTMENTS.water_influx, water_outflow=INIT_COMPARTMENTS.water_outflow, cohorts=None, cohorts_replications=None,
                 leaf_Lmax=INIT_COMPARTMENTS.leaf_Lmax, leaf_is_growing=INIT_COMPARTMENTS.leaf_is_growing, lamina_Lmax=INIT_COMPARTMENTS.lamina_Lmax, sheath_Lmax = INIT_COMPARTMENTS.sheath_Lmax):

        super(HiddenZone, self).__init__(label)

        if cohorts is None:
            cohorts = []
        self.cohorts = cohorts  #: list of cohort values - Hack to treat tillering cases : TEMPORARY. Devrait être porté à l'échelle de la plante uniquement mais je ne vois pas comment faire mieux
        self.cohorts_replications = cohorts_replications  #: dictionary of number of replications per cohort rank

        self.label = label

        # state parameters
        self.amino_acids = amino_acids           #: :math:`\mu mol N
        self.proteins = proteins                 #: :math:`\mu mol N
        self.sucrose = sucrose                   #: :math:`\mu mol C
        self.temperature = temperature           #: °C
        self.leaf_pseudo_age = leaf_pseudo_age   #: °Cd
        self.leaf_L = leaf_L                     #: m
        self.leaf_Lmax = leaf_Lmax               #: m
        self.lamina_Lmax = lamina_Lmax          #: m
        self.sheath_Lmax = sheath_Lmax          #: m
        self.leaf_is_growing = leaf_is_growing   #: -
        self.mstruct = mstruct                   #: g
        self.age = age                          #: °Cd
        self.Tr = Tr                            #: #: mmol H20 m-2 s-1
        self.green_area = green_area            #: m²

        # state variables
        self.water_content = water_content      #: g H2O

        # fluxes from xylem
        self.water_influx = water_influx  #: current flow of water from xylem to hiddenzone integrated over delta t (g H2O)
        self.water_outflow = water_outflow  #: current flow of water from hiddenzone to emerged lamina if any integrated over delta t (g H2O)

        # other fluxes
        # self.initial_volume = None #: m3
        self.osmotic_water_potential = osmotic_water_potential  #: MPa
        self.total_water_potential = total_water_potential    #: MPa
        self.resistance = None  #: resistance of water flux between two organs (MPa s g-1)
        self.extensibility = None    #: MPa-1
        self.turgor_water_potential = turgor_water_potential    #: MPa

        # intermediate variables
        self.length = min(leaf_L, leaf_pseudostem_length)       #: m
        self.leaf_pseudostem_length = leaf_pseudostem_length    #: m
        self.width = width                                      #: m
        self.thickness = thickness                              #: m

    @property
    def nb_replications(self):
        return sum(int(v <= self.index) * self.cohorts_replications.get(v, 0) for v in self.cohorts) + 1

    def calculate_aggregated_variables(self):
        """
        :return:
        """
        self.Growth = self.calculate_delta_water_content(self.water_influx, self.water_outflow)

    #:  Model equations for water flux
    @staticmethod
    def calculate_initial_volume(mstruct):
        """ Hidden zone initial volume calculated from mstruct. This calculation is only performed at t = previous leaf emergence

        :param float mstruct: (g)

        :return: volume (m3), water content (g)
        :rtype: (float, float)
        """
        dry_mass = mstruct / HiddenZone.PARAMETERS.RATIO_MSTRUCT_DM  #: total dry mass (g)
        volume = dry_mass * HiddenZone.PARAMETERS.SLOPE_MASS_VOLUME + HiddenZone.PARAMETERS.OFFSET_MASS_VOLUME  #: m3
        water_content = volume * parameters.RHO_WATER  #: g
        return volume, water_content

    @staticmethod
    def calculate_volume(water_content):
        """ Hidden zone volume, assumed to be proportional to water content.

        :param float water_content: g H2O

        :return: volume (m3)
        :rtype: float
        """
        return water_content / parameters.RHO_WATER

    @staticmethod
    def calculate_osmotic_water_potential(sucrose, amino_acids, proteins, volume, temperature, age):
        """ Osmotic water potential of the organ calculated according to metabolites

        :param float sucrose: µmol C under the form of sucrose
        :param float amino_acids: µmol N under the form of amino acids
        :param float proteins:µmol N under the form of proteins
        :param float volume: (m3)
        :param float temperature: degree Celsius
        :param float age: hiddenzone temperature-compensated age (s)

        :return: Osmotic water potential (MPa)
        :rtype: float
        """
        temperature_K = temperature + parameters.CELSIUS_2_KELVIN

        sucrose = ((sucrose * 1E-6) / parameters.NB_C_SUCROSE) * parameters.VANT_HOFF_SUCROSE
        amino_acids = ((amino_acids * 1E-6) / parameters.AMINO_ACIDS_N_RATIO) * parameters.VANT_HOFF_AMINO_ACIDS
        proteins_actual = proteins * min(1., (age * 5E-6 + 0.1))  # TODO: temp hack to account for the fact that N is mainly Nstruct in young hz
        proteins = ((proteins_actual * 1E-6) / parameters.AMINO_ACIDS_N_RATIO) * parameters.VANT_HOFF_AMINO_ACIDS

        osmotic_water_potential = - parameters.R * temperature_K * ((sucrose + amino_acids + proteins) / volume ) * 1E-6

        return osmotic_water_potential

    @staticmethod
    def calculate_water_potential(turgor_water_potential, osmotic_water_potential):
        """ Total water potential of the organ

        :param float turgor_water_potential: MPa
        :param float osmotic_water_potential: MPa

        :return: Total water potential (MPa)
        :rtype: float
        """
        return turgor_water_potential + osmotic_water_potential

    @staticmethod
    def calculate_hiddenzone_length(leaf_L, leaf_pseudostem_length):
        """ Length of the hidden zone

        :param float leaf_L: Total leaf length (m)
        :param float leaf_pseudostem_length: Length of the pseudostem (m)

        :return: Length of the hidden zone (m)
        :rtype: float
        """
        return min(leaf_L, leaf_pseudostem_length)

    @staticmethod
    def calculate_resistance(hiddenzone_dimensions):
        """ Resistance of water flow between the hiddenzone and the xylem.
        Relations were set proportional to the length and inversely proportional to the area of organ's cross section.

        :param dict hiddenzone_dimensions: dict of hidden zone dimensions at time t. Keys = ['length', 'thickness', 'width] (m)

        :return: resistance (MPa s g-1)
        :rtype: float
        """
        resistance = 0.5 * Xylem.PARAMETERS.R_xylem_hz * (hiddenzone_dimensions['length'] / (hiddenzone_dimensions['width'] * hiddenzone_dimensions['thickness']))

        return resistance

    @staticmethod
    def calculate_water_flux(organ_water_potential, xylem_water_potential, resistance, delta_t):
        """ Water flow into the organ according to water potential gradient with the xylem.

        :param float organ_water_potential: water potential of the current organ (MPa)
        :param float xylem_water_potential: water potential of the xylem (MPa)
        :param float resistance: transport resistance between organ and xylem (MPa s g-1)
        :param float delta_t: time step of the simulation (s)

        :return: Water influx into the current organ integrated over delta_t (g)
        :rtype: float
        """

        water_influx = ((xylem_water_potential - organ_water_potential) / resistance) * delta_t

        return water_influx

    @staticmethod
    def calculate_delta_water_content(water_influx, water_outflow):
        """ delta of water flow for the hidden zone.

        :param float water_flux: Water influx integrated over delta_t (g H2O)
        :param float water_outflow: Water loss through the emerged lamina or sheath if any (g H2O)

        :return: Delta of water flow into the organ (g)
        :rtype: float
        """

        return water_influx - water_outflow

    @staticmethod
    def calculate_extensibility_init(age, delta_t):

        """ Hidden zone extensibility in each dimension in relation to non-reversible dimensional changes.

        :param float age: hidden zone age (°Cd)
        :param float delta_t: time step of the simulation (s)

        :return: Extensibility z and x (MPa-1): {'z': float, 'x': float}
        :rtype: dict
        """
        if age <= HiddenZone.PARAMETERS.tend:
            beta_function_norm = (1 - (1 + (HiddenZone.PARAMETERS.tend - age) / (HiddenZone.PARAMETERS.tend - HiddenZone.PARAMETERS.tmax))
                                  * ((age - HiddenZone.PARAMETERS.tbase) / (HiddenZone.PARAMETERS.tend - HiddenZone.PARAMETERS.tbase)) **
                                  ((HiddenZone.PARAMETERS.tend - HiddenZone.PARAMETERS.tbase) / (HiddenZone.PARAMETERS.tend - HiddenZone.PARAMETERS.tmax))) + HiddenZone.PARAMETERS.OFFSET_LEAF
        else:
            beta_function_norm = 0

        phi = {}

        for phi_init_dimensions, phi_init_value in HiddenZone.PARAMETERS.phi_initial.items():
            phi[phi_init_dimensions] = phi_init_value * beta_function_norm * delta_t
        return phi

    @staticmethod
    def calculate_extensibility_dimensions(age, delta_t):

        """ Hidden zone extensibility in each dimension in relation to non-reversible dimensional changes.

        :param float age: hidden zone age (°Cd)
        :param float delta_t: time step of the simulation (s)

        :return: Extensibility z and x (MPa-1): {'z': float, 'x': float}
        :rtype: dict
        """

        #: TO DO parameters into parameters.py

        phi = {}

        for phi_init_dimensions, phi_init_value in HiddenZone.PARAMETERS.phi_initial.items():

            # FUNCTION 1

            # if phi_init_dimensions == 'x': #: width
            #     if age <= HiddenZone.PARAMETERS.tend:
            #         beta_function_norm = 1 / (1 + exp(3*(age/1E06)))
            #     else:
            #         beta_function_norm = 0
            #     phi[phi_init_dimensions] = phi_init_value * beta_function_norm * delta_t
            # if phi_init_dimensions == 'y': #: thickness
            #     if age <= HiddenZone.PARAMETERS.tend:
            #         beta_function_norm = 1 / (1 + exp(3*(age/1E06)))
            #     else:
            #         beta_function_norm = 0
            #     phi[phi_init_dimensions] = phi_init_value * beta_function_norm * delta_t

            # FUNCTION 2

            if phi_init_dimensions == 'x': #: width
                if age <= HiddenZone.PARAMETERS.tend:
                # if age <= 2E06:
                    beta_function_norm = 1 / (1 + exp(1/1E06 * age))
                else:
                    beta_function_norm = 0
                phi[phi_init_dimensions] = phi_init_value * beta_function_norm * delta_t
            if phi_init_dimensions == 'y': #: thickness
                if age <= HiddenZone.PARAMETERS.tend:
                # if age <= 2E06:
                    beta_function_norm = 1 / (1 + exp(1/1E06 * age))
                else:
                    beta_function_norm = 0
                phi[phi_init_dimensions] = phi_init_value * beta_function_norm * delta_t
            if phi_init_dimensions == 'z': #: length
                if age <= HiddenZone.PARAMETERS.tend:
                # if age <= 2E06:
                    beta_function_norm = (1 - (1 + (HiddenZone.PARAMETERS.tend - age) / (HiddenZone.PARAMETERS.tend - HiddenZone.PARAMETERS.tmax))
                                          * ((age - HiddenZone.PARAMETERS.tbase) / (HiddenZone.PARAMETERS.tend - HiddenZone.PARAMETERS.tbase)) **
                                          ((HiddenZone.PARAMETERS.tend - HiddenZone.PARAMETERS.tbase) / (HiddenZone.PARAMETERS.tend - HiddenZone.PARAMETERS.tmax)))
                                          # + HiddenZone.PARAMETERS.OFFSET_LEAF     #: No offset
                else:
                    beta_function_norm = 0
                phi[phi_init_dimensions] = phi_init_value * beta_function_norm * delta_t

            # # FUNCTION 3
            #
            # if phi_init_dimensions == 'x':  #: width
            #     if age <= 1E09:
            #         beta_function_norm = 1 / (1 + exp(3 / 1E09 * age))
            #     else:
            #         beta_function_norm = 0
            #     phi[phi_init_dimensions] = phi_init_value * beta_function_norm * delta_t
            # if phi_init_dimensions == 'y':  #: thickness
            #     if age <= 1E09:
            #         beta_function_norm = 1 / (1 + exp(3 / 1E09 * age))
            #     else:
            #         beta_function_norm = 0
            #     phi[phi_init_dimensions] = phi_init_value * beta_function_norm * delta_t
            # if phi_init_dimensions == 'z':  #: length
            #     if age <= 1E09:
            #         beta_function_norm = 1 / (1 + exp(1 / 1E09 * age))
            #     else:
            #         beta_function_norm = 0
            #     phi[phi_init_dimensions] = phi_init_value * beta_function_norm * delta_t

        return phi

    @staticmethod
    def calculate_end_elongation(leaf_Lmax, leaf_L):
        """ End of elongation of the leaf in length dimension.

        :param float leaf_Lmax: maximum leaf length (m)
        :param float leaf_L: actual leaf length (m)

        :return: end_elongation
        :rtype: float
        """
        end_elongation = leaf_Lmax - leaf_L
        return end_elongation

    @staticmethod
    def calculate_organ_volume(hiddenzone_dimensions):
        """ HiddenZone volume, assumed to be equal to a box dimensions.

        :param float length: (m)
        :param float width: (m)
        :param float thickness: (m)

        :return: volume (m3)
        :rtype: float
        """
        organ_volume = hiddenzone_dimensions['length'] * hiddenzone_dimensions['width'] * hiddenzone_dimensions['thickness']

        return organ_volume

    @staticmethod
    def calculate_delta_turgor_water_potential(phi, turgor_water_potential, volume, delta_water_content, water_content):
        """ Delta of turgor water potential according to hidden zone water content, turgor water potential, dimensions, elasticity and plasticity.
        Hidden zone geometry is supposed to be a rectangular prism.

        :param dict [str, float] phi: float phi: dict of cell wall extensibility (MPa). Keys = ['x', 'y', 'z]
        :param float turgor_water_potential: MPa
        :param float volume: m3
        :param float delta_water_content: delta water content integrated over delta t (g)
        :param float water_content: water content (g)

        :return: Delta of turgor water potential (MPa)
        :rtype: float
        """
        epsilon_x, epsilon_y, epsilon_z = HiddenZone.PARAMETERS.epsilon['x'], HiddenZone.PARAMETERS.epsilon['y'], HiddenZone.PARAMETERS.epsilon['z']
        elastic_component = (epsilon_x * epsilon_y * epsilon_z) / (epsilon_z * epsilon_x + epsilon_z * epsilon_y + epsilon_x * epsilon_y)   #: Elastic reversible growth (MPa)
        plastic_component = (phi['x'] + phi['y'] + phi['z'])   #: Plastic irreversible growth
        delta_turgor_water_potential = ((1 / (parameters.RHO_WATER * volume)) * delta_water_content - plastic_component * (max(turgor_water_potential, HiddenZone.PARAMETERS.GAMMA) - HiddenZone.PARAMETERS.GAMMA)) * elastic_component  #: (MPa)

        return delta_turgor_water_potential

    @staticmethod
    def calculate_delta_organ_dimensions(delta_turgor_water_potential, turgor_water_potential, phi, organ_dimensions):
        """ Delta of organ dimensions according to turgor water potential, dimensions, extensibility and plasticity.
        Hidden zone geometry is supposed to be a rectangular prism.

        :param float delta_turgor_water_potential: delta of turgor water potential integrated over delta t (MPa)
        :param float turgor_water_potential: MPa
        :param dict phi: dict of cell wall extensibility (MPa). Keys = ['x', 'y', 'z]
        :param dict organ_dimensions: dict of organ dimensions at time t. Keys = [length', 'width', 'thickness'] (m)
        :param float water_content: water content (g)

        :return: Delta of organ specific-dimensions (m). Keys = [leaf_pseudostem_length', 'width', 'thickness']
        :rtype: dict
        """
        delta_organ_dimensions = {}
        epsilon_dict = HiddenZone.PARAMETERS.epsilon
        mapping_dimensions = {'x': 'width', 'y': 'thickness', 'z': 'length'}

        for epsilon_dimension, epsilon_value in epsilon_dict.items():
            if epsilon_dimension != 'z':
                delta_organ_dimensions[mapping_dimensions[epsilon_dimension]] = 0

            else:
                delta_organ_dimensions[mapping_dimensions[epsilon_dimension]] = ((1 / epsilon_value) * delta_turgor_water_potential +
                                                                             phi[epsilon_dimension] * (max(turgor_water_potential, HiddenZone.PARAMETERS.GAMMA) - HiddenZone.PARAMETERS.GAMMA)) *\
                                                                            organ_dimensions[mapping_dimensions[epsilon_dimension]]

            # delta_organ_dimensions[mapping_dimensions[epsilon_dimension]] = ((1 / epsilon_value) * delta_turgor_water_potential +
            #                                                                  phi[epsilon_dimension] * (max(turgor_water_potential, HiddenZone.PARAMETERS.GAMMA) - HiddenZone.PARAMETERS.GAMMA)) *\
            #                                                                 organ_dimensions[mapping_dimensions[epsilon_dimension]]

        return delta_organ_dimensions


class PhotosyntheticOrgan(Organ):
    """
    The class :class:`PhotosyntheticOrgan` defines the water flow in a photosynthetic organ.

    A :class:`photosynthetic organ <PhotosyntheticOrgan>` must have at least 1
    :class:`photosynthetic organ element <PhotosyntheticOrganElement>`:
    :class:`lamina element <LaminaElement>`, :class:`internode element <InternodeElement>`, or :class:`sheath element <SheathElement>`.

    :class:`PhotosyntheticOrgan` is the base class of all photosynthetic organs. DO NOT INSTANTIATE IT.
    """

    def __init__(self, label, exposed_element, enclosed_element):
        """
        :param str label: Photosynthetic organ label
        :param LaminaElement or InternodeElement or SheathElement exposed_element: the exposed element
        :param LaminaElement or InternodeElement or SheathElement enclosed_element: the enclosed element
        """
        super(PhotosyntheticOrgan, self).__init__(label)

        self.exposed_element = exposed_element
        self.enclosed_element = enclosed_element
        self.Tr = None  #: mmol H20 m-2 s-1
        self.green_area = None  #: m2

    def calculate_aggregated_variables(self):
        self.Total_Transpiration = 0
        self.total_water_influx = 0
        for element in (self.exposed_element, self.enclosed_element):
            if element is not None:
                element.calculate_aggregated_variables()
                self.Total_Transpiration += element.Total_Transpiration
                self.total_water_influx += element.total_water_influx


class Lamina(PhotosyntheticOrgan):
    """
    The class :class:`Lamina`.
    """

    def __init__(self, label='lamina', exposed_element=None, enclosed_element=None):
        """
        :param str label: lamina label
        :param LaminaElement exposed_element: the exposed lamina object
        :param LaminaElement enclosed_element: the enclosed lamina object
        """
        super(Lamina, self).__init__(label, exposed_element, enclosed_element)


class Internode(PhotosyntheticOrgan):
    """
    The class :class:`Internode`.
    """

    def __init__(self, label=None, exposed_element=None, enclosed_element=None):
        """
        :param str label: Internode label
        :param InternodeElement exposed_element: the exposed internode object
        :param InternodeElement enclosed_element: the enclosed internode object
        """
        super(Internode, self).__init__(label, exposed_element, enclosed_element)


class Sheath(PhotosyntheticOrgan):
    """
    The class :class:`Sheath`.
    """

    def __init__(self, label=None, exposed_element=None, enclosed_element=None):
        """
        :param str label: Sheath label
        :param SheathElement exposed_element: the exposed sheath object
        :param SheathElement enclosed_element: the enclosed sheath object
        """
        super(Sheath, self).__init__(label, exposed_element, enclosed_element)


class PhotosyntheticOrganElement(object):
    """
    The class :class:`PhotosyntheticOrganElement` defines the water flow in a photosynthetic organ element.

    An element must belong to an organ of the same type (e.g. a class:`LaminaElement` must belong to a class:`Lamina`).

    A :class:`photosynthetic organ element <PhotosyntheticOrganElement>` must have at least 1
    :class:`lamina element<LaminaElement>`,
    :class:`internode element <InternodeElement>`,
    or :class:`sheath element <SheathElement>`.

    :class:`PhotosyntheticOrganElement` is the base class of all photosynthetic organ elements. DO NOT INSTANTIATE IT.
    """

    PARAMETERS = parameters.PHOTOSYNTHETIC_ORGAN_ELEMENT_PARAMETERS  #: the internal parameters of the photosynthetic organs elements
    INIT_COMPARTMENTS = parameters.PHOTOSYNTHETIC_ORGAN_ELEMENT_INIT_COMPARTMENTS  #: the initial values of compartments and state parameters

    def __init__(self, label=None, temperature = INIT_COMPARTMENTS.temperature, age = INIT_COMPARTMENTS.age, green_area = INIT_COMPARTMENTS.green_area, mstruct = INIT_COMPARTMENTS.mstruct, Ts = INIT_COMPARTMENTS.Ts,
                 Tr = INIT_COMPARTMENTS.Tr, length=INIT_COMPARTMENTS.length, sucrose = INIT_COMPARTMENTS.sucrose, amino_acids = INIT_COMPARTMENTS.amino_acids, proteins = INIT_COMPARTMENTS.proteins,
                 osmotic_water_potential = INIT_COMPARTMENTS.osmotic_water_potential, total_water_potential = INIT_COMPARTMENTS.total_water_potential,
                 turgor_water_potential = INIT_COMPARTMENTS.turgor_water_potential, water_content = INIT_COMPARTMENTS.water_content,
                 water_influx = INIT_COMPARTMENTS.water_influx, thickness= INIT_COMPARTMENTS.thickness, width = INIT_COMPARTMENTS.width, cohorts=None, cohorts_replications=None):

        self.label = label                                      #: the label of the element
        if cohorts is None:  #: list of cohort values - Hack to treat tillering cases : TEMPORARY. Devrait être porté à l'échelle de la plante uniquement mais je ne vois pas comment faire mieux
            cohorts = []
        self.cohorts = cohorts  #: list of cohort values - Hack to treat tillering cases : TEMPORARY. Devrait être porté à l'échelle de la plante uniquement mais je ne vois pas comment faire mieux
        self.cohorts_replications = cohorts_replications  #: dictionary of number of replications per cohort rank

        # state parameters
        self.age = age                                          #: °Cd
        self.amino_acids = amino_acids                          #: :math:`\mu mol N
        self.green_area = green_area                            #: m2
        self.mstruct = mstruct                                  #: g
        self.proteins = proteins                                #: :math:`\mu mol N
        self.sucrose = sucrose                                  #: :math:`\mu mol C
        self.Ts = Ts                                            #: °C
        self.temperature = temperature                          #: °C
        self.Tr = Tr                                            #: mmol H20 m-2 s-1
        self.thickness = thickness  #: m
        self.width = width          #: m

        # intermediate variables
        self.turgor_water_potential = turgor_water_potential    #: MPa
        self.osmotic_water_potential = osmotic_water_potential  #: MPa
        self.total_water_potential = total_water_potential #: MPa
        self.length = length  #: m

        # state variables
        self.water_content = water_content                     #: g H2O

        # fluxes to xylem
        self.water_influx = water_influx  #: current flow of water from xylem to organ integrated over delta t (g H2O)

        # other fluxes
        self.resistance = None  #: resistance of water flux between two organs (MPa s g-1)

        # Integrated variables
        self.delta_t = 3600     #: the delta t of the simulation (in seconds)

    @property
    def nb_replications(self):
        return sum(int(v <= self.index) * self.cohorts_replications.get(v, 0) for v in self.cohorts) + 1

    def calculate_aggregated_variables(self):
        """Calculate the integrative variables of the element.
        """

        self.Total_Transpiration = self.calculate_Total_Transpiration(self.Tr, self.green_area, self.delta_t)
        self.total_water_influx = self.calculate_total_water_influx(self.water_influx)

    # VARIABLES
    @staticmethod
    def calculate_total_water_influx(water_influx):
        """
        Water influx from xylem to organ
        :param float water_influx: Water influx (g H2O)
        :return: Total water influx (g H2O)
        """

        total_water_influx = water_influx

        return total_water_influx

    @staticmethod
    def calculate_Total_Transpiration(Tr, green_area, delta_t):
        """Surfacic transpiration rate of an element

        :param float Tr: Transpiration rate (mmol H2O m-2 s-1)
        :param float green_area: Green area (m2)
        :param float delta_t: time step of the simulation (s)

        :return: Total transpiration (g H2O)
        :rtype: float
        """

        conversion_ratio = 1E-3 * parameters.WATER_MOLAR_MASS * delta_t # en g H2O

        # return Tr * green_area
        return Tr * green_area * conversion_ratio

    # FLUXES
    #: Water flow equations common to all photosynthetic organ elements

    @staticmethod
    def calculate_volume(water_content):
        """ Photosynthetic element volume, assumed to be proportional to water content.

        :param float water_content: (g H2O)

        :return: volume (m3)
        :rtype: float
        """
        return water_content / parameters.RHO_WATER

    @staticmethod
    def calculate_organ_volume(organ_dimensions):
        """ Photosynthetic element volume, assumed to be equal to a box dimensions.

        :param float length: (m)
        :param float width: (m)
        :param float thickness: (m)

        :return: volume (m3)
        :rtype: float
        """
        organ_volume = organ_dimensions['length'] * organ_dimensions['width'] * organ_dimensions['thickness']
        return organ_volume

    def calculate_initial_water_content(self, hiddenzone_osmotic_water_potential, sucrose, amino_acids, proteins, temperature):
        """ Initial water content of photosynthetic element at emergence
        calculated in order that the resulting osmotic pressure of the element will be similar to that of the hidden zone

        :param float hiddenzone_osmotic_water_potential: the osmotic pressure of the hidden zone (MPa)
        :param float sucrose: amount of sucrose in the element (µmol C)
        :param float amino_acids: amount of amino acids in the element (µmol N)
        :param float proteins: amount of proteins in the element (µmol N)
        :param float temperature: element temperature (°C)

        :return: element initial water content (g H2O)
        :rtype: float
        """
        temperature_K = temperature + parameters.CELSIUS_2_KELVIN

        sucrose = ((sucrose * 1E-6) / parameters.NB_C_SUCROSE) * parameters.VANT_HOFF_SUCROSE
        amino_acids = ((amino_acids * 1E-6) / parameters.AMINO_ACIDS_N_RATIO) * parameters.VANT_HOFF_AMINO_ACIDS
        proteins = ((proteins * 1E-6) / parameters.AMINO_ACIDS_N_RATIO) * parameters.VANT_HOFF_AMINO_ACIDS

        initial_water_content = (- parameters.R * temperature_K * (sucrose + amino_acids + proteins)) / hiddenzone_osmotic_water_potential

        return initial_water_content

    @staticmethod
    def calculate_osmotic_water_potential(sucrose, amino_acids, proteins, volume, temperature):
        """ Osmotic water potential of the organ calculated according to metabolites

        :param float sucrose: µmol C under the form of sucrose
        :param float amino_acids: µmol N under the form of amino acids
        :param float proteins:µmol N under the form of proteins
        :param float volume: (m3)
        :param float temperature: degree Celsius

        :return: Osmotic water potential (MPa)
        :rtype: float
        """
        temperature_K = temperature + parameters.CELSIUS_2_KELVIN

        sucrose = ((sucrose * 1E-6) / parameters.NB_C_SUCROSE) * parameters.VANT_HOFF_SUCROSE
        amino_acids = ((amino_acids * 1E-6) / parameters.AMINO_ACIDS_N_RATIO) * parameters.VANT_HOFF_AMINO_ACIDS
        non_chloroplastic_proteins = proteins / 1
        proteins = ((non_chloroplastic_proteins * 1E-6) / parameters.AMINO_ACIDS_N_RATIO) * parameters.VANT_HOFF_AMINO_ACIDS

        osmotic_water_potential = - parameters.R * temperature_K * ((sucrose + amino_acids + proteins) / volume) * 1E-6
        return osmotic_water_potential

    @staticmethod
    def calculate_resistance(organ_dimensions, R_xylem):
        """
        Resistance of water flow between the lamina and xylem
        Relations were set proportional to the length and inversely proportional to the area of organ's cross section.

        :param dict lamina_dimensions: dict of lamina dimensions (m). Keys = ['length', 'width', 'thickness']

        :return: resistance (MPa s g-1)
        :rtype: float
        """
        resistance = 0.5 * Xylem.PARAMETERS.R_xylem_organ * organ_dimensions['length'] / (organ_dimensions['width'] * organ_dimensions['thickness'])

        return resistance

    @staticmethod
    def calculate_water_potential(turgor_water_potential, osmotic_water_potential):
        """ Total water potential of the organ

        :param float turgor_water_potential: MPa
        :param float osmotic_water_potential: MPa

        :return: Total water potential (MPa)
        :rtype: float
        """
        return turgor_water_potential + osmotic_water_potential

    @staticmethod
    def calculate_water_flux(total_water_potential, xylem_water_potential, resistance, delta_t):
        """ Water flow into the organ according to water potential gradient with the xylem.

        :param float organ_water_potential: water potential of the current organ (MPa)
        :param float xylem_water_potential: water potential of the xylem (MPa)
        :param float resistance: transport resistance between organ and xylem (MPa s g-1)
        :param float delta_t: time step of the simulation (s)

        :return: Water influx into the current organ integrated over delta_t (g H2O)
        :rtype: float
        """
        return ((xylem_water_potential - total_water_potential) / resistance) * delta_t

    @staticmethod
    def calculate_delta_water_content(water_influx, Total_Transpiration, delta_t):
        """ Delta of water flow for the lamina.

        :param float water_flux: Water influx from xylem integrated over delta_t (g)
        :param float Tr: Lamina surfacic transpiration rate (mmol H20 m-2 s-1)
        :param float green_area: Lamina sgreen area (m2)
        :param float delta_t: Time step of the simulation (s)

        :return: Delta of water flow into the organ (g H2O)
        :rtype: float
        """
        # return water_influx - (Total_Transpiration * parameters.WATER_MOLAR_MASS * 1E-3 * delta_t) # transpiration en g
        return water_influx - Total_Transpiration # transpiration en g H2O

    @staticmethod
    def calculate_delta_turgor_water_potential(organ_dimensions, volume, delta_water_content):
        """ Delta of turgor water potential according to organ water content, turgor water potential, dimensions and elasticity.
        Extensibility (psi) is supposed to be 0 as this tissue is mature (growth completed).

        :param dict organ_dimensions: dict of organ dimensions at time t. Keys = ['length', 'thickness', 'width'] (m)
        :param float delta_water_content: delta water content integrated over delta t (g)

        :return: Delta of turgor water potential (MPa)
        :rtype: float
        """
        epsilon_z, epsilon_x, epsilon_y = PhotosyntheticOrganElement.PARAMETERS.epsilon['z'], PhotosyntheticOrganElement.PARAMETERS.epsilon['x'], PhotosyntheticOrganElement.PARAMETERS.epsilon['y']
        elastic_component = (epsilon_z * epsilon_x * epsilon_y) / (epsilon_z * epsilon_x + epsilon_z * epsilon_y + epsilon_x * epsilon_y)  #: Elastic reversible growth (MPa)
        plastic_component = 0
        # organ_volume = organ_dimensions['length'] * organ_dimensions['width'] * organ_dimensions['thickness']  #: (m3)
        delta_turgor_water_potential = ((1 / (parameters.RHO_WATER * volume)) * delta_water_content - plastic_component) * elastic_component  #: (MPa)

        return delta_turgor_water_potential

    @staticmethod
    def calculate_delta_organ_dimensions(delta_turgor_water_potential, organ_dimensions):
        """Delta of lamina dimensions according to turgor water potential, dimensions, and elasticity

        :param float delta_turgor_water_potential: delta of turgor water potential integrated over delta t (MPa)
        :param dict organ_dimensions: dict of organ dimensions at time t. Keys = ['length', 'thickness', 'width'] (m)
        :return: Delta of organ specific-dimensions (m). Keys = ['length', 'width', 'thickness']
        :rtype: dict
        """
        delta_organ_dimensions = {}
        epsilon_dict = PhotosyntheticOrganElement.PARAMETERS.epsilon.items()
        mapping_dimensions = {'x': 'width', 'y': 'thickness', 'z': 'length'}

        for epsilon_dimension, epsilon_value in epsilon_dict:
            delta_organ_dimensions[mapping_dimensions[epsilon_dimension]] = ((1 / epsilon_value) * delta_turgor_water_potential) * organ_dimensions[mapping_dimensions[epsilon_dimension]]

        return delta_organ_dimensions


class LaminaElement(PhotosyntheticOrganElement):
    """
    The class :class:`LaminaElement`.
    """

    PARAMETERS = parameters.LAMINA_ELEMENT_PARAMETERS                   #: the internal parameters of the lamina
    INIT_COMPARTMENTS = parameters.LAMINA_ELEMENT_INIT_COMPARTMENTS     #: the initial values of compartments and state parameters

    def __init__(self, thickness=INIT_COMPARTMENTS.thickness, width=INIT_COMPARTMENTS.width, **kwargs):
        super(LaminaElement, self).__init__(**kwargs)

        self.thickness = thickness  #: m
        self.width = width          #: m


class InternodeElement(PhotosyntheticOrganElement):
    """
    The class :class:`InternodeElement`.
    """

    PARAMETERS = parameters.INTERNODE_ELEMENT_PARAMETERS                           #: the internal parameters of the internode
    INIT_COMPARTMENTS = parameters.INTERNODE_ELEMENT_INIT_COMPARTMENTS             #: the initial values of compartments and state parameters

    def __init__(self, thickness=INIT_COMPARTMENTS.thickness, width=INIT_COMPARTMENTS.width, **kwargs):
        super(InternodeElement, self).__init__(**kwargs)

        self.thickness = thickness  #: m
        self.width = width          #: m


class SheathElement(PhotosyntheticOrganElement):
    """
    The class :class:`SheathElement`.
    """

    PARAMETERS = parameters.SHEATH_ELEMENT_PARAMETERS                   #: the internal parameters of the sheath
    INIT_COMPARTMENTS = parameters.SHEATH_ELEMENT_INIT_COMPARTMENTS     #: the initial values of compartments and state parameters

    def __init__(self, thickness=INIT_COMPARTMENTS.thickness, width=INIT_COMPARTMENTS.width, **kwargs):
        super(SheathElement, self).__init__(**kwargs)

        self.thickness = thickness  #: m
        self.width = width          #: m
