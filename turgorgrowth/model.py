# -*- coding: latin-1 -*-
"""
    turgorgrowth.model
    ~~~~~~~~~~~~~~~~~~

    The module :mod:`turgorgrowth.model` defines the equations of wtaer flow, turgor pressure and growth.

    :license: CeCILL-C, see LICENSE for details.

"""

from __future__ import division  # use "//" to do integer division
import numpy as np
from math import exp, log10, log
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
    INIT_COMPARTMENTS = parameters.AXIS_INIT_COMPARTMENTS  #: the initial values of compartments and state parameters

    def __init__(self, label=None, roots=None, xylem=None, phytomers=None, mstruct=INIT_COMPARTMENTS.mstruct):
        """
        :param str label: the label of the axis
        :param Root roots: Root object
        :param Xylem xylem: Xylem object
        :param list [Phytomer] phytomers: list of Phytomer objects
        :param float mstruct: structural mass (g) of the axis
        """
        self.label = label
        self.roots = roots
        self.xylem = xylem
        if phytomers is None:
            phytomers = []
        self.phytomers = phytomers  #: the list of phytomers

        # state parameters
        self.mstruct = mstruct

        # integrative variables
        self.Total_Transpiration_turgor = None  #: the total transpiration (mmol s-1)
        self.Growth = None  #: the hiddenzone growth (g H2O)
        self.total_water_influx = None  #: water influx between xylem and photosynthetic organ element (g H2O)
        self.xylem_water_potential = None  #: xylem water potential (Mpa)
        self.plant_water_content = None #: plant water content (g H2O)
        self.plant_WC_DM = None #: relative plant water content per structural mass of the shoot axis (%)
        self.green_area_rep = None  #: gree area (m2)
        self.green_area = None  #: gree area (m2)
        self.LAI_turgor = None  #: leaf area index
        self.delta_plant_water_content = None   #: plant delta water content (g H2O)

    def calculate_aggregated_variables(self):
        """Calculate the integrative variables of the axis recursively.
        """
        self.Total_Transpiration_turgor = 0
        self.Growth = 0
        self.total_water_influx = 0
        self.xylem_water_potential = 0
        self.plant_water_content = 0
        self.plant_WC_DM = 0
        self.green_area_rep = 0
        self.green_area = 0
        self.LAI_turgor = 0
        self.delta_plant_water_content = 0

        for phytomer in self.phytomers:
            phytomer.calculate_aggregated_variables()
            self.Total_Transpiration_turgor += phytomer.Total_Transpiration_turgor * phytomer.nb_replications
            self.green_area_rep += phytomer.green_area * phytomer.nb_replications
            self.green_area += phytomer.green_area * phytomer.nb_replications
            self.LAI_turgor = self.green_area_rep * parameters.PLANT_PARAMETERS.plant_density
            self.plant_water_content += phytomer.water_content * phytomer.nb_replications
            if phytomer != phytomer.hiddenzone:
                self.total_water_influx += phytomer.total_water_influx
            if phytomer.hiddenzone is not None:
                self.Growth += phytomer.hiddenzone.water_influx * phytomer.hiddenzone.nb_replications

        self.xylem_water_potential = self.xylem.total_water_potential
        self.delta_plant_water_content += (self.total_water_influx + self.Growth) - self.Total_Transpiration_turgor
        self.plant_WC_DM = self.plant_water_content / (self.mstruct + self.plant_water_content) * 100

    @staticmethod
    def calculate_ratio_WC_mstruct(plant_water_content, mstruct):
        """
        Ratio between water content and structural mass of the axis

        :param float plant_water_content: g
        :param float mstruct: g

        :return: Water content : Structural mass ratio (%)
        :rtype: float
        """
        plant_WC_DM = plant_water_content / mstruct * 100

        return plant_WC_DM

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
        self.index = index  #: the index of the phytomer
        self.lamina = lamina  #: the lamina
        self.internode = internode  #: the internode
        self.sheath = sheath  #: the sheath
        self.hiddenzone = hiddenzone  #: the hidden zone
        if cohorts is None:
            cohorts = []
        self.cohorts = cohorts  #: list of cohort values - Hack to treat tillering cases : TEMPORARY. Devrait être porté à l'échelle de la plante uniquement mais je ne vois pas comment faire mieux
        self.cohorts_replications = cohorts_replications  #: dictionary of number of replications per cohort rank

        # integrative variables
        self.green_area = None  #: m2
        self.Total_Transpiration_turgor = None #: g H20
        self.water_influx = None    #: g H20
        self.Growth = None  #: g H20
        self.total_water_influx = None  #: g H20
        self.water_content = None   #: g H2O

    def calculate_aggregated_variables(self):
        """Calculate the integrative variables of the phytomer recursively.
        """
        self.Total_Transpiration_turgor = 0
        self.water_influx = 0
        self.total_water_influx = 0
        self.green_area = 0
        self.water_content = 0
        for organ_ in (self.lamina, self.internode, self.sheath, self.hiddenzone):
            if organ_ is not None:
                organ_.calculate_aggregated_variables()
                self.green_area += organ_.green_area
                self.water_content += organ_.water_content
                if hasattr(organ_, 'Total_Transpiration_turgor'):
                    self.Total_Transpiration_turgor += organ_.Total_Transpiration_turgor
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

    def __init__(self, label='xylem', total_water_potential=INIT_COMPARTMENTS.total_water_potential, soil_water_potential=INIT_COMPARTMENTS.soil_water_potential):

        super(Xylem, self).__init__(label)

        # state parameters
        self.soil_water_potential = soil_water_potential  #: MPa

        # other fluxes
        self.total_water_potential = total_water_potential    #: MPa

        # integrative variables
        self.delta_t = 3600     #: the delta t of the simulation (in seconds)


    def calculate_aggregated_variables(self):
        """Calculate the integrative variables of the phytomer recursively.
        """
        self.xylem_water_potential = self.calculate_xylem_water_potential(self.soil_water_potential, self.total_water_influx, self.Growth, self.Total_Transpiration_turgor, self.delta_t)

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
        # soil_water_potential = - exp((-SRWC + 39.765) / 18.902)

        # Test calibration
        soil_water_potential = - exp((-SRWC + (39.765 / 2.5)) / 18.902)

        return soil_water_potential

    @staticmethod
    def calculate_soil_conductivity(LAI_turgor):
        """Axial conductivity between soil and xylem as function of LAI

        :param float LAI_turgor

        :return: Axial conductivity between soil and xylem g s-1 Mpa -1
        :rtype: float
        """
        # Implementation of a soil conductivity as function of LAI
        Ksoil = 0.0015 * LAI_turgor + 0.032

        return Ksoil

    @staticmethod
    def calculate_xylem_water_potential(soil_water_potential, total_water_influx, Growth, Ksoil, delta_t):
        """Total water potential of the xylem

        :param float soil_water_potential: MPa
        :param float total_water_influx: g H2O
        :param float Growth: g H2O
        :param float delta_t: time step of the simulation (s)

        :return: Total water potential (MPa)
        :rtype: float
        """

        #: TODO: doc et paramètres fonctions

        #: Axial resistance between soil and xylem is a fixed parameter : R_soil

        # xylem_water_potential = soil_water_potential - ((Total_Transpiration * 1E-3 * 18) + Growth) * Xylem.PARAMETERS.R_soil * delta_t
        total_water_potential = soil_water_potential - ((Growth + total_water_influx) * Xylem.PARAMETERS.R_soil * delta_t)

        # Martre et al. (2001)
        Ksoil = ((21.2 + 7.2) * 18 * 1E-03) * 0.2
        # total_water_potential = soil_water_potential - ((Growth + total_water_influx) * (1E-04 / Ksoil * delta_t))
        # total_water_potential = soil_water_potential - ((Growth + total_water_influx) * (1E-02 / Ksoil * delta_t))

        return total_water_potential


class HiddenZone(Organ):
    """
    The class :class:`HiddenZone`.
    """

    PARAMETERS = parameters.HIDDEN_ZONE_PARAMETERS                #: the internal parameters of the hidden zone
    INIT_COMPARTMENTS = parameters.HIDDEN_ZONE_INIT_COMPARTMENTS  #: the initial values of compartments and state parameters

    def __init__(self, label='hiddenzone', fructan=INIT_COMPARTMENTS.fructan, leaf_enclosed_mstruct=INIT_COMPARTMENTS.leaf_enclosed_mstruct, leaf_pseudo_age=INIT_COMPARTMENTS.leaf_pseudo_age, hiddenzone_age=INIT_COMPARTMENTS.hiddenzone_age, amino_acids=INIT_COMPARTMENTS.amino_acids, proteins=INIT_COMPARTMENTS.proteins, sucrose=INIT_COMPARTMENTS.sucrose,
                 width_prev = INIT_COMPARTMENTS.width_prev, thickness_prev = INIT_COMPARTMENTS.thickness_prev,
                 temperature=INIT_COMPARTMENTS.temperature, mstruct=INIT_COMPARTMENTS.mstruct, osmotic_water_potential=INIT_COMPARTMENTS.osmotic_water_potential,
                 total_water_potential=INIT_COMPARTMENTS.total_water_potential, leaf_pseudostem_length=INIT_COMPARTMENTS.leaf_pseudostem_length,
                 leaf_L=INIT_COMPARTMENTS.leaf_L, thickness=INIT_COMPARTMENTS.thickness, width=INIT_COMPARTMENTS.width, omega=INIT_COMPARTMENTS.omega,
                 turgor_water_potential=INIT_COMPARTMENTS.turgor_water_potential, water_content=INIT_COMPARTMENTS.water_content,
                 Tr=INIT_COMPARTMENTS.Tr, green_area=INIT_COMPARTMENTS.green_area, water_influx=INIT_COMPARTMENTS.water_influx, water_outflow=INIT_COMPARTMENTS.water_outflow, cohorts=None, cohorts_replications=None, leaf_Wmax = INIT_COMPARTMENTS.leaf_Wmax, lamina_Lmax = INIT_COMPARTMENTS.lamina_Lmax,
                 leaf_is_growing=INIT_COMPARTMENTS.leaf_is_growing, delta_teq = INIT_COMPARTMENTS.delta_teq):

        super(HiddenZone, self).__init__(label)

        if cohorts is None:
            cohorts = []
        self.cohorts = cohorts  #: list of cohort values - Hack to treat tillering cases : TEMPORARY. Devrait être porté à l'échelle de la plante uniquement mais je ne vois pas comment faire mieux
        self.cohorts_replications = cohorts_replications  #: dictionary of number of replications per cohort rank

        self.label = label

        # state parameters
        self.fructan = fructan                   #: :math:`\mu mol C
        self.amino_acids = amino_acids           #: :math:`\mu mol N
        self.proteins = proteins                 #: :math:`\mu mol N
        self.sucrose = sucrose                   #: :math:`\mu mol C
        self.temperature = temperature           #: °C
        self.leaf_pseudo_age = leaf_pseudo_age   #: °Cd
        self.leaf_L = leaf_L                     #: m
        self.leaf_is_growing = leaf_is_growing   #: -
        self.mstruct = mstruct                   #: g
        self.leaf_enclosed_mstruct = leaf_enclosed_mstruct                   #: g
        self.hiddenzone_age = hiddenzone_age                          #: °Cd
        self.Tr = Tr                            #: #: mmol H20 m-2 s-1
        self.green_area = green_area            #: m²
        self.delta_teq = delta_teq
        self.length = min(leaf_L, leaf_pseudostem_length)       #: m
        self.leaf_pseudostem_length = leaf_pseudostem_length    #: m
        self.width = width                                      #: m
        self.thickness = thickness                              #: m
        self.width_prev = width_prev   #: m
        self.thickness_prev = thickness_prev   #: m
        self.water_content = water_content      #: g H2O

        # fluxes from xylem
        self.water_influx = water_influx  #: current flow of water from xylem to hiddenzone integrated over delta t (g H2O)
        self.water_outflow = water_outflow  #: current flow of water from hiddenzone to emerged lamina if any integrated over delta t (g H2O)

        # other fluxes
        self.initial_volume = None #: m3
        self.osmotic_water_potential = osmotic_water_potential  #: MPa
        self.total_water_potential = total_water_potential    #: MPa
        self.resistance = None  #: resistance of water flux between two organs (MPa s g-1)
        self.extensibility = None    #: MPa-1
        self.turgor_water_potential = turgor_water_potential    #: MPa
        self.leaf_Lmax = None   #: m

        # intermediate variables
        self.omega = omega    #: -


    @property
    def nb_replications(self):
        return sum(int(v <= self.index) * self.cohorts_replications.get(v, 0) for v in self.cohorts) + 1

    def calculate_aggregated_variables(self):
        """
        :return:
        """
        self.Growth = self.calculate_delta_water_content(self.water_influx, self.water_outflow)

    #:  Model equations for water flux
    # @staticmethod
    # def calculate_vstorage(leaf_pseudo_age, vstorage):
    #
    #     """ Storage portion of the Hidden zone volume. The left over is assumed to be xylem vessels.
    #
    #     :param float leaf_pseudo_age: Hidden zone leaf pseudo age (°Cd)
    #
    #     :return: vstorage
    #     :rtype: float
    #     """
    #
    #     if leaf_pseudo_age <= HiddenZone.PARAMETERS.tend:
    #         vstorage = 2 / (1 + exp((0.045 * leaf_pseudo_age / 216000))) #: v2
    #     else:
    #         vstorage = vstorage
    #
    #     return vstorage


    @staticmethod
    def calculate_initial_volume(mstruct, leaf_L, thickness, width):
        """ Hidden zone initial volume calculated from mstruct. This calculation is only performed at t = previous leaf emergence

        :param float mstruct: (g)

        :return: volume (m3), water content (g)
        :rtype: (float, float)
        """
        dry_mass = mstruct / HiddenZone.PARAMETERS.RATIO_MSTRUCT_DM  #: total dry mass (g)
        volume = dry_mass * HiddenZone.PARAMETERS.SLOPE_MASS_VOLUME + HiddenZone.PARAMETERS.OFFSET_MASS_VOLUME  #: m3

        # volume = leaf_L * thickness * width

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
    def calculate_solutes_contribution(fructan, sucrose, amino_acids, volume):
        """
        #: TODO: doc
        """

        solutes = (sucrose + fructan + amino_acids) / (volume * parameters.RHO_WATER)
        # omega = solutes * 0.00043 + 0.006       # Po = -0.7 MPa
        # omega = solutes * 0.0004 + 0.0064       # Po = -0.65 MPa
        # omega = solutes * 0.0003 + 0.0052       # Po = -0.8 MPa
        # omega = solutes * 0.0004 + 0.0056       # Po = -0.75 MPa

        # omega = solutes * 0.000375 + 0.015       # Po = -0.8 MPa   bis

        # SIGMOID
        # omega = 0.033 / (0.033 + exp(- solutes / 220))
        omega = 0.03 / (0.03 + exp(- solutes / 172))

        return omega


    @staticmethod
    def calculate_osmotic_water_potential(fructan, sucrose, amino_acids, volume, temperature, omega):
        """ Osmotic water potential of the organ calculated according to metabolites

        :param float fructan: µmol C under the form of fructan
        :param float sucrose: µmol C under the form of sucrose
        :param float amino_acids: µmol N under the form of amino acids
        :param float volume: (m3)
        :param float temperature: air temperature (°C)
        :param float omega: TODO doc

        :return: Osmotic water potential (MPa)
        :rtype: float
        """
        temperature_K = temperature + parameters.CELSIUS_2_KELVIN

        # sucrose = (sucrose * 1E-6) / (parameters.NB_C_SUCROSE * DP) * parameters.VANT_HOFF_SUCROSE
        # amino_acids = ((amino_acids * 1E-6) / parameters.AMINO_ACIDS_N_RATIO) * parameters.VANT_HOFF_AMINO_ACIDS
        # fructan = (fructan * 1E-6) / (parameters.NB_C_SUCROSE * DP) * parameters.VANT_HOFF_SUCROSE

        #: update contribution of solutes - 07.2024
        sucrose = (sucrose * 1E-6) / (parameters.NB_C_SUCROSE)
        amino_acids = ((amino_acids * 1E-6) / parameters.AMINO_ACIDS_N_RATIO)
        fructan = (fructan * 1E-6) / (parameters.NB_C_SUCROSE)
        osmotic_water_potential = - parameters.R * temperature_K * (fructan + sucrose + amino_acids) / (omega * volume * parameters.RHO_WATER * parameters.VSTORAGE)

        # osmotic_water_potential = -0.8

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

        #: TODO : doc
        """
        #: Coussement et al. (2018)
        resistance = 0.5 * Xylem.PARAMETERS.R_xylem_hz * (hiddenzone_dimensions['length'] / (hiddenzone_dimensions['width'] * hiddenzone_dimensions['thickness']))

        #: Martre et al. (2001)
        Kr = 21.2 * 18 * 1E-03
        # resistance = 0.5 * (1E-04 / Kr) * hiddenzone_dimensions['length'] / (hiddenzone_dimensions['width'] * hiddenzone_dimensions['thickness'])
        # resistance = 0.5 * (1/ Kr) * hiddenzone_dimensions['length'] / (hiddenzone_dimensions['width'] * hiddenzone_dimensions['thickness'])

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

        :return: Extensibility z, y and x (MPa-1): {'z': float, 'y': float, 'x': float}
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

        :return: Extensibility z, y and x (MPa-1): {'z': float, 'y': float, 'x': float}
        :rtype: dict
        """
        #: TODO parameters into parameters.py
        phi = {}

        for phi_init_dimensions, phi_init_value in HiddenZone.PARAMETERS.phi_initial.items():
            if phi_init_dimensions == 'x': #: width
                # Function of Liu et al. (2007)
                if age <= HiddenZone.PARAMETERS.tend:
                    beta_function_norm = 2 / (1 + exp(3/1E06 * age))
                else:
                    beta_function_norm = 0
            if phi_init_dimensions == 'y': #: thickness
                # Function of Liu et al. (2007)
                if age <= HiddenZone.PARAMETERS.tend:
                    beta_function_norm = 2 / (1 + exp(3/1E06 * age))
                else:
                    beta_function_norm = 0
            if phi_init_dimensions == 'z': #: length
                # Function of Coussement et al. (2018)
                if age <= HiddenZone.PARAMETERS.tend:
                    beta_function_norm = (1 - (1 + (HiddenZone.PARAMETERS.tend - age) / (HiddenZone.PARAMETERS.tend - HiddenZone.PARAMETERS.tmax))
                                          * ((age - HiddenZone.PARAMETERS.tbase) / (HiddenZone.PARAMETERS.tend - HiddenZone.PARAMETERS.tbase)) **
                                          ((HiddenZone.PARAMETERS.tend - HiddenZone.PARAMETERS.tbase) / (HiddenZone.PARAMETERS.tend - HiddenZone.PARAMETERS.tmax)))
                else:
                    beta_function_norm = 0

            phi[phi_init_dimensions] = phi_init_value * beta_function_norm * delta_t

        return phi

    @staticmethod
    def calculate_extensibility_temperature(age, delta_teq, delta_t):

        """ Hidden zone extensibility in each dimension in relation to non-reversible dimensional changes.

        :param float age: hidden zone age (°Cd)
        :param float delta_teq: temperature-compensated time (s)
        :param float delta_t: time step of the simulation (s)

        :return: Extensibility z, y and x (MPa-1): {'z': float, 'y': float, 'x': float}
        :rtype: dict
        """
        #: TODO parameters into parameters.py

        phi = {}

        for phi_init_dimensions, phi_init_value in HiddenZone.PARAMETERS.phi_initial.items():
            if phi_init_dimensions == 'x': #: width
                # Function of Liu et al. (2007)
                beta_function_norm = 2 / (1 + exp(3 / 1E06 * age))
            if phi_init_dimensions == 'y': #: thickness
                # Function of Liu et al. (2007)
                beta_function_norm = 2 / (1 + exp(3 / 1E06 * age))
            if phi_init_dimensions == 'z': #: length
                if age <= HiddenZone.PARAMETERS.tend:
                    # Function of Coussement et al. (2018)
                    beta_function_norm = (1 - (1 + (HiddenZone.PARAMETERS.tend - age) / (HiddenZone.PARAMETERS.tend - HiddenZone.PARAMETERS.tmax))
                                          * ((age - HiddenZone.PARAMETERS.tbase) / (HiddenZone.PARAMETERS.tend - HiddenZone.PARAMETERS.tbase)) **
                                          ((HiddenZone.PARAMETERS.tend - HiddenZone.PARAMETERS.tbase) / (HiddenZone.PARAMETERS.tend - HiddenZone.PARAMETERS.tmax)))
                else:
                    beta_function_norm = 0

            phi[phi_init_dimensions] = phi_init_value * beta_function_norm * delta_t * (delta_teq / delta_t)

        return phi


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
    def calculate_delta_turgor_water_potential(phi, turgor_water_potential, volume, delta_water_content):
        """ Delta of turgor water potential of hidden zone.

        :param dict [str, float] phi: float phi: dict of cell wall extensibility (MPa). Keys = ['x', 'y', 'z]
        :param float turgor_water_potential: MPa
        :param float volume: m3
        :param float delta_water_content: delta water content integrated over delta t (g)

        :return: Delta of turgor water potential (MPa)
        :rtype: float
        """

        epsilon_x, epsilon_y, epsilon_z = HiddenZone.PARAMETERS.epsilon['x'], HiddenZone.PARAMETERS.epsilon['y'], HiddenZone.PARAMETERS.epsilon['z']
        elastic_component = (epsilon_x * epsilon_y * epsilon_z) / (epsilon_z * epsilon_x + epsilon_z * epsilon_y + epsilon_x * epsilon_y)   #: Elastic reversible growth (MPa)
        plastic_component = (phi['x'] + phi['y'] + phi['z'])   #: Plastic irreversible growth
        organ_volume = volume

        delta_turgor_water_potential = ((1 / (parameters.RHO_WATER * organ_volume * parameters.VSTORAGE)) * delta_water_content - plastic_component * (max(turgor_water_potential, HiddenZone.PARAMETERS.GAMMA) - HiddenZone.PARAMETERS.GAMMA)) * elastic_component  #: (MPa)

        return delta_turgor_water_potential

    @staticmethod
    def calculate_delta_organ_dimensions(delta_turgor_water_potential, turgor_water_potential, phi, organ_dimensions):
        """ Delta of organ dimensions according to turgor water potential, dimensions, extensibility and plasticity.
        Hidden zone geometry is supposed to be a rectangular prism.

        :param float delta_turgor_water_potential: delta of turgor water potential integrated over delta t (MPa)
        :param float turgor_water_potential: MPa
        :param dict phi: dict of cell wall extensibility (MPa). Keys = ['x', 'y', 'z]
        :param dict organ_dimensions: dict of organ dimensions at time t. Keys = [length', 'width', 'thickness'] (m)

        :return: Delta of organ specific-dimensions (m). Keys = ['leaf_L', 'width', 'thickness']
        :rtype: dict
        """
        delta_organ_dimensions = {}
        epsilon_dict = HiddenZone.PARAMETERS.epsilon
        mapping_dimensions = {'x': 'width', 'y': 'thickness', 'z': 'length'}

        for epsilon_dimension, epsilon_value in epsilon_dict.items():
            delta_organ_dimensions[mapping_dimensions[epsilon_dimension]] = ((1 / epsilon_value) * delta_turgor_water_potential +
                                                                             phi[epsilon_dimension] * (max(turgor_water_potential, HiddenZone.PARAMETERS.gamma) - HiddenZone.PARAMETERS.gamma)) *\
                                                                            organ_dimensions[mapping_dimensions[epsilon_dimension]]

        return delta_organ_dimensions

    @staticmethod
    def calculate_delta_organ_dimensions_plastic(turgor_water_potential, phi, organ_dimensions):
        """ Irreversible delta of organ dimensions according to turgor water potential, dimensions and plasticity.
        Hidden zone geometry is supposed to be a rectangular prism.

        :param float turgor_water_potential: MPa
        :param dict phi: dict of cell wall extensibility (MPa). Keys = ['x', 'y', 'z]
        :param dict organ_dimensions: dict of organ dimensions at time t. Keys = [length', 'width', 'thickness'] (m)

        :return: Delta of organ specific-dimensions (m). Keys = ['Leaf_L'', 'width', 'thickness']
        :rtype: dict
        """
        delta_organ_dimensions_plastic = {}
        epsilon_dict = HiddenZone.PARAMETERS.epsilon
        mapping_dimensions = {'x': 'width', 'y': 'thickness', 'z': 'length'}

        for epsilon_dimension, epsilon_value in epsilon_dict.items():
            delta_organ_dimensions_plastic[mapping_dimensions[epsilon_dimension]] = (phi[epsilon_dimension] * (max(turgor_water_potential, HiddenZone.PARAMETERS.GAMMA) - HiddenZone.PARAMETERS.GAMMA)) *\
                                                                            organ_dimensions[mapping_dimensions[epsilon_dimension]]

        return delta_organ_dimensions_plastic

    @staticmethod
    def calculate_delta_organ_dimensions_elastic(delta_turgor_water_potential, organ_dimensions):
        """ Reversible delta of organ dimensions according to turgor water potential, dimensions and extensibility.
        Hidden zone geometry is supposed to be a rectangular prism.

        :param float delta_turgor_water_potential: delta of turgor water potential integrated over delta t (MPa)
        :param dict organ_dimensions: dict of organ dimensions at time t. Keys = [length', 'width', 'thickness'] (m)

        :return: Delta of organ specific-dimensions (m). Keys = ['leaf_L'', 'width', 'thickness']
        :rtype: dict
        """
        delta_organ_dimensions_elastic = {}
        epsilon_dict = HiddenZone.PARAMETERS.epsilon
        mapping_dimensions = {'x': 'width', 'y': 'thickness', 'z': 'length'}

        for epsilon_dimension, epsilon_value in epsilon_dict.items():
            delta_organ_dimensions_elastic[mapping_dimensions[epsilon_dimension]] = ((1 / epsilon_value) * delta_turgor_water_potential) *\
                                                                            organ_dimensions[mapping_dimensions[epsilon_dimension]]

        return delta_organ_dimensions_elastic


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
        self.green_area = None  #: m2
        self.water_content = None   #: g H2O

    def calculate_aggregated_variables(self):
        self.Total_Transpiration_turgor = 0
        self.total_water_influx = 0
        self.green_area = 0
        self.water_content = 0
        for element in (self.exposed_element, self.enclosed_element):
            if element is not None:
                element.calculate_aggregated_variables()
                self.Total_Transpiration_turgor += element.Total_Transpiration_turgor
                self.total_water_influx += element.total_water_influx
                self.green_area += element.green_area
                self.water_content += element.water_content


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

    def __init__(self, label=None, is_growing = INIT_COMPARTMENTS.is_growing, temperature = INIT_COMPARTMENTS.temperature, age = INIT_COMPARTMENTS.age, green_area = INIT_COMPARTMENTS.green_area, mstruct = INIT_COMPARTMENTS.mstruct, Ts = INIT_COMPARTMENTS.Ts,
                 Tr = INIT_COMPARTMENTS.Tr, sucrose = INIT_COMPARTMENTS.sucrose, amino_acids = INIT_COMPARTMENTS.amino_acids, proteins = INIT_COMPARTMENTS.proteins, fructan = INIT_COMPARTMENTS.fructan,
                 osmotic_water_potential = INIT_COMPARTMENTS.osmotic_water_potential, total_water_potential = INIT_COMPARTMENTS.total_water_potential,
                 turgor_water_potential = INIT_COMPARTMENTS.turgor_water_potential, water_influx = INIT_COMPARTMENTS.water_influx, Wmax = INIT_COMPARTMENTS.Wmax,
                 length=INIT_COMPARTMENTS.length, thickness= INIT_COMPARTMENTS.thickness, width = INIT_COMPARTMENTS.width, water_content = INIT_COMPARTMENTS.water_content, cohorts=None, cohorts_replications=None):

        self.label = label                                      #: the label of the element
        if cohorts is None:  #: list of cohort values - Hack to treat tillering cases : TEMPORARY. Devrait être porté à l'échelle de la plante uniquement mais je ne vois pas comment faire mieux
            cohorts = []
        self.cohorts = cohorts  #: list of cohort values - Hack to treat tillering cases : TEMPORARY. Devrait être porté à l'échelle de la plante uniquement mais je ne vois pas comment faire mieux
        self.cohorts_replications = cohorts_replications  #: dictionary of number of replications per cohort rank

        # state parameters
        self.is_growing = is_growing                                        #: -
        self.age = age                                          #: °Cd
        self.Wmax = Wmax                                          #: m
        self.amino_acids = amino_acids                          #: :math:`\mu mol N
        self.green_area = green_area                            #: m2
        self.mstruct = mstruct                                  #: g
        self.proteins = proteins                                #: :math:`\mu mol N
        self.sucrose = sucrose                                  #: :math:`\mu mol C
        self.fructan = fructan                                  #: :math:`\mu mol C
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
        self.Total_Transpiration_turgor = self.calculate_Total_Transpiration(self.Tr, self.green_area, self.delta_t)
        self.total_water_influx = self.calculate_total_water_influx(self.water_influx)

    # VARIABLES
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

    @staticmethod
    def calculate_resistance(organ_dimensions):
        """
        Resistance of water flow between the lamina and xylem
        Relations were set proportional to the length and inversely proportional to the area of organ's cross section.
        :param dict lamina_dimensions: dict of lamina dimensions (m). Keys = ['length', 'width', 'thickness']
        :return: resistance (MPa s g-1)
        :rtype: float

        TODO : doc
        """
        #: Coussement et al. (2018)
        resistance = 0.5 * Xylem.PARAMETERS.R_xylem_organ * organ_dimensions['length'] / (organ_dimensions['width'] * organ_dimensions['thickness'])

        #: Martre et al. (2001)
        Kr = 7.2 * 18 * 1E-03
        # resistance = 0.5 * (1E-04 / Kr) * organ_dimensions['length'] / (organ_dimensions['width'] * organ_dimensions['thickness'])
        # resistance = 0.5 * (1/ Kr) * organ_dimensions['length'] / (organ_dimensions['width'] * organ_dimensions['thickness'])

        return resistance

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

        conversion_ratio = 18 * 1E-03 * delta_t # gH2O

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
    def calculate_osmotic_water_potential(sucrose, amino_acids, volume, temperature, fructan):
        """ Osmotic water potential of the hiddenzone calculated according to metabolites

        :param float sucrose: µmol C under the form of sucrose
        :param float amino_acids: µmol N under the form of amino acids
        :param float volume: (m3)
        :param float temperature: air temperature (°C)
        :param float fructan: µmol C under the form of fructan

        :return: Osmotic water potential (MPa)
        :rtype: float

        TODO : contribution dans parameters.py
        """
        temperature_K = temperature + parameters.CELSIUS_2_KELVIN

        sucrose = ((sucrose * 1E-6) / parameters.NB_C_SUCROSE)
        amino_acids = ((amino_acids * 1E-6) / parameters.AMINO_ACIDS_N_RATIO)
        fructan = (fructan * 1E-6) / (parameters.NB_C_SUCROSE)

        # Contribution des solutés organiques au potentiel osmotique (Baca Cabrera et al., 2020; Bajji et al., 2000; Thomas, 1991)
        contribution = 0.3

        osmotic_water_potential = - parameters.R * temperature_K * (fructan + sucrose + amino_acids) / (volume * parameters.RHO_WATER * parameters.VSTORAGE * contribution)

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
    def calculate_delta_water_content(water_influx, Total_Transpiration_turgor,):
        """ Delta of water flow for the lamina.

        :param float water_flux: Water influx from xylem integrated over delta_t (g)
        :param float Total_Transpiration_turgor: Element transpiration (g H2O)

        :return: Delta of water flow into the organ (g H2O)
        :rtype: float
        """
        return water_influx - Total_Transpiration_turgor


    @staticmethod
    def calculate_delta_turgor_water_potential(volume, delta_water_content):
        """ Delta of turgor water potential according to organ volume and elasticity.
        Extensibility (phi) is supposed to be 0 as this tissue is mature (growth completed).

        :param dict organ_dimensions: dict of organ dimensions at time t. Keys = ['length', 'thickness', 'width'] (m)
        :param float delta_water_content: delta water content integrated over delta t (g)

        :return: Delta of turgor water potential (MPa)
        :rtype: float
        """
        epsilon_z, epsilon_x, epsilon_y = PhotosyntheticOrganElement.PARAMETERS.epsilon['z'], PhotosyntheticOrganElement.PARAMETERS.epsilon['x'], PhotosyntheticOrganElement.PARAMETERS.epsilon['y']
        elastic_component = (epsilon_z * epsilon_x * epsilon_y) / (epsilon_z * epsilon_x + epsilon_z * epsilon_y + epsilon_x * epsilon_y)  #: Elastic reversible growth (MPa)
        plastic_component = 0   #: Plastic irreversible growth (MPa)
        organ_volume = volume
        delta_turgor_water_potential = ((1 / (parameters.RHO_WATER * organ_volume * parameters.VSTORAGE)) * delta_water_content - plastic_component) * elastic_component  #: (MPa)

        return delta_turgor_water_potential

    @staticmethod
    def calculate_delta_organ_dimensions(delta_turgor_water_potential, organ_dimensions):
        """Delta of lamina dimensions according to turgor water potential, dimensions, and elasticity.

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


class InternodeElement(PhotosyntheticOrganElement):
    """
    The class :class:`InternodeElement`.
    """

    PARAMETERS = parameters.INTERNODE_ELEMENT_PARAMETERS                           #: the internal parameters of the internode
    INIT_COMPARTMENTS = parameters.INTERNODE_ELEMENT_INIT_COMPARTMENTS             #: the initial values of compartments and state parameters


class SheathElement(PhotosyntheticOrganElement):
    """
    The class :class:`SheathElement`.
    """

    PARAMETERS = parameters.SHEATH_ELEMENT_PARAMETERS                   #: the internal parameters of the sheath
    INIT_COMPARTMENTS = parameters.SHEATH_ELEMENT_INIT_COMPARTMENTS     #: the initial values of compartments and state parameters
