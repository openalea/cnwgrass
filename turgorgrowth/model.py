# -*- coding: latin-1 -*-

from __future__ import division  # use "//" to do integer division
import math
import parameters


"""
    turgorgrowth.model
    ~~~~~~~~~~~~~

    The module :mod:`turgorgrowth.model` defines the equations of wtaer flow, turgor pressure and growth.

    :license: CeCILL-C, see LICENSE for details.

"""


class Population(object):
    """
    The class :class:`Population`.

    A :class:`population <Population>` must have at least one :class:`plant <Plant>`.
    """

    PARAMETERS = parameters.POPULATION_PARAMETERS  #: the internal parameters of the population

    def __init__(self, plants=None):
        if plants is None:
            plants = []
        self.plants = plants  #: the list of plants

    # def calculate_aggregated_variables(self):
    #     """Calculate the integrative variables of the population recursively.
    #     """
    #     for plant in self.plants:
    #         plant.calculate_aggregated_variables()


class Plant(object):
    """
    The class :class:`Plant` defines the CN exchanges at plant scale.

    A :class:`plant <Plant>` must have at least one :class:`axis <Axis>`.
    """

    PARAMETERS = parameters.PLANT_PARAMETERS  #: the internal parameters of the plants

    def __init__(self, index=None, axes=None):
        self.index = index  #: the index of the plant
        if axes is None:
            axes = []
        self.axes = axes  #: the list of axes

    def calculate_aggregated_variables(self):
        """Calculate the integrative variables of the plant recursively.
        """
        for axis in self.axes:
            axis.calculate_aggregated_variables()


class Axis(object):
    """
    The class :class:`Axis`.

    An :class:`axis <Axis>` must have:
        * one :class:`set of roots <Roots>`,
        * at least one :class:`phytomer<Phytomer>`.
    """

    PARAMETERS = parameters.AXIS_PARAMETERS  #: the internal parameters of the axes

    def __init__(self, label=None, roots=None, phytomers=None):
        self.label = label  #: the label of the axis
        self.roots = roots  #: the roots
        if phytomers is None:
            phytomers = []
        self.phytomers = phytomers  #: the list of phytomers

    # def calculate_aggregated_variables(self):
    #     """Calculate the integrative variables of the axis recursively.
    #     """
    #     if self.roots is not None:
    #         self.roots.calculate_aggregated_variables()
    #     for phytomer in self.phytomers:
    #         phytomer.calculate_aggregated_variables()


class Roots(object):
    """
    The class :class:`Roots`.
    """

    def __init__(self, label='roots', total_water_potential=-0.1):
        self.label = label  #: the index of the phytomer
        self.total_water_potential = total_water_potential  #: MPa


class Phytomer(object):
    """
    The class :class:`Phytomer`.

    A :class:`phytomer <Phytomer>` must have at least:
        * 1 photosynthetic organ: :class:`lamina <Lamina>`, :class:`internode <Internode>`,
                                  or :class:`sheath <Sheath>`.
        * or 1 :class:`hiddenzone <HiddenZone>`.
    """

    PARAMETERS = parameters.PHYTOMER_PARAMETERS  #: the internal parameters of the phytomers

    def __init__(self, index=None, lamina=None, internode=None, sheath=None, hiddenzone=None):
        self.index = index  #: the index of the phytomer
        self.lamina = lamina  #: the lamina
        self.internode = internode  #: the internode
        self.sheath = sheath  #: the sheath
        self.hiddenzone = hiddenzone  #: the hidden zone

    # def calculate_aggregated_variables(self):
    #     """Calculate the integrative variables of the phytomer recursively.
    #     """
    #     self.mstruct = 0
    #     for organ_ in (self.chaff, self.peduncle, self.lamina, self.internode, self.sheath, self.hiddenzone):
    #         if organ_ is not None:
    #             organ_.calculate_aggregated_variables()
    #             self.mstruct += organ_.mstruct


class Organ(object):
    """
    The class :class:`Organ`.

    :class:`Organ` is the base class of all organs. DO NOT INSTANTIATE IT.
    """
    PARAMETERS = parameters.ORGAN_PARAMETERS  #: the internal parameters of the organ

    def __init__(self, label):
        self.label = label  #: the label of the organ

    def initialize(self):
        """Initialize the derived attributes of the organ.
        """
        pass

    def calculate_aggregated_variables(self):
        """Calculate the integrative variables of the organ recursively.
        """
        pass


class HiddenZone(Organ):
    """
    The class :class:`HiddenZone`.
    """

    PARAMETERS = parameters.HIDDEN_ZONE_PARAMETERS                #: the internal parameters of the hidden zone
    INIT_COMPARTMENTS = parameters.HIDDEN_ZONE_INIT_COMPARTMENTS  #: the initial values of compartments and state parameters

    def __init__(self, label='hiddenzone', age=INIT_COMPARTMENTS.age, amino_acids=INIT_COMPARTMENTS.amino_acids, length=INIT_COMPARTMENTS.length, proteins=INIT_COMPARTMENTS.proteins,
                 radius=INIT_COMPARTMENTS.radius, sucrose=INIT_COMPARTMENTS.sucrose, temperature=INIT_COMPARTMENTS.temperature, turgor_water_potential=INIT_COMPARTMENTS.turgor_water_potential,
                 water_content=INIT_COMPARTMENTS.water_content):

        super(HiddenZone, self).__init__(label)
        self.label = label

        # state parameters
        self.amino_acids = amino_acids  #: µmol N
        self.proteins = proteins  #: µmol N
        self.sucrose = sucrose  #: :math:`\mu mol C
        self.temperature = temperature  #: °C

        # state variables
        self.age = age                                          #: °Cd
        self.length = length                                    #: m
        self.radius = radius                                    #: m
        self.turgor_water_potential = turgor_water_potential    #: MPa
        self.water_content = water_content                      #: g

        # fluxes
        self.water_influx = None  #: water influx into organ (g)

        # Intermediate variables
        self.osmotic_water_potential = None  #: MPa
        self.resistance = None               #: resistance of water flux between two organs (MPa s g-1)
        self.total_water_potential = None    #: MPa

    #:  Model equations for water flux
    def calculate_volume(self, length, radius):
        """ Hidden zone volume, assumed to be a cylinder.

        :Parameters:.
            - `length` (:class:`float`) - (m)
            - `radius` (:class:`float`) - (m)

        :Returns:
            volume (m3)
        :Returns Type:
            :class:`float`
        """
        return length * math.pi * radius ** 2

    def calculate_osmotic_water_potential(self, sucrose, amino_acids, proteins, volume, temperature):
        """ Osmotic water potential of the organ calculated according to metabolites

        :Parameters:
            - `sucrose` (:class:`float`) - µmol C under the form of sucrose
            - `amino_acids` (:class:`float`) - µmol N under the form of amino acids
            - `proteins` (:class:`float`) - µmol N under the form of proteins
            - `volume` (:class:`float`) - m3
            - `temperature` (:class:`float`) - degree Celsius

        :Returns:
            Osmotic water potential (MPa)
        :Returns Type:
            :class:`float`
        """
        # TODO: (i) add other metabolites, (ii) convert from µmol C to µmol metabolites, (iii) estimates of organ volume
        temperature_K = temperature + parameters.celsius_2_kelvin

        sucrose = ((sucrose * 1E-6) / parameters.NB_C_SUCROSE) * parameters.VANT_HOFF_SUCROSE
        amino_acids = ((amino_acids * 1E-6) / parameters.AMINO_ACIDS_N_RATIO) * parameters.VANT_HOFF_AMINO_ACIDS
        proteins = ((proteins * 1E-6) / parameters.AMINO_ACIDS_N_RATIO) * parameters.VANT_HOFF_AMINO_ACIDS

        osmotic_water_potential = - parameters.R * temperature_K * ((sucrose + amino_acids + proteins) / (volume * Organ.PARAMETERS.vstorage)) * 1E-6
        return osmotic_water_potential

    def calculate_water_potential(self, turgor_water_potential, osmotic_water_potential):
        """ Total water potential of the organ

        :Parameters:
            - `turgor_water_potential` (:class:`float`) - MPa
            - `osmotic_water_potential` (:class:`float`) - MPa

        :Returns:
            Total water potential (MPa)
        :Returns Type:
            :class:`float`
        """
        return turgor_water_potential + osmotic_water_potential

    def calculate_resistance(self, hiddenzone_dimensions, previous_organ):
        """ Resistance of water flow between the hiddenzone and its predecessor (can be root or the last elongated internode if any).
            Relations were set proportional to the length and inversely proportional to the area of organ's cross section.

        :Parameters:.
            - `hiddenzone_dimensions` (:class:`dict`) - dict of hidden zone dimensions at time t. Keys = ['length', 'radius'] (m)
            - `previous organ` (:class:`object`) - predecessor of the hiddenzone, either root or internode element with length and radius attributes

        :Returns:
            resistance (MPa s g-1)
        :Returns Type:
            :class:`float`
        """
        if previous_organ.label == 'roots':
            resistance = 0.5 * HiddenZone.PARAMETERS.R_xylem * (hiddenzone_dimensions['length'] / hiddenzone_dimensions['radius'] ** 2)  # TODO: check this with Tom
        else:
            resistance = 0.5 * HiddenZone.PARAMETERS.R_xylem * (hiddenzone_dimensions['length'] / hiddenzone_dimensions['radius'] ** 2 +
                                                                previous_organ.length / previous_organ.radius ** 2)
        return resistance

    def calculate_water_influx(self, organ_water_potential, previous_organ_water_potential, resistance, delta_t):
        """ Sap flow into the organ according to water potential gradient with the previous organ
        For lamina, the previous organ is the hiddenzone if leaf not mature, else it is the sheath
        For hiddenzones and sheaths, the previous organ is the root compartment if located above short internodes, else it is the previous internode
        For internodes the previous organ is the previous one else the organ is neglected.

        :Parameters:
            - `organ_water_potential` (:class:`float`) - water potential of the current organ (MPa)
            - `previous_organ_water_potential` (:class:`float`) - water potential of the previous organ (MPa)
            - `resistance` (:class:`float`) - transport resistance between organs (MPa s g-1)
            - `delta_t` (:class:`float`) - Time step of the simulation (s)

        :Returns:
            Water influx into the current organ integrated over delta_t (g)
        :Returns Type:
            :class:`float`
        """
        return ((previous_organ_water_potential - organ_water_potential) / resistance) * delta_t

    def calculate_delta_water_content(self, water_influx, water_outflow=0):
        """ delta of sap flow for the hidden zone.

        :Parameters:
            - `water_influx` (:class:`float`) - Water influx integrated over delta_t (g)
            - `water_outflow` (:class:`float`) - Water loss through the emerged lamina or sheath if any (g).

        :Returns:
            Delta of sap flow into the organ (g)
        :Returns Type:
            :class:`float`
        """
        return water_influx - water_outflow

    def calculate_extensibility(self, age, delta_t):

        """ Hidden zone extensibility in each dimension in relation to non-reversible dimensional changes.

        :Parameters:
            - `age` (:class:`float`) - hidden zone age (°Cd)
            - `delta_t` (:class:`float`) - Time step of the simulation (s)

        :Returns:
            Dict of extensibility (MPa-1)). Keys = ['z', 'x'].
        :Returns Type:
            :class:`dict`
        """
        if age <= HiddenZone.PARAMETERS.tend:
            beta_function_norm = (1 - (1 + (HiddenZone.PARAMETERS.tend - age) / (HiddenZone.PARAMETERS.tend - HiddenZone.PARAMETERS.tmax)) * (age / HiddenZone.PARAMETERS.tend) **
                                  (HiddenZone.PARAMETERS.tend / (HiddenZone.PARAMETERS.tend - HiddenZone.PARAMETERS.tmax)))
        else:
            beta_function_norm = 0

        phi = {}
        for phi_init_dimensions, phi_init_value in HiddenZone.PARAMETERS.phi_initial.items():
            phi[phi_init_dimensions] = phi_init_value * beta_function_norm * delta_t
        return phi

    def calculate_delta_turgor_water_potential(self, phi, turgor_water_potential, organ_dimensions, delta_water_content):
        """ Delta of turgor water potential according to organ water content, turgor water potential, dimensions, extensibility and plasticity.

        :Parameters:
            - `phi` (:class:`dict`) - dict of cell wall extensibility (MPa). Keys = ['z', 'x']
            - `turgor_water_potential` (:class:`float`) - MPa
            - `organ_dimensions` (:class:`dict`) - dict of organ dimensions at time t. Keys = ['length', 'radius'] (m)
            - `delta_water_content` (:class:`float`) - delta water content integrated over delta t (g)

        :Returns:
            Delta of turgor water potential (MPa)
        :Returns Type:
            :class:`float`
        """
        epsilon_z, epsilon_x = HiddenZone.PARAMETERS.epsilon['z'], HiddenZone.PARAMETERS.epsilon['x']
        elastic_component = (epsilon_z * epsilon_x) / (2 * epsilon_z + epsilon_x)   #: Elastic reversible growth (MPa)
        plastic_component = (phi['z'] + 2 * phi['x']) * (max(turgor_water_potential, HiddenZone.PARAMETERS.gamma) - HiddenZone.PARAMETERS.gamma)  #: Plastic irreversible growth
        organ_volume = math.pi * organ_dimensions['radius'] ** 2 * organ_dimensions['length']  #: (m3)
        delta_turgor_water_potential = ((1 / (parameters.rho_water * organ_volume)) * delta_water_content * 1E-3 - plastic_component) * elastic_component  #: (MPa)

        return delta_turgor_water_potential

    def calculate_delta_organ_dimensions(self, delta_turgor_water_potential, turgor_water_potential, phi, organ_dimensions):
        """ Delta of organ dimensions according to turgor water potential, dimensions, extensibility and plasticity.

        :Parameters:
            - `delta_turgor_water_potential` (:class:`float`) - delta of turgor water potential integrated over delta t (MPa)
            - `turgor_water_potential` (:class:`float`) - MPa
            - `phi` (:class:`dict`) - dict of cell wall extensibility (MPa). Keys = ['z', 'x']
            - `organ_dimensions` (:class:`dict`) - dict of organ dimensions at time t. Keys = ['length', 'radius'] (m)

        :Returns:
            Delta of organ specific-dimensions (m). Keys = ['length', 'radius']
        :Returns Type:
            :class:`dict`
        """
        delta_organ_dimensions = {}
        epsilon_dict = HiddenZone.PARAMETERS.epsilon.items()
        mapping_dimensions = {'x': 'radius', 'z': 'length'}

        for epsilon_dimension, epsilon_value in epsilon_dict:
            delta_organ_dimensions[mapping_dimensions[epsilon_dimension]] = ((1 / epsilon_value) * delta_turgor_water_potential +
                                                                             phi[epsilon_dimension] * (max(turgor_water_potential, LaminaElement.PARAMETERS.gamma) - LaminaElement.PARAMETERS.gamma)) *\
                                                                            organ_dimensions[mapping_dimensions[epsilon_dimension]]
        return delta_organ_dimensions

    def calculate_growing_degree_days(self, delta_t, air_temperature=20):
        """ Organ thermal age. Tbase = 0°C.

        :Parameters:
            - `air_temperature` (:class:`float`) - Mean air temperature during delta t (°C)
            - `delta_t` (:class:`float`) - Time step of the simulation (s)

        :Returns:
            Growing_degree_days (°C).
        :Returns Type:
            :class:`float`
        """
        return max(0., (air_temperature / (parameters.HOUR_TO_SECOND_CONVERSION_FACTOR * 24)) * delta_t)


class PhotosyntheticOrgan(Organ):
    """
    The class :class:`PhotosyntheticOrgan` defines the water flow in a photosynthetic organ.

    A :class:`photosynthetic organ <PhotosyntheticOrgan>` must have at least 1
    :class:`lamina <Lamina>`,
    :class:`internode  <Internodet>`,
    or :class:`sheath  <Sheath>`.

    :class:`PhotosyntheticOrgan` is the base class of all photosynthetic organs. DO NOT INSTANTIATE IT.
    """

    def __init__(self, label, exposed_element, enclosed_element):

        super(PhotosyntheticOrgan, self).__init__(label)

        self.exposed_element = exposed_element    #: the exposed element
        self.enclosed_element = enclosed_element  #: the enclosed element


class Lamina(PhotosyntheticOrgan):
    """
    The class :class:`Lamina`.
    """

    def __init__(self, label=None, exposed_element=None, enclosed_element=None):
        super(Lamina, self).__init__(label, exposed_element, enclosed_element)


class Internode(PhotosyntheticOrgan):
    """
    The class :class:`Internode`.
    """

    def __init__(self, label=None, exposed_element=None, enclosed_element=None):
        super(Internode, self).__init__(label, exposed_element, enclosed_element)


class Sheath(PhotosyntheticOrgan):
    """
    The class :class:`Sheath`.
    """

    def __init__(self, label=None, exposed_element=None, enclosed_element=None):
        super(Sheath, self).__init__(label, exposed_element, enclosed_element)


class PhotosyntheticOrganElement(object):
    """
    The class :class:`PhotosyntheticOrganElement` defines the water flow in a photosynthetic organ element.

    A :class:`photosynthetic organ element <PhotosyntheticOrganElement>` must have at least 1
    :class:`lamina element<LaminaElement>`,
    :class:`internode element <InternodeElement>`,
    or :class:`sheath element <SheathElement>`.

    :class:`PhotosyntheticOrganElement` is the base class of all photosynthetic organ elements. DO NOT INSTANTIATE IT.
    """

    INIT_COMPARTMENTS = parameters.PHOTOSYNTHETIC_ORGAN_ELEMENT_INIT_COMPARTMENTS  #: the initial values of compartments and state parameters

    def __init__(self, age=INIT_COMPARTMENTS.age, amino_acids=INIT_COMPARTMENTS.amino_acids, label=None, length=INIT_COMPARTMENTS.length,
                 osmotic_water_potential=INIT_COMPARTMENTS.osmotic_water_potential, proteins=INIT_COMPARTMENTS.proteins, sucrose=INIT_COMPARTMENTS.sucrose, thickness=INIT_COMPARTMENTS.thickness,
                 temperature=INIT_COMPARTMENTS.temperature, total_water_potential=INIT_COMPARTMENTS.total_water_potential, transpiration=INIT_COMPARTMENTS.transpiration,
                 turgor_water_potential=INIT_COMPARTMENTS.turgor_water_potential, water_content=INIT_COMPARTMENTS.water_content):

        self.label = label                                      #: the label of the element

        # state parameters
        self.amino_acids = amino_acids  #: µmol N
        self.proteins = proteins  #: µmol N
        self.sucrose = sucrose                                  #: :math:`\mu mol C
        self.temperature = temperature                          #: °C
        self.transpiration = transpiration                      #: mmol H20 s-1

        # state variables
        self.age = age                                          #: °Cd
        self.length = length                                    #: m
        self.osmotic_water_potential = osmotic_water_potential  #: MPa
        self.thickness = thickness                              #: m
        self.total_water_potential = total_water_potential      #: MPa
        self.turgor_water_potential = turgor_water_potential    #: MPa
        self.water_content = water_content                      #: g

        # fluxes
        self.water_influx = 0                                   #: water influx into organ (g)
        self.resistance = None                                  #: resistance of water flux between two organs (MPa s g-1)

    #: Water flow equations common to all photosynthetic organ elements

    def calculate_water_potential(self, turgor_water_potential, osmotic_water_potential):
        """ Total water potential of the organ

        :Parameters:
            - `turgor_water_potential` (:class:`float`) - MPa
            - `osmotic_water_potential` (:class:`float`) - MPa

        :Returns:
            Total water potential (MPa)
        :Returns Type:
            :class:`float`
        """
        return turgor_water_potential + osmotic_water_potential

    def calculate_water_influx(self, organ_water_potential, previous_organ_water_potential, resistance, delta_t):
        """ Sap flow into the organ according to water potential gradient with the previous organ
        For lamina, the previous organ is the hiddenzone if leaf not mature, else it is the sheath
        For hiddenzones and sheaths, the previous organ is the root compartment if located above short internodes, else it is the previous internode
        For internodes the previous organ is the previous one else the organ is neglected.

        :Parameters:
            - `organ_water_potential` (:class:`float`) - water potential of the current organ (MPa)
            - `previous_organ_water_potential` (:class:`float`) - water potential of the previous organ (MPa)
            - `resistance` (:class:`float`) - transport resistance between organs (MPa s g-1)
            - `delta_t` (:class:`float`) - Time step of the simulation (s)

        :Returns:
            Water influx into the current organ integrated over delta_t (g)
        :Returns Type:
            :class:`float`
        """
        return ((previous_organ_water_potential - organ_water_potential) / resistance) * delta_t

    def calculate_growing_degree_days(self, delta_t, air_temperature=20):
        """ Organ thermal age. Tbase = 0°C.

        :Parameters:
            - `air_temperature` (:class:`float`) - Mean air temperature during delta t (°C)
            - `delta_t` (:class:`float`) - Time step of the simulation (s)

        :Returns:
            Growing_degree_days (°C).
        :Returns Type:
            :class:`float`
        """
        return max(0., (air_temperature / (parameters.HOUR_TO_SECOND_CONVERSION_FACTOR * 24)) * delta_t)


class LaminaElement(PhotosyntheticOrganElement):
    """
    The class :class:`LaminaElement`.
    """

    PARAMETERS = parameters.LAMINA_PARAMETERS                   #: the internal parameters of the lamina

    def __init__(self, thickness, width, **kwargs):
        super(LaminaElement, self).__init__(**kwargs)

        self.thickness = thickness  #: m
        self.width = width          #: m

    # Model equations for water flux
    def calculate_volume(self, length, width, thickness):
        """ Lamina volume, assumed to be a parallelepiped.

        :Parameters:.
            - `length` (:class:`float`) - (m)
            - `width` (:class:`float`) - (m)
            - `thickness` (:class:`float`) - (m)

        :Returns:
            volume (m3)
        :Returns Type:
            :class:`float`
        """
        return length * width * thickness

    def calculate_osmotic_water_potential(self, sucrose, amino_acids, proteins, organ_volume, organ_temperature):
        """ Osmotic water potential of the organ calculated according to metabolites

        :Parameters:
            - `sucrose` (:class:`float`) - µmol C under the form of sucrose
            - `amino_acids` (:class:`float`) - µmol N under the form of amino acids
            - `proteins` (:class:`float`) - µmol N under the form of proteins
            - `organ_volume` (:class:`float`) - m3
            - `organ_temperature` (:class:`float`) - degree Celsius

        :Returns:
            Osmotic water potential (MPa)
        :Returns Type:
            :class:`float`
        """
        # TODO: (i) add other metabolites, (ii) convert from µmol C to µmol metabolites, (iii) estimates of organ volume
        temperature_K = organ_temperature + parameters.celsius_2_kelvin

        sucrose = ((sucrose * 1E-6) / parameters.NB_C_SUCROSE) * parameters.VANT_HOFF_SUCROSE
        amino_acids = ((amino_acids * 1E-6) / parameters.AMINO_ACIDS_N_RATIO) * parameters.VANT_HOFF_AMINO_ACIDS
        proteins = ((proteins * 1E-6) / parameters.AMINO_ACIDS_N_RATIO) * parameters.VANT_HOFF_AMINO_ACIDS

        osmotic_water_potential = - parameters.R * temperature_K * ((sucrose + amino_acids + proteins) / (organ_volume * LaminaElement.PARAMETERS.vstorage)) * 1E-6
        return osmotic_water_potential

    def calculate_resistance(self, lamina_dimensions, previous_organ):
        """ Resistance of water flow between the lamina and hidden zone (if leaf is growing) or the sheath (if lamina is mature).
            Relations were set proportional to the length and inversely proportional to the area of organ's cross section.

        :Parameters:.
            - `lamina_dimensions` (:class:`dict`) - dict of lamina dimensions (m). Keys = ['length', 'width', 'thickness']
            - `previous_organ` (:class:`object`) - predecessor of the lamina, either sheath or hiddenzone with length and radius attributes

        :Returns:
            resistance (MPa s g-1)
        :Returns Type:
            :class:`float`
        """
        resistance = 0.5 * LaminaElement.PARAMETERS.R_xylem * (previous_organ.length / previous_organ.radius ** 2 +
                                                               lamina_dimensions['length'] / (lamina_dimensions['width'] * lamina_dimensions['thickness']))
        return resistance

    def calculate_delta_water_content(self, water_influx, transpiration, delta_t):
        """ delta of sap flow for the lamina.

        :Parameters:
            - `water_influx` (:class:`float`) - Water influx integrated over delta_t (g)
            - `transpiration` (:class:`float`) - Lamina transpiration rate (mmol H20 s-1)
            - `delta_t` (:class:`float`) - Time step of the simulation (s).

        :Returns:
            Delta of sap flow into the organ (g)
        :Returns Type:
            :class:`float`
        """
        return water_influx - (transpiration * parameters.WATER_MOLAR_MASS * 1E-3 * delta_t)

    def calculate_extensibility(self, age, delta_t):

        """ Lamina extensibility in each dimension in relation to non-reversible dimensional changes.

        :Parameters:
            - `age` (:class:`float`) - Lamina age (°Cd)
            - `delta_t` (:class:`float`) - Time step of the simulation (s).

        :Returns:
            Dict of extensibility (MPa-1)). Keys = ['z', 'x', 'y'].
        :Returns Type:
            :class:`dict`
        """
        if age <= LaminaElement.PARAMETERS.tend:
            beta_function_norm = (1 - (1 + (LaminaElement.PARAMETERS.tend - age) / (LaminaElement.PARAMETERS.tend - LaminaElement.PARAMETERS.tmax)) * (age / LaminaElement.PARAMETERS.tend) **
                                  (LaminaElement.PARAMETERS.tend / (LaminaElement.PARAMETERS.tend - LaminaElement.PARAMETERS.tmax)))
        else:
            beta_function_norm = 0

        phi = {}
        for phi_init_dimensions, phi_init_value in LaminaElement.PARAMETERS.phi_initial.items():
            phi[phi_init_dimensions] = phi_init_value * beta_function_norm * delta_t
        return phi

    def calculate_delta_turgor_water_potential(self, phi, turgor_water_potential, organ_dimensions, delta_water_content):
        """ Delta of turgor water potential according to organ water content, turgor water potential, dimensions, extensibility and plasticity.

        :Parameters:
            - `phi` (:class:`dict`) - dict of cell wall extensibility (MPa). Keys = ['z', 'x', 'y']
            - `turgor_water_potential` (:class:`float`) - MPa
            - `organ_dimensions` (:class:`dict`) - dict of organ dimensions at time t. Keys = ['length', 'thickness', 'width'] (m)
            - `delta_water_content` (:class:`float`) - delta water content integrated over delta t (g)

        :Returns:
            Delta of turgor water potential (MPa)
        :Returns Type:
            :class:`float`
        """
        epsilon_z, epsilon_x, epsilon_y = LaminaElement.PARAMETERS.epsilon['z'], LaminaElement.PARAMETERS.epsilon['x'], LaminaElement.PARAMETERS.epsilon['y']
        elastic_component = (epsilon_z * epsilon_x * epsilon_y) / (epsilon_z * epsilon_x + epsilon_z * epsilon_y + epsilon_x * epsilon_y)  #: Elastic reversible growth (MPa)
        plastic_component = (phi['z'] + phi['x'] + phi['y']) * (max(turgor_water_potential, LaminaElement.PARAMETERS.gamma) - LaminaElement.PARAMETERS.gamma)  #: Plastic irreversible growth
        organ_volume = organ_dimensions['length'] * organ_dimensions['width'] * organ_dimensions['thickness']  #: (m3)
        delta_turgor_water_potential = ((1 / (parameters.rho_water * organ_volume)) * delta_water_content - plastic_component) * elastic_component  #: (MPa)

        return delta_turgor_water_potential

    def calculate_delta_organ_dimensions(self, delta_turgor_water_potential, turgor_water_potential, phi, organ_dimensions):
        """ Delta of organ dimensions according to turgor water potential, dimensions, extensibility and plasticity.

        :Parameters:
            - `delta_turgor_water_potential` (:class:`float`) - delta of turgor water potential integrated over delta t (MPa)
            - `turgor_water_potential` (:class:`float`) - MPa
            - `phi` (:class:`dict`) - dict of cell wall extensibility (MPa). Keys = ['z', 'x', 'y]
            - `organ_dimensions` (:class:`dict`) - dict of organ dimensions at time t. Keys = ['length', 'thickness', 'width'] (m)

        :Returns:
            Delta of organ specific-dimensions (m). Keys = ['length', 'width', 'thickness']
        :Returns Type:
            :class:`dict`
        """
        delta_organ_dimensions = {}
        epsilon_dict = LaminaElement.PARAMETERS.epsilon.items()
        mapping_dimensions = {'x': 'width', 'y': 'thickness', 'z': 'length'}

        for epsilon_dimension, epsilon_value in epsilon_dict:
            delta_organ_dimensions[mapping_dimensions[epsilon_dimension]] = ((1 / epsilon_value) * delta_turgor_water_potential +
                                                                             phi[epsilon_dimension] * (max(turgor_water_potential, LaminaElement.PARAMETERS.gamma) - LaminaElement.PARAMETERS.gamma)) *\
                                                                             organ_dimensions[mapping_dimensions[epsilon_dimension]]
        return delta_organ_dimensions


class InternodeElement(Internode):
    """
    The class :class:`InternodeElement`.
    """

    PARAMETERS = parameters.INTERNODE_PARAMETERS                #: the internal parameters of the lamina

    def __init__(self, radius, **kwargs):
        super(InternodeElement, self).__init__(**kwargs)

        self.radius = radius  #: m

    # Model equations for water flux
    def calculate_volume(self, length, radius):
        """ Internode volume, assumed to be a cylinder.

        :Parameters:.
            - `length` (:class:`float`) - (m)
            - `radius` (:class:`float`) - (m)

        :Returns:
            volume (m3)
        :Returns Type:
            :class:`float`
        """
        return length * math.pi * radius ** 2

    def calculate_osmotic_water_potential(self, sucrose, amino_acids, proteins, organ_volume, organ_temperature):
        """ Osmotic water potential of the organ calculated according to metabolites

        :Parameters:
            - `sucrose` (:class:`float`) - µmol C under the form of sucrose
            - `amino_acids` (:class:`float`) - µmol N under the form of amino acids
            - `proteins` (:class:`float`) - µmol N under the form of proteins
            - `organ_volume` (:class:`float`) - m3
            - `organ_temperature` (:class:`float`) - degree Celsius

        :Returns:
            Osmotic water potential (MPa)
        :Returns Type:
            :class:`float`
        """
        # TODO: (i) add other metabolites, (ii) convert from µmol C to µmol metabolites, (iii) estimates of organ volume
        temperature_K = organ_temperature + parameters.celsius_2_kelvin

        sucrose = ((sucrose * 1E-6) / parameters.NB_C_SUCROSE) * parameters.VANT_HOFF_SUCROSE
        amino_acids = ((amino_acids * 1E-6) / parameters.AMINO_ACIDS_N_RATIO) * parameters.VANT_HOFF_AMINO_ACIDS
        proteins = ((proteins * 1E-6) / parameters.AMINO_ACIDS_N_RATIO) * parameters.VANT_HOFF_AMINO_ACIDS

        osmotic_water_potential = - parameters.R * temperature_K * ((sucrose + amino_acids + proteins) / (organ_volume * InternodeElement.PARAMETERS.vstorage)) * 1E-6
        return osmotic_water_potential

    def calculate_resistance(self, internode_dimensions, previous_organ):
        """ Resistance of water flow between the internode and the root compartment or the previous internode if elongated.
            Relations were set proportional to the length and inversely proportional to the area of organ's cross section.

        :Parameters:.
            - `internode_dimensions` (:class:`dict`) - dict of current internode dimensions. Keys = ['length', 'radius'] (m)
            - `previous_organ` (:class:`object`) - predecessor of the internode, either root or internode element with length and radius attributes

        :Returns:
            resistance (MPa s g-1)
        :Returns Type:
            :class:`float`
        """
        if previous_organ.label == 'roots':
            resistance = 0.5 * HiddenZone.PARAMETERS.R_xylem * (internode_dimensions['length'] / internode_dimensions['radius'] ** 2)
        else:
            resistance = 0.5 * InternodeElement.PARAMETERS.R_xylem * (internode_dimensions['length'] / internode_dimensions['radius'] ** 2 +
                                                                      previous_organ.length / previous_organ.radius ** 2)
        return resistance

    def calculate_delta_water_content(self, water_influx, water_influx_next_internode, transpiration, delta_t):
        """ delta of sap flow for the internode.

        :Parameters:
            - `water_influx` (:class:`float`) - Water influx integrated of the current internode over delta_t (g)
            - `water_influx_next_internode` (:class:`float`) - Water influx of the next organ (if any) integrated over delta_t (g)
            - `transpiration` (:class:`float`) - Internode transpiration rate (if any, mmol H20 s-1)
            - `delta_t` (:class:`float`) - Time step of the simulation (s)

        :Returns:
            Delta of sap flow into the organ (g)
        :Returns Type:
            :class:`float`
        """
        return water_influx - water_influx_next_internode - (transpiration * parameters.WATER_MOLAR_MASS * 1E-3 * delta_t)

    def calculate_extensibility(self, age, delta_t):

        """ Internode extensibility in each dimension in relation to non-reversible dimensional changes.

        :Parameters:
            - `age` (:class:`float`) - internode age (°Cd)
            - `delta_t` (:class:`float`) - Time step of the simulation (s).

        :Returns:
            Dict of extensibility (MPa-1). Keys = ['z', 'x'].
        :Returns Type:
            :class:`dict`
        """
        if age <= InternodeElement.PARAMETERS.tend:
            beta_function_norm = (1 - (1 + (InternodeElement.PARAMETERS.tend - age) / (InternodeElement.PARAMETERS.tend - InternodeElement.PARAMETERS.tmax)) * (age / InternodeElement.PARAMETERS.tend)
                                  ** (InternodeElement.PARAMETERS.tend / (InternodeElement.PARAMETERS.tend - InternodeElement.PARAMETERS.tmax)))
        else:
            beta_function_norm = 0

        phi = {}
        for phi_init_dimensions, phi_init_value in InternodeElement.PARAMETERS.phi_initial.items():
            phi[phi_init_dimensions] = phi_init_value * beta_function_norm * delta_t
        return phi

    def calculate_delta_turgor_water_potential(self, phi, turgor_water_potential, organ_dimensions, delta_water_content):
        """ Delta of turgor water potential according to organ water content, turgor water potential, dimensions, extensibility and plasticity.

        :Parameters:
            - `phi` (:class:`dict`) - dict of cell wall extensibility (MPa). Keys = ['z', 'x']
            - `turgor_water_potential` (:class:`float`) - MPa
            - `organ_dimensions` (:class:`dict`) - dict of organ dimensions at time t. Keys = ['length', 'radius'] (m)
            - `delta_water_content` (:class:`float`) - delta water content integrated over delta t (g)

        :Returns:
            Delta of turgor water potential (MPa)
        :Returns Type:
            :class:`float`
        """
        epsilon_z, epsilon_x = InternodeElement.PARAMETERS.epsilon['z'], InternodeElement.PARAMETERS.epsilon['x']
        elastic_component = (epsilon_z * epsilon_x) / (2 * epsilon_z + epsilon_x)   #: Elastic reversible growth (MPa)
        plastic_component = (phi['z'] + 2 * phi['x']) * (max(turgor_water_potential, InternodeElement.PARAMETERS.gamma) - InternodeElement.PARAMETERS.gamma)  #: Plastic irreversible growth
        organ_volume = math.pi * organ_dimensions['radius'] ** 2 * organ_dimensions['length']  #: (m3)
        delta_turgor_water_potential = ((1 / (parameters.rho_water * organ_volume)) * delta_water_content * 1E-3 - plastic_component) * elastic_component  #: (MPa)

        return delta_turgor_water_potential

    def calculate_delta_organ_dimensions(self, delta_turgor_water_potential, turgor_water_potential, phi, organ_dimensions):
        """ Delta of organ dimensions according to turgor water potential, dimensions, extensibility and plasticity.

        :Parameters:
            - `delta_turgor_water_potential` (:class:`float`) - delta of turgor water potential integrated over delta t (MPa)
            - `turgor_water_potential` (:class:`float`) - MPa
            - `phi` (:class:`dict`) - dict of cell wall extensibility at time t  (MPa). Keys = ['z', 'x']
            - `organ_dimensions` (:class:`dict`) - dict of organ dimensions at time t. Keys = ['length', 'radius'] (m)

        :Returns:
            Delta of organ specific-dimensions (m). Keys = ['length', 'radius']
        :Returns Type:
            :class:`dict`
        """
        delta_organ_dimensions = {}
        epsilon_dict = InternodeElement.PARAMETERS.epsilon.items()
        mapping_dimensions = {'x': 'radius', 'z': 'length'}

        for epsilon_dimension, epsilon_value in epsilon_dict:
            delta_organ_dimensions[mapping_dimensions[epsilon_dimension]] = ((1 / epsilon_value) * delta_turgor_water_potential +
                                                                             phi[epsilon_dimension] * (max(turgor_water_potential, LaminaElement.PARAMETERS.gamma) - LaminaElement.PARAMETERS.gamma)) *\
                                                                            organ_dimensions[mapping_dimensions[epsilon_dimension]]
        return delta_organ_dimensions


class SheathElement(PhotosyntheticOrganElement):
    """
    The class :class:`SheathElement`.
    """

    PARAMETERS = parameters.SHEATH_PARAMETERS                #: the internal parameters of the sheath

    def __init__(self, radius, **kwargs):
        super(SheathElement, self).__init__(**kwargs)

        self.radius = radius  #: m

    # Model equations for water flux
    def calculate_volume(self, length, radius):
        """ Sheath volume, assumed to be a cylinder.

        :Parameters:.
            - `length` (:class:`float`) - (m)
            - `radius` (:class:`float`) - (m)

        :Returns:
            volume (m3)
        :Returns Type:
            :class:`float`
        """
        return length * math.pi * radius ** 2

    def calculate_osmotic_water_potential(self, sucrose, amino_acids, proteins, organ_volume, organ_temperature):
        """ Osmotic water potential of the organ calculated according to metabolites

        :Parameters:
            - `sucrose` (:class:`float`) - µmol C under the form of sucrose
            - `amino_acids` (:class:`float`) - µmol N under the form of amino acids
            - `proteins` (:class:`float`) - µmol N under the form of proteins
            - `organ_volume` (:class:`float`) - m3
            - `organ_temperature` (:class:`float`) - degree Celsius

        :Returns:
            Osmotic water potential (MPa)
        :Returns Type:
            :class:`float`
        """
        # TODO: (i) add other metabolites, (ii) convert from µmol C to µmol metabolites, (iii) estimates of organ volume
        temperature_K = organ_temperature + parameters.celsius_2_kelvin

        sucrose = ((sucrose * 1E-6) / parameters.NB_C_SUCROSE) * parameters.VANT_HOFF_SUCROSE
        amino_acids = ((amino_acids * 1E-6) / parameters.AMINO_ACIDS_N_RATIO) * parameters.VANT_HOFF_AMINO_ACIDS
        proteins = ((proteins * 1E-6) / parameters.AMINO_ACIDS_N_RATIO) * parameters.VANT_HOFF_AMINO_ACIDS

        osmotic_water_potential = - parameters.R * temperature_K * ((sucrose + amino_acids + proteins) / (organ_volume * SheathElement.PARAMETERS.vstorage)) * 1E-6
        return osmotic_water_potential

    def calculate_resistance(self, sheath_dimensions, previous_organ):
        """ Resistance of water flow between the sheath and the previous sheath element if any, the hiddenzone if any , the previous elongated internode if any or the root compartment.
            Relations were set proportional to the length and inversely proportional to the area of organ's cross section.

        :Parameters:.
            - `sheath_dimensions` (:class:`dict`) - dict of current hidden zone dimensions. Keys = ['length', 'radius'] (m)
            - `previous_organ` (:class:`object`) - predecessor of the sheath, either sheath, internode, hiddenzone or root with length and radius attributes

        :Returns:
            resistance (MPa s g-1)
        :Returns Type:
            :class:`float`
        """
        if previous_organ.label == 'roots':
            resistance = 0.5 * HiddenZone.PARAMETERS.R_xylem * (sheath_dimensions['length'] / sheath_dimensions['radius'] ** 2)
        else:
            resistance = 0.5 * SheathElement.PARAMETERS.R_xylem * (sheath_dimensions['length'] / sheath_dimensions['radius'] ** 2 +
                                                                   previous_organ.length / previous_organ.radius ** 2)
        return resistance

    def calculate_delta_water_content(self, water_influx, water_influx_next_internode, transpiration, delta_t):
        """ delta of sap flow for the sheath.

        :Parameters:
            - `water_influx` (:class:`float`) - Water influx integrated of the current sheath over delta_t (g)
            - `water_influx_next_internode` (:class:`float`) - Water influx of the next organ (if any) integrated over delta_t (g)
            - `transpiration` (:class:`float`) - Internode transpiration rate (if any, mmol H20 s-1)
            - `delta_t` (:class:`float`) - Time step of the simulation (s)

        :Returns:
            Delta of sap flow into the organ (g)
        :Returns Type:
            :class:`float`
        """
        return water_influx - water_influx_next_internode - (transpiration * parameters.WATER_MOLAR_MASS * 1E-3 * delta_t)

    def calculate_extensibility(self, age, delta_t):
        """ Sheath extensibility in each dimension in relation to non-reversible dimensional changes.

        :Parameters:
            - `age` (:class:`float`) - sheath age (°Cd)
            - `delta_t` (:class:`float`) - Time step of the simulation (s).

        :Returns:
            Dict of extensibility (MPa-1)). Keys = ['z', 'x'].
        :Returns Type:
            :class:`dict`
        """
        if age <= SheathElement.PARAMETERS.tend:
            beta_function_norm = (1 - (1 + (SheathElement.PARAMETERS.tend - age) / (SheathElement.PARAMETERS.tend - SheathElement.PARAMETERS.tmax)) * (age / SheathElement.PARAMETERS.tend) **
                                  (SheathElement.PARAMETERS.tend / (SheathElement.PARAMETERS.tend - SheathElement.PARAMETERS.tmax)))
        else:
            beta_function_norm = 0

        phi = {}
        for phi_init_dimensions, phi_init_value in SheathElement.PARAMETERS.phi_initial.items():
            phi[phi_init_dimensions] = phi_init_value * beta_function_norm * delta_t
        return phi

    def calculate_delta_turgor_water_potential(self, phi, turgor_water_potential, organ_dimensions, delta_water_content):
        """ Delta of turgor water potential according to organ water content, turgor water potential, dimensions, extensibility and plasticity.

        :Parameters:
            - `phi` (:class:`dict`) - dict of cell wall extensibility (MPa). Keys = ['z', 'x']
            - `turgor_water_potential` (:class:`float`) - MPa
            - `organ_dimensions` (:class:`dict`) - dict of organ dimensions at time t. Keys = ['length', 'radius'] (m)
            - `delta_water_content` (:class:`float`) - delta water content integrated over delta t (g)

        :Returns:
            Delta of turgor water potential (MPa)
        :Returns Type:
            :class:`float`
        """
        epsilon_z, epsilon_x = SheathElement.PARAMETERS.epsilon['z'], SheathElement.PARAMETERS.epsilon['x']
        elastic_component = (epsilon_z * epsilon_x) / (2 * epsilon_z + epsilon_x)   #: Elastic reversible growth (MPa)
        plastic_component = (phi['z'] + 2 * phi['x']) * (max(turgor_water_potential, SheathElement.PARAMETERS.gamma) - SheathElement.PARAMETERS.gamma)  #: Plastic irreversible growth
        organ_volume = math.pi * organ_dimensions['radius'] ** 2 * organ_dimensions['length']  #: (m3)
        delta_turgor_water_potential = ((1 / (parameters.rho_water * organ_volume)) * delta_water_content * 1E-3 - plastic_component) * elastic_component  #: (MPa)

        return delta_turgor_water_potential

    def calculate_delta_organ_dimensions(self, delta_turgor_water_potential, turgor_water_potential, phi, organ_dimensions):
        """ Delta of organ dimensions according to turgor water potential, dimensions, extensibility and plasticity.

        :Parameters:
            - `delta_turgor_water_potential` (:class:`float`) - delta of turgor water potential integrated over delta t (MPa)
            - `turgor_water_potential` (:class:`float`) - MPa
            - `phi` (:class:`dict`) - dict of cell wall extensibility (MPa). Keys = ['z', 'x']
            - `organ_dimensions` (:class:`dict`) - dict of organ dimensions at time t. Keys = ['length', 'radius'] (m)

        :Returns:
            Delta of organ specific-dimensions (m). Keys = ['length', 'radius']
        :Returns Type:
            :class:`dict`
        """
        delta_organ_dimensions = {}
        epsilon_dict = SheathElement.PARAMETERS.epsilon.items()
        mapping_dimensions = {'x': 'radius', 'z': 'length'}

        for epsilon_dimension, epsilon_value in epsilon_dict:
            delta_organ_dimensions[mapping_dimensions[epsilon_dimension]] = ((1 / epsilon_value) * delta_turgor_water_potential +
                                                                             phi[epsilon_dimension] * (max(turgor_water_potential, LaminaElement.PARAMETERS.gamma) - LaminaElement.PARAMETERS.gamma)) *\
                                                                            organ_dimensions[mapping_dimensions[epsilon_dimension]]
        return delta_organ_dimensions
