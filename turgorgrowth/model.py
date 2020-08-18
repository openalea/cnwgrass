# -*- coding: latin-1 -*-
"""
    turgorgrowth.model
    ~~~~~~~~~~~~~~~~~~

    The module :mod:`turgorgrowth.model` defines the equations of wtaer flow, turgor pressure and growth.

    :license: CeCILL-C, see LICENSE for details.

"""

from __future__ import division  # use "//" to do integer division
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


class Axis(object):
    """
    The class :class:`Axis`.

    An :class:`axis <Axis>` must have:
        * one :class:`root compartment <Roots>`,
        * at least one :class:`phytomer<Phytomer>`.
    """

    PARAMETERS = parameters.AXIS_PARAMETERS  #: the internal parameters of the axes

    def __init__(self, label=None, roots=None, phytomers=None):
        """
        :param str label: the label of the axis
        :param Roots roots: Root object
        :param list [Phytomer] phytomers: list of Phytomer objects
        """
        self.label = label
        self.roots = roots
        if phytomers is None:
            phytomers = []
        self.phytomers = phytomers  #: the list of phytomers


class Roots(object):
    """
    The class :class:`Roots`.
    """

    def __init__(self, label='roots', total_water_potential=-0.1):
        """
        :param str label: root label
        :param float total_water_potential: Water potential of roots (MPa)
        """
        self.label = label
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


class HiddenZone(Organ):
    """
    The class :class:`HiddenZone`.
    """

    PARAMETERS = parameters.HIDDEN_ZONE_PARAMETERS                #: the internal parameters of the hidden zone
    INIT_COMPARTMENTS = parameters.HIDDEN_ZONE_INIT_COMPARTMENTS  #: the initial values of compartments and state parameters

    def __init__(self, label='hiddenzone', leaf_pseudo_age=INIT_COMPARTMENTS.age, amino_acids=INIT_COMPARTMENTS.amino_acids, leaf_L=INIT_COMPARTMENTS.leaf_L,
                 leaf_pseudostem_length=INIT_COMPARTMENTS.leaf_pseudostem_length, proteins=INIT_COMPARTMENTS.proteins, sucrose=INIT_COMPARTMENTS.sucrose, temperature=INIT_COMPARTMENTS.temperature,
                 thickness=INIT_COMPARTMENTS.thickness, turgor_water_potential=INIT_COMPARTMENTS.turgor_water_potential, water_content=INIT_COMPARTMENTS.water_content,
                 mstruct=INIT_COMPARTMENTS.mstruct, osmotic_water_potential=INIT_COMPARTMENTS.osmotic_water_potential, total_water_potential=INIT_COMPARTMENTS.total_water_potential,
                 width=INIT_COMPARTMENTS.width):

        super(HiddenZone, self).__init__(label)
        self.label = label

        # state parameters
        self.amino_acids = amino_acids           #: µmol N
        self.proteins = proteins                 #: µmol N
        self.sucrose = sucrose                   #: :math:`\mu mol C
        self.temperature = temperature           #: °C
        self.leaf_pseudo_age = leaf_pseudo_age   #: °Cd
        self.leaf_L = leaf_L                     #: m
        self.mstruct = mstruct                   #: g

        # state variables
        self.leaf_pseudostem_length = leaf_pseudostem_length    #: m
        self.width = width                                      #: m
        self.thickness = thickness                              #: m
        self.turgor_water_potential = turgor_water_potential    #: MPa
        self.water_content = water_content                      #: g

        # fluxes
        self.water_influx = None  #: water influx into organ (g)
        self.resistance = None               #: resistance of water flux between two organs (MPa s g-1)

        # Intermediate variables
        self.length = min(leaf_L, leaf_pseudostem_length)
        self.osmotic_water_potential = osmotic_water_potential  #: MPa
        self.total_water_potential = total_water_potential    #: MPa

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
        water_content = volume * parameters.rho_water  #: g
        return volume, water_content

    @staticmethod
    def calculate_volume(water_content):
        """ Hidden zone volume, assumed to be proportional to water content.

        :param float water_content: (g)

        :return: volume (m3)
        :rtype: float
        """
        return water_content / parameters.rho_water

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
        temperature_K = temperature + parameters.celsius_2_kelvin

        sucrose = ((sucrose * 1E-6) / parameters.NB_C_SUCROSE) * parameters.VANT_HOFF_SUCROSE
        amino_acids = ((amino_acids * 1E-6) / parameters.AMINO_ACIDS_N_RATIO) * parameters.VANT_HOFF_AMINO_ACIDS
        proteins_actual = proteins * min(1., (age * 5E-6 + 0.1))  # TODO: temp hack to account for the fact that N is mainly Nstruct in young hz
        proteins = ((proteins_actual * 1E-6) / parameters.AMINO_ACIDS_N_RATIO) * parameters.VANT_HOFF_AMINO_ACIDS

        osmotic_water_potential = - parameters.R * temperature_K * ((sucrose + amino_acids + proteins) / (volume * Organ.PARAMETERS.vstorage)) * 1E-6
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
    def calculate_resistance(hiddenzone_dimensions, previous_organ):
        """ Resistance of water flow between the hiddenzone and its predecessor (can be root or the last elongated internode if any).
        Relations were set proportional to the length and inversely proportional to the area of organ's cross section.

        :param dict hiddenzone_dimensions: dict of hidden zone dimensions at time t. Keys = ['length', 'thickness', 'width] (m)
        :param Roots or InternodeElement previous_organ: predecessor of the hiddenzone, either root or internode element with length and radius attributes

        :return: resistance (MPa s g-1)
        :rtype: float
        """
        if previous_organ.label == 'roots':
            resistance = 0.5 * HiddenZone.PARAMETERS.R_xylem * (hiddenzone_dimensions['length'] / (hiddenzone_dimensions['width'] * hiddenzone_dimensions['thickness']))
        else:
            resistance = 0.5 * HiddenZone.PARAMETERS.R_xylem * (previous_organ.length / (previous_organ.width * previous_organ.thickness) +
                                                                hiddenzone_dimensions['length'] / (hiddenzone_dimensions['width'] * hiddenzone_dimensions['thickness']))
        return resistance

    @staticmethod
    def calculate_water_influx(organ_water_potential, previous_organ_water_potential, resistance, delta_t):
        """ Sap flow into the organ according to water potential gradient with the previous organ
        For hiddenzones the previous organ is the root compartment if located above short internodes, else it is the previous internode

        :param float organ_water_potential: water potential of the current organ (MPa)
        :param float previous_organ_water_potential: water potential of the previous organ (MPa)
        :param float resistance: transport resistance between organs (MPa s g-1)
        :param float delta_t: time step of the simulation (s)

        :return: Water influx into the current organ integrated over delta_t (g)
        :rtype: float
        """
        return ((previous_organ_water_potential - organ_water_potential) / resistance) * delta_t

    @staticmethod
    def calculate_delta_water_content(water_influx, water_outflow=0):
        """ delta of sap flow for the hidden zone.

        :param float water_influx: Water influx integrated over delta_t (g)
        :param float water_outflow: Water loss through the emerged lamina or sheath if any (g)

        :return: Delta of sap flow into the organ (g)
        :rtype: float
        """
        return water_influx - water_outflow

    @staticmethod
    def calculate_extensibility(age, delta_t):

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
        epsilon_x, epsilon_y, epsilon_z = [1000*water_content]*3  # #HiddenZone.PARAMETERS.epsilon['x'], HiddenZone.PARAMETERS.epsilon['y'], HiddenZone.PARAMETERS.epsilon['z']
        elastic_component = (epsilon_x * epsilon_y * epsilon_z) / (epsilon_z * epsilon_x + epsilon_z * epsilon_y + epsilon_x * epsilon_y)   #: Elastic reversible growth (MPa)
        plastic_component = (phi['x'] + phi['y'] + phi['z']) * (max(turgor_water_potential, HiddenZone.PARAMETERS.gamma) - HiddenZone.PARAMETERS.gamma)  #: Plastic irreversible growth
        delta_turgor_water_potential = ((1 / (parameters.rho_water * volume)) * delta_water_content - plastic_component) * elastic_component  #: (MPa)

        return delta_turgor_water_potential

    @staticmethod
    def calculate_delta_organ_dimensions(delta_turgor_water_potential, turgor_water_potential, phi, organ_dimensions, water_content):
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
        # epsilon_dict = HiddenZone.PARAMETERS.epsilon
        epsilon_dict = {'x': water_content*1000, 'y': water_content*1000, 'z': water_content*1000}
        mapping_dimensions = {'x': 'width', 'y': 'thickness', 'z': 'length'}

        for epsilon_dimension, epsilon_value in epsilon_dict.items():
            delta_organ_dimensions[mapping_dimensions[epsilon_dimension]] = ((1 / epsilon_value) * delta_turgor_water_potential +
                                                                             phi[epsilon_dimension] * (max(turgor_water_potential, HiddenZone.PARAMETERS.gamma) - HiddenZone.PARAMETERS.gamma)) *\
                                                                            organ_dimensions[mapping_dimensions[epsilon_dimension]]
        return delta_organ_dimensions


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
        """
        :param str label: Photosynthetic organ label
        :param LaminaElement or InternodeElement or SheathElement exposed_element: the exposed element
        :param LaminaElement or InternodeElement or SheathElement enclosed_element: the enclosed element
        """
        super(PhotosyntheticOrgan, self).__init__(label)

        self.exposed_element = exposed_element
        self.enclosed_element = enclosed_element


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

    A :class:`photosynthetic organ element <PhotosyntheticOrganElement>` must have at least 1
    :class:`lamina element<LaminaElement>`,
    :class:`internode element <InternodeElement>`,
    or :class:`sheath element <SheathElement>`.

    :class:`PhotosyntheticOrganElement` is the base class of all photosynthetic organ elements. DO NOT INSTANTIATE IT.
    """

    INIT_COMPARTMENTS = parameters.PHOTOSYNTHETIC_ORGAN_ELEMENT_INIT_COMPARTMENTS  #: the initial values of compartments and state parameters

    def __init__(self, label=None, age=INIT_COMPARTMENTS.age, amino_acids=INIT_COMPARTMENTS.amino_acids, green_area=INIT_COMPARTMENTS.green_area, length=INIT_COMPARTMENTS.length,
                 mstruct=INIT_COMPARTMENTS.mstruct, osmotic_water_potential=INIT_COMPARTMENTS.osmotic_water_potential, proteins=INIT_COMPARTMENTS.proteins, sucrose=INIT_COMPARTMENTS.sucrose,
                 Ts=INIT_COMPARTMENTS.temperature, total_water_potential=INIT_COMPARTMENTS.total_water_potential, Tr=INIT_COMPARTMENTS.Tr,
                 turgor_water_potential=INIT_COMPARTMENTS.turgor_water_potential, water_content=INIT_COMPARTMENTS.water_content):

        self.label = label                                      #: the label of the element

        # state parameters
        self.age = age                                          #: °Cd
        self.amino_acids = amino_acids                          #: µmol N
        self.green_area = green_area                            #: m2
        self.mstruct = mstruct                                  #: g
        self.proteins = proteins                                #: µmol N
        self.sucrose = sucrose                                  #: :math:`\mu mol C
        self.temperature = Ts                                   #: °C
        self.Tr = Tr                                            #: mmol H20 m-2 s-1

        # state variables
        self.length = length                                    #: m
        self.turgor_water_potential = turgor_water_potential    #: MPa
        self.water_content = water_content                      #: g

        # fluxes
        self.water_influx = None                                #: water influx into organ (g)
        self.resistance = None                                  #: resistance of water flux between two organs (MPa s g-1)

        # Intermediate variables
        self.osmotic_water_potential = osmotic_water_potential  #: MPa
        self.total_water_potential = total_water_potential      #: MPa

    #: Water flow equations common to all photosynthetic organ elements

    @staticmethod
    def calculate_volume(water_content):
        """ Photosynthetic element volume, assumed to be proportional to water content.

        :param float water_content: (g)

        :return: volume (m3)
        :rtype: float
        """
        return water_content / parameters.rho_water

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
    def calculate_water_influx(organ_water_potential, previous_organ_water_potential, resistance, delta_t):
        """ Sap flow into the organ according to water potential gradient with the previous organ
        For lamina, the previous organ is the hiddenzone if leaf not mature, else it is the sheath
        For sheaths, the previous organ is the root compartment if located above short internodes, else it is the previous internode
        For internodes the previous organ is the previous one else the organ is neglected.

        :param float organ_water_potential: water potential of the current organ (MPa)
        :param float previous_organ_water_potential: water potential of the previous organ (MPa)
        :param float resistance: transport resistance between organs (MPa s g-1)
        :param float delta_t: time step of the simulation (s)

        :return: Water influx into the current organ integrated over delta_t (g)
        :rtype: float
        """
        return ((previous_organ_water_potential - organ_water_potential) / resistance) * delta_t


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

    # Model equations for water flux

    @staticmethod
    def calculate_initial_water_content(hiddenzone_osmotic_water_potential, sucrose, amino_acids, proteins, temperature):
        """ Initial water content of photosynthetic element at emergence
        calculated in order that the resulting osmotic pressure of the element will be similar to that of the hidden zone

        :param float hiddenzone_osmotic_water_potential: the osmotic pressure of the hidden zone (MPa)
        :param float sucrose: amount of sucrose in the element (µmol C)
        :param float amino_acids: amount of amino acids in the element (µmol N)
        :param float proteins: amount of proteins in the element (µmol N)
        :param float temperature: element temperature (°C)

        :return: element initial water content (g)
        :rtype: float
        """
        temperature_K = temperature + parameters.celsius_2_kelvin

        sucrose = ((sucrose * 1E-6) / parameters.NB_C_SUCROSE) * parameters.VANT_HOFF_SUCROSE
        amino_acids = ((amino_acids * 1E-6) / parameters.AMINO_ACIDS_N_RATIO) * parameters.VANT_HOFF_AMINO_ACIDS
        proteins = ((proteins * 1E-6) / parameters.AMINO_ACIDS_N_RATIO) * parameters.VANT_HOFF_AMINO_ACIDS

        initial_water_content = (- parameters.R * temperature_K * (sucrose + amino_acids + proteins)) / (hiddenzone_osmotic_water_potential * LaminaElement.PARAMETERS.vstorage)
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
        temperature_K = temperature + parameters.celsius_2_kelvin

        sucrose = ((sucrose * 1E-6) / parameters.NB_C_SUCROSE) * parameters.VANT_HOFF_SUCROSE
        amino_acids = ((amino_acids * 1E-6) / parameters.AMINO_ACIDS_N_RATIO) * parameters.VANT_HOFF_AMINO_ACIDS
        non_chloroplastic_proteins = proteins / 1
        proteins = ((non_chloroplastic_proteins * 1E-6) / parameters.AMINO_ACIDS_N_RATIO) * parameters.VANT_HOFF_AMINO_ACIDS

        osmotic_water_potential = - parameters.R * temperature_K * ((sucrose + amino_acids + proteins) / (volume * LaminaElement.PARAMETERS.vstorage)) * 1E-6
        return osmotic_water_potential

    @staticmethod
    def calculate_resistance(lamina_dimensions, previous_element):
        """ Resistance of water flow between the lamina and hidden zone (if leaf is growing) or the sheath (if lamina is mature).
        Relations were set proportional to the length and inversely proportional to the area of organ's cross section.

        :param dict lamina_dimensions: dict of lamina dimensions (m). Keys = ['length', 'width', 'thickness']
        :param HiddenZone or SheathElement previous_element: predecessor of the lamina, either sheath or hiddenzone

        :return: resistance (MPa s g-1)
        :rtype: float
        """
        resistance = 0.5 * LaminaElement.PARAMETERS.R_xylem * (previous_element.length / (previous_element.width * previous_element.thickness) +
                                                               lamina_dimensions['length'] / (lamina_dimensions['width'] * lamina_dimensions['thickness']))
        return resistance

    @staticmethod
    def calculate_delta_water_content(water_influx, Tr, green_area, delta_t):
        """ Delta of sap flow for the lamina.

        :param float water_influx: Water influx integrated over delta_t (g)
        :param float Tr: Lamina surfacic transpiration rate (mmol H20 m-2 s-1)
        :param float green_area: Lamina sgreen area (m2)
        :param float delta_t: Time step of the simulation (s)

        :return: Delta of sap flow into the organ (g)
        :rtype: float
        """
        total_transpiration = Tr * green_area  # mmol H20 s-1
        return water_influx - (total_transpiration * parameters.WATER_MOLAR_MASS * 1E-3 * delta_t)

    @staticmethod
    def calculate_delta_turgor_water_potential(organ_dimensions, delta_water_content):
        """ Delta of turgor water potential according to organ water content, turgor water potential, dimensions and elasticity.
        Extensibility (psi) is supposed to be 0 as this tissue is mature (growth completed).

        :param dict organ_dimensions: dict of organ dimensions at time t. Keys = ['length', 'thickness', 'width'] (m)
        :param float delta_water_content: delta water content integrated over delta t (g)

        :return: Delta of turgor water potential (MPa)
        :rtype: float
        """
        epsilon_z, epsilon_x, epsilon_y = LaminaElement.PARAMETERS.epsilon['z'], LaminaElement.PARAMETERS.epsilon['x'], LaminaElement.PARAMETERS.epsilon['y']
        elastic_component = (epsilon_z * epsilon_x * epsilon_y) / (epsilon_z * epsilon_x + epsilon_z * epsilon_y + epsilon_x * epsilon_y)  #: Elastic reversible growth (MPa)
        organ_volume = organ_dimensions['length'] * organ_dimensions['width'] * organ_dimensions['thickness']  #: (m3)
        delta_turgor_water_potential = ((1 / (parameters.rho_water * organ_volume)) * delta_water_content) * elastic_component  #: (MPa)

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
        epsilon_dict = LaminaElement.PARAMETERS.epsilon.items()
        mapping_dimensions = {'x': 'width', 'y': 'thickness', 'z': 'length'}

        for epsilon_dimension, epsilon_value in epsilon_dict:
            delta_organ_dimensions[mapping_dimensions[epsilon_dimension]] = ((1 / epsilon_value) * delta_turgor_water_potential) * organ_dimensions[mapping_dimensions[epsilon_dimension]]
        return delta_organ_dimensions


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

    # Model equations for water flux

    @staticmethod
    def calculate_initial_water_content(hiddenzone_osmotic_water_potential, sucrose, amino_acids, proteins, temperature):
        """ Initial water content of photosynthetic element at emergence
        calculated in order that the resulting osmotic pressure of the element will be similar to that of the hidden zone

        :param float hiddenzone_osmotic_water_potential: the osmotic pressure of the hidden zone (MPa)
        :param float sucrose: amount of sucrose in the element (µmol C)
        :param float amino_acids: amount of amino acids in the element (µmol N)
        :param float proteins: amount of proteins in the element (µmol N)
        :param float temperature: element temperature (°C)

        :return: element initial water content (g)
        :rtype: float
        """
        temperature_K = temperature + parameters.celsius_2_kelvin

        sucrose = ((sucrose * 1E-6) / parameters.NB_C_SUCROSE) * parameters.VANT_HOFF_SUCROSE
        amino_acids = ((amino_acids * 1E-6) / parameters.AMINO_ACIDS_N_RATIO) * parameters.VANT_HOFF_AMINO_ACIDS
        proteins = ((proteins * 1E-6) / parameters.AMINO_ACIDS_N_RATIO) * parameters.VANT_HOFF_AMINO_ACIDS

        initial_water_content = (- parameters.R * temperature_K * (sucrose + amino_acids + proteins)) / (hiddenzone_osmotic_water_potential * InternodeElement.PARAMETERS.vstorage)
        return initial_water_content

    @staticmethod
    def calculate_initial_water_potential(resistance_dict):
        """ Initial water potential of emerging internode calculated in order that the influx from the hiddenzone matches the outflux towards the next organ

        :param dict resistance_dict: a dictionary with the resistances of each boundary element

        :return: element initial water potential (MPa)
        :rtype: float
        """
        numerator, denominator = 0, 0
        for element in resistance_dict.keys():
            other_element_resistances = sum([resistance_dict[other_element] for other_element in resistance_dict.keys() if other_element != element])
            numerator += element.total_water_potential * other_element_resistances
            denominator += other_element_resistances

        initial_water_potential = numerator / denominator
        return initial_water_potential

    @staticmethod
    def calculate_osmotic_water_potential(sucrose, amino_acids, proteins, volume, temperature):
        """ Osmotic water potential of the organ calculated according to metabolites

        :param float sucrose: µmol C under the form of sucrose
        :param float amino_acids: µmol N under the form of amino acids
        :param float proteins: µmol N under the form of proteins
        :param float volume: (m3)
        :param float temperature: degree Celsius

        :return: Osmotic water potential (MPa)
        :rtype: float
        """
        temperature_K = temperature + parameters.celsius_2_kelvin

        sucrose = ((sucrose * 1E-6) / parameters.NB_C_SUCROSE) * parameters.VANT_HOFF_SUCROSE
        amino_acids = ((amino_acids * 1E-6) / parameters.AMINO_ACIDS_N_RATIO) * parameters.VANT_HOFF_AMINO_ACIDS
        non_chloroplastic_proteins = proteins / 1
        proteins = ((non_chloroplastic_proteins * 1E-6) / parameters.AMINO_ACIDS_N_RATIO) * parameters.VANT_HOFF_AMINO_ACIDS

        osmotic_water_potential = - parameters.R * temperature_K * ((sucrose + amino_acids + proteins) / (volume * InternodeElement.PARAMETERS.vstorage)) * 1E-6
        return osmotic_water_potential

    @staticmethod
    def calculate_resistance(internode_dimensions, previous_organ):
        """ Resistance of water flow between the internode and the root compartment or the previous internode if elongated.
            Relations were set proportional to the length and inversely proportional to the area of organ's cross section.

        :param dict internode_dimensions: dict of current internode dimensions. Keys = ['length', 'width', 'thickness'] (m)
        :param Roots or InternodeElement previous_organ: predecessor dimensions of the internode, either root or internode element

        :return: resistance (MPa s g-1)
        :rtype: float
        """
        if previous_organ.label == 'roots':
            resistance = 0.5 * InternodeElement.PARAMETERS.R_xylem * (internode_dimensions['length'] / (internode_dimensions['width'] * internode_dimensions['thickness']))
        else:
            resistance = 0.5 * InternodeElement.PARAMETERS.R_xylem * (internode_dimensions['length'] / (internode_dimensions['width'] * internode_dimensions['thickness']) +
                                                                      previous_organ.length / (previous_organ.width * previous_organ.thickness))
        return resistance

    @staticmethod
    def calculate_delta_water_content(water_influx, water_influx_next_internode, Tr, green_area, delta_t):
        """ delta of sap flow for the internode.

        :param float water_influx: Water influx integrated of the current internode over delta_t (g)
        :param float water_influx_next_internode: Water influx of the next organ (if any) integrated over delta_t (g)
        :param float Tr: Internode surfacic transpiration rate (mmol H20 m-2 s-1)
        :param float green_area: Internode green area (m2)
        :param float delta_t: Time step of the simulation (s)

        :return: Delta of sap flow into the organ (g)
        :rtype: float
        """
        total_transpiration = Tr * green_area  # m2
        return water_influx - water_influx_next_internode - (total_transpiration * parameters.WATER_MOLAR_MASS * 1E-3 * delta_t)

    @staticmethod
    def calculate_delta_turgor_water_potential(organ_dimensions, delta_water_content):
        """ Delta of turgor water potential according to organ water content, turgor water potential, dimensions and elasticity.
        Extensibility (psi) is supposed to be 0 as this tissue is mature (growth completed).

        :param dict organ_dimensions: dict of organ dimensions at time t. Keys = ['length', 'width', 'thickness'] (m)
        :param float delta_water_content: delta water content integrated over delta t (g)

        :return: Delta of turgor water potential (MPa)
        :rtype: float
        """
        epsilon_z, epsilon_x, epsilon_y = InternodeElement.PARAMETERS.epsilon['z'], InternodeElement.PARAMETERS.epsilon['x'], InternodeElement.PARAMETERS.epsilon['y']
        elastic_component = (epsilon_z * epsilon_x * epsilon_y) / (epsilon_z * epsilon_x + epsilon_z * epsilon_y + epsilon_x * epsilon_y)  #: Elastic reversible growth (MPa)
        organ_volume = organ_dimensions['length'] * organ_dimensions['width'] * organ_dimensions['thickness']  #: (m3)
        delta_turgor_water_potential = ((1 / (parameters.rho_water * organ_volume)) * delta_water_content) * elastic_component  #: (MPa)

        return delta_turgor_water_potential

    @staticmethod
    def calculate_delta_organ_dimensions(delta_turgor_water_potential, organ_dimensions):
        """Delta of internode dimensions according to turgor water potential, dimensions, and elasticity

        :param float delta_turgor_water_potential: delta of turgor water potential integrated over delta t (MPa)
        :param dict organ_dimensions: dict of organ dimensions at time t. Keys = ['length', 'width', 'thickness'] (m)

        :return: Delta of organ specific-dimensions (m). Keys = ['length', 'width', 'thickness']
        :rtype: dict
        """

        delta_organ_dimensions = {}
        epsilon_dict = InternodeElement.PARAMETERS.epsilon.items()
        mapping_dimensions = {'x': 'width', 'y': 'thickness', 'z': 'length'}

        for epsilon_dimension, epsilon_value in epsilon_dict:
            delta_organ_dimensions[mapping_dimensions[epsilon_dimension]] = ((1 / epsilon_value) * delta_turgor_water_potential) * organ_dimensions[mapping_dimensions[epsilon_dimension]]
        return delta_organ_dimensions


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

    # Model equations for water flux

    @staticmethod
    def calculate_initial_water_content(hiddenzone_osmotic_water_potential, sucrose, amino_acids, proteins, temperature):
        """ Initial water content of photosynthetic element at emergence
        calculated in order that the resulting osmotic pressure of the element will be similar to that of the hidden zone

        :param float hiddenzone_osmotic_water_potential: the osmotic pressure of the hidden zone (MPa)
        :param float sucrose: amount of sucrose in the element (µmol C)
        :param float amino_acids: amount of amino acids in the element (µmol N)
        :param float proteins: amount of proteins in the element (µmol N)
        :param float temperature: element temperature (°C)

        :return: element initial water content (g)
        :rtype: float
        """
        temperature_K = temperature + parameters.celsius_2_kelvin

        sucrose = ((sucrose * 1E-6) / parameters.NB_C_SUCROSE) * parameters.VANT_HOFF_SUCROSE
        amino_acids = ((amino_acids * 1E-6) / parameters.AMINO_ACIDS_N_RATIO) * parameters.VANT_HOFF_AMINO_ACIDS
        proteins = ((proteins * 1E-6) / parameters.AMINO_ACIDS_N_RATIO) * parameters.VANT_HOFF_AMINO_ACIDS

        initial_water_content = (- parameters.R * temperature_K * (sucrose + amino_acids + proteins)) / (hiddenzone_osmotic_water_potential * SheathElement.PARAMETERS.vstorage)
        return initial_water_content

    @staticmethod
    def calculate_initial_water_potential(resistance_dict):
        """ Initial water potential of emerging sheath calculated in order that the influx from the hiddenzone matches the outflux towards the next organ

        :param dict resistance_dict: a dictionary with the resistances of each boundary element

        :return: element initial water potential (MPa)
        :rtype: float
        """
        nominator, denominator = 0, 0
        for element in resistance_dict.keys():
            other_element_resistances = sum([resistance_dict[other_element] for other_element in resistance_dict.keys() if other_element != element])
            nominator += element.total_water_potential * other_element_resistances
            denominator += other_element_resistances

        initial_water_potential = nominator / denominator
        return initial_water_potential

    @staticmethod
    def calculate_osmotic_water_potential(sucrose, amino_acids, proteins, volume, temperature):
        """ Osmotic water potential of the organ calculated according to metabolites

        :param float sucrose: µmol C under the form of sucrose
        :param float amino_acids: µmol N under the form of amino acids
        :param float proteins: µmol N under the form of proteins
        :param float volume: (m3)
        :param float temperature: degree Celsius

        :return: Osmotic water potential (MPa)
        :rtype: float
        """
        temperature_K = temperature + parameters.celsius_2_kelvin

        sucrose = ((sucrose * 1E-6) / parameters.NB_C_SUCROSE) * parameters.VANT_HOFF_SUCROSE
        amino_acids = ((amino_acids * 1E-6) / parameters.AMINO_ACIDS_N_RATIO) * parameters.VANT_HOFF_AMINO_ACIDS
        non_chloroplastic_proteins = proteins / 1
        proteins = ((non_chloroplastic_proteins * 1E-6) / parameters.AMINO_ACIDS_N_RATIO) * parameters.VANT_HOFF_AMINO_ACIDS

        osmotic_water_potential = (- parameters.R * temperature_K * (sucrose + amino_acids + proteins)) / (volume * SheathElement.PARAMETERS.vstorage) * 1E-6
        return osmotic_water_potential

    @staticmethod
    def calculate_resistance(sheath_dimensions, previous_organ):
        """ Resistance of water flow between the sheath and the previous element which could be one of these: enclosed sheath element, hiddenzone, previous elongated internode or root compartment.
            Relations were set proportional to the length and inversely proportional to the area of organ's cross section.

        :param dict sheath_dimensions: dict of current hidden zone dimensions. Keys = ['length', 'width', 'thickness'] (m)
        :param SheathElement or HiddenZone or InternodeElement or Roots previous_organ: predecessor dimensions of the sheath, either sheath, internode, hiddenzone or root

        :return: resistance (MPa s g-1).
        :rtype: float
        """
        if isinstance(previous_organ, Roots):
            resistance = 0.5 * HiddenZone.PARAMETERS.R_xylem * (sheath_dimensions['length'] / (sheath_dimensions['width'] * sheath_dimensions['thickness']))
        else:
            resistance = 0.5 * SheathElement.PARAMETERS.R_xylem * (sheath_dimensions['length'] / (sheath_dimensions['width'] * sheath_dimensions['thickness']) +
                                                                   previous_organ.length / (previous_organ.width * previous_organ.thickness))
        return resistance

    @staticmethod
    def calculate_delta_water_content(water_influx, water_influx_next_element, Tr, green_area, delta_t):
        """ Delta of sap flow for the sheath.

        :param float water_influx: Water influx integrated of the current sheath over delta_t (g)
        :param float water_influx_next_element: Water influx of the next element (if any) integrated over delta_t (g)
        :param float Tr: Sheath surfacic transpiration rate (mmol H20 m-2 s-1)
        :param float green_area: Sheath green area (m2)
        :param float delta_t: Time step of the simulation (s)

        :return: Delta of sap flow into the organ (g)
        :rtype: float
        """
        total_transpiration = Tr * green_area  # m2
        return water_influx - water_influx_next_element - (total_transpiration * parameters.WATER_MOLAR_MASS * 1E-3 * delta_t)

    @staticmethod
    def calculate_delta_turgor_water_potential(organ_dimensions, delta_water_content):
        """ Delta of turgor water potential according to organ water content, turgor water potential, dimensions and elasticity.
        Extensibility (psi) is supposed to be 0 as this tissue is mature (growth completed).

        :param dict organ_dimensions: dict of organ dimensions at time t. Keys = ['length', 'width', 'thickness'] (m)
        :param float delta_water_content: delta water content integrated over delta t (g)

        :return: Delta of turgor water potential (MPa)
        :rtype: float
        """
        epsilon_z, epsilon_x, epsilon_y = SheathElement.PARAMETERS.epsilon['z'], SheathElement.PARAMETERS.epsilon['x'], SheathElement.PARAMETERS.epsilon['y']
        elastic_component = (epsilon_z * epsilon_x * epsilon_y) / (epsilon_z * epsilon_x + epsilon_z * epsilon_y + epsilon_x * epsilon_y)  #: Elastic reversible growth (MPa)
        organ_volume = organ_dimensions['length'] * organ_dimensions['width'] * organ_dimensions['thickness']  #: (m3)
        delta_turgor_water_potential = ((1 / (parameters.rho_water * organ_volume)) * delta_water_content) * elastic_component  #: (MPa)

        return delta_turgor_water_potential

    @staticmethod
    def calculate_delta_organ_dimensions(delta_turgor_water_potential, organ_dimensions):
        """Delta of sheath dimensions according to turgor water potential, dimensions, and elasticity

        :param float delta_turgor_water_potential: delta of turgor water potential integrated over delta t (MPa)
        :param dict organ_dimensions: dict of organ dimensions at time t. Keys = ['length', 'width', 'thickness'] (m)

        :return: Delta of organ specific-dimensions (m). Keys = ['length', 'width', 'thickness']
        :rtype: dict
        """
        delta_organ_dimensions = {}
        epsilon_dict = SheathElement.PARAMETERS.epsilon.items()
        mapping_dimensions = {'x': 'width', 'y': 'thickness', 'z': 'length'}

        for epsilon_dimension, epsilon_value in epsilon_dict:
            delta_organ_dimensions[mapping_dimensions[epsilon_dimension]] = ((1 / epsilon_value) * delta_turgor_water_potential) * organ_dimensions[mapping_dimensions[epsilon_dimension]]
        return delta_organ_dimensions
