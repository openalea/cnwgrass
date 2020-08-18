"""
    turgorgrowth.parameters
    ~~~~~~~~~~~~~~~~~~

    The module :mod:`turgorgrowth.parameters` defines the parameters of the model.

    :license: CeCILL-C, see LICENSE for details.
"""

HOUR_TO_SECOND_CONVERSION_FACTOR = 3600.  #: Number of seconds in 1 hour

celsius_2_kelvin = 273.15  #: conversion factor from degree Celsius to Kelvin
R = 8.31  #: Perfect gas constant (J mol-1 K-1)
rho_water = 1E6  #: Water density (g m-3)

NB_C_SUCROSE = 12  #: Number of C in 1 mol of sucrose
SUCROSE_MOLAR_MASS = 342  #: g mol-1
AMINO_ACIDS_N_RATIO = 1.17  #: Mean number of mol of N in 1 mol of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly)
WATER_MOLAR_MASS = 18  #: g mol-1

VANT_HOFF_SUCROSE = 1  #: Van't Hoff coefficient of sucrose (dimensionless)
VANT_HOFF_AMINO_ACIDS = 1.25  #: Van't Hoff coefficient estimated for amino acids (dimensionless)


class PopulationParameters(object):
    """
    Internal parameters of populations.
    """
    def __init__(self):
        pass


#: The instance of class :class:`turgorgrowth.parameters.PopulationParameters` for current process
POPULATION_PARAMETERS = PopulationParameters()


class PlantParameters(object):
    """
    Internal parameters of plants.
    """
    def __init__(self):
        pass


#: The instance of class :class:`turgorgrowth.parameters.PlantParameters` for current process
PLANT_PARAMETERS = PlantParameters()


class AxisParameters(object):
    """
    Internal parameters of axes.
    """
    def __init__(self):
        pass


#: The instance of class :class:`turgorgrowth.parameters.AxisParameters` for current process
AXIS_PARAMETERS = AxisParameters()


class PhytomerParameters(object):
    """
    Internal parameters of phytomers.
    """
    def __init__(self):
        pass


#: The instance of class :class:`turgorgrowth.parameters.PhytomerParameters` for current process
PHYTOMER_PARAMETERS = PhytomerParameters()


class OrganParameters(object):
    """
    Internal parameters of organs.
    """
    def __init__(self):
        self.R_xylem = 2                                     #: Elemental xylem flow resistence (MPa s g-1 m)
        self.vstorage = 0.8                                     #: Storage portion of the organ volume (0.2 attributed to the xylem)
        self.gamma = 0.3                                        #: Critical value for the pressure component which must be exceeded for irreversible volume changes (MPa)


#: The instance of class :class:`turgorgrowth.parameters.PhytomerParameters` for current process
ORGAN_PARAMETERS = OrganParameters()


class HiddenZoneParameters(OrganParameters):
    """
    Internal parameters of hidden growing zones.
    """
    def __init__(self):
        super(HiddenZoneParameters, self).__init__()

        self.epsilon = {'x': 50, 'y': 50, 'z': 50}  #: 0.9 Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        self.phi_initial = {'x': 1E-09, 'y': 1E-09, 'z': 0.4E-05}  #: Initial dimension-specific cell wall extensibility in relation to non-reversible dimensional changes (MPa-1 s-1)
        self.tend = 2500000  # 300 * 3600 * 24 / 12   #: Lamina age when extensibility reaches 0 (s at 12°c). Calculated from elongwheat parameter for phase 2
        self.tmax = 2000000  # 190 * 3600 * 24 / 12   #: Lamina age when organ extensibility is reduced by half of the initial value (s at 12°c). Calculated from elongwheat parameter for phase 2
        self.tbase = -25 * 3600 * 24 / 12  #: beginning of leaf elongation in automate growth (s at 12°c); fitted from adapted data from Fournier 2005
        L0 = abs((1 + (self.tend / (self.tend - self.tmax))) * (min(1.0, float(-self.tbase) / float(self.tend - self.tbase)) ** ((self.tend - self.tbase) / (self.tend - self.tmax))))  #: Leaf length at t=0 in automate growth (beta function) (m)
        FITTED_L0 = 0.01557936  #: Fitted value of leaf length at t=0 after rescaling the beta function with L0 (m); Fournier 2005 sur courbe corrigee
        self.OFFSET_LEAF = FITTED_L0 - L0  #: Offset used for the final fitting of the beta function (m)
        self.RATIO_MSTRUCT_DM = 0.8  #: Ratio mstruct/dry matter (dimensionless). From growthwheat model.
        self.SLOPE_MASS_VOLUME = 3.23337E-06  #: Slope of the relation between leaf dry mass and its volume at the time of the previous leaf emergence (m3 g-1). Found from Williams 1960, Fig 11.
        self.OFFSET_MASS_VOLUME = 1.82312E-13  #: Offset of the relation between leaf dry mass and its volume at the time of the previous leaf emergence (m3). Found from Williams 1960, Fig 11.


#: The instance of class :class:`turgorgrowth.parameters.HiddenZoneParameters` for current process
HIDDEN_ZONE_PARAMETERS = HiddenZoneParameters()


class HiddenZoneInitCompartments(object):
    """
    Initial values for compartments of hidden zones.
    """
    def __init__(self):
        # state parameters
        self.age = None  #: °Cd
        self.amino_acids = 0.000075  #: :math:`\mu mol N
        self.proteins = 0.0011  #: :math:`\mu mol N
        self.sucrose = 0.000384  #: :math:`\mu mol C
        self.temperature = 20  #: °C
        self.mstruct = 2.65E-08    #: g TODO: viens d'ou?

        # intermediate variables
        self.osmotic_water_potential = -0.8  #: MPa
        self.total_water_potential = -0.1    #: MPa

        # state variables
        self.leaf_pseudostem_length = 4E-6   #: m
        self.leaf_L = 4E-6                   #: m
        self.thickness = 1E-3                #: m
        self.width = 0.01                    #: m
        self.turgor_water_potential = 0.7    #: MPa
        self.water_content = -R * (self.temperature + celsius_2_kelvin) * (
                ((self.sucrose * 1E-6 / NB_C_SUCROSE) * VANT_HOFF_SUCROSE) +
                ((self.amino_acids * 1E-6 / AMINO_ACIDS_N_RATIO) * VANT_HOFF_AMINO_ACIDS) +
                ((self.proteins * 1E-6 / AMINO_ACIDS_N_RATIO) * VANT_HOFF_AMINO_ACIDS)) / self.osmotic_water_potential  #: g


#: The instance of class :class:`turgorgrowth.parameters.HiddenZoneInitCompartments` for current process
HIDDEN_ZONE_INIT_COMPARTMENTS = HiddenZoneInitCompartments()


class PhotosyntheticOrganElementInitCompartments(object):
    """
    Initial values for compartments of photosynthetic organ elements.
    """
    def __init__(self):
        # state parameters
        self.age = 0                         #: °Cd
        self.sucrose = 15                    #: :math:`\mu mol C
        self.amino_acids = 2                 #: :math:`\mu mol N
        self.proteins = 25                   #: :math:`\mu mol N
        self.temperature = 20                #: °C
        self.Tr = 0                          #: mmol H20 m-2 s-1
        self.mstruct = 0                     #: g
        self.green_area = 0                  #: m2

        # intermediate variables
        self.osmotic_water_potential = -0.8  #: MPa
        self.total_water_potential = -0.1    #: MPa

        # state variables
        self.length = 4E-5                   #: m
        self.turgor_water_potential = 0.7    #: MPa
        self.water_content = -R * (self.temperature + celsius_2_kelvin) * (
                ((self.sucrose * 1E-6 / NB_C_SUCROSE) * VANT_HOFF_SUCROSE) +
                ((self.amino_acids * 1E-6 / AMINO_ACIDS_N_RATIO) * VANT_HOFF_AMINO_ACIDS) +
                ((self.proteins * 1E-6 / AMINO_ACIDS_N_RATIO) * VANT_HOFF_AMINO_ACIDS)) / self.osmotic_water_potential  #: g


#: The instance of class :class:`turgorgrowth.parameters.LaminaInitCompartments` for current process
PHOTOSYNTHETIC_ORGAN_ELEMENT_INIT_COMPARTMENTS = PhotosyntheticOrganElementInitCompartments()


class LaminaElementParameters(OrganParameters):
    """
    Internal parameters of lamina.
    """
    def __init__(self):
        super(LaminaElementParameters, self).__init__()

        self.epsilon = {'x': 10, 'y': 10, 'z': 10}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.


#: The instance of class :class:`turgorgrowth.parameters.LaminaParameters` for current process
LAMINA_ELEMENT_PARAMETERS = LaminaElementParameters()


class LaminaElementInitCompartments(object):
    """
    Initial values for compartments of lamina elements.
    """
    def __init__(self):
        # state variables
        self.width = 0.003       #: m
        self.thickness = 0.0005  #: m


#: The instance of class :class:`turgorgrowth.parameters.LaminaElementInitCompartments` for current process
LAMINA_ELEMENT_INIT_COMPARTMENTS = LaminaElementInitCompartments()


class InternodeElementParameters(OrganParameters):
    """
    Internal parameters of internodes.
    """
    def __init__(self):
        super(InternodeElementParameters, self).__init__()

        self.epsilon = {'x': 10, 'y': 10, 'z': 10}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        self.vstorage = 0.1  #: Storage portion of the organ volume


#: The instance of class :class:`turgorgrowth.parameters.InternodeParameters` for current process
INTERNODE_ELEMENT_PARAMETERS = InternodeElementParameters()


class InternodeElementInitCompartments(object):
    """
    Initial values for compartments of internode elements.
    """
    def __init__(self):
        # state variables
        self.width = 0.003       #: m
        self.thickness = 0.0005  #: m


#: The instance of class :class:`turgorgrowth.parameters.InternodeElementInitCompartments` for current process
INTERNODE_ELEMENT_INIT_COMPARTMENTS = InternodeElementInitCompartments()


class SheathElementParameters(OrganParameters):
    """
    Internal parameters of sheaths.
    """
    def __init__(self):
        super(SheathElementParameters, self).__init__()

        self.epsilon = {'x': 10, 'y': 10, 'z': 10}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        self.vstorage = 0.1  #: Storage portion of the organ volume


#: The instance of class :class:`turgorgrowth.parameters.InternodeParameters` for current process
SHEATH_ELEMENT_PARAMETERS = SheathElementParameters()


class SheathElementInitCompartments(object):
    """
    Initial values for compartments of sheath elements.
    """
    def __init__(self):
        # state variables
        self.width = 0.003       #: m
        self.thickness = 0.0005  #: m


#: The instance of class :class:`turgorgrowth.parameters.SheathElementInitCompartments` for current process
SHEATH_ELEMENT_INIT_COMPARTMENTS = SheathElementInitCompartments()
