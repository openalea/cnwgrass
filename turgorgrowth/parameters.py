# -*- coding: latin-1 -*-

"""
    turgorgrowth.parameters
    ~~~~~~~~~~~~~~~~~~

    The module :mod:`turgorgrowth.parameters` defines the parameters of the model.

    :license: CeCILL-C, see LICENSE for details.
"""

HOUR_TO_SECOND_CONVERSION_FACTOR = 3600.  #: Number of seconds in 1 hour

celsius_2_kelvin = 273.15  #: conversion factor from degree Celsius to Kelvin
R = 8.31  #: Perfect gas constant (J mol-1 K-1)
rho_water = 1E6  #: Water density (kg m-3)

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


#: The instance of class :class:`cnwheat.parameters.PopulationParameters` for current process
POPULATION_PARAMETERS = PopulationParameters()


class PlantParameters(object):
    """
    Internal parameters of plants.
    """
    def __init__(self):
        pass


#: The instance of class :class:`cnwheat.parameters.PlantParameters` for current process
PLANT_PARAMETERS = PlantParameters()


class AxisParameters(object):
    """
    Internal parameters of axes.
    """
    def __init__(self):
        pass


#: The instance of class :class:`cnwheat.parameters.AxisParameters` for current process
AXIS_PARAMETERS = AxisParameters()


class PhytomerParameters(object):
    """
    Internal parameters of phytomers.
    """
    def __init__(self):
        pass


#: The instance of class :class:`cnwheat.parameters.PhytomerParameters` for current process
PHYTOMER_PARAMETERS = PhytomerParameters()


class OrganParameters(object):
    """
    Internal parameters of organs.
    """
    def __init__(self):
        self.R_xylem = 1E-5 * HOUR_TO_SECOND_CONVERSION_FACTOR  #: Elemental xylem flow resistence (MPa s g-1 m, from Coussement et al., 2018)
        self.vstorage = 0.8                                     #: Storage portion of the organ volume (0.2 attributed to the xylem)
        self.gamma = 0.3                                        #: Critical value for the pressure component which must be exceeded for irreversible volume changes (MPa)


#: The instance of class :class:`cnwheat.parameters.PhytomerParameters` for current process
ORGAN_PARAMETERS = OrganParameters()


class HiddenZoneParameters(OrganParameters):
    """
    Internal parameters of hidden growing zones.
    """
    def __init__(self):
        super(HiddenZoneParameters, self).__init__()

        self.epsilon = {'x': 0.2, 'z': 0.2}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        self.phi_initial = {'x': 0.25 / HOUR_TO_SECOND_CONVERSION_FACTOR,  #: Initial dimension-specific cell wall extensibility in relation to non-reversible dimensional changes (MPa-1 s-1)
                            'z': 0.45 / HOUR_TO_SECOND_CONVERSION_FACTOR}
        self.tend = 200.  #: Hidden zone age when extensibility reaches 0 (°Cd)
        self.tmax = 100.  #: Hidden zone age when organ extensibility is reduced by half of the initial value (°Cd)


#: The instance of class :class:`turgorgrowth.parameters.HiddenZoneParameters` for current process
HIDDEN_ZONE_PARAMETERS = HiddenZoneParameters()


class HiddenZoneInitCompartments(object):
    """
    Initial values for compartments of hidden zones.
    """
    def __init__(self):
        # state parameters
        self.age = 0  #: °Cd
        self.amino_acids = 0.000075  #: :math:`\mu mol N
        self.proteins = 0.0011  #: :math:`\mu mol N
        self.sucrose = 0.000384  #: :math:`\mu mol C
        self.temperature = 20  #: °C

        # state variables
        self.length = 4E-5                   #: m
        self.radius = 4E-5                   #: m
        self.turgor_water_potential = 0.7    #: MPa
        self.water_content = 0               #: g


#: The instance of class :class:`turgorgrowth.parameters.HiddenZoneInitCompartments` for current process
HIDDEN_ZONE_INIT_COMPARTMENTS = HiddenZoneInitCompartments()


class PhotosyntheticOrganElementInitCompartments(object):
    """
    Initial values for compartments of photosynthetic organ elements.
    """
    def __init__(self):
        # state parameters
        self.age = 0  #: °Cd
        self.sucrose = 15  #: :math:`\mu mol C
        self.amino_acids = 2  #: :math:`\mu mol N
        self.proteins = 25  #: :math:`\mu mol N
        self.temperature = 20  #: °C
        self.transpiration = 0  #: mmol H20

        # state variables
        self.length = 4E-5                   #: m
        self.osmotic_water_potential = -0.8  #: MPa
        self.thickness = 0.25E-3             #: m
        self.total_water_potential = -0.1    #: MPa
        self.turgor_water_potential = 0.7    #: MPa
        self.water_content = 0               #: g
        self.width = 2E-3                    #: m


#: The instance of class :class:`turgorgrowth.parameters.LaminaInitCompartments` for current process
PHOTOSYNTHETIC_ORGAN_ELEMENT_INIT_COMPARTMENTS = PhotosyntheticOrganElementInitCompartments()


class LaminaParameters(OrganParameters):
    """
    Internal parameters of lamina.
    """
    def __init__(self):
        super(LaminaParameters, self).__init__()

        self.epsilon = {'x': 0.5, 'y': 0.5, 'z': 1.}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        self.phi_initial = {'x': 0.40 / HOUR_TO_SECOND_CONVERSION_FACTOR,  #: Initial dimension-specific cell wall extensibility in relation to non-reversible dimensional changes (MPa-1 s-1)
                            'y': 0.01 / HOUR_TO_SECOND_CONVERSION_FACTOR,
                            'z': 0.45 / HOUR_TO_SECOND_CONVERSION_FACTOR}
        self.tend = 200  #: Lamina age when extensibility reaches 0 (°Cd)
        self.tmax = 100  #: Lamina age when organ extensibility is reduced by half of the initial value (°Cd)


#: The instance of class :class:`turgorgrowth.parameters.LaminaParameters` for current process
LAMINA_PARAMETERS = LaminaParameters()


class InternodeParameters(OrganParameters):
    """
    Internal parameters of internodes.
    """
    def __init__(self):
        super(InternodeParameters, self).__init__()

        self.epsilon = {'x': 0.2, 'z': 0.2}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: radius, z: length.
        self.phi_initial = {'x': 0.25 / HOUR_TO_SECOND_CONVERSION_FACTOR,  #: Initial dimension-specific cell wall extensibility in relation to non-reversible dimensional changes (MPa-1 s-1)
                            'z': 0.45 / HOUR_TO_SECOND_CONVERSION_FACTOR}
        self.tend = 200  #: Hidden zone age when extensibility reaches 0 (°Cd)
        self.tmax = 100  #: Hidden zone age when organ extensibility is reduced by half of the initial value (°Cd)

        self.vstorage = 0.1  #: Storage portion of the organ volume


#: The instance of class :class:`turgorgrowth.parameters.InternodeParameters` for current process
INTERNODE_PARAMETERS = InternodeParameters()


class SheathParameters(OrganParameters):
    """
    Internal parameters of sheaths.
    """
    def __init__(self):
        super(SheathParameters, self).__init__()

        self.epsilon = {'x': 1, 'z': 0.2}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: radius, z: length.
        self.phi_initial = {'x': 0.25 / HOUR_TO_SECOND_CONVERSION_FACTOR,  #: Initial dimension-specific cell wall extensibility in relation to non-reversible dimensional changes (MPa-1 s-1)
                            'z': 0.45 / HOUR_TO_SECOND_CONVERSION_FACTOR}
        self.tend = 200  #: Hidden zone age when extensibility reaches 0 (°Cd)
        self.tmax = 100  #: Hidden zone age when organ extensibility is reduced by half of the initial value (°Cd)

        self.vstorage = 0.05  #: Storage portion of the organ volume


#: The instance of class :class:`turgorgrowth.parameters.InternodeParameters` for current process
SHEATH_PARAMETERS = SheathParameters()
