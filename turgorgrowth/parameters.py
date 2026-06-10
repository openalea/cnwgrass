# -*- coding: latin-1 -*-
"""
    turgorgrowth.parameters
    ~~~~~~~~~~~~~~~~~~

    The module :mod:`turgorgrowth.parameters` defines the parameters of the model.

    :license: CeCILL-C, see LICENCE for details.
"""

from math import exp

HOUR_TO_SECOND_CONVERSION_FACTOR = 3600.  #: Number of seconds in 1 hour

CELSIUS_2_KELVIN = 273.15  #: conversion factor from degree Celsius to Kelvin
R = 8.31  #: Perfect gas constant (J mol-1 K-1)
RHO_WATER = 1E6  #: Water density (g m-3)
PI = 3.141592653    #: Pi (?)

NB_C_SUCROSE = 12  #: Number of C in 1 mol of sucrose
SUCROSE_MOLAR_MASS = 342  #: g mol-1
AMINO_ACIDS_N_RATIO = 1.17  #: Mean number of mol of N in 1 mol of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly)
WATER_MOLAR_MASS = 18  #: g mol-1
VSTORAGE = 0.8  #: Storage portion of the volume of the organ (-)

VANT_HOFF_SUCROSE = 1  #: Van't Hoff coefficient of sucrose (dimensionless)
VANT_HOFF_AMINO_ACIDS = 1.25  #: Van't Hoff coefficient estimated for amino acids (dimensionless)

RATIO_MSTRUCT_DM = 0.8  #: Ratio mstruct/dry matter (dimensionless). From growthwheat model.
SLOPE_MASS_VOLUME = 3.23337E-06  #: Slope of the relation between leaf dry mass and its volume at the time of the previous leaf emergence (m3 g-1). Found from Williams 1960, Fig 11.
OFFSET_MASS_VOLUME = 1.82312E-13  #: Offset of the relation between leaf dry mass and its volume at the time of the previous leaf emergence (m3). Found from Williams 1960, Fig 11.


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
        super(PlantParameters, self).__init__()

        self.plant_density = 250


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

class AxisInitCompartments(object):
    """
    Initial values for compartments of axis.
    """
    def __init__(self):
        # state parameters
        self.SAM_temperature = 12  #: initial temperature of shoot apical meristem (°C)

#: The instance of class :class:`turgorgrowth.parameters.PhytomerInitCompartments` for current process
AXIS_INIT_COMPARTMENTS = AxisInitCompartments()

class PhytomerParameters(object):
    """
    Internal parameters of phytomers.
    """
    def __init__(self):
        pass


#: The instance of class :class:`turgorgrowth.parameters.PhytomerParameters` for current process
PHYTOMER_PARAMETERS = PhytomerParameters()

class PhytomerInitCompartments(object):
    """
    Initial values for compartments of hidden zones.
    """
    def __init__(self):
        # state parameters
        self.Tr = 0         #: mmol H20 m-2 s-1
        self.green_area = 0         #: m2

#: The instance of class :class:`turgorgrowth.parameters.PhytomerInitCompartments` for current process
PHYTOMER_INIT_COMPARTMENTS = PhytomerInitCompartments()

class OrganParameters(object):
    """
    Internal parameters of organs.
    """
    def __init__(self):
        super(OrganParameters, self).__init__()


#: The instance of class :class:`turgorgrowth.parameters.PhytomerParameters` for current process
ORGAN_PARAMETERS = OrganParameters()


class HiddenZoneParameters(OrganParameters):
    """
    Internal parameters of hidden growing zones.
    """
    def __init__(self):
        super(HiddenZoneParameters, self).__init__()

        # Elasticity
        self.epsilon = {'x': 50, 'y': 40, 'z': 50}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.

        # Extensibility
        self.phi_initial = {'x': 13E-09, 'y': 10E-09, 'z': 34E-06}  #: Initial dimension-specific cell wall extensibility in relation to non-reversible dimensional changes (MPa-1 s-1)

        # Length
        self.tend = 2160000  #: end of leaf elongation in automate growth (s at 12°c); fitted from adapted data from Fournier 2005
        self.tmax = 1473120  #: time at which leaf elongation rate is maximal in automate growth (s at 12°c); fitted from adapted data from Fournier 2005
        self.tbase = -822960  #: beginning of leaf elongation in automate growth (s at 12°c); fitted from adapted data from Fournier 2005

        # Parameters for temperature responses
        self.Temp_Tref = 12  # Arbitrary reference temperature (°C)
        self.Temp_Ea_R = 8900  # Parameter Ea/R in Eyring equation from Johnson and Lewin (1946) - Parameter value fitted from Kemp and Blacklow (1982) (K)
        self.Temp_DS_R = 68.432  # Parameter deltaS/R in Eyring equation from Johnson and Lewin (1946) - Parameter value fitted from Kemp and Blacklow (1982) (dimensionless)
        self.Temp_DH_R = 20735.5  # Parameter deltaH/R in Eyring equation from Johnson and Lewin (1946) - Parameter value fitted from Kemp and Blacklow (1982) (K)
        self.Temp_Ttransition = 9  # Below this temperature f = linear function of temperature instead of Arrhenius-like(°C)

        # Maximal leaf length
        self.leaf_Lmax_MAX = 0.45  #: Maximum leaf_Lmax (m) (Gauthier et al., 2021)

        # Maximal leaf width
        self.leaf_Wmax_Marion = {1: 0.0030, 2: 0.0033, 3: 0.0040, 4: 0.0045, 5: 0.0056, 6: 0.0075, 7: 0.010, 8: 0.012, 9: 0.013, 10: 0.014, 11: 0.018}  #: m

        L0 = abs((1 + (self.tend / (self.tend - self.tmax))) * (min(1.0, float(-self.tbase) / float(self.tend - self.tbase)) ** ((self.tend - self.tbase) / (self.tend - self.tmax))))  #: Leaf length at t=0 in automate growth (beta function) (m)
        FITTED_L0 = 0.01557936  #: Fitted value of leaf length at t=0 after rescaling the beta function with L0 (m); Fournier 2005 sur courbe corrigee
        self.OFFSET_LEAF = FITTED_L0 - L0  #: Offset used for the final fitting of the beta function (m)
        self.RATIO_MSTRUCT_DM = 0.8     #: Ratio mstruct/dry matter (dimensionless). From growthwheat model.
        self.SLOPE_MASS_VOLUME = 3.23337E-06  #: Slope of the relation between leaf dry mass and its volume at the time of the previous leaf emergence (m3 g-1). Found from Williams 1960, Fig 11.
        self.OFFSET_MASS_VOLUME = 1.82312E-13  #: Offset of the relation between leaf dry mass and its volume at the time of the previous leaf emergence (m3). Found from Williams 1960, Fig 11.
        self.GAMMA = 0.2    #: Critical value for the pressure component which must be exceeded for irreversible volume changes (MPa). Found from Coussement et al., 2018 : 0.3 Mpa for soybean.

        self.Sa = 260     #: (mol m-3)    Parameter of a sigmoidal function of equivalent solutes concentration used in osmotic water potential
        self.Sb = 0.9      #: (-)    Parameter of a sigmoidal function of equivalent solutes concentration used in osmotic water potential
        self.Sc = -8.5     #: (-)    Parameter of a sigmoidal function of equivalent solutes concentration used in osmotic water potential
        self.Sd = 1000    #: (mol m-3)  Parameter of a sigmoidal function of equivalent solutes concentration used in osmotic water potential

        self.WL_ratio = 0.675    #: -
        self.TL_ratio = 0.14    #: -

#: The instance of class :class:`turgorgrowth.parameters.HiddenZoneParameters` for current process
HIDDEN_ZONE_PARAMETERS = HiddenZoneParameters()


class HiddenZoneInitCompartments(object):
    """
    Initial values for compartments of hidden zones.
    """
    def __init__(self):
        self.Tr = 0     #: mmol H20 m-2 s-1
        self.green_area = 0     #: m2
        self.temperature = 8  #: °C
        self.hiddenzone_age = 0  #: °Cd
        self.leaf_pseudo_age = -1  #: °Cd
        self.amino_acids = 0.000075  #: :math:`\mu mol N
        self.proteins = 0.0011  #: :math:`\mu mol N
        self.sucrose = 0.000384  #: :math:`\mu mol C
        self.fructan = 0    #: :math:`\mu mol C
        self.mstruct = 1.26E-07    #: g
        self.leaf_enclosed_mstruct = 1.26E-07    #: g
        self.delta_teq = 3600   #: s     #: time equivalent to a reference temperature i.e. temperature-compensated time (Parent, 2010)

        self.leaf_pseudostem_length = 4E-5   #: m
        self.leaf_L = 5E-5                   #: m
        self.length_hz_En = None                #: m
        self.lamina_Lmax = None                 #: m
        self.leaf_Wmax = None                 #: m
        self.leaf_is_growing = True                #: -
        self.width = 0.003                 #: m
        self.thickness = 0.0005     #: m
        self.volume = self.leaf_L * self.width * self.thickness

        self.water_content = self.volume * RHO_WATER
        self.water_influx = 0                 #: g H2O
        self.water_outflow = 0                 #: g H2O
        self.SRWC = 80  #: %
        self.osmotic_water_potential = -0.8 #: Mpa
        self.water_potential = - exp((-self.SRWC + 39.765) / 18.902)  #: MPa
        self.turgor_water_potential = self.water_potential - self.osmotic_water_potential  #: MPa


#: The instance of class :class:`turgorgrowth.parameters.HiddenZoneInitCompartments` for current process
HIDDEN_ZONE_INIT_COMPARTMENTS = HiddenZoneInitCompartments()

class RootsParameters(object):
    """
    Internal parameters of roots.
    """
    def __init__(self):
        super(RootsParameters, self).__init__()

#: The instance of class :class:`cnwheat.parameters.XylemParameters` for current process
ROOTS_PARAMETERS = RootsParameters()


class RootsInitCompartments(object):
    """
    Initial values for compartments of roots
    """

#: The instance of class :class:`cnwheat.parameters.RootsInitCompartments` for current process
ROOTS_INIT_COMPARTMENTS = RootsInitCompartments()


class XylemParameters(object):
    """
    Internal parameters of xylem.
    """
    def __init__(self):
        super(XylemParameters, self).__init__()

        self.R_xylem_hz = 1     #: Flow resistance between xylem and shoot organs (Mpa s g-1 m) : 1
        self.R_soil = 1E-05    #: Flow resistance between soil and xylem (Mpa s g-1 m) : 1E-05
        self.R_xylem_organ = 0.25      #: Flow resistance between xylem and shoot organs (Mpa s g-1 m) : 0.25

#: The instance of class :class:`cnwheat.parameters.XylemParameters` for current process
XYLEM_PARAMETERS = XylemParameters()


class XylemInitCompartments(object):
    """
    Initial values for compartments of xylem.
    """
    def __init__(self):
        # state parameters

        # intermediate variables
        self.SRWC = 80  #: %
        self.soil_water_potential = - exp((-self.SRWC + 39.765) / 18.902)  #: MPa
        self.water_potential = - exp((-self.SRWC + 39.765) / 18.902)  #: MPa


#: The instance of class :class:`cnwheat.parameters.XylemInitCompartments` for current process
XYLEM_INIT_COMPARTMENTS = XylemInitCompartments()


class PhotosyntheticOrganElementParameters(object):
    """
    Internal parameters of Photosynthetic Organ Element.
    """
    def __init__(self):
        super(PhotosyntheticOrganElementParameters, self).__init__()

        self.epsilon = {'x': 150, 'y': 120, 'z': 150}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.

        self.Sa = 260     #: (mol m-3)    Parameter of a sigmoidal function of equivalent solutes concentration used in osmotic water potential
        self.Sb = 0.9      #: (-)    Parameter of a sigmoidal function of equivalent solutes concentration used in osmotic water potential
        self.Sc = -8.5     #: (-)    Parameter of a sigmoidal function of equivalent solutes concentration used in osmotic water potential
        self.Sd = 1000    #: (mol m-3)  Parameter of a sigmoidal function of equivalent solutes concentration used in osmotic water potential

#: The instance of class :class:`cnwheat.parameters.PhotosyntheticOrganElementParameters` for current process
PHOTOSYNTHETIC_ORGAN_ELEMENT_PARAMETERS = PhotosyntheticOrganElementParameters()

class PhotosyntheticOrganElementInitCompartments(object):
    """
    Initial values for compartments of photosynthetic organ elements.
    """
    def __init__(self):
        # state parameters
        self.is_growing = None                       #: -
        self.age = None                       #: °Cd
        self.Wmax = None                       #: °Cd
        self.amino_acids = 0                 #: :math:`\mu mol N
        self.proteins = 0                   #: :math:`\mu mol N
        self.sucrose = 0                    #: :math:`\mu mol C
        self.fructan = 0    #: :math:`\mu mol C
        self.temperature = 0                #: °C
        self.Ts = 12                         #: °C
        self.green_area = 1E-4               #: initial value of green_area (m2)
        self.Tr = 0                          #: mmol H20 m-2 s-1
        self.mstruct = 0                     #: g
        self.width = 0.003                   #: m init
        self.thickness = 0.0005              #: m init
        self.length = 4E-5                   #: m init

        self.volume = self.length * self.width * self.thickness
        self.water_content = self.volume * RHO_WATER
        self.water_influx = 0   #: g H2O
        self.water_outflow = 0   #: g H2O
        self.SRWC = 80  #: %
        self.osmotic_water_potential = -0.8   #: MPa
        self.water_potential = - exp((-self.SRWC + 39.765) / 18.902)  #: MPa
        self.turgor_water_potential = self.water_potential - self.osmotic_water_potential   #: MPa



#: The instance of class :class:`turgorgrowth.parameters.LaminaInitCompartments` for current process
PHOTOSYNTHETIC_ORGAN_ELEMENT_INIT_COMPARTMENTS = PhotosyntheticOrganElementInitCompartments()


class LaminaElementParameters(OrganParameters):
    """
    Internal parameters of lamina.
    """
    def __init__(self):
        super(LaminaElementParameters, self).__init__()

        self.epsilon = {'x': 150, 'y': 120, 'z': 150}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.

#: The instance of class :class:`turgorgrowth.parameters.LaminaParameters` for current process
LAMINA_ELEMENT_PARAMETERS = LaminaElementParameters()


class LaminaElementInitCompartments(object):
    """
    Initial values for compartments of lamina elements.
    """
    def __init__(self):
        # state parameters
        self.width = 0.003                   #: m
        self.thickness = 0.0005              #: m

#: The instance of class :class:`turgorgrowth.parameters.LaminaElementInitCompartments` for current process
LAMINA_ELEMENT_INIT_COMPARTMENTS = LaminaElementInitCompartments()


class InternodeElementParameters(OrganParameters):
    """
    Internal parameters of internodes.
    """
    def __init__(self):
        super(InternodeElementParameters, self).__init__()

        self.epsilon = {'x': 10, 'y': 10, 'z': 10}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.

#: The instance of class :class:`turgorgrowth.parameters.InternodeParameters` for current process
INTERNODE_ELEMENT_PARAMETERS = InternodeElementParameters()


class InternodeElementInitCompartments(object):
    """
    Initial values for compartments of internode elements.
    """
    def __init__(self):
        # state parameters
        self.width = 0.003                   #: m
        self.thickness = 0.0005              #: m

#: The instance of class :class:`turgorgrowth.parameters.InternodeElementInitCompartments` for current process
INTERNODE_ELEMENT_INIT_COMPARTMENTS = InternodeElementInitCompartments()


class SheathElementParameters(OrganParameters):
    """
    Internal parameters of sheaths.
    """
    def __init__(self):
        super(SheathElementParameters, self).__init__()

        self.epsilon = {'x': 150, 'y': 120, 'z': 150}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.


#: The instance of class :class:`turgorgrowth.parameters.InternodeParameters` for current process
SHEATH_ELEMENT_PARAMETERS = SheathElementParameters()


class SheathElementInitCompartments(object):
    """
    Initial values for compartments of sheath elements.
    """
    def __init__(self):
        # state parameters
        self.width = 0.003                   #: m
        self.thickness = 0.0005              #: m

#: The instance of class :class:`turgorgrowth.parameters.SheathElementInitCompartments` for current process
SHEATH_ELEMENT_INIT_COMPARTMENTS = SheathElementInitCompartments()


class SoilParameters(object):
    """
    Internal parameters of soil.
    """
    def __init__(self):
        self.AWC = 50  # Available Water Capacity (g)
        self.Soil_a = 15.906  #: Mpa - Parameter for soil water function (adapté pour sol limono-argileux profond, Grignon)
        self.Soil_b = 18.902  #: % - Parameter for soil water function (adapté pour sol limono-argileux profond, Grignon)

#: The instance of class :class:`cnwheat.parameters.SoilParameters` for current process
SOIL_PARAMETERS = SoilParameters()
