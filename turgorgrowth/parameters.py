# -*- coding: latin-1 -*-
"""
    turgorgrowth.parameters
    ~~~~~~~~~~~~~~~~~~

    The module :mod:`turgorgrowth.parameters` defines the parameters of the model.

    :license: CeCILL-C, see LICENSE for details.
"""

from math import exp

HOUR_TO_SECOND_CONVERSION_FACTOR = 3600.  #: Number of seconds in 1 hour

CELSIUS_2_KELVIN = 273.15  #: conversion factor from degree Celsius to Kelvin
R = 8.31  #: Perfect gas constant (J mol-1 K-1)
RHO_WATER = 1E6  #: Water density (g m-3)

NB_C_SUCROSE = 12  #: Number of C in 1 mol of sucrose
SUCROSE_MOLAR_MASS = 342  #: g mol-1
AMINO_ACIDS_N_RATIO = 1.17  #: Mean number of mol of N in 1 mol of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly)
WATER_MOLAR_MASS = 18  #: g mol-1

VANT_HOFF_SUCROSE = 1  #: Van't Hoff coefficient of sucrose (dimensionless)
VANT_HOFF_AMINO_ACIDS = 1.25  #: Van't Hoff coefficient estimated for amino acids (dimensionless)

# TEST
RATIO_MSTRUCT_DM = 0.8  #: Ratio mstruct/dry matter (dimensionless). From growthwheat model.
SLOPE_MASS_VOLUME = 3.23337E-06  #: Slope of the relation between leaf dry mass and its volume at the time of the previous leaf emergence (m3 g-1). Found from Williams 1960, Fig 11.
OFFSET_MASS_VOLUME = 1.82312E-13  #: Offset of the relation between leaf dry mass and its volume at the time of the previous leaf emergence (m3). Found from Williams 1960, Fig 11.

def from_dataframe(object_, dataframe_):
    """Set attributes of *object_* from data in *dataframe_*.

    :Parameters:
        - `object_` (:class:`object`) - The object to set.
        - `dataframe_` (:class:`pandas.DataFrame`) - The dataframe used to set the attribute(s)
          of *object_*.
          *dataframe_* must have only 2 rows:

              * one row is for the header and contains the name of each attribute,
              * and one row contains the value of each attribute.
    """
    object_.__dict__.update(dataframe_.to_dict(orient='index')[dataframe_.first_valid_index()])


def to_dataframe(object_):
    """Create and return a dataframe from attributes of *object_*.

    :Parameters:
        - `object_` (:class:`object`) - The object used to create the dataframe.

    :Returns:
        A dataframe which contains the attributes of *object_*, with only 2 rows:

          * one row is for the header and contains the name of each attribute,
          * and one row contains the value of each attribute.

    :Returns Type:
        :class:`pandas.DataFrame`
    """
    return pd.DataFrame(object_.__dict__, index=[0]).sort_index(axis=1)


class PopulationParameters(object):
    """
    Internal parameters of populations.
    """
    def __init__(self):
        pass


#: The instance of class :class:`turgorgrowth.parameters.PopulationParameters` for current process
POPULATION_PARAMETERS = PopulationParameters()

class SoilParameters(object):
    """
    Internal parameters of soil.
    """
    def __init__(self):
        pass


#: The instance of class :class:`turgorgrowth.parameters.SoilParameters` for current process
SOIL_PARAMETERS = SoilParameters()

class SoilInitCompartments(object):
    """
    Initial values for compartments of soil.
    """
    def __init__(self):
        # state parameters
        self.SRWC = 80  #: %
        self.soil_water_potential = - exp((-self.SRWC + 39.765) / 18.902)  #: MPa

#: The instance of class :class:`turgorgrowth.parameters.SoilInitCompartments` for current process
SOIL_INIT_COMPARTMENTS = SoilInitCompartments()


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

class PhytomerInitCompartments(object):
    """
    Initial values for compartments of hidden zones.
    """
    def __init__(self):
        # state parameters
        self.Tr = 0
        self.green_area = 0

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

        #INITAL
        # self.epsilon = {'x': 50, 'y': 50, 'z': 50}  #: 0.9 Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        # self.epsilon = {'x': 3, 'y': 2, 'z': 3}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        # self.phi_initial = {'x': 1E-09, 'y': 1E-09, 'z': 0.4E-05}  #: Initial dimension-specific cell wall extensibility in relation to non-reversible dimensional changes (MPa-1 s-1)

        # volume work package
        # self.phi_initial = {'x': 1E-03, 'y': 1E-03, 'z': 4E-02}  #: Initial dimension-specific cell wall extensibility in relation to non-reversible dimensional changes (MPa-1 s-1)
        # self.phi_initial = {'x': 1E-07, 'y': 1E-07, 'z': 1E-05}  #: Initial dimension-specific cell wall extensibility in relation to non-reversible dimensional changes (MPa-1 s-1)
        # self.epsilon = {'x': 1, 'y': 1, 'z': 1}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        self.epsilon = {'x': 1.5, 'y': 1.5, 'z': 3.0}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.

        # phi work package
        self.phi_initial = {'x': 3.5E-05, 'y': 2.5E-05, 'z': 2E-05}  #: Initial dimension-specific cell wall extensibility in relation to non-reversible dimensional changes (MPa-1 s-1)

        self.tend = 2500000  # 300 * 3600 * 24 / 12   #: Lamina age when extensibility reaches 0 (s at 12°C). Calculated from elongwheat parameter for phase 2
        self.tmax = 2000000  # 190 * 3600 * 24 / 12   #: Lamina age when organ extensibility is reduced by half of the initial value (s at 12°C). Calculated from elongwheat parameter for phase 2
        self.tbase = -25 * 3600 * 24 / 12 #: - 180 000 #: beginning of leaf elongation in automate growth (s at 12°C); fitted from adapted data from Fournier 2005
        L0 = abs((1 + (self.tend / (self.tend - self.tmax))) * (min(1.0, float(-self.tbase) / float(self.tend - self.tbase)) ** ((self.tend - self.tbase) / (self.tend - self.tmax))))  #: Leaf length at t=0 in automate growth (beta function) (m)
        FITTED_L0 = 0.01557936  #: Fitted value of leaf length at t=0 after rescaling the beta function with L0 (m); Fournier 2005 sur courbe corrigee
        self.OFFSET_LEAF = FITTED_L0 - L0  #: Offset used for the final fitting of the beta function (m)
        self.RATIO_MSTRUCT_DM = 0.8  #: Ratio mstruct/dry matter (dimensionless). From growthwheat model.
        self.SLOPE_MASS_VOLUME = 3.23337E-06  #: Slope of the relation between leaf dry mass and its volume at the time of the previous leaf emergence (m3 g-1). Found from Williams 1960, Fig 11.
        self.OFFSET_MASS_VOLUME = 1.82312E-13  #: Offset of the relation between leaf dry mass and its volume at the time of the previous leaf emergence (m3). Found from Williams 1960, Fig 11.
        self.GAMMA = 0.3


#: The instance of class :class:`turgorgrowth.parameters.HiddenZoneParameters` for current process
HIDDEN_ZONE_PARAMETERS = HiddenZoneParameters()


class HiddenZoneInitCompartments(object):
    """
    Initial values for compartments of hidden zones.
    """
    def __init__(self):
        self.Tr = 0
        self.green_area = 0
        self.temperature = 8  #: °C
        self.age = 0  #: °Cd
        self.amino_acids = 0.000075  #: :math:`\mu mol N
        self.proteins = 0.0011  #: :math:`\mu mol N
        self.sucrose = 0.000384  #: :math:`\mu mol C
        # self.mstruct = 2.65E-08    #: g TODO: viens d'ou?
        self.mstruct = 2.65E-10    #: g TODO: viens d'ou?

        self. volume = self.mstruct / RATIO_MSTRUCT_DM * SLOPE_MASS_VOLUME + OFFSET_MASS_VOLUME  #: m3
        self.water_content = ((self.mstruct / RATIO_MSTRUCT_DM) * SLOPE_MASS_VOLUME + OFFSET_MASS_VOLUME) * RHO_WATER  #: g
        # WC = ((mstruct / 0.8) * 3.23337E-06 + 1.82312E-13 ) * 1E06

        # self.osmotic_water_potential = - R * (self.temperature + CELSIUS_2_KELVIN) * (
        #          ((self.sucrose * 1E-6 / NB_C_SUCROSE) * VANT_HOFF_SUCROSE) +
        #          ((self.amino_acids * 1E-6 / AMINO_ACIDS_N_RATIO) * VANT_HOFF_AMINO_ACIDS) +
        #          ((self.proteins * 1E-6 / AMINO_ACIDS_N_RATIO) * VANT_HOFF_AMINO_ACIDS) ) / (self.volume * RHO_WATER)

        self.osmotic_water_potential = -0.8

        # self.osmotic_water_potential = - R * (self.temperature + CELSIUS_2_KELVIN) * (
        #          ((self.sucrose * 1E-6 / NB_C_SUCROSE) * VANT_HOFF_SUCROSE) +
        #          ((self.amino_acids * 1E-6 / AMINO_ACIDS_N_RATIO) * VANT_HOFF_AMINO_ACIDS) ) / (self.volume * RHO_WATER)

        self.SRWC = 80  #: %
        self.total_water_potential = - exp((-self.SRWC + 39.765) / 18.902)  #: MPa
        self.turgor_water_potential = self.total_water_potential - self.osmotic_water_potential   #: MPa

        self.leaf_pseudostem_length = 4E-6   #: m
        self.leaf_L = 4E-6                   #: m
        self.leaf_Lmax = 0.1331                    #: m
        self.lamina_Lmax = 0.0991                  #: m
        self.leaf_is_growing = True                #:
        self.width = 4E-6

        # self.thickness = 4E-6
        self.thickness = self.volume / (self.width * self.leaf_L)                     #: m
        # TO DO: vérifier que length = leaf_L car pas d'émergence

        self.water_influx = 0
        self.water_outflow = 0

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
    # def __init__(self):

#: The instance of class :class:`cnwheat.parameters.RootsInitCompartments` for current process
ROOTS_INIT_COMPARTMENTS = RootsInitCompartments()


class XylemParameters(object):
    """
    Internal parameters of xylem.
    """
    def __init__(self):
        super(XylemParameters, self).__init__()

        # self.R_xylem_blade = 2    #: Flow resistance between xylem and shoot organs (Mpa s g-1 m) # change after Tom's discussion 10/2023
        # self.R_xylem_hz = 0.2    #: Flow resistance between xylem and shoot organs (Mpa s g-1 m) # change after Tom's discussion 10/2023
        self.R_soil = 0.001    #: Flow resistance between soil and xylem (Mpa s g-1 m)

        # volume
        self.R_xylem_blade = 0.5  #: Flow resistance between xylem and shoot organs (Mpa s g-1 m)
        self.R_xylem_hz = 0.2  #: Flow resistance between xylem and shoot organs (Mpa s g-1 m)

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
        self.total_water_potential = - exp((-self.SRWC + 39.765) / 18.902)  #: MPa

#: The instance of class :class:`cnwheat.parameters.XylemInitCompartments` for current process
XYLEM_INIT_COMPARTMENTS = XylemInitCompartments()


class PhotosyntheticOrganElementParameters(object):
    """
    Internal parameters of Photosynthetic Organ Element.
    """
    def __init__(self):
        super(PhotosyntheticOrganElementParameters, self).__init__()

        # self.epsilon = {'x': 50, 'y': 50, 'z': 50}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        # self.epsilon = {'x': 4, 'y': 4, 'z': 0.7}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        # self.epsilon = {'x': 0.4, 'y': 0.4, 'z': 0.7}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        # self.epsilon = {'x': 15, 'y': 15, 'z': 30}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        self.epsilon = {'x': 1.5, 'y': 1.5, 'z': 3.0}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        # self.epsilon = {'x': 0.15, 'y': 0.15, 'z': 0.30}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        # self.epsilon = {'x': 1.5, 'y': 1.5, 'z': 30}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.

        # https://acsess.onlinelibrary.wiley.com/doi/epdf/10.2134/agronj1979.00021962007100010008x
        # entre 10 et 30 MPa pour wheat et 18 wheatgrass

#: The instance of class :class:`cnwheat.parameters.PhotosyntheticOrganElementParameters` for current process
PHOTOSYNTHETIC_ORGAN_ELEMENT_PARAMETERS = PhotosyntheticOrganElementParameters()

class PhotosyntheticOrganElementInitCompartments(object):
    """
    Initial values for compartments of photosynthetic organ elements.
    """
    def __init__(self):
        # state parameters
        self.age = None                       #: °Cd
        #: TO DO : self.age = 0
        self.sucrose = 15                    #: :math:`\mu mol C
        self.amino_acids = 2                 #: :math:`\mu mol N
        self.proteins = 25                   #: :math:`\mu mol N
        self.temperature = 8                #: °C
        self.Ts = 12                         #: °C
        self.green_area = 1E-4               #: initial value of green_area (m2)
        self.Tr = 0                          #: mmol H20 m-2 s-1
        self.mstruct = 0                     #: g
        # self.width = 0.003                   #: m init
        # self.thickness = 0.0005              #: m init
        # self.length = 4E-5                   #: m init
        self.width = 0.001                   #: m
        self.thickness = 0.0001              #: m
        self.length = 0.0001                   #: m
        # self.sucrose = 0.00974303  #: :math:`\mu mol C
        # self.amino_acids = 0.00009251  #: :math:`\mu mol N
        # self.proteins = 0.01017986  #: :math:`\mu mol N
        self.volume = self.width * self.thickness * self.length
        self.water_content = self.volume * RHO_WATER

        # self.osmotic_water_potential = - R * (self.Ts + CELSIUS_2_KELVIN) * (
        #         (((self.sucrose * 1E-6) / NB_C_SUCROSE) * VANT_HOFF_SUCROSE) +
        #         (((self.amino_acids * 1E-6) / AMINO_ACIDS_N_RATIO) * VANT_HOFF_AMINO_ACIDS) +
        #         (((self.proteins * 1E-6) / AMINO_ACIDS_N_RATIO) * VANT_HOFF_AMINO_ACIDS) ) / (self.volume * RHO_WATER)

        # self.osmotic_water_potential = - R * (self.Ts + CELSIUS_2_KELVIN) * (
        #         (((self.sucrose * 1E-6) / NB_C_SUCROSE) * VANT_HOFF_SUCROSE) +
        #         (((self.amino_acids * 1E-6) / AMINO_ACIDS_N_RATIO) * VANT_HOFF_AMINO_ACIDS) ) / (self.volume * RHO_WATER)
        self.osmotic_water_potential = -0.8

        self.SRWC = 80  #: %
        self.total_water_potential = - exp((-self.SRWC + 39.765) / 18.902)  #: MPa
        self.turgor_water_potential = self.total_water_potential - self.osmotic_water_potential   #: MPa

        self.water_influx = 0
        self.water_outflow = 0


#: The instance of class :class:`turgorgrowth.parameters.LaminaInitCompartments` for current process
PHOTOSYNTHETIC_ORGAN_ELEMENT_INIT_COMPARTMENTS = PhotosyntheticOrganElementInitCompartments()


class LaminaElementParameters(OrganParameters):
    """
    Internal parameters of lamina.
    """
    def __init__(self):
        super(LaminaElementParameters, self).__init__()

        # self.epsilon = {'x': 50, 'y': 50, 'z': 50}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        # self.epsilon = {'x': 4, 'y': 4, 'z': 0.7}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        # self.epsilon = {'x': 0.4, 'y': 0.4, 'z': 0.7}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        # self.epsilon = {'x': 15, 'y': 15, 'z': 30}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        self.epsilon = {'x': 1.5, 'y': 1.5, 'z': 3.0}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        # self.epsilon = {'x': 0.15, 'y': 0.15, 'z': 0.30}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        # self.epsilon = {'x': 1.5, 'y': 1.5, 'z': 30}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.

#: The instance of class :class:`turgorgrowth.parameters.LaminaParameters` for current process
LAMINA_ELEMENT_PARAMETERS = LaminaElementParameters()


class LaminaElementInitCompartments(object):
    """
    Initial values for compartments of lamina elements.
    """
    def __init__(self):
        # state parameters
        # self.width = 0.003       #: m
        # self.thickness = 0.0005  #: m
        self.width = 0.001                   #: m init
        self.thickness = 0.0001              #: m init

#: The instance of class :class:`turgorgrowth.parameters.LaminaElementInitCompartments` for current process
LAMINA_ELEMENT_INIT_COMPARTMENTS = LaminaElementInitCompartments()


class InternodeElementParameters(OrganParameters):
    """
    Internal parameters of internodes.
    """
    def __init__(self):
        super(InternodeElementParameters, self).__init__()

        # self.epsilon = {'x': 50, 'y': 50, 'z': 50}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        # self.epsilon = {'x': 4, 'y': 4, 'z': 0.7}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        # self.epsilon = {'x': 0.4, 'y': 0.4, 'z': 0.7}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        # self.epsilon = {'x': 15, 'y': 15, 'z': 30}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        self.epsilon = {'x': 1.5, 'y': 1.5, 'z': 3.0}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        # self.epsilon = {'x': 0.15, 'y': 0.15, 'z': 0.30}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        # self.epsilon = {'x': 1.5, 'y': 1.5, 'z': 30}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.

#: The instance of class :class:`turgorgrowth.parameters.InternodeParameters` for current process
INTERNODE_ELEMENT_PARAMETERS = InternodeElementParameters()


class InternodeElementInitCompartments(object):
    """
    Initial values for compartments of internode elements.
    """
    def __init__(self):
        # state parameters
        # self.width = 0.003       #: m
        # self.thickness = 0.0005  #: m
        self.width = 0.001                   #: m init
        self.thickness = 0.0001              #: m init

#: The instance of class :class:`turgorgrowth.parameters.InternodeElementInitCompartments` for current process
INTERNODE_ELEMENT_INIT_COMPARTMENTS = InternodeElementInitCompartments()


class SheathElementParameters(OrganParameters):
    """
    Internal parameters of sheaths.
    """
    def __init__(self):
        super(SheathElementParameters, self).__init__()

        # self.epsilon = {'x': 50, 'y': 50, 'z': 50}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        # self.epsilon = {'x': 4, 'y': 4, 'z': 0.7}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        # self.epsilon = {'x': 0.4, 'y': 0.4, 'z': 0.7}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        # self.epsilon = {'x': 15, 'y': 15, 'z': 30}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        self.epsilon = {'x': 1.5, 'y': 1.5, 'z': 3.0}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        # self.epsilon = {'x': 0.15, 'y': 0.15, 'z': 0.30}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        # self.epsilon = {'x': 1.5, 'y': 1.5, 'z': 30}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.

#: The instance of class :class:`turgorgrowth.parameters.InternodeParameters` for current process
SHEATH_ELEMENT_PARAMETERS = SheathElementParameters()


class SheathElementInitCompartments(object):
    """
    Initial values for compartments of sheath elements.
    """
    def __init__(self):
        # state parameters
        # self.width = 0.003       #: m
        # self.thickness = 0.0005  #: m
        self.width = 0.001                   #: m init
        self.thickness = 0.0001              #: m init

#: The instance of class :class:`turgorgrowth.parameters.SheathElementInitCompartments` for current process
SHEATH_ELEMENT_INIT_COMPARTMENTS = SheathElementInitCompartments()
