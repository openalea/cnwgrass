# -*- coding: latin-1 -*-
"""
    turgorgrowth.parameters
    ~~~~~~~~~~~~~~~~~~

    The module :mod:`turgorgrowth.parameters` defines the parameters of the model.

    :license: CeCILL-C, see LICENSE for details.
"""

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
        self.SRWC = 90 #: %
        self.soil_water_potential = -0.1 #: Mpa

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

        self.VSTORAGE = 0.8  #: Storage portion of the organ volume (0.2 for xylem)



#: The instance of class :class:`turgorgrowth.parameters.PhytomerParameters` for current process
ORGAN_PARAMETERS = OrganParameters()


class HiddenZoneParameters(OrganParameters):
    """
    Internal parameters of hidden growing zones.
    """
    def __init__(self):
        super(HiddenZoneParameters, self).__init__()

        #INITAL
        # self.epsilon = {'x': 5, 'y': 5, 'z': 5}  #: 0.9 Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        #V1
        self.epsilon = {'x': 3, 'y': 2, 'z': 3}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.

        #INITIAL
        # self.phi_initial = {'x': 1E-09, 'y': 1E-09, 'z': 0.4E-05}  #: Initial dimension-specific cell wall extensibility in relation to non-reversible dimensional changes (MPa-1 s-1)
        #V1
        # self.phi_initial = {'x': 0.4E-05, 'y': 0.01E-05, 'z': 0.45E-05}  #: Initial dimension-specific cell wall extensibility in relation to non-reversible dimensional changes (MPa-1 s-1)
        #V2
        # self.phi_initial = {'x': 0.4E-04, 'y': 0.01E-04, 'z': 0.45E-04}  #: Initial dimension-specific cell wall extensibility in relation to non-reversible dimensional changes (MPa-1 s-1)
        #V3
        # self.phi_initial = {'x': 0.4E-06, 'y': 0.1E-06, 'z': 0.1E-04}  #: Initial dimension-specific cell wall extensibility in relation to non-reversible dimensional changes (MPa-1 s-1)
        #V4
        # self.phi_initial = {'x': 0.4E-06, 'y': 0.1E-06, 'z': 0.2E-04}  #: Initial dimension-specific cell wall extensibility in relation to non-reversible dimensional changes (MPa-1 s-1)

        self.phi_initial = {'x': 1E-09, 'y': 1E-09, 'z': 1E-05}  #: Initial dimension-specific cell wall extensibility in relation to non-reversible dimensional changes (MPa-1 s-1)

        self.tend = 2500000  # 300 * 3600 * 24 / 12   #: Lamina age when extensibility reaches 0 (s at 12°C). Calculated from elongwheat parameter for phase 2
        self.tmax = 2000000  # 190 * 3600 * 24 / 12   #: Lamina age when organ extensibility is reduced by half of the initial value (s at 12°C). Calculated from elongwheat parameter for phase 2
        self.tbase = -25 * 3600 * 24 / 12  #: beginning of leaf elongation in automate growth (s at 12°C); fitted from adapted data from Fournier 2005
        L0 = abs((1 + (self.tend / (self.tend - self.tmax))) * (min(1.0, float(-self.tbase) / float(self.tend - self.tbase)) ** ((self.tend - self.tbase) / (self.tend - self.tmax))))  #: Leaf length at t=0 in automate growth (beta function) (m)
        FITTED_L0 = 0.01557936  #: Fitted value of leaf length at t=0 after rescaling the beta function with L0 (m); Fournier 2005 sur courbe corrigee
        self.OFFSET_LEAF = FITTED_L0 - L0  #: Offset used for the final fitting of the beta function (m)
        self.RATIO_MSTRUCT_DM = 0.8  #: Ratio mstruct/dry matter (dimensionless). From growthwheat model.
        self.SLOPE_MASS_VOLUME = 3.23337E-06  #: Slope of the relation between leaf dry mass and its volume at the time of the previous leaf emergence (m3 g-1). Found from Williams 1960, Fig 11.
        self.OFFSET_MASS_VOLUME = 1.82312E-13  #: Offset of the relation between leaf dry mass and its volume at the time of the previous leaf emergence (m3). Found from Williams 1960, Fig 11.
        self.GAMMA = 0.3
        self.VSTORAGE = 0.8  #: Storage portion of the organ volume (0.2 for xylem)


#: The instance of class :class:`turgorgrowth.parameters.HiddenZoneParameters` for current process
HIDDEN_ZONE_PARAMETERS = HiddenZoneParameters()


class HiddenZoneInitCompartments(object):
    """
    Initial values for compartments of hidden zones.
    """
    def __init__(self):
        # state parameters
        self.age = 200  #: °Cd
        self.amino_acids = 0.000075  #: :math:`\mu mol N
        self.proteins = 0.0011  #: :math:`\mu mol N
        self.sucrose = 0.000384  #: :math:`\mu mol C
        self.temperature = 20  #: °C
        self.mstruct = 2.65E-08    #: g TODO: viens d'ou?
        self.Tr = 0
        self.green_area = 0

        # intermediate variables
        self.osmotic_water_potential = -0.8  #: MPa
        self.total_water_potential = -0.1    #: MPa
        self.leaf_pseudostem_length = 4E-6   #: m
        self.leaf_L = 4E-6                   #: m
        self.leaf_Lmax = 0.1331              #: m
        self.lamina_Lmax = 0.0991            #: m
        self.leaf_is_growing = True          #:
        self.thickness = 1E-3                #: m
        self.width = 0.01                    #: m
        self.turgor_water_potential = 0.7    #: MPa

        # state variables
        self.water_content = -R * (self.temperature + CELSIUS_2_KELVIN) * (
                 ((self.sucrose * 1E-6 / NB_C_SUCROSE) * VANT_HOFF_SUCROSE) +
                 ((self.amino_acids * 1E-6 / AMINO_ACIDS_N_RATIO) * VANT_HOFF_AMINO_ACIDS) +
                 ((self.proteins * 1E-6 / AMINO_ACIDS_N_RATIO) * VANT_HOFF_AMINO_ACIDS)) / self.osmotic_water_potential  #: g

        # fluxes with ylem
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

        self.epsilon = {'x': 10, 'y': 10, 'z': 10}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.

#: The instance of class :class:`cnwheat.parameters.XylemParameters` for current process
ROOTS_PARAMETERS = RootsParameters()


class RootsInitCompartments(object):
    """
    Initial values for compartments of roots
    """
    def __init__(self):
        # state parameters
        self.age = None  #: °Cd
        self.amino_acids = 0.000075  #: :math:`\mu mol N
        self.proteins = 0.0011  #: :math:`\mu mol N
        self.sucrose = 0.000384  #: :math:`\mu mol C
        self.temperature = 20  #: °C
        self.mstruct = 2.65E-08    #: g TODO: viens d'où?

        # intermediate variables
        self.osmotic_water_potential = -0.8  #: MPa
        self.total_water_potential = -0.1    #: MPa
        self.turgor_water_potential = 0.7    #: MPa
        self.xylem_water_potential = -0.1  #: MPa
        self.soil_water_potential = -0.1  #: MPa
        self.length = 0.05 #: m
        self.width = 1E-3  #: m
        self.thickness = 1E-3 #: m

        # state variables
        self.water_content = -R * (self.temperature + CELSIUS_2_KELVIN) * (
                 ((self.sucrose * 1E-6 / NB_C_SUCROSE) * VANT_HOFF_SUCROSE) +
                 ((self.amino_acids * 1E-6 / AMINO_ACIDS_N_RATIO) * VANT_HOFF_AMINO_ACIDS) +
                 ((self.proteins * 1E-6 / AMINO_ACIDS_N_RATIO) * VANT_HOFF_AMINO_ACIDS)) / self.osmotic_water_potential  #: g

#: The instance of class :class:`cnwheat.parameters.RootsInitCompartments` for current process
ROOTS_INIT_COMPARTMENTS = RootsInitCompartments()


class XylemParameters(object):
    """
    Internal parameters of xylem.
    """
    def __init__(self):
        super(XylemParameters, self).__init__()

        self.VSTORAGE = 0.2  #: Storage portion of the xylem volume
        ##self.R_xylem = 0.00001  #: Elemental xylem flow resistance (Mpa h g-1 m)
        # self.R_xylem = 0.5    #: Flow resistance between xylem and shoot organs (Mpa s g-1 m)
        # self.R_soil = 0.01    #: Flow resistance between soil and xylem (Mpa s g-1 m)
        self.R_xylem = 2    #: Flow resistance between xylem and shoot organs (Mpa s g-1 m)
        self.R_soil = 0.001    #: Flow resistance between soil and xylem (Mpa s g-1 m)

#: The instance of class :class:`cnwheat.parameters.XylemParameters` for current process
XYLEM_PARAMETERS = XylemParameters()


class XylemInitCompartments(object):
    """
    Initial values for compartments of xylem.
    """
    def __init__(self):
        # state parameters

        # intermediate variables
        self.total_water_potential = -0.1    #: MPa
        self.soil_water_potential = -0.1  #: MPa
        self.SRWC = 80  #: %

#: The instance of class :class:`cnwheat.parameters.XylemInitCompartments` for current process
XYLEM_INIT_COMPARTMENTS = XylemInitCompartments()


class PhotosyntheticOrganElementParameters(object):
    """
    Internal parameters of Photosynthetic Organ Element.
    """
    def __init__(self):
        super(PhotosyntheticOrganElementParameters, self).__init__()

        # V1
        ##self.epsilon = {'x': 50, 'y': 50, 'z': 50}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        # V2
        ##self.epsilon = {'x': 0.4E-04, 'y': 0.1E-06, 'z': 0.4E-02}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        # V3
        ##self.epsilon = {'x': 0.4, 'y': 0.1, 'z': 0.4}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.

        self.epsilon = {'x': 4, 'y': 4, 'z': 0.7}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.

        ##self.epsilon = {'x': 50, 'y': 50, 'z': 50}  #: Dimension-specific elasticity in relation to reversible dimensional changes (MPa). x: width, y: thickness, z: length.
        ##https://acsess.onlinelibrary.wiley.com/doi/epdf/10.2134/agronj1979.00021962007100010008x
        ##entre 10 et 30 MPa pour wheat et 18 wheatgrass

        self.vstorage = 0.8  #: Storage portion of the organ volume


#: The instance of class :class:`cnwheat.parameters.PhotosyntheticOrganElementParameters` for current process
PHOTOSYNTHETIC_ORGAN_ELEMENT_PARAMETERS = PhotosyntheticOrganElementParameters()

class PhotosyntheticOrganElementInitCompartments(object):
    """
    Initial values for compartments of photosynthetic organ elements.
    """
    def __init__(self):
        # state parameters
        self.age = None                         #: °Cd
        self.sucrose = 15                    #: :math:`\mu mol C
        self.amino_acids = 2                 #: :math:`\mu mol N
        self.proteins = 25                   #: :math:`\mu mol N
        self.temperature = 20                #: °C
        self.Ts = 12                         #: °C
        self.mstruct = 0                     #: g
        self.width = 0.003                   #: m
        self.thickness = 0.0005              #: m
        self.green_area = 1E-4   #: initial value of green_area (m2)
        self.Tr = 0                          #: mmol H20 m-2 s-1


        # intermediate variables
        self.osmotic_water_potential = -0.8  #: MPa
        self.total_water_potential = -0.1    #: MPa
        self.turgor_water_potential = 0.7    #: MPa
        self.length = 4E-5                   #: m

        # state variables
        self.water_content = -R * (self.temperature + CELSIUS_2_KELVIN) * (
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

#: The instance of class :class:`turgorgrowth.parameters.LaminaParameters` for current process
LAMINA_ELEMENT_PARAMETERS = LaminaElementParameters()


class LaminaElementInitCompartments(object):
    """
    Initial values for compartments of lamina elements.
    """
    def __init__(self):
        # state parameters
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

        self.vstorage = 0.8  #: Storage portion of the organ volume

#: The instance of class :class:`turgorgrowth.parameters.InternodeParameters` for current process
INTERNODE_ELEMENT_PARAMETERS = InternodeElementParameters()


class InternodeElementInitCompartments(object):
    """
    Initial values for compartments of internode elements.
    """
    def __init__(self):
        # state parameters
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

        self.vstorage = 0.8  #: Storage portion of the organ volume


#: The instance of class :class:`turgorgrowth.parameters.InternodeParameters` for current process
SHEATH_ELEMENT_PARAMETERS = SheathElementParameters()


class SheathElementInitCompartments(object):
    """
    Initial values for compartments of sheath elements.
    """
    def __init__(self):
        # state parameters
        self.width = 0.003       #: m
        self.thickness = 0.0005  #: m

#: The instance of class :class:`turgorgrowth.parameters.SheathElementInitCompartments` for current process
SHEATH_ELEMENT_INIT_COMPARTMENTS = SheathElementInitCompartments()
