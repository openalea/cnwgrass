# -*- coding: latin-1 -*-

"""
    farquharwheat.paramters
    ~~~~~~~~~~~~~~~~~~~

    The module :mod:`farquharwheat.parameters` defines the constant parameters.

    :copyright: Copyright 2014-2015 INRA-ECOSYS, see AUTHORS.
    :license: see LICENSE for details.

"""
# -- Version of the photosynthesis sub model
SurfacicProteins = True      #: If True, surfacic proteins used to regulate photosynthesis ; if not total N
NSC_Retroinhibition = True   #: If True, NSC (Non-Structural Carbohydrates) downregulate photosynthesis
prim_scale = False           #: If True, photosynthesis calculated at primitive scale, if not at organ scale

if not SurfacicProteins:
    # Used in Barillot et al. (2016) and Gauthier et al. (2020)
    # Dependence to surfacic_nitrogen including structural nitrogen

    # -- Nitrogen dependance of photosynthetic parameters (derived from Braune et al. (2009) and Evers et al. (2010):
    #     * S_surfacic_nitrogen: slope of the relation between surfacic_nitrogen and the parameter
    #         * alpha: mol e- m2 mol-1 photon g-1 N
    #         * Vc_max25: µmol CO2 g-1 N s-1
    #         * Jmax25: µmol e- g-1 N s-1
    #         * TPU25: µmol CO2 g-1 N s-1
    #         * Rdark25: µmol CO2 g-1 N s-1
    #     * surfacic_nitrogen_min: minimum amount of nitrogen below which photosynthesis rate is zero (g (N) m-2)
    #     * beta: intercept parameter of the relation between alpha and surfacic_nitrogen (mol e- mol-1 photons)
    #     * delta1 and delta2: parameters of m (scaling factor of gs) dependance to surfacic_nitrogen (m2 g-1 and dimensionless, respectively)
    PARAM_N = {'S_surfacic_nitrogen': {'Vc_max25': 84.965, 'Jmax25': 117.6, 'alpha': 0.0413, 'TPU25': 9.25, 'Rdark25': 0.493},
               'surfacic_nitrogen_min': {'Vc_max25': 0.17, 'Jmax25': 0.17, 'TPU25': 0.229, 'Rdark25': 0.118}, 'beta': 0.2101, 'delta1': 14.7, 'delta2': -0.548}
    NA_0 = 2  # Initial value of surfacic_nitrogen (g m-2), used if no surfacic_nitrogen is provided by user
    Psurf_to_SLNnonstruct = 1.06  #: Conversion factor from surfacic protein content to non-structural SLN (estimation from NEMA and Ljutovac simulations)

else:
    # Used in Gauthier et al. (2021)
    # Dependence to surfacic_nitrogen without structural nitrogen

    # -- Nitrogen dependance of photosynthetic parameters (derived from Braune et al. (2009) and Evers et al. (2010):
    #     * S_surfacic_nitrogen: slope of the relation between surfacic_nitrogen and the parameter
    #         * alpha: mol e- m2 mol-1 photon g-1 N
    #         * Vc_max25: µmol CO2 g-1 N s-1
    #         * Jmax25: µmol e- g-1 N s-1
    #         * TPU25: µmol CO2 g-1 N s-1
    #         * Rdark25: µmol CO2 g-1 N s-1
    #     * surfacic_nitrogen_min: minimum amount of nitrogen below which photosynthesis rate is zero (g (N) m-2)
    #     * beta: intercept parameter of the relation between alpha and surfacic_nitrogen (mol e- mol-1 photons)
    #     * delta1 and delta2: parameters of m (scaling factor of gs) dependance to surfacic_nitrogen (m2 g-1 and dimensionless, respectively)
    PARAM_N = {'S_surfacic_nitrogen': {'Vc_max25': 84.965, 'Jmax25': 117.6, 'alpha': 0.0413, 'TPU25': 9.25, 'Rdark25': 0.493},
               'surfacic_nitrogen_min': {'Vc_max25': 0., 'Jmax25': 0., 'TPU25': 0., 'Rdark25': 0.}, 'beta': 0.2101 + 0.0083, 'delta1': 14.7, 'delta2': -0.548}
    NA_0 = 2  # Initial value of surfacic_nitrogen (g m-2), used if no surfacic_nitrogen is provided by user
    Psurf_to_SLNnonstruct = 1.06  #: Conversion factor from surfacic protein content to non-structural SLN (estimation from NEMA and Ljutovac simulations)

# -- Molecular weights
MM_WATER = 18  # Molar mass of water (g mol-1)
N_MOLAR_MASS = 14  # Molar mass of nitrogen (g mol-1)

# -- Physical parameter
A = 2.5  # Attenuation coefficient of wind within a wheat canopy. From Campbell and Norman (1998), 2nd edition. Can also be estimated by: A = sqrt((0.2*LAI*h)/sqrt((4*width*h)/(pi*LAI))
GAMMA = 66E-3  # Psychrometric constant (KPa K-1). Mean value
I0 = 1370  # Extraterrestrial solar radiation (W m-2)
K = 0.40  # Von Kármán's constant (dimensionless)
LAMBDA = 2260E3  # Latent heat for vaporisation of water (J kg-1)
RHOCP = 1256  # Volumetric heat capacity of air (J m-3 K-1)
SIGMA = 5.6704E-8  # Stefan-Bolzmann constant (W-2 K-4)
ZR = 2  # Height above canopy at which reference wind (Ur) is measured (m)

s_A = 239  # Factor in the calculation of the Saturation vapour pressure
s_B = 17.4  # Factor in the calculation of the Saturation vapour pressure
s_C = 0.611  # Factor in the calculation of the Saturation vapour pressure
s_expo = 2  # exponent in the calculation of the Saturation vapour pressure

R = 8.3144  # Gas constant (J mol-1 K-1)
KELVIN_DEGREE = 273.15  #: Conversion factor from degree C to Kelvin
PATM = 1.01325E5  # Atmospheric pressure (Pa)

Ur_min = 0.1  # Minimum wind speed (m s-1) TODO: Temporary patch to avoid div 0 error
Zh_d = 0.7  # estimation factor canopy heigh into Zero plane displacement height
Zh_Zo = 0.1  # estimation factor canopy heigh into roughness length
rhb_blade_A = 154  # factor in Boundary layer resistance to heat calculation for blade
rhb_other_A = 1.2E-5  # factor A in Boundary layer resistance to heat calculation for other organs than blades
rhb_other_B = 1.5E-5  # factor B in Boundary layer resistance to heat calculation for other organs than blades
rhb_other_C = 0.47  # factor C in Boundary layer resistance to heat calculation for other organs than blades
ra_expo = 2  # exponent in calculation of Turbulence resistance to heat

PARa_to_RGa = 1.53  # Used to convert PAR absorbed into RG absorbed (see details in notice entitiled "Notes sur le calcul du rayonnement net à partir du PAR absorbé")
Watt_to_PPFD = 4.55  # It is assumed that 1 W m-2 of PAR is equivalent to 4.55 µmol m-2 s-1 of PAR (Goudriaan and Laar, 1994)
rbh_rbw = 0.96  # estimation factor Boundary layer resistance for water from Boundary layer resistance to heat

# -- Photosynthetic parameter
O2 = 21000  # Intercellular O2 concentration, µmol mol(air)-1 or Pa, from Bernacchi et al. (2001)
KC25 = 404  # Affinity constant of RuBisCO for C, µmol mol-1 or Pa, from Bernacchi et al. (2001) (estimation in Braune et al. (2009) not enough accurate)
KO25 = 278.4E3  # Affinity constant of RuBisCO for O2, µmol mol-1 or Pa, from Bernacchi et al. (2001) (estimation in Braune et al. (2009) not enough accurate)
GAMMA25 = 39  # CO2 compensation point, µmol(CO2) mol-1 (air), from Braune et al. (2009)
THETA = 0.72  # curvature parameter of J, dimensionless
EFFICENCY_STEM = 0.78
J_expo = 2  # exponent in calculation of the Electron transport rate (Muller et al. (2005), Evers et al. (2010))
J_A = 4  # factor A in calculation of the Electron transport rate (Muller et al. (2005), Evers et al. (2010))
J_B = 2  # factor B in calculation of the Electron transport rate (Muller et al. (2005), Evers et al. (2010))
Aj_A = 4  # parameter in the calculation of the Rate of assimilation under RuBP regeneration limitation
Aj_B = 8  # parameter in the calculation of the Rate of assimilation under RuBP regeneration limitation
Vomax_A = 0.5  # parameter in the calculation of the Maximum rate of Vo
Ap_A = 3  # parameter in the calculation of the Rate of assimilation under TPU limitation
Rd_A = 0.33  # parameter A in the calculation of the Rd
Rd_B = 0.5  # parameter B in the calculation of the Rd
Rd_C = 15  # parameter C in the calculation of the Rd
Ci_init_ratio = 0.7

# -- Stomatal conductance parameter
GSMIN = 0.05  # Minimum gsw, measured in the dark (mol m-2 s-1). Braune et al. (2009).
GB = 3.5  # Boundary layer conductance to water vapour (mol m-2 s-1). Muller et al., (2005)
K_Cs = 1.37  # factor in CO2 concentration at organ surface. From Prieto et al. (2012).
gsw_gs_CO2 = 1.6  # conversion factor from gsw into gs_CO2
Ci_A = 1.37  # factor in Ci calculation. Comes from (1.6)^(2/3)

# -- Temperature dependance of photosynthetic parameters (parameter values derived from Braune et al. (2009) except for Kc, Ko, and Rdark (Bernacchi et al., 2001))
#     * deltaHa, deltaHd: enthalpie of activation and deactivation respectively (kJ mol-1)
#     * deltaS: entropy term (kJ mol-1 K-1)
#     * Tref: reference temperature (K)

PARAM_TEMP = {'deltaHa': {'Vc_max': 89.7, 'Jmax': 48.9, 'TPU': 47., 'Kc': 79.43, 'Ko': 36.38, 'Gamma': 35., 'Rdark': 46.39},
              'deltaHd': {'Vc_max': 149.3, 'Jmax': 152.3, 'TPU': 152.3},
              'deltaS': {'Vc_max': 0.486, 'Jmax': 0.495, 'TPU': 0.495}, 'Tref': 298.15}

DELTA_CONVERGENCE = 0.01  #: The relative delta for Ci and Ts convergence.

# -- Inhibition of the photosynthesis by carbohydrates (from Azcon-Bieto 1983)
WSC_min = 100000  # Surfacic WSC content above which inhibition of the photosynthesis by WSC occures (µmol C m-2)
Inhibition_max = 1  # Maximum inhibition ratio
K_Inhibition = 938000  # 'Affinity' coefficient for the inhibition of the photosynthesis by WSC (µmol C m-2)


class ElementDefaultProperties(object):
    """
    Properties by default for the elements. Used in FarquharWheat facade.
    """
    def __init__(self):
        self.PARa = 0
        self.PARa_prim = []  #: For calculations at primitive scale
        self.area_prim = 0   #: For calculations at primitive scale
        self.nitrates = 0
        self.proteins = 0
        self.Nstruct = 0
        self.sucrose = 0
        self.fructan = 0
        self.starch = 0
        self.green_area = 0


class AxisDefaultProperties(object):
    """
    Properties by default for the axis. Used in FarquharWheat facade.
    """
    def __init__(self):
        self.height_canopy = 0.78
