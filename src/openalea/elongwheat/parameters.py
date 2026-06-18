# -*- coding: latin-1 -*-

from __future__ import division  # use "//" to do integer division

"""
    elongwheat.parameters
    ~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`elongwheat.parameters` defines the constant parameters.

"""


class Parameters(object):
    def __init__(self):
        # -------------------------------------------------------------------------------------------------------------------
        # --- SAM
        # -------------------------------------------------------------------------------------------------------------------
        # Parameters for temperature responses
        self.Temp_Tref = 12  # Arbitrary reference temperature (°C)
        self.Temp_Ea_R = 8900  # Parameter Ea/R in Eyring equation from Johnson and Lewin (1946) - Parameter value fitted from Kemp and Blacklow (1982) (K)
        self.Temp_DS_R = 68.432  # Parameter deltaS/R in Eyring equation from Johnson and Lewin (1946) - Parameter value fitted from Kemp and Blacklow (1982) (dimensionless)
        self.Temp_DH_R = 20735.5  # Parameter deltaH/R in Eyring equation from Johnson and Lewin (1946) - Parameter value fitted from Kemp and Blacklow (1982) (K)
        self.Temp_Ttransition = 9  # Below this temperature f = linear function of temperature instead of Arrhenius-like(°C)

        # Development
        self.PLASTOCHRONE = 70 / self.Temp_Tref * 24 * 3600  #76.1/12*24*3600    #: Leaf plastochron (s at 12°C) calculated from Ljutovac 2002 with primordia of 5E-5 m (76 dd) ; Malvoisin 35dd associated with init 3E-5 m
        self._delta_TT_GA = None
        self.max_nb_leaves = 11  #: Max number of leaves per axis

        # -------------------------------------------------------------------------------------------------------------------
        # --- Leaves
        # -------------------------------------------------------------------------------------------------------------------
        self.conc_sucrose_offset = 240  #: sucrose concentration under which no leaf and internode growth occurs (µmol of C g-1 mstruct) - same for root growth, see growthwheat

        # Exponential elongation
        self.RERmax_Ljutovac_fit = {2: 2.55e-06, 3: 1.85e-06, 4: 1.7e-06, 5: 1.6e-06, 6: 1.5e-06, 7: 1.5e-06,
                                    8: 1.5e-06, 9: 1.5e-06, 10: 1.4e-06,
                                    11: 1.4e-06}  # Optimal RERmax (s-1 at 12°C) calibrated from Ljutovac (2002)
        # { 5 : 0.00000279 , 6 : 0.00000176 , 7 : 0.00000162 , 8 : 0.00000144 , 9 : 0.00000144 , 10 : 0.00000144 , 11 : 0.00000142 } # Observed RER on data of Ljutovac 2002 RER (s-1 at 12°C)
        #RERmax = {2: 2.1e-06, 3: 2.5e-06, 4: 2.68e-06, 5: 2.65e-06, 6: 2.65e-06, 7: 2.65e-06, 8: 2.5e-06, 9: 2.3e-06, 10: 2.2e-06, 11: 2.2e-06}  # RERmax (s-1 at 12°C) fitted for simulations accounting for metabolic regulation
        # self.RERmax = {2: 2.e-06, 3: 2.33e-06, 4: 2.19e-06, 5: 2.19e-06, 6: 2.25e-06, 7: 2.15e-06, 8: 2.1e-06, 9: 1.9e-06, 10: 1.75e-06, 11: 1.7e-06}  # RERmax (s-1 at 12°C) fitted for simulations accounting for metabolic regulation
        self.RERmax = {2: 1.85e-06, 3: 2.57e-06, 4: 2.33e-06, 5: 2.22e-06, 6: 2.18e-06, 7: 2.05e-06, 8: 1.9e-06,
                       9: 1.75e-06, 10: 1.7e-06,
                       11: 1.7e-06}  # RERmax (s-1 at 12°C) fitted for simulations accounting for metabolic regulation

        self.RER_Kc = 100  #: affinity coefficient of RER to C (µmol g-1)
        self.RER_Kn = 15  #: affinity coefficient of RER to N (µmol g-1)

        # Automate elongation
        #- Parameters for leaf 1 and 2
        self.te_init = 250 * 3600 * 24 / 12
        self.tm_init = 150 * 3600 * 24 / 12
        self.tb_init = -35 * 3600 * 24 / 12
        # - Following leaves
        self.te = 300 * 3600 * 24 / self.Temp_Tref  #: end of leaf elongation in automate growth (s at 12°c); fitted from adapted data from Fournier 2005
        self.tm = 204.6 * 3600 * 24 / self.Temp_Tref  #: time at which leaf elongation rate is maximal in automate growth (s at 12°c); fitted from adapted data from Fournier 2005
        self.tb = -114.3 * 3600 * 24 / self.Temp_Tref  #: beginning of leaf elongation in automate growth (s at 12°c); fitted from adapted data from Fournier 2005
        # NB : Previous fit on adapted data from Fournier 2005 in phyllochronic time te = 271, tm=176, tb=-25
        self.leaf_Lmax_MAX = 0.45  #: Maximum leaf_Lmax (m)

        self.leaf_pseudo_age_Vmax = 1.2  #: Maximal regulation of leaf length after emergence by CN status (dimensionless)
        self.leaf_pseudo_age_Kc = 150  #: affinity coefficient to C (µmol g-1)
        self.leaf_pseudo_age_Kn = 4  #: affinity coefficient to N (µmol g-1)

        # Leaf maximal width
        #- If optimal growth used
        self.leaf_Wmax_dict = {1: 0.0030, 2: 0.0033, 3: 0.0040, 4: 0.0048, 5: 0.00545, 6: 0.00734, 7: 0.00913,
                               8: 0.0111, 9: 0.0163, 10: 0.014, 11: 0.018}  #: m (Ljutovac 2002)
        self.leaf_Wmax_Marion = {1: 0.0030, 2: 0.0033, 3: 0.0040, 4: 0.0045, 5: 0.0056, 6: 0.0075, 7: 0.010, 8: 0.012,
                                 9: 0.013, 10: 0.014, 11: 0.018}  #: m (Ljutovac 2002)
        # - If optimal growth not used
        self.leaf_W_L_MIN = 0.025  #: Minimal leaf width
        self.leaf_W_L_a = -0.02033728
        self.leaf_W_L_b = -0.00005445836
        self.leaf_W_L_c = 0.000459551

        # Structural Specific Lamina Weight
        self.leaf_SSLW = {1: 21, 2: 21, 3: 21, 4: 22, 5: 22, 6: 22, 7: 22, 8: 24, 9: 25, 10: 28, 11: 31}  # SSLW (g m-2)
        self.leaf_SSLW_Marion = {1: 21, 2: 21, 3: 22, 4: 17, 5: 19.58, 6: 20.9, 7: 23.55, 8: 25.63, 9: 28.65, 10: 28,
                                 11: 31}  # SSLW (g m-2)
        self.leaf_SSLW_NEMA = {1: 15, 2: 23, 3: 25, 4: 24, 5: 21, 6: 18, 7: 16, 8: 18, 9: 21, 10: 26,
                               11: 33}  # Manip NEMA 05/06 traitments N+ (from data of J. Bertheloot, 2004) sauf pour F7/F8
        self.leaf_SSLW_MIN = 5.
        self.leaf_SSLW_MAX = 45.
        self.leaf_SSLW_a = 47.50516  # g m-2
        self.leaf_SSLW_b = 2927.944  # µmol C g-1

        self.leaf_LSSW_dict = {1: 0.06, 2: 0.08, 3: 0.08, 4: 0.18, 5: 0.17, 6: 0.21, 7: 0.24, 8: 0.4, 9: 0.5, 10: 0.55,
                               11: 0.65}  # Manip NEMA 05/06 Soissons N+ (from data of J. Bertheloot, 2004)
        self.leaf_LSSW_dict_Marion = {1: 0.06, 2: 0.08, 3: 0.09, 4: 0.16, 5: 0.21, 6: 0.26, 7: 0.32, 8: 0.37, 9: 0.43,
                                      10: 0.55, 11: 0.65}
        self.leaf_LSSW_a = 0.00005
        self.leaf_LSSW_integral_min = 1700
        self.leaf_LSSW_nominal_A = 0.0403
        self.leaf_LSSW_nominal_B = -0.0099
        self.leaf_LSSW_MIN = 0.05
        self.leaf_LSSW_MAX = 0.8

        # Sheath: lamina ratio. Parameters of a polynomial function, from Dornbush 2011
        self.SL_ratio_a = -0.0021
        self.SL_ratio_b = 0.037
        self.SL_ratio_c = - 0.1527
        self.SL_ratio_d = 0.4962

        # Ratio of leaf length composed by the division zone.
        # Parameters are used for an inverse beta function representing the ratio of the leaf composed by the division zone according to its relative length in log.
        # The model was fitted on literature data on wheat (Fournier 2005, Beemster and Masle 1996, Schuppler 1998).
        self.ratio_DZ_l_init = 0.065  #: normalized log of leaf length at which the leaf is fully composed by the division zone (dimensionless).
        self.ratio_DZ_l_mid = 0.075  #: intermediate point of the beta function (dimensionless).
        self.ratio_DZ_l_end = 0.7  #: normalized log of leaf length at which the leaf has no more division zone (dimensionless).

        # -------------------------------------------------------------------------------------------------------------------
        # --- Internodes
        # -------------------------------------------------------------------------------------------------------------------
        # Exponential elongation
        self.RERmax_dict_IN = {1: 2.4E-06, 2: 2.4E-06, 3: 1.85E-06, 4: 3.E-08, 5: 3.7E-06, 6: 3.3E-06, 7: 3.E-06,
                               8: 2.E-06, 9: 2.E-06, 10: 2.5E-06, 11: 2.5E-06, 12: 1.5E-06}  #: s-1 at 12°C FIT jan 20
        # RERmax_dict_IN = {1: 2.4E-06, 2: 2.4E-06, 3: 2.4E-06, 4: 2.4E-06, 5: 2.6E-06, 6: 2.6E-06, 7: 2.3E-06, 8: 1.9E-06, 9: 1.91E-06, 10: 1.9E-06, 11: 1.76E-06, 12: 1.5E-06}  #: s-1 at 12°C FIT jan 20

        # { 3 : 2.48E-06 ,4 : 2.48E-06 ,5 : 2.48E-06 , 6 : 2.48E-06 , 7 : 2.48E-06 , 8 : 2.48E-06 , 9 : 2.48E-06 , 10 : 1.9E-06 , 11 : 1.6E-06 }#: s-1 at 12°C
        # estimate from Ljutovac 2002 over the period until leaf ligulation i.e. wider than in the model.
        # Because i) not enough data if we consider only up to previous leaf ligulation, ii) same exponential like period

        # SCALING_FACTOR_INT = 59    #: Scaling factor of the internode in automate growth (dimensionless), Malvoisin 1984 II
        # SCALING_FACTOR_INT = 53    #: Scaling factor of the internode in automate growth (dimensionless), Ljutovac 2002, 250pl.m-2

        # Initiation of internode elongation
        self.nb_PLASTO_internode_init = 5  #: Delay between leaf initiation and internode initiation expressed as a number of plastochron. From Malvoisin 1984b, associated with primordia of 5.10-4 m
        self.internode_L_init = 5E-5       #: Initial internode length (m)

        # Automate elongation
        # te_IN = 210 * 3600 * 24 /12  #: end of internode elongation in automate growth; Malvoisin 1984 II
        # tm_IN = 156 * 3600 * 24 /12  #: time at which internode elongation rate is maximal in automate growth (s); Malvoisin 1984 II
        # tb_IN = -70 * 3600 * 24 /12  #: beginning of internode elongation in automate growth (s); Malvoisin 1984 II
        self.te_IN = 331.7538 * 3600 * 24 / 12  #: end of internode elongation in automate growth; Ljutovac 2002, 250pl.m-2
        self.tm_IN = 252.7798 * 3600 * 24 / 12  #: time at which internode elongation rate is maximal in automate growth (s);Ljutovac 2002, 250pl.m-2
        self.tb_IN = -374.2918 * 3600 * 24 / 12  #: beginning of internode elongation in automate growth (s);Ljutovac 2002, 250pl.m-2

        self.ratio_LSIW_LSSW = 2.5  #: ratio lineic structural internode mass / lineic structural sheath mass  of the specific structural dry masses (from data of J. Bertheloot, 2004)
        self.internode_LSIW_dict = {1: 2.8, 2: 2.8, 3: 2.8, 4: 2.8, 5: 2.8, 6: 2.8, 7: 2.8, 8: 2.8, 9: 2.3, 10: 1.7,
                                    11: 1.6, 12: 1.4,
                                    13: 0.7}  #: experiment of M.Gauthier 2017/18, consistent with R.Barillot 2014

        # -------------------------------------------------------------------------------------------------------------------
        # --- Coleoptiles
        # -------------------------------------------------------------------------------------------------------------------
        self.delta_L_emergence = 0.003  #: length of visible coleoptile ie emerged from soil (m). Current value set allows coleo to match length of sheath 1 in Sreten.


    @property
    def delta_TT_GA(self):
        if self._delta_TT_GA is not None:
            return self._delta_TT_GA
        else:
            return self.PLASTOCHRONE * 5  #: Thermal time between floral transition of SAM and Gibberelin production expressed as a function of plastochron (s at 12°C) ; Malvoisin's data give 7 plastochrons

    @delta_TT_GA.setter
    def delta_TT_GA(self, value):
        """ Allows the value to be modified externally from simulation.py"""
        self._delta_TT_GA = value


class ParametersHydraulics(Parameters):
    def __init__(self):
        super().__init__()
        self.PLASTOCHRONE = 76.1 / self.Temp_Tref * 24 * 3600
        self.RERmax = {1: 2.4E-06, 2: 2.4E-06, 3: 1.85E-06, 4: 3.E-08, 5: 3.45e-06, 6: 2.3e-06, 7: 2.2e-06, 8: 2.1e-06, 9: 2.0e-06, 10: 1.85e-06, 11: 1.76e-06}  # RERmax (s-1 at 12°C) fitted for simulations accounting for metabolic regulation
        self.n = 4  #: Shape parameter of the xylem water potential function in RER (-)
        self.psi_ref = -0.185  #: Critical leaf water potential threshold at which RER is reduced by 50% (Mpa) (adapted from Boyer, 1968)
        self.lamina_Lmax_dict = {3: 0.09896253, 4: 0.10463127, 5: 0.10867923, 6: 0.1352669, 7: 0.16634665,
                                 8: 0.19288286,
                                 9: 0.25080917, 10: 0.28677033, 11: 0.3}  #: Maximum lamina length (m)
        self.GAMMA = 0.2  #: Critical value for the pressure component which must be exceeded for irreversible volume changes (MPa). Found from Coussement et al., 2018 : 0.3 Mpa for soybean.


def parameters_factory(hydraulics=False):
    """
    Factory for parameter instantiation according to the version of the model used

    :param bool hydraulics: if True the model will assume the coupling to the turgor-driven growth model

    :return: Instance of Parameters oo ParametersHydraulics
    """

    if not hydraulics:
        return Parameters()
    else:
        return ParametersHydraulics()

class HiddenZoneInit(object):
    """
    Initial values for hidden zones
    """

    def __init__(self):
        self.leaf_is_growing = True
        self.internode_is_growing = False
        self.leaf_pseudostem_length = 4E-5  #: m
        self.delta_leaf_pseudostem_length = 0  #: m
        self.internode_distance_to_emerge = 0  #: m
        self.delta_internode_distance_to_emerge = 0  #: m, needed for growthwheat
        self.leaf_L = 5E-5  #: m, should be consistent with PLASTOCHRONE
        self.delta_leaf_L = 0  #: m, needed for growthwheat
        self.internode_L = 0  #: m
        self.delta_internode_L = 0  #: m, needed for growthwheat
        self.leaf_Lmax = None  #: m, no calculation before emergence Ln-1
        self.leaf_Lmax_em = None  #: m, no calculation before emergence Ln-1
        self.lamina_Lmax = None  #: m, no calculation before emergence Ln-1
        self.sheath_Lmax = None  #: m, no calculation before emergence Ln-1
        self.leaf_Wmax = None  #: m, no calculation before emergence Ln-1
        self.SSLW = None  #: g m-2, no calculation before emergence Ln-1
        self.LSSW = None  #: g m-1, no calculation before emergence Ln-1 (about 2)
        self.leaf_is_emerged = False
        self.internode_Lmax = None  #: m, no calculation before ligulation Ln
        self.internode_Lmax_lig = None  #: m, no calculation before ligulation Ln
        self.LSIW = None  #: g m-1, no calculation before ligulation Ln
        self.internode_is_visible = False
        self.leaf_pseudo_age = 0
        self.internode_pseudo_age = 0
        self.delta_leaf_pseudo_age = 0
        self.delta_internode_pseudo_age = 0
        self.hiddenzone_age = 0
        self.is_over = False
        self.leaf_is_remobilizing = False
        self.internode_is_remobilizing = False
        self.ratio_DZ = 1.0

        # Default values used for RER calculation in elongwheat
        self.sucrose = 5E-6  #: µmol C
        self.amino_acids = 4E-6  #: µmol N
        self.fructan = 0  #: µmol C - about 10% DM
        self.leaf_enclosed_mstruct = 1.26E-07  #: g
        self.internode_enclosed_mstruct = 0  #: g
        self.mstruct = self.leaf_enclosed_mstruct + self.internode_enclosed_mstruct  #: g
        self.leaf_enclosed_Nstruct = self.leaf_enclosed_mstruct * 0.005  #: g, parameter value in growth wheat
        self.internode_enclosed_Nstruct = self.internode_enclosed_mstruct * 0.0322  #: g, parameter value in growth wheat
        self.Nstruct = self.leaf_enclosed_Nstruct + self.internode_enclosed_Nstruct  #: g
        self.proteins = 2.6E-03  #: µmol N - about 9% N
        self.conc_cytokinins = 150  #: AU / g mstruct
        self.mean_conc_sucrose = 0  #: µmol C / g mstruct


class HiddenZoneInitHydraulics(HiddenZoneInit):
    """
    Initial values for hidden zones
    """

    def __init__(self):
        super().__init__()
        self.sheath_is_emerged = False
        self.init_leaf_L = 0  #: m, needed for growth-wheat
        self.leaf_pseudo_age = -1  #: s, no calculation before emergence Ln-1


class ElementInit(object):
    """
    Initial values for emerged and growing elements
    """

    def __init__(self):
        self.is_growing = True
        self.is_over = False
        self.length = 0  #: m
        self.Wmax = None  # :m maximum lamina width (could store sheath and internode diameter as well)
        self.senesced_length_element = 0  #: m
        self.green_area = 0  #: m2
        self.age = 0  #: Thermal Time
        self.age_teq = 0  #: s
        self.max_proteins = 0  #: µmol N
        self.Nresidual = 0  #: g
        self.sucrose = 0  #: µmol C
        self.amino_acids = 0  #: µmol N
        self.nitrates = 0  #: µmol N
        self.fructan = 0  #: µmol C
        self.starch = 0  #: µmol C
        self.proteins = 0  #: µmol N
        self.mstruct = 0  #: g
        self.senesced_mstruct = 0  #: g
        self.max_mstruct = 0  #: g
        self.Nstruct = 0  #: g
        self.cytokinins = 0  #: g
