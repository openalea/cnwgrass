# -*- coding: latin-1 -*-

from __future__ import division  # use "//" to do integer division

"""
    elongwheat.parameters
    ~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`elongwheat.parameters` defines the constant parameters.

    :copyright: Copyright 2015 INRA-ECOSYS, see AUTHORS.
    :license: see LICENSE for details.

"""

# -------------------------------------------------------------------------------------------------------------------
# --- SAM
# -------------------------------------------------------------------------------------------------------------------
PLASTOCHRONE = 76.1/12*24*3600    #: Leaf plastochron (s at 12°C) calculated from Ljutovac 2002 with primordia of 5E-5 m (76 dd) ; Malvoisin 35dd associated with init 3E-5 m
max_nb_leaves = 11   #: Max number of leaves per axis
delta_TT_GA = PLASTOCHRONE * 5  #: Thermal time between floral transition of SAM and Gibberelin production expressed as a function of plastochron (s at 12°C) ; Malvoisin's data give 7 plastochrons

sowing_depth = 0.05  #: Sowing depth (m) used to define plant emergence

# Parameters for temperature responses
Temp_Tref = 12        # Arbitrary reference temperature (°C)
Temp_Ea_R = 8900      # Parameter Ea/R in Eyring equation from Johnson and Lewin (1946) - Parameter value fitted from Kemp and Blacklow (1982) (K)
Temp_DS_R = 68.432    # Parameter deltaS/R in Eyring equation from Johnson and Lewin (1946) - Parameter value fitted from Kemp and Blacklow (1982) (dimensionless)
Temp_DH_R = 20735.5   # Parameter deltaH/R in Eyring equation from Johnson and Lewin (1946) - Parameter value fitted from Kemp and Blacklow (1982) (K)
Temp_Ttransition = 9  # Below this temperature f = linear function of temperature instead of Arrhenius-like(°C)

# -------------------------------------------------------------------------------------------------------------------
# --- Leaves
# -------------------------------------------------------------------------------------------------------------------
conc_sucrose_offset = 240  #: sucrose concentration under which no leaf and internode growth occurs (µmol of C g-1 mstruct) - same for root growth, see growthwheat

# Exponential elongation
RERmax_Ljutovac_fit = {5: 0.000003, 6: 0.00000175, 7: 0.00000164, 8: 0.00000154, 9: 0.00000151, 10: 0.00000134, 11: 0.00000129} # Optimal RERmax (s-1 at 12°C) allowing to simulate leaf dimensions of Ljutovac (2002)
# { 5 : 0.00000279 , 6 : 0.00000176 , 7 : 0.00000162 , 8 : 0.00000144 , 9 : 0.00000144 , 10 : 0.00000144 , 11 : 0.00000142 } # Observed RER on data of Ljutovac 2002 RER (s-1 at 12°C)
RERmax = {5: 3.35e-06, 6: 2.1e-06, 7: 2.e-06, 8: 1.83e-06, 9: 1.8e-06, 10: 1.65e-06, 11: 1.56e-06}  # RERmax (s-1 at 12°C) fitted for simulations accounting for metabolic regulation
RER_Kc = 100  #: affinity coefficient of RER to C (µmol g-1)
RER_Kn = 15   #: affinity coefficient of RER to N (µmol g-1)

# Automate elongation
te = 300 * 3600 * 24 / 12     #: end of leaf elongation in automate growth (s at 12°c); fitted from adapted data from Fournier 2005
tm = 204.6 * 3600 * 24 / 12   #: time at which leaf elongation rate is maximal in automate growth (s at 12°c); fitted from adapted data from Fournier 2005
tb = -114.3 * 3600 * 24 / 12  #: beginning of leaf elongation in automate growth (s at 12°c); fitted from adapted data from Fournier 2005
# NB : Previous fit on adapted data from Fournier 2005 in phyllochronic time te = 271, tm=176, tb=-25
leaf_Lmax_MAX = 0.45           #: Maximum leaf_Lmax (m)

leaf_pseudo_age_Vmax = 1.2    #: Maximal regulation of leaf length after emergence by CN status (dimensionless)
leaf_pseudo_age_Kc = 150      #: affinity coefficient to C (µmol g-1)
leaf_pseudo_age_Kn = 4        #: affinity coefficient to N (µmol g-1)

# Leaf maximal width  TODO doc
leaf_Wmax_dict = {3: 0.0040, 4: 0.0045, 5: 0.0056, 6: 0.0075, 7: 0.010, 8: 0.012, 9: 0.013, 10: 0.014, 11: 0.018}  #: m (Ljutovac 2002)
leaf_W_L_MIN = 0.025  #: Minimal leaf width
leaf_W_L_a = -0.02033728
leaf_W_L_b = -0.00005445836
leaf_W_L_c = 0.000459551

# Structural Specific Lamina Weight
leaf_SSLW = {1: 22, 2: 22, 3: 22, 4: 22, 5: 22, 6: 22, 7: 22, 8: 24, 9: 25, 10: 28, 11: 31}  # SSLW (g m-2)
leaf_SSLW_NEMA = {1: 15, 2: 23, 3: 25, 4: 24, 5: 21, 6: 18, 7: 16, 8: 18, 9: 21, 10: 26, 11: 33}  # Manip NEMA 05/06 traitments N+ (from data of J. Bertheloot, 2004) sauf pour F7/F8
leaf_SSLW_MIN = 5.
leaf_SSLW_MAX = 45.
leaf_SSLW_a = 47.50516
leaf_SSLW_b = 2927.944

leaf_LSSW_dict = {1: 0.08, 2: 0.09, 3: 0.11, 4: 0.18, 5: 0.17, 6: 0.21, 7: 0.24, 8: 0.4, 9: 0.5, 10: 0.55, 11: 0.65}  # Manip NEMA 05/06 Soissons N+ (from data of J. Bertheloot, 2004)
leaf_LSSW_a = 0.00005
leaf_LSSW_integral_min = 1700
leaf_LSSW_nominal_A = 0.0403
leaf_LSSW_nominal_B = -0.0099
leaf_LSSW_MIN = 0.05
leaf_LSSW_MAX = 0.8

# Sheath: lamina ratio. Parameters of a polynomial function, from Dornbush 2011
SL_ratio_a = -0.0021
SL_ratio_b = 0.037
SL_ratio_c = - 0.1527
SL_ratio_d = 0.4962

# Ratio of leaf length composed by the division zone.
# Parameters are used for an inverse beta function representing the ratio of the leaf composed by the division zone accoring to its relative length in log.
# The model was fitted on litterature data on wheat (Fournier 2005, Beemster and Masle 1996, Schuppler 1998).
ratio_DZ_l_init = 0.065     #: normalized log of leaf length at which the the leaf is fully composed by the division zone (dimensionless).
ratio_DZ_l_mid = 0.075      #: intermediate point of the beta function (dimensionless).
ratio_DZ_l_end = 0.7        #: normalized log of leaf length at which the the leaf has no more division zone (dimensionless).

# -------------------------------------------------------------------------------------------------------------------
# --- Internodes
# -------------------------------------------------------------------------------------------------------------------
# Exponential elongation
RERmax_dict_IN = {3: 2.4E-06, 4: 2.4E-06, 5: 2.4E-06, 6: 2.4E-06, 7: 2.16E-06, 8: 1.8E-06, 9: 1.91E-06, 10: 1.9E-06, 11: 1.76E-06, 12: 1.5E-06}  #: s-1 at 12°C FIT jan 20
# { 3 : 2.48E-06 ,4 : 2.48E-06 ,5 : 2.48E-06 , 6 : 2.48E-06 , 7 : 2.48E-06 , 8 : 2.48E-06 , 9 : 2.48E-06 , 10 : 1.9E-06 , 11 : 1.6E-06 }#: s-1 at 12°C
# estimate from Ljutovac 2002 over the period until leaf ligulation i.e. wider than in the model.
# Because i) not enought data if we consider only up to previous leaf ligulation, ii) same exponential like period

# SCALING_FACTOR_INT = 59  #: Scaling factor of the internode in automate growth (dimensionless), Malvoisin 1984 II
# SCALING_FACTOR_INT = 53    #: Scaling factor of the internode in automate growth (dimensionless), Ljutovac 2002, 250pl.m-2

# Initiation of internode elongation
nb_PLASTO_internode_init = 5  #: Delay between leaf initiation and internode initiation expressed as a number of plastochron. From Malvoisin 1984b, associated with priodia of 5.10-4 m
internode_L_init = 5E-5       #: Initial internode length (m)

# Automate elongation
# te_IN = 210 * 3600 * 24 /12  #: end of internode elongation in automate growth; Malvoisin 1984 II
# tm_IN = 156 * 3600 * 24 /12  #: time at which internode elongation rate is maximal in automate growth (s); Malvoisin 1984 II
# tb_IN = -70 * 3600 * 24 /12  #: beginning of internode elongation in automate growth (s); Malvoisin 1984 II
te_IN = 331.7538 * 3600 * 24 / 12  #: end of internode elongation in automate growth; Ljutovac 2002, 250pl.m-2
tm_IN = 252.7798 * 3600 * 24 / 12  #: time at which internode elongation rate is maximal in automate growth (s);Ljutovac 2002, 250pl.m-2
tb_IN = -374.2918 * 3600 * 24 / 12  #: beginning of internode elongation in automate growth (s);Ljutovac 2002, 250pl.m-2

ratio_LSIW_LSSW = 2.5  #: ratio lineic structural internode mass / lineic structural sheath mass  of the specific structural dry masses (from data of J. Bertheloot, 2004)
internode_LSIW_dict = {1: 2.8, 2: 2.8, 3: 2.8, 4: 2.8, 5: 2.8, 6: 2.8, 7: 2.8, 8: 2.8, 9: 2.3, 10: 1.7, 11: 1.6, 12: 1.4, 13: 0.7}  #: experiment of M.Gauthier 2017/18, consistent with that of R.Barillot 2014


class HiddenZoneInit(object):
    """
    Initial values for hidden zones
    """
    def __init__(self):
        self.leaf_is_growing = True
        self.internode_is_growing = False
        self.leaf_pseudostem_length = 4E-5           #: m TODO: ok comme valeur?
        self.delta_leaf_pseudostem_length = 0        #: m
        self.internode_distance_to_emerge = 0        #: m
        self.delta_internode_distance_to_emerge = 0  #: m, needded for growthwheat
        self.leaf_L = 5E-5                       #: m, should be consistent with PLASTOCHRONE
        self.delta_leaf_L = 0                    #: m, needded for growthwheat
        self.internode_L = 0                     #: m
        self.delta_internode_L = 0               #: m, needded for growthwheat
        self.leaf_Lmax = None                    #: m, no calculation before emergence Ln-1
        self.leaf_Lmax_em = None                 #: m, no calculation before emergence Ln-1
        self.lamina_Lmax = None                  #: m, no calculation before emergence Ln-1
        self.sheath_Lmax = None                  #: m, no calculation before emergence Ln-1
        self.leaf_Wmax = None                    #: m, no calculation before emergence Ln-1
        self.SSLW = None                         #: g m-2, no calculation before emergence Ln-1
        self.LSSW = None                         #: g m-1, no calculation before emergence Ln-1 (about 2)
        self.leaf_is_emerged = False
        self.internode_Lmax = None               #: m, no calculation before ligulation Ln
        self.internode_Lmax_lig = None           #: m, no calculation before ligulation Ln
        self.LSIW = None                         #: g m-1, no calculation before ligulation Ln
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

        # Default values used for RER calculation in elong wheat
        self.sucrose = 5E-6                      #: µmol C
        self.amino_acids = 4E-6                  #: µmol N
        self.fructan = 0                         #: µmol C - about 10% DM
        self.leaf_enclosed_mstruct = 1.26E-07    #: g
        self.internode_enclosed_mstruct = 0      #: g
        self.mstruct = self.leaf_enclosed_mstruct + self.internode_enclosed_mstruct  #: g
        self.leaf_enclosed_Nstruct = self.leaf_enclosed_mstruct * 0.005  #: g, parameter value in growth wheat
        self.internode_enclosed_Nstruct = self.internode_enclosed_mstruct * 0.0322   #: g, parameter value in growth wheat
        self.Nstruct = self.leaf_enclosed_Nstruct + self.internode_enclosed_Nstruct  #: g
        self.proteins = 2.6E-03                   #: µmol N - about 9% N
        self.conc_cytokinins = 150                #: AU / g mstruct
        self.mean_conc_sucrose = 0  #: µmol C / g mstruct


class ElementInit(object):
    """
    Initial values for emerged and growing elements
    """
    def __init__(self):
        self.is_growing = True
        self.is_over = False
        self.length = 0               #: m
        self.Wmax = None                  # :m maximum lamina width (could store sheath and internode diameter as well)
        self.senesced_length_element = 0  #: m
        self.green_area = 0  #: m2
        self.age = 0                  #: Thermal Time
        self.max_proteins = 0         #: µmol N
        self.Nresidual = 0            #: g
        self.sucrose = 0              #: µmol C
        self.amino_acids = 0          #: µmol N
        self.nitrates = 0             #: µmol N
        self.fructan = 0              #: µmol C
        self.starch = 0              #: µmol C
        self.proteins = 0             #: µmol N
        self.mstruct = 0              #: g
        self.senesced_mstruct = 0  #: g
        self.max_mstruct = 0          #: g
        self.Nstruct = 0              #: g
        self.cytokinins = 0           #: g
