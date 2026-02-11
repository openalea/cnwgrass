# -*- coding: latin-1 -*-

from __future__ import division  # use "//" to do integer division

"""
    growthwheat.parameters
    ~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`growthwheat.parameters` defines the constant parameters.

    :copyright: Copyright 2015 INRA-ECOSYS, see AUTHORS.
    :license: see LICENSE for details.

"""

C_MOLAR_MASS = 12                     #: Carbon molar mass (g mol-1)
N_MOLAR_MASS = 14                     #: Nitrogen molar mass (g mol-1)

CONVERSION_FACTOR_20_TO_12 = 0.45     #: modified_Arrhenius_equation(12)/modified_Arrhenius_equation(20)

HEXOSE_MOLAR_MASS_C_RATIO = 0.42        #: Contribution of C in hexose mass
PROTEINS_MOLAR_MASS_N_RATIO = 0.151     #: Mean contribution of N in protein mass (Penning De Vries 1989)
AMINO_ACIDS_MOLAR_MASS_N_RATIO = 0.135  #: Mean contribution of N in amino acids mass of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly)

# Shoot
ALPHA = 0.106  # 1.537e-02            #: Parameter of the relation between leaf mass and leaf length (g m^(-BETA))
BETA = 1.28                           #: Parameter of the relation between leaf mass and leaf length (dimensionless)
RATIO_SUCROSE_MSTRUCT = 0.444         #: Mass of C (under carbohydrate form, g) in 1 g of mstruct : 0.384 from (Penning de Vries, Witlage and Kremer, 1978), 0.444 for cellulose (about 45% C in
# senescent leaves)
RATIO_AMINO_ACIDS_MSTRUCT = 0.005     #: Mass of N (under amino acid/protein form, g) in 1 g of mstruct (Penning de Vries, Witlage and Kremer, 1978)
AMINO_ACIDS_C_RATIO = 4.15            #: Mean number of mol of C in 1 mol of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly) from (Penning de Vries, Witlage and Kremer, 1978)
AMINO_ACIDS_N_RATIO = 1.25            #: Mean number of mol of N in 1 mol of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly) from (Penning de Vries, Witlage and Kremer, 1978)
RATIO_MSTRUCT_DM = 0.8                #: Ratio mstruct/dry matter (dimensionless)
RATIO_ENCLOSED_LEAF_INTERNODE = 5     #: We use ratio sheath:lamina of the specific structural dry masses (from data of J. Bertheloot, 2004)
INIT_CYTOKININS_EMERGED_TISSUE = 200  #: Initial amount of cytokinins allocated in the mstruct of a newly emerged tissue
# Leaf Automate elongation
te = 300 * 3600 * 24 / 12  #: end of leaf elongation in automate growth (s at 12°c); fitted from adapted data from Fournier 2005
FITTED_L0 = 0.01557936     #: Fitted value of leaf length at t=0 after rescaling the beta function with L0 (m); Fournier 2005 sur courbe corrigee

# Internode Automate elongation
FITTED_L0_IN = 1/59.0         #: Scaling factor of the internode in automate growth (dimensionless), fitted from Malvoisin 1984 II
te_IN = 210 * 3600 * 24 / 12  #: end of internode elongation in automate growth (s at 12°c) ; fitted from Malvoisin 1984 II
internode_rapid_growth_t = 288000  # Pseudo age at which the internode starts the rapid growth stage (s)

# Roots
conc_sucrose_offset = 240  #: sucrose concentration under which no root growth occurs (µmol of C g-1 mstruct) - same for leaf and internode growth, see elongwheat

VMAX_ROOTS_GROWTH_POSTFLO = 0.015 * CONVERSION_FACTOR_20_TO_12  #: Maximal rate of root structural dry matter growth (µmol C s-1 g-1 MS) post flo at 12°C
VMAX_ROOTS_GROWTH_PREFLO = 0.0855 * CONVERSION_FACTOR_20_TO_12  #: Maximal rate of root structural dry matter growth (µmol C s-1 g-1 MS) pre flo at 12°C
K_ROOTS_GROWTH = 1250                 #: Affinity coefficient of root structural dry matter growth (µmol C g-1 MS) post flo
N_ROOTS_GROWTH = 1.8

RATIO_C_MSTRUCT_ROOTS = 0.444         #: Mean contribution of carbon to root structural dry mass (g C g-1 Mstruct) : same as shoot
RATIO_N_MSTRUCT_ROOTS_ = 0.005        #: Mean contribution of nitrogen to root structural dry mass (g N g-1 Mstruct) : same as shoot

CONVERSION_MMOL_C_G_MSTRUCT_ROOTS = (1E-6 * C_MOLAR_MASS) / RATIO_C_MSTRUCT_ROOTS

MINERAL_LIVING_TISSUE = 0.05  #: Mineral content of the structural mass (g g-1 mstruct) (Thornley and Cannell, 2000)
MINERAL_SENESCED_TISSUE = 0.05/2  #: Mineral content of the structural mass of a senescent tissue (g g-1 mstruct) (Thornley and Cannell, 2000)


class OrganInit:
    """
    Initial values for organs
    """
    def __init__(self):
        self.is_growing = True
        self.senesced_length_element = 0      #: m
        self.senesced_mstruct = 0     #: g
        self.green_area = 0           #: m2
        self.sucrose = 0              #: µmol C
        self.amino_acids = 0          #: µmol N
        self.fructan = 0              #: µmol C
        self.proteins = 0             #: µmol N
        self.nitrates = 0             #: µmol N
        self.mstruct = 0              #: g
        self.Nstruct = 0              #: g
        self.Nresidual = 0            #: g
        self.cytokinins = 0           #: g
        self.conc_cytokinins = 200.0  #: AU / g mstruct # Not used
