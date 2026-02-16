# -*- coding: latin-1 -*-

from __future__ import division  # use "//" to do integer division

import numpy as np

from openalea.growthwheat import parameters

"""
    growthwheat.model
    ~~~~~~~~~~~~~

    The module :mod:`growthwheat.model` defines the equations of the kinetic of leaf growth (mass flows) according to leaf elongation. Also includes root growth.

    :copyright: Copyright 2014-2015 INRA-ECOSYS, see AUTHORS.
    :license: see LICENSE for details.

"""


def calculate_ratio_mstruct_DM(mstruct, sucrose, fructans, amino_acids, proteins):
    """
    Ratio mstruct/dry matter (dimensionless)

    :param float mstruct: Strutural mass (g)
    :param float sucrose: Sucrose amount (�mol C)
    :param float fructans: Fructans amount (�mol C)
    :param float amino_acids: Amino acids amount (�mol N)
    :param float proteins: proteins amount (�mol N)

    :return: Ratio mstruct/dry matter (dimensionless)
    :rtype: float
    """
    C_MOLAR_MASS = parameters.C_MOLAR_MASS
    N_MOLAR_MASS = parameters.N_MOLAR_MASS
    dry_mass = ((sucrose * 1E-6 * C_MOLAR_MASS) / parameters.HEXOSE_MOLAR_MASS_C_RATIO +
                (fructans * 1E-6 * C_MOLAR_MASS) / parameters.HEXOSE_MOLAR_MASS_C_RATIO +
                (amino_acids * 1E-6 * N_MOLAR_MASS) / parameters.AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                (proteins * 1E-6 * N_MOLAR_MASS) / parameters.AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                mstruct)

    return mstruct / dry_mass


def calculate_delta_leaf_enclosed_mstruct(leaf_L, delta_leaf_L, ratio_mstruct_DM):
    """ Relation between length and mstruct for the leaf segment located in the hidden zone during the exponential-like growth phase.
    Parameters alpha_mass_growth and beta_mass_growth estimated from Williams (1960) and expressed in g of dry mass)
    The actual ratio_mstruct_DM is then used to convert in g of structural dry mass.

    :param float leaf_L: Total leaf length (m)
    :param float delta_leaf_L: delta of leaf length (m)
    :param float ratio_mstruct_DM: Ratio mstruct/dry matter (dimensionless)

    :return: delta_leaf_enclosed_mstruct (g)
    :rtype: float
    """
    return parameters.ALPHA * parameters.BETA * leaf_L ** (parameters.BETA - 1) * delta_leaf_L * ratio_mstruct_DM


def calculate_delta_leaf_enclosed_mstruct_postE(delta_leaf_pseudo_age, leaf_pseudo_age, leaf_pseudostem_L, enclosed_mstruct, LSSW):
    """ mstruct of the enclosed leaf from the emergence of the leaf to the end of elongation.
    Final mstruct of the enclosed leaf matches sheath mstruct calculation when it is mature.
    #TODO : Hiddenzone mstruct calculation is not correct for sheath shorten than previous one.

    :param float delta_leaf_pseudo_age: Delta of Pseudo age of the leaf since beginning of automate elongation (s)
    :param float leaf_pseudo_age: Pseudo age of the leaf since beginning of automate elongation (s)
    :param float leaf_pseudostem_L: Pseudostem length (m)
    :param float enclosed_mstruct: mstruct of the enclosed leaf (g)
    :param float LSSW: Lineic Structural Sheath Weight (g m-1).

    :return: delta_leaf_enclosed_mstruct (g)
    :rtype: float
    """
    enclosed_mstruct_max = leaf_pseudostem_L * LSSW

    if leaf_pseudo_age < parameters.te:
        delta_enclosed_mstruct = (enclosed_mstruct_max - enclosed_mstruct) / (parameters.te - leaf_pseudo_age) * delta_leaf_pseudo_age
    else:
        delta_enclosed_mstruct = 0

    return max(0., delta_enclosed_mstruct)


def calculate_delta_internode_enclosed_mstruct(internode_L, delta_internode_L, ratio_mstruct_DM):
    """ Relation between length and mstruct for the internode segment located in the hidden zone.
    Same relationship than for enclosed leaf corrected by RATIO_ENCLOSED_LEAF_INTERNODE.
    Parameters alpha_mass_growth and beta_mass_growth estimated from Williams (1975) and expressed in g of dry mass.
    The actual ratio_mstruct_DM is then used to convert in g of structural dry mass.

    :param float internode_L: Enclosed internode length (m)
    :param float delta_internode_L: delta of enclosed internode length (m)
    :param float ratio_mstruct_DM: Ratio mstruct/dry matter (dimensionless)

    :return: delta_enclosed_internode_mstruct (g)
    :rtype: float
    """
    return parameters.RATIO_ENCLOSED_LEAF_INTERNODE * parameters.ALPHA * parameters.BETA * internode_L ** (parameters.BETA - 1) * delta_internode_L * ratio_mstruct_DM


def calculate_delta_internode_enclosed_mstruct_postL(delta_internode_pseudo_age, internode_pseudo_age, internode_L, internode_pseudostem_L, internode_Lmax, LSIW, enclosed_mstruct):
    """ mstruct of the enclosed internode from the ligulation of the leaf to the end of elongation.
    Final mstruct of the enclosed internode matches internode mstruct calculation when it is mature.

    :param float delta_internode_pseudo_age: Delta of Pseudo age of the internode since beginning of automate elongation (s)
    :param float internode_pseudo_age: Pseudo age of the internode since beginning of automate elongation (s)
    :param float internode_L: Current lenght of the internode (m)
    :param float internode_pseudostem_L: Pseudostem length of the internode (m)
    :param float internode_Lmax: Final length of the internode (m)
    :param float LSIW: Lineic Structural Internode Weight (g m-1).
    :param float enclosed_mstruct: mstruct of the enclosed leaf (g)


    :return: delta_internode_enclosed_mstruct (g)
    :rtype: float
    """
    if not internode_Lmax:
        enclosed_mstruct_max = internode_L * LSIW
    else:
        enclosed_mstruct_max = min(internode_pseudostem_L, internode_Lmax) * LSIW

    if internode_pseudo_age < parameters.te_IN:
        delta_enclosed_mstruct = (enclosed_mstruct_max - enclosed_mstruct) / (parameters.te_IN - internode_pseudo_age) * delta_internode_pseudo_age
    else:
        delta_enclosed_mstruct = 0

    return max(0., delta_enclosed_mstruct)


def calculate_delta_emerged_tissue_mstruct(SW, previous_mstruct, metric):
    """ delta mstruct of emerged tissue (lamina, sheath and internode). Calculated from tissue area.

    :param float SW: For Lamina : Structural Specific Weight (g m-2); For sheath and internode : Lineic Structural Weight (g m-1)
    :param float previous_mstruct: mstruct at the previous time step i.e. not yet updated (g)
    :param float metric: For Lamina : Area at the current time step, as updated by the geometrical model (m2); For sheath and internode : Length at the current time step (m)

    :return: delta mstruct (g)
    :rtype: float
    """
    updated_mstruct = SW * metric
    delta_mstruct = updated_mstruct - previous_mstruct
    return max(0., delta_mstruct)


def calculate_delta_Nstruct(delta_mstruct):
    """ delta Nstruct of hidden zone and emerged tissue (lamina and sheath).

    :param float delta_mstruct: delta of mstruct (g)

    :return: delta Nstruct (g)
    :rtype: float
    """
    return delta_mstruct * parameters.RATIO_AMINO_ACIDS_MSTRUCT


def calculate_export(delta_mstruct, metabolite, hiddenzone_mstruct):
    """Export of metabolite from the hidden zone towards the emerged part of the leaf integrated over delta_t.

    :param float delta_mstruct: Delta of structural dry mass of the emerged part of the leaf (g)
    :param float metabolite: Metabolite amount in the hidden zone (�mol C or N)
    :param float hiddenzone_mstruct: Structural mass of the hidden zone (g)

    :return: metabolite export (�mol N)
    :rtype: float
    """
    return delta_mstruct * max(0., (metabolite / hiddenzone_mstruct))


def calculate_init_cytokinins_emerged_tissue(delta_mstruct):
    """Initial amount of cytokinins allocated in the mstruct of a newly emerged tissue.

    :param float delta_mstruct: Delta of structural dry mass of the emerged part of the leaf (g)

    :return: cytokinins addition (AU)
    :rtype: float
    """
    return delta_mstruct * parameters.INIT_CYTOKININS_EMERGED_TISSUE  # TODO: Set according to protein concentration ?


def calculate_s_Nstruct_amino_acids(delta_hiddenzone_Nstruct, delta_lamina_Nstruct, delta_sheath_Nstruct, delta_internode_Nstruct):
    """Consumption of amino acids for the calculated mstruct growth (�mol N consumed by mstruct growth)

    :param float delta_hiddenzone_Nstruct: Nstruct growth of the hidden zone (g)
    :param float delta_lamina_Nstruct: Nstruct growth of the lamina (g)
    :param float delta_sheath_Nstruct: Nstruct growth of the sheath (g)
    :param float delta_internode_Nstruct: Nstruct growth of the internode (g)

    :return: Amino acid consumption (�mol N)
    :rtype: float
    """
    return (delta_hiddenzone_Nstruct + delta_lamina_Nstruct + delta_sheath_Nstruct + delta_internode_Nstruct) / parameters.N_MOLAR_MASS * 1E6


def calculate_s_mstruct_sucrose(delta_hiddenzone_mstruct, delta_lamina_mstruct, delta_sheath_mstruct, s_Nstruct_amino_acids_N):
    """Consumption of sucrose for the calculated mstruct growth (�mol C consumed by mstruct growth)

    :param float delta_hiddenzone_mstruct: mstruct growth of the hidden zone (g)
    :param float delta_lamina_mstruct: mstruct growth of the lamina (g)
    :param float delta_sheath_mstruct: mstruct growth of the sheath (g)
    :param float s_Nstruct_amino_acids_N: Total amino acid consumption (�mol N) due to Nstruct (�mol N)

    :return: Sucrose consumption (�mol C)
    :rtype: float
    """
    s_Nstruct_amino_acids = s_Nstruct_amino_acids_N / parameters.AMINO_ACIDS_N_RATIO  #: �mol of AA
    s_mstruct_amino_acids_C = s_Nstruct_amino_acids * parameters.AMINO_ACIDS_C_RATIO  #: �mol of C coming from AA
    s_mstruct_C = (delta_hiddenzone_mstruct + delta_lamina_mstruct + delta_sheath_mstruct) * parameters.RATIO_SUCROSE_MSTRUCT / parameters.C_MOLAR_MASS * 1E6  #: Total C used for mstruct growth (�mol C)
    s_mstruct_sucrose_C = s_mstruct_C - s_mstruct_amino_acids_C  #: �mol of C coming from sucrose

    return s_mstruct_sucrose_C


def calculate_sheath_mstruct(sheath_L, LSSW):
    """ mstruct of the sheath.
      Final mstruct of the enclosed leaf matches sheath mstruct calculation when it is mature.

      :param float sheath_L: Sheath length (m)
      :param float LSSW: Lineic Structural Sheath Weight (g m-1).

      :return: Structural mass of the sheath (g)
      :rtype: float
      """
    return sheath_L * LSSW


def calculate_roots_mstruct_growth(sucrose, amino_acids, mstruct, delta_teq, postflowering_stages):
    """Root structural dry mass growth integrated over delta_t

    :param float sucrose: Amount of sucrose in roots (�mol C)
    :param float amino_acids: Amount of amino acids in roots (�mol N)
    :param float mstruct: Root structural mass (g)
    :param float delta_teq: Time compensated for the effect of temperature - Time equivalent at Tref (s)
    :param bool postflowering_stages: Option : True to run a simulation with postflo parameter

    :return: mstruct_C_growth (�mol C), mstruct_growth (g), Nstruct_growth (g), Nstruct_N_growth (�mol N)
    :rtype: (float, float, float, float)
    """
    conc_sucrose_effective = max(0., sucrose / mstruct - parameters.conc_sucrose_offset)

    if postflowering_stages:
        Vmax = parameters.VMAX_ROOTS_GROWTH_POSTFLO
    else:
        Vmax = parameters.VMAX_ROOTS_GROWTH_PREFLO
    N = parameters.N_ROOTS_GROWTH

    if conc_sucrose_effective > 0.:
        mstruct_C_growth = max(0., ((conc_sucrose_effective ** N) * Vmax) / ((conc_sucrose_effective ** N) + (parameters.K_ROOTS_GROWTH ** N)) * delta_teq * mstruct)  #: root growth in C (�mol of C)
    else:
        mstruct_C_growth = 0.
    mstruct_growth = mstruct_C_growth * parameters.CONVERSION_MMOL_C_G_MSTRUCT_ROOTS  #: root growth (g of structural dry mass)

    Nstruct_growth = mstruct_growth * parameters.RATIO_N_MSTRUCT_ROOTS_  #: root growth in N (g of structural dry mass)
    Nstruct_N_growth = min(amino_acids, (Nstruct_growth / parameters.N_MOLAR_MASS) * 1E6)  #: root growth in nitrogen (�mol N)

    return mstruct_C_growth, mstruct_growth, Nstruct_growth, Nstruct_N_growth


def calculate_roots_s_mstruct_sucrose(delta_roots_mstruct, s_Nstruct_amino_acids_N):
    """Consumption of sucrose for the calculated mstruct growth (�mol C consumed by mstruct growth)

    :param float delta_roots_mstruct: mstruct growth of the roots (g)
    :param float s_Nstruct_amino_acids_N: Total amino acid consumption (�mol N) due to Nstruct (�mol N)

    :return: Sucrose consumption (�mol C)
    :rtype: float
    """
    s_Nstruct_amino_acids = s_Nstruct_amino_acids_N / parameters.AMINO_ACIDS_N_RATIO  #: �mol of AA
    s_mstruct_amino_acids_C = s_Nstruct_amino_acids * parameters.AMINO_ACIDS_C_RATIO  #: �mol of C coming from AA
    s_mstruct_C = delta_roots_mstruct * parameters.RATIO_SUCROSE_MSTRUCT / parameters.C_MOLAR_MASS * 1E6  #: Total C used for mstruct growth (�mol C)
    s_mstruct_sucrose_C = s_mstruct_C - s_mstruct_amino_acids_C  #: �mol of coming from sucrose

    return s_mstruct_sucrose_C


def calculate_mineral_plant(mstruct, senesced_mstruct):
    """ Mineral mass.

    :param float mstruct: structural mass of the plant (g)
    :param float senesced_mstruct: senesced structural mass of the plant (g)

    :return: Mineral mass of the plant (g)
    :rtype: float
    """
    mineral_plant = (mstruct * parameters.MINERAL_LIVING_TISSUE) + (senesced_mstruct * parameters.MINERAL_SENESCED_TISSUE)
    return mineral_plant
