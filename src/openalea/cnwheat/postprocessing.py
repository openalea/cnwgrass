# -*- coding: latin-1 -*-

from __future__ import division  # use "//" to do integer division
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from openalea.cnwheat import simulation as cnwheat_simulation, model as cnwheat_model, parameters as cnwheat_parameters, tools as cnwheat_tools
from openalea.respiwheat import model as respiwheat_model

"""
    cnwheat.postprocessing
    ~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`cnwheat.postprocessing` defines post-processing to apply
    on CN-Wheat outputs, and provides a front-end to automatize the generation of graphs
    for validation of the outputs.

    Please use front-ends :func:`postprocessing` and :func:`generate_graphs`.

    :copyright: Copyright 2014-2017 INRA-ECOSYS, see AUTHORS.
    :license: CeCILL-C, see LICENSE for details.

    **Acknowledgments**: The research leading these results has received funding through the
    Investment for the Future programme managed by the Research National Agency
    (BreedWheat project ANR-10-BTBR-03).

    .. seealso:: Barillot et al. 2016.
"""

#: the time index
T_INDEX = cnwheat_simulation.Simulation.T_INDEX

#: the index to locate the plants in the modeled system
PLANTS_INDEXES = cnwheat_simulation.Simulation.PLANTS_INDEXES
#: concatenation of :attr:`T_INDEX` and :attr:`PLANTS_INDEXES`
PLANTS_T_INDEXES = cnwheat_simulation.Simulation.PLANTS_T_INDEXES
#: plants post-processing variables
PLANTS_POSTPROCESSING_VARIABLES = []
#: concatenation of :attr:`PLANTS_T_INDEXES`, :attr:`PLANTS_RUN_VARIABLES <cnwheat.simulation.Simulation.PLANTS_RUN_VARIABLES>` and :attr:`PLANTS_POSTPROCESSING_VARIABLES`
PLANTS_RUN_POSTPROCESSING_VARIABLES = set(PLANTS_T_INDEXES + cnwheat_simulation.Simulation.PLANTS_RUN_VARIABLES + PLANTS_POSTPROCESSING_VARIABLES)

#: the indexes to locate the axes in the modeled system
AXES_INDEXES = cnwheat_simulation.Simulation.AXES_INDEXES
#: concatenation of :attr:`T_INDEX` and :attr:`AXES_INDEXES`
AXES_T_INDEXES = cnwheat_simulation.Simulation.AXES_T_INDEXES
#: axes post-processing variables
AXES_POSTPROCESSING_VARIABLES = ['C_N_ratio', 'C_N_ratio_shoot', 'N_content', 'N_content_shoot', 'N_content_roots', 'N_content_mstruct', 'N_content_mstruct_shoot', 'N_content_total_DM_shoot',
                                 'N_content_mstruct_roots', 'sum_N_g', 'sum_N_g_shoot',
                                 'sum_dry_mass', 'sum_dry_mass_shoot', 'sum_dry_mass_laminae', 'sum_dry_mass_stem', 'sum_dry_mass_roots', 'dry_mass_phloem', 'shoot_roots_ratio',
                                 'shoot_roots_mstruct_ratio', 'Total_Photosynthesis', 'Tillers_Photosynthesis', 'Tillers_Photosynthesis_An', 'NNI',
                                 'NS', 'NS_shoot', 'NS_stem', 'NS_laminae', 'NS_roots', 'mstruct_shoot', 'mstruct_laminae', 'mstruct_stem',
                                 'C_respired_shoot', 'C_respired_roots', 'Cont_WSC_DM', 'Cont_WSC_DM_shoot', 'Cont_WSC_DM_roots', 'Cont_WSC_DM_laminae', 'Cont_WSC_DM_stem',
                                 'sum_C_g', 'sum_C_g_shoot', 'sum_NSC_g', 'C_exudated']

#: concatenation of :attr:`AXES_T_INDEXES`, :attr:`AXES_RUN_VARIABLES <cnwheat.simulation.Simulation.AXES_RUN_VARIABLES>` and :attr:`AXES_POSTPROCESSING_VARIABLES`
AXES_RUN_POSTPROCESSING_VARIABLES = set(AXES_T_INDEXES + cnwheat_simulation.Simulation.AXES_RUN_VARIABLES + AXES_POSTPROCESSING_VARIABLES)

#: the indexes to locate the phytomers in the modeled system
PHYTOMERS_INDEXES = cnwheat_simulation.Simulation.PHYTOMERS_INDEXES
#: concatenation of :attr:`T_INDEX` and :attr:`PHYTOMERS_INDEXES`
PHYTOMERS_T_INDEXES = cnwheat_simulation.Simulation.PHYTOMERS_T_INDEXES
#: phytomers post-processing variables
PHYTOMERS_POSTPROCESSING_VARIABLES = []
#: concatenation of :attr:`PHYTOMERS_T_INDEXES`, :attr:`PHYTOMERS_RUN_VARIABLES <cnwheat.simulation.Simulation.PHYTOMERS_RUN_VARIABLES>` and :attr:`PHYTOMERS_POSTPROCESSING_VARIABLES`
PHYTOMERS_RUN_POSTPROCESSING_VARIABLES = set(PHYTOMERS_T_INDEXES + cnwheat_simulation.Simulation.PHYTOMERS_RUN_VARIABLES + PHYTOMERS_POSTPROCESSING_VARIABLES)

#: the indexes to locate the organs in the modeled system
ORGANS_INDEXES = cnwheat_simulation.Simulation.ORGANS_INDEXES
#: concatenation of :attr:`T_INDEX` and :attr:`ORGANS_INDEXES`
ORGANS_T_INDEXES = cnwheat_simulation.Simulation.ORGANS_T_INDEXES
#: organs post-processing variables
ORGANS_POSTPROCESSING_VARIABLES = ['Conc_Amino_Acids', 'Conc_Nitrates', 'Conc_Sucrose', 'Conc_cytokinins', 'Dry_Mass', 'Proteins_N_Mass',  'N_tot', 'WSC_g']
ORGANS_RUN_VARIABLES_ADDITIONAL = ['sucrose_consumption_mstruct', 'AA_consumption_mstruct', 'synthetized_mstruct']
#: concatenation of :attr:`ORGANS_T_INDEXES`, :attr:`ORGANS_RUN_VARIABLES <cnwheat.simulation.Simulation.ORGANS_RUN_VARIABLES>` and :attr:`ORGANS_POSTPROCESSING_VARIABLES`
ORGANS_RUN_POSTPROCESSING_VARIABLES = set(ORGANS_T_INDEXES + cnwheat_simulation.Simulation.ORGANS_RUN_VARIABLES + ORGANS_POSTPROCESSING_VARIABLES + ORGANS_RUN_VARIABLES_ADDITIONAL)

#: the indexes to locate the hidden zones in the modeled system
HIDDENZONE_INDEXES = cnwheat_simulation.Simulation.HIDDENZONE_INDEXES
#: concatenation of :attr:`T_INDEX` and :attr:`HIDDENZONE_INDEXES`
HIDDENZONE_T_INDEXES = cnwheat_simulation.Simulation.HIDDENZONE_T_INDEXES
#: hidden zones post-processing variables
HIDDENZONE_POSTPROCESSING_VARIABLES = ['Conc_Amino_Acids', 'Conc_Fructan', 'Conc_Proteins', 'Conc_Sucrose', 'RER', 'nb_replications', 'Cont_Fructan_DM', 'Cont_Proteins_DM', 'sum_dry_mass',
                                       'N_content', 'Cont_WSC_DM']
HIDDENZONE_RUN_VARIABLES_ADDITIONAL = ['leaf_L', 'delta_leaf_L', 'internode_L', 'leaf_pseudostem_length', 'leaf_is_emerged', 'Respi_growth', 'sucrose_consumption_mstruct', 'AA_consumption_mstruct']
#: concatenation of :attr:`HIDDENZONE_T_INDEXES`, :attr:`HIDDENZONE_RUN_VARIABLES <cnwheat.simulation.Simulation.HIDDENZONE_RUN_VARIABLES>` and :attr:`HIDDENZONE_POSTPROCESSING_VARIABLES`
HIDDENZONE_RUN_POSTPROCESSING_VARIABLES = set(HIDDENZONE_T_INDEXES + cnwheat_simulation.Simulation.HIDDENZONE_RUN_VARIABLES + HIDDENZONE_RUN_VARIABLES_ADDITIONAL + HIDDENZONE_POSTPROCESSING_VARIABLES)

#: the indexes to locate the elements in the modeled system
ELEMENTS_INDEXES = cnwheat_simulation.Simulation.ELEMENTS_INDEXES
#: concatenation of :attr:`T_INDEX` and :attr:`ELEMENTS_INDEXES`
ELEMENTS_T_INDEXES = cnwheat_simulation.Simulation.ELEMENTS_T_INDEXES
#: elements post-processing variables
ELEMENTS_POSTPROCESSING_VARIABLES = ['Conc_Amino_Acids', 'Conc_Fructan', 'Conc_Nitrates', 'Conc_Proteins', 'Conc_Starch', 'Conc_Sucrose', 'Conc_TriosesP', 'Cont_Fructan_DM',
                                     'Conc_cytokinins', 'Surfacic_NS', 'NS', 'N_content', 'N_content_total_DM', 'N_tot', 'nb_replications', 'SLA', 'SLN',
                                     'SLN_nonstruct', 'sum_dry_mass', 'Photosynthetic_efficiency', 'Cont_WSC_DM']
ELEMENTS_RUN_VARIABLES_ADDITIONAL = ['length', 'PARa']
#: concatenation of :attr:`ELEMENTS_T_INDEXES`, :attr:`ELEMENTS_RUN_VARIABLES <cnwheat.simulation.Simulation.ELEMENTS_RUN_VARIABLES>` and :attr:`ELEMENTS_POSTPROCESSING_VARIABLES`
ELEMENTS_RUN_POSTPROCESSING_VARIABLES = set(ELEMENTS_T_INDEXES + cnwheat_simulation.Simulation.ELEMENTS_RUN_VARIABLES + ELEMENTS_RUN_VARIABLES_ADDITIONAL + ELEMENTS_POSTPROCESSING_VARIABLES)

#: the indexes to locate the soils in the modeled system
SOILS_INDEXES = cnwheat_simulation.Simulation.SOILS_INDEXES
#: concatenation of :attr:`T_INDEX` and :attr:`SOILS_INDEXES`
SOILS_T_INDEXES = cnwheat_simulation.Simulation.SOILS_T_INDEXES
#: soils post-processing variables
SOILS_POSTPROCESSING_VARIABLES = []
#: concatenation of :attr:`SOILS_T_INDEXES`, :attr:`SOILS_RUN_VARIABLES <cnwheat.simulation.Simulation.SOILS_RUN_VARIABLES>` and :attr:`SOILS_POSTPROCESSING_VARIABLES`
SOILS_RUN_POSTPROCESSING_VARIABLES = SOILS_T_INDEXES + cnwheat_simulation.Simulation.SOILS_RUN_VARIABLES + SOILS_POSTPROCESSING_VARIABLES


# -------------------------------------------------
# ----------- POST-PROCESSING FUNCTIONS -----------
# ----------- DO NOT USE THEM DIRECTLY ------------
# -------------------------------------------------

class Roots:
    """
    Post-processing to apply on Roots outputs.
    """

    def __init__(self):
        pass

    @staticmethod
    def calculate_WSC_g(sucrose):
        """Mass of Water Soluble Carbohydrates

        :param float sucrose: Amount of sucrose (�mol` C)

        :return: Water Soluble Carbohydrates (g)
        :rtype: float
        """
        C_MOLAR_MASS = cnwheat_model.EcophysiologicalConstants.C_MOLAR_MASS

        WSC_g = (sucrose * 1E-6 * C_MOLAR_MASS) / cnwheat_model.EcophysiologicalConstants.HEXOSE_MOLAR_MASS_C_RATIO

        return WSC_g

    @staticmethod
    def calculate_Conc_Nitrates(nitrates, mstruct):
        """Nitrate concentration.

        :param float nitrates: Amount of nitrates (�mol` N)
        :param float mstruct: Structural mass (g)

        :return: Nitrate concentration (�mol` nitrates g-1 mstruct)
        :rtype: float
        """
        return nitrates / mstruct

    @staticmethod
    def calculate_Conc_Amino_Acids(amino_acids, mstruct):
        """Amino_acid concentration.

        :param float amino_acids: Amount of amino acids (�mol` N)
        :param float mstruct: Structural mass (g)

        :return: Amino_acid concentration (�mol` amino_acids g-1 mstruct)
        :rtype: float
        """
        return (amino_acids / cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_N_RATIO) / mstruct

    @staticmethod
    def calculate_conc_sucrose(sucrose, mstruct):
        """Sucrose concentration.

        :param float sucrose: Amount of sucrose (�mol` C)
        :param float mstruct: Structural mass (g)

        :return: Sucrose concentration (�mol` sucrose g-1 mstruct)
        :rtype: float
        """
        return (sucrose / mstruct) / cnwheat_model.EcophysiologicalConstants.NB_C_SUCROSE

    @staticmethod
    def calculate_conc_cytokinins(cytokinins, mstruct):
        """Cytokinin concentration.

        :param float cytokinins: Amount of cytokinins (AU)
        :param float mstruct: Structural mass (g)

        :return: cytokinins concentration (AU cytokinins g-1 mstruct)
        :rtype: float
        """
        return cytokinins / mstruct


class Phloem:
    """
    Post-processing to apply on Phloem outputs.
    """

    def __init__(self):
        pass

    @staticmethod
    def calculate_WSC_g(sucrose):
        """Mass of Water Soluble Carbohydrates

        :param float sucrose: Amount of sucrose (�mol` C)

        :return: Water Soluble Carbohydrates (g)
        :rtype: float
        """
        C_MOLAR_MASS = cnwheat_model.EcophysiologicalConstants.C_MOLAR_MASS

        WSC_g = (sucrose * 1E-6 * C_MOLAR_MASS) / cnwheat_model.EcophysiologicalConstants.HEXOSE_MOLAR_MASS_C_RATIO

        return WSC_g

    @staticmethod
    def calculate_conc_amino_acids(amino_acids, mstruct_axis):
        """Amino_acids concentration. Related to the structural dry mass of the culm.

        :param float amino_acids: Amount of amino_acids in phloem (�mol` N)
        :param float mstruct_axis: The structural dry mass of the axis (g)

        :return: Amino_acids concentration (�mol` amino_acids g-1 mstruct)
        :rtype: float
        """
        return (amino_acids / cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_N_RATIO) / mstruct_axis

    @staticmethod
    def calculate_conc_sucrose(sucrose, mstruct_axis):
        """Sucrose concentration. Related to the structural dry mass of the culm

        :param float sucrose: Amount of sucrose in phloem (�mol` C)
        :param float mstruct_axis: The structural dry mass of the axis (g)

        :return: Sucrose concentration (�mol` sucrose g-1 mstruct)
        :rtype: float
        """
        return (sucrose / cnwheat_model.EcophysiologicalConstants.NB_C_SUCROSE) / mstruct_axis


class Grains:
    """
    Post-processing to apply on Grains outputs.
    """

    def __init__(self):
        pass

    @staticmethod
    def calculate_dry_mass(structure, starch, proteins):
        """Grain total dry mass.

        :param float structure: Grain structural C mass (�mol` C)
        :param float starch:  Grain starch content (�mol` C)
        :param float proteins: Grain protein content (�mol` N)

        :return: Grain total dry mass (g)
        :rtype: float
        """
        #: Carbohydrates mass, grain carbohydrates supposed to be mainly starch i.e. glucose polymers (C6 H12 O6)
        C_mass = ((structure + starch) * 1E-6 * cnwheat_model.EcophysiologicalConstants.C_MOLAR_MASS) / cnwheat_model.EcophysiologicalConstants.HEXOSE_MOLAR_MASS_C_RATIO

        #: N mass, grain proteins were supposed to be gluten mainly composed of Glu, Gln and Pro
        N_mass = (proteins * 1E-6 * cnwheat_model.EcophysiologicalConstants.N_MOLAR_MASS) / cnwheat_model.Grains.AMINO_ACIDS_MOLAR_MASS_N_RATIO

        return C_mass + N_mass

    @staticmethod
    def calculate_WSC_g(sucrose, starch):
        """ Mass of Water Soluble Carbohydrates

        :param float sucrose: Amount of sucrose (�mol` C)
        :param float starch: Amount of sucrose (�mol` C)

        :return: Water Soluble Carbohydrates (g)
        :rtype: float
        """
        C_MOLAR_MASS = cnwheat_model.EcophysiologicalConstants.C_MOLAR_MASS
        WSC_g = ((sucrose * 1E-6 * C_MOLAR_MASS) / cnwheat_model.EcophysiologicalConstants.HEXOSE_MOLAR_MASS_C_RATIO +
                 (starch * 1E-6 * C_MOLAR_MASS) / cnwheat_model.EcophysiologicalConstants.HEXOSE_MOLAR_MASS_C_RATIO)

        return WSC_g

    @staticmethod
    def calculate_protein_N_mass(proteins):
        """Grain total protein mass.

        :param float proteins: Grain protein content (�mol` N)

        :return: Grain total protein mass (g)
        :rtype: float
        """
        mass_N_proteins = proteins * 1E-6 * cnwheat_model.EcophysiologicalConstants.N_MOLAR_MASS  #: Mass of nitrogen in proteins (g)
        # masS_proteins = mass_N_proteins / EcophysiologicalConstants.AMINO_ACIDS_MOLAR_MASS_N_RATIO     #: Total mass of proteins (g)
        return mass_N_proteins


class HiddenZone:
    """
    Post-processing to apply on HiddenZone outputs.
    """

    def __init__(self):
        pass

    @staticmethod
    def calculate_dry_mass(sucrose, starch, fructan, amino_acids, proteins, mstruct):
        """Dry mass

        :param float sucrose: Amount of sucrose (�mol` C)
        :param float starch: Amount of starch (�mol` C)
        :param float fructan: Amount of fructan (�mol` C)
        :param float amino_acids: Amount of amino acids (�mol` N)
        :param float proteins: Amount of proteins (�mol` N)
        :param float mstruct: strcural mass (g)

        :return: Dry mass (g)
        :rtype: float
        """
        C_MOLAR_MASS = cnwheat_model.EcophysiologicalConstants.C_MOLAR_MASS
        N_MOLAR_MASS = cnwheat_model.EcophysiologicalConstants.N_MOLAR_MASS
        dry_mass = ((sucrose * 1E-6 * C_MOLAR_MASS) / cnwheat_model.EcophysiologicalConstants.HEXOSE_MOLAR_MASS_C_RATIO +
                    (starch * 1E-6 * C_MOLAR_MASS) / cnwheat_model.EcophysiologicalConstants.HEXOSE_MOLAR_MASS_C_RATIO +
                    (fructan * 1E-6 * C_MOLAR_MASS) / cnwheat_model.EcophysiologicalConstants.HEXOSE_MOLAR_MASS_C_RATIO +
                    (amino_acids * 1E-6 * N_MOLAR_MASS) / cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                    (proteins * 1E-6 * N_MOLAR_MASS) / cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                    mstruct)

        return dry_mass

    @staticmethod
    def calculate_WSC_g(sucrose, fructan):
        """Mass of Water Soluble Carbohydrates

        :param float sucrose: Amount of sucrose (�mol` C)
        :param float fructan: Amount of fructan (�mol` C)

        :return: Water Soluble Carbohydrates (g)
        :rtype: float
        """
        C_MOLAR_MASS = cnwheat_model.EcophysiologicalConstants.C_MOLAR_MASS

        WSC_g = ((sucrose * 1E-6 * C_MOLAR_MASS) / cnwheat_model.EcophysiologicalConstants.HEXOSE_MOLAR_MASS_C_RATIO +
                 (fructan * 1E-6 * C_MOLAR_MASS) / cnwheat_model.EcophysiologicalConstants.HEXOSE_MOLAR_MASS_C_RATIO
                 )

        return WSC_g

    @staticmethod
    def calculate_C_g(sucrose, starch, fructan, amino_acids, proteins, mstruct):
        """Mass of carbon metabolites

        :param float sucrose: Amount of sucrose (�mol` C)
        :param float starch: Amount of starch (�mol` C)
        :param float fructan: Amount of fructan (�mol` C)
        :param float amino_acids: Amount of amino acids (�mol` N)
        :param float proteins: Amount of proteins (�mol` N)
        :param float mstruct: strcural mass (g)

        :return: Dry mass (g)
        :rtype: float
        """
        C_MOLAR_MASS = cnwheat_model.EcophysiologicalConstants.C_MOLAR_MASS
        N_MOLAR_MASS = cnwheat_model.EcophysiologicalConstants.N_MOLAR_MASS
        C_mass = ((sucrose * 1E-6 * C_MOLAR_MASS) +
                  (starch * 1E-6 * C_MOLAR_MASS) +
                  (fructan * 1E-6 * C_MOLAR_MASS) +
                  (
                              amino_acids * 1E-6 * N_MOLAR_MASS) * cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_MOLAR_MASS_C_RATIO / cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                  (proteins * 1E-6 * N_MOLAR_MASS) * cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_MOLAR_MASS_C_RATIO / cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                  mstruct * cnwheat_model.EcophysiologicalConstants.RATIO_C_mstruct)

        return C_mass

    @staticmethod
    def calculate_N_g(amino_acids, proteins, Nstruct):
        """Mass of N metabolites

        :param float amino_acids: Amount of amino acids (�mol` N)
        :param float proteins: Amount of proteins (�mol` N)
        :param float Nstruct: N structural mass(g)

        :return: Dry mass (g)
        :rtype: float
        """
        N_MOLAR_MASS = cnwheat_model.EcophysiologicalConstants.N_MOLAR_MASS
        N_mass = ((amino_acids * 1E-6 * N_MOLAR_MASS) +
                  (proteins * 1E-6 * N_MOLAR_MASS) +
                  Nstruct)

        return N_mass

    @staticmethod
    def calculate_fructan_g(fructan):
        """Mass of fructans

        :param float fructan: Amount of fructans (�mol` C)

        :return: Dry mass (g)
        :rtype: float
        """
        C_MOLAR_MASS = cnwheat_model.EcophysiologicalConstants.C_MOLAR_MASS
        return (fructan * 1E-6 * C_MOLAR_MASS) / cnwheat_model.EcophysiologicalConstants.HEXOSE_MOLAR_MASS_C_RATIO

    @staticmethod
    def calculate_proteins_g(proteins):
        """Mass of proteins

        :param float proteins: Amount of proteins (�mol` N)

        :return: Dry mass (g)
        :rtype: float
        """
        N_MOLAR_MASS = cnwheat_model.EcophysiologicalConstants.N_MOLAR_MASS
        return (proteins * 1E-6 * N_MOLAR_MASS) / cnwheat_model.EcophysiologicalConstants.PROTEINS_MOLAR_MASS_N_RATIO

    @staticmethod
    def calculate_Conc_Amino_Acids(amino_acids, mstruct):
        """Amino acid concentration.

        :param float amino_acids: N amino acids (�mol` N)
        :param float mstruct: Structural mass

        :return: Amino_acid concentration (�mol` amino acids g-1 mstruct)
        :rtype: float
        """
        return (amino_acids / cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_N_RATIO) / mstruct

    @staticmethod
    def calculate_conc_sucrose(sucrose, mstruct):
        """Sucrose concentration.

        :param float sucrose: C sucrose (�mol` C)
        :param float mstruct: Structural mass

        :return: Sucrose concentration (�mol` sucrose g-1 mstruct)
        :rtype: float
        """
        return (sucrose / mstruct) / cnwheat_model.EcophysiologicalConstants.NB_C_SUCROSE

    @staticmethod
    def calculate_conc_fructan(fructan, mstruct):
        """Fructan concentration.

        :param float fructan: C fructan (�mol` C)
        :param float mstruct: Structural mass

        :return: Fructan concentration (�mol` fructan g-1 mstruct, eq. glucose).
        :rtype: float
        """
        return (fructan / mstruct) / cnwheat_model.EcophysiologicalConstants.NB_C_HEXOSES

    @staticmethod
    def calculate_conc_protein(proteins, mstruct):
        """Proteins concentration.

        :param float proteins: N proteins (�mol` N)
        :param float mstruct: Structural mass

        :return: Protein concentration (g proteins g-1 mstruct)
        :rtype: float
        """
        mass_N_proteins = proteins * 1E-6 * cnwheat_model.EcophysiologicalConstants.N_MOLAR_MASS  #: Mass of N in proteins (g)
        masS_proteins = mass_N_proteins / cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_MOLAR_MASS_N_RATIO  #: Total mass of proteins (g)
        return masS_proteins / mstruct

    @staticmethod
    def calculate_RER(delta_leaf_L, leaf_L, delta_t):
        """Relative Extension Rate.

        :param float delta_leaf_L: delta of leaf length between t and t-1 (m)
        :param float leaf_L: leaf length (m)
        :param float delta_t: delta_t (s)

        :return: Relative Extension Rate (s-1)
        :rtype: float
        """
        return (delta_leaf_L / delta_t) / leaf_L


class Element:
    """
    Post-processing to apply on Element outputs.
    """

    def __init__(self):
        pass

    @staticmethod
    def calculate_dry_mass(triosesP, sucrose, starch, fructan, nitrates, amino_acids, proteins, mstruct):
        """Dry mass

        :param triosesP: Amount of triose phosphates (�mol` C)
        :param float sucrose: Amount of sucrose (�mol` C)
        :param float starch: Amount of sucrose (�mol` C)
        :param float fructan: Amount of sucrose (�mol` C)
        :param float nitrates: Amount of nitrates (�mol` N)
        :param float amino_acids: Amount of sucrose (�mol` N)
        :param float proteins: Amount of sucrose (�mol` N)
        :param float mstruct: strcural mass (g)

        :return: Dry mass (g)
        :rtype: float
        """
        C_MOLAR_MASS = cnwheat_model.EcophysiologicalConstants.C_MOLAR_MASS
        N_MOLAR_MASS = cnwheat_model.EcophysiologicalConstants.N_MOLAR_MASS
        dry_mass = ((triosesP * 1E-6 * C_MOLAR_MASS) / cnwheat_model.EcophysiologicalConstants.TRIOSESP_MOLAR_MASS_C_RATIO +
                    (sucrose * 1E-6 * C_MOLAR_MASS) / cnwheat_model.EcophysiologicalConstants.HEXOSE_MOLAR_MASS_C_RATIO +
                    (starch * 1E-6 * C_MOLAR_MASS) / cnwheat_model.EcophysiologicalConstants.HEXOSE_MOLAR_MASS_C_RATIO +
                    (fructan * 1E-6 * C_MOLAR_MASS) / cnwheat_model.EcophysiologicalConstants.HEXOSE_MOLAR_MASS_C_RATIO +
                    (nitrates * 1E-6 * N_MOLAR_MASS) / cnwheat_model.EcophysiologicalConstants.NITRATES_MOLAR_MASS_N_RATIO +
                    (amino_acids * 1E-6 * N_MOLAR_MASS) / cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                    (proteins * 1E-6 * N_MOLAR_MASS) / cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                    mstruct)

        return dry_mass

    @staticmethod
    def calculate_WSC_g(triosesP, sucrose, starch, fructan):
        """ Mass of Water Soluble Carbohydrates

        :param triosesP: Amount of triose phosphates (�mol` C)
        :param float sucrose: Amount of sucrose (�mol` C)
        :param float starch: Amount of sucrose (�mol` C)
        :param float fructan: Amount of sucrose (�mol` C)

        :return: Water Soluble Carbohydrates (g)
        :rtype: float
        """
        C_MOLAR_MASS = cnwheat_model.EcophysiologicalConstants.C_MOLAR_MASS
        WSC_g = ((triosesP * 1E-6 * C_MOLAR_MASS) / cnwheat_model.EcophysiologicalConstants.TRIOSESP_MOLAR_MASS_C_RATIO +
                 (sucrose * 1E-6 * C_MOLAR_MASS) / cnwheat_model.EcophysiologicalConstants.HEXOSE_MOLAR_MASS_C_RATIO +
                 (starch * 1E-6 * C_MOLAR_MASS) / cnwheat_model.EcophysiologicalConstants.HEXOSE_MOLAR_MASS_C_RATIO +
                 (fructan * 1E-6 * C_MOLAR_MASS) / cnwheat_model.EcophysiologicalConstants.HEXOSE_MOLAR_MASS_C_RATIO
                 )

        return WSC_g

    @staticmethod
    def calculate_fructan_g(fructan):
        """Mass of fructans

        :param float fructan: Amount of fructans (�mol` C)

        :return: Dry mass (g)
        :rtype: float
        """
        C_MOLAR_MASS = cnwheat_model.EcophysiologicalConstants.C_MOLAR_MASS
        return (fructan * 1E-6 * C_MOLAR_MASS) / cnwheat_model.EcophysiologicalConstants.HEXOSE_MOLAR_MASS_C_RATIO

    @staticmethod
    def calculate_C_g(triosesP, sucrose, starch, fructan, amino_acids, proteins, mstruct):
        """Mass of carbon metabolites

        :param float triosesP: Amount of triose phosphates (�mol` C)
        :param float sucrose: Amount of sucrose (�mol` C)
        :param float starch: Amount of starch (�mol` C)
        :param float fructan: Amount of fructan (�mol` C)
        :param float amino_acids: Amount of amino acids (�mol` N)
        :param float proteins: Amount of proteins (�mol` N)
        :param float mstruct: strcural mass (g)

        :return: Dry mass (g)
        :rtype: float
        """
        C_MOLAR_MASS = cnwheat_model.EcophysiologicalConstants.C_MOLAR_MASS
        N_MOLAR_MASS = cnwheat_model.EcophysiologicalConstants.N_MOLAR_MASS
        C_mass = ((triosesP * 1E-6 * C_MOLAR_MASS) +
                  (sucrose * 1E-6 * C_MOLAR_MASS) +
                  (starch * 1E-6 * C_MOLAR_MASS) +
                  (fructan * 1E-6 * C_MOLAR_MASS) +
                  (
                              amino_acids * 1E-6 * N_MOLAR_MASS) * cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_MOLAR_MASS_C_RATIO / cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                  (proteins * 1E-6 * N_MOLAR_MASS) * cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_MOLAR_MASS_C_RATIO / cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                  mstruct * cnwheat_model.EcophysiologicalConstants.RATIO_C_mstruct)

        return C_mass

    @staticmethod
    def calculate_N_g(nitrates, amino_acids, proteins, Nstruct):
        """Mass of N metabolites

        :param float nitrates: Amount of nitrates (�mol` N)
        :param float amino_acids: Amount of amino acids (�mol` N)
        :param float proteins: Amount of proteins (�mol` N)
        :param float Nstruct: N structural mass(g)

        :return: Dry mass (g)
        :rtype: float
        """
        N_MOLAR_MASS = cnwheat_model.EcophysiologicalConstants.N_MOLAR_MASS
        N_mass = ((nitrates * 1E-6 * N_MOLAR_MASS) +
                  (amino_acids * 1E-6 * N_MOLAR_MASS) +
                  (proteins * 1E-6 * N_MOLAR_MASS) +
                  Nstruct)

        return N_mass

    @staticmethod
    def calculate_conc_triosesP(triosesP, mstruct):
        """Triose Phosphates concentration.

        :param float triosesP: Amount of triose phosphates (�mol` C)
        :param float mstruct: Structural mass

        :return: Triose phosphates concentration (�mol` triosesP g-1 mstruct)
        :rtype: float
        """
        return (triosesP / mstruct) / cnwheat_model.EcophysiologicalConstants.NB_C_TRIOSEP

    @staticmethod
    def calculate_conc_sucrose(sucrose, mstruct):
        """Sucrose concentration.

        :param float sucrose: Amount of sucrose (�mol` C)
        :param float mstruct: Structural mass

        :return: Sucrose concentration (�mol` sucrose g-1 mstruct)
        :rtype: float
        """
        return (sucrose / mstruct) / cnwheat_model.EcophysiologicalConstants.NB_C_SUCROSE

    @staticmethod
    def calculate_conc_starch(starch, mstruct):
        """Starch concentration.

        :param float starch: Amount of sucrose (�mol` C)
        :param float mstruct: Structural mass

        :return: Starch concentration (�mol` starch g-1 mstruct)
        :rtype: float
        """
        return (starch / mstruct) / cnwheat_model.EcophysiologicalConstants.NB_C_HEXOSES

    @staticmethod
    def calculate_conc_fructan(fructan, mstruct):
        """Fructan concentration.

        :param float fructan: Amount of fructan (�mol` C)
        :param float mstruct: Structural mass

        :return: Fructan concentration (�mol` fructan g-1 mstruct, eq. glucose).
        :rtype: float
        """
        return (fructan / mstruct) / cnwheat_model.EcophysiologicalConstants.NB_C_HEXOSES

    @staticmethod
    def calculate_Conc_Nitrates(nitrates, mstruct):
        """Nitrate concentration.

        :param float nitrates: Amount of nitrates (�mol` N)
        :param float mstruct: Structural mass

        :return: Nitrate concentration (�mol` nitrates g-1 mstruct)
        :rtype: float
        """
        return nitrates / mstruct

    @staticmethod
    def calculate_Conc_Amino_Acids(amino_acids, mstruct):
        """Amino_acid concentration.

        :param float amino_acids: Amount of amino acids (�mol` N)
        :param float mstruct: Structural mass

        :return: Amino_acid concentration (�mol` amino acids g-1 mstruct)
        :rtype: float
        """
        return (amino_acids / cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_N_RATIO) / mstruct

    @staticmethod
    def calculate_conc_proteins(proteins, mstruct):
        """Protein concentration.

        :param float proteins: Amount of proteins (�mol` N)
        :param float mstruct: Structural mass

        :return: Protein concentration (g proteins g-1 mstruct)
        :rtype: float
        """
        mass_N_proteins = proteins * 1E-6 * cnwheat_model.EcophysiologicalConstants.N_MOLAR_MASS  #: Mass of N in proteins (g)
        mass_proteins = mass_N_proteins / cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_MOLAR_MASS_N_RATIO  #: Total mass of proteins (g)
        return mass_proteins / mstruct

    @staticmethod
    def calculate_conc_cytokinins(cytokinins, mstruct):
        """Cytokinin concentration.

        :param float cytokinins: Amount of cytokinins (AU)
        :param float mstruct: Structural mass

        :return: Cytokinin concentration (AU g-1 mstruct)
        :rtype: float
        """
        return cytokinins / mstruct

    @staticmethod
    def calculate_SLN(nitrates, amino_acids, proteins, Nstruct, green_area):
        """ Surfacic Leaf Nitrogen (SLN, g.m-2)

        :param float nitrates: Amount of nitrates (�mol` N)
        :param float amino_acids: Amount of amino_acids (�mol` N)
        :param float proteins: Amount of proteins (�mol` N)
        :param float Nstruct: Structural N (g)
        :param float green_area: Green area (m-2)

        :return: Surfacic Leaf Nitrogen (SLN, g.m-2)
        :rtype: float
        """
        mass_N_tot = (nitrates + amino_acids + proteins) * 1E-6 * cnwheat_model.EcophysiologicalConstants.N_MOLAR_MASS + Nstruct
        return mass_N_tot / green_area

    @staticmethod
    def calculate_SLN_nonstruct(nitrates, amino_acids, proteins, green_area):
        """ Surfacic Leaf Nitrogen (SLN, g.m-2)

        :param float nitrates: Amount of nitrates (�mol` N)
        :param float amino_acids: Amount of amino_acids (�mol` N)
        :param float proteins: Amount of proteins (�mol` N)
        :param float green_area: Green area (m-2)

        :return: Surfacic Leaf Nitrogen (SLN, g.m-2)
        :rtype: float
        """
        mass_N_tot = (nitrates + amino_acids + proteins) * 1E-6 * cnwheat_model.EcophysiologicalConstants.N_MOLAR_MASS
        return mass_N_tot / green_area

    @staticmethod
    def calculate_SLA(dry_mass, green_area):
        """Specific Leaf Area (SLA, m2.kg-1)

        :param float dry_mass: Dry mass (g)
        :param float green_area: Green area (m-2)

        :return: Specific Leaf Area (SLA, m2 kg-1)
        :rtype: float
        """
        return green_area / (dry_mass * 10 ** -3)

    @staticmethod
    def calculate_surfacic_non_structural(dry_mass, mstruct, green_area):
        """Surfacic content of non structural mass

        :param float dry_mass: Dry mass (g)
        :param float mstruct: Structural mass (g)
        :param float green_area: Green area (m-2)

        :return: Surfacic non structural mass (g m-2)
        :rtype: float
        """
        return (dry_mass - mstruct) / green_area

    @staticmethod
    def calculate_ratio_non_structural(dry_mass, mstruct):
        """Surfacic content of nitrogen

        :param float dry_mass: Dry mass (g)
        :param float mstruct: Structural mass (g)

        :return: Surfacic non structural mass (g m-2)
        :rtype: float
        """
        return (1 - mstruct / dry_mass) * 100


# -----------------------------------------------------------------------------
# --------------------- POST-PROCESSING FRONT-END -----------------------------
# PLEASE USE THIS FUNCTION TO APPLY POST-PROCESSING ON THE OUTPUT OF CN-WHEAT -
# -----------------------------------------------------------------------------


def postprocessing(plants_df=None, axes_df=None, metamers_df=None, hiddenzones_df=None, organs_df=None, elements_df=None, soils_df=None, delta_t=1):
    """
    Compute post-processing from CN-Wheat outputs, and format the post-processing to :class:`dataframes <pandas.DataFrame>`.

    For each post-processing output dataframe:

        * compute post-processing from CN-Wheat outputs,
        * concatenate CN-Wheat outputs and post-processing and place the results in a jointed dataframe,
        * reorder the columns of the dataframes according to :attr:`PLANTS_RUN_POSTPROCESSING_VARIABLES`,
          :attr:`AXES_RUN_POSTPROCESSING_VARIABLES`, :attr:`PHYTOMERS_RUN_POSTPROCESSING_VARIABLES`,
          :attr:`ORGANS_RUN_POSTPROCESSING_VARIABLES`, :attr:`HIDDENZONE_RUN_POSTPROCESSING_VARIABLES`,
          :attr:`ELEMENTS_RUN_POSTPROCESSING_VARIABLES` and :attr:`SOILS_RUN_POSTPROCESSING_VARIABLES`,
        * and convert the indexes of plants and metamers to integers (if relevant).

    :param pandas.DataFrame plants_df: CN-Wheat outputs at plant scale (see :attr:`simulation.Simulation.PLANTS_RUN_VARIABLES`)
    :param pandas.DataFrame axes_df: CN-Wheat outputs at axis scale (see :attr:`simulation.Simulation.AXES_RUN_VARIABLES`)
    :param pandas.DataFrame metamers_df: CN-Wheat outputs at phytomer scale (see :attr:`simulation.Simulation.PHYTOMERS_RUN_VARIABLES`)
    :param pandas.DataFrame hiddenzones_df: CN-Wheat outputs at hidden zone scale (see :attr:`simulation.Simulation.HIDDENZONE_RUN_VARIABLES`)
    :param pandas.DataFrame organs_df: CN-Wheat outputs at organ scale (see :attr:`simulation.Simulation.ORGANS_RUN_VARIABLES`)
    :param pandas.DataFrame elements_df: CN-Wheat outputs at element scale (see :attr:`simulation.Simulation.ELEMENTS_RUN_VARIABLES`)
    :param pandas.DataFrame soils_df: CN-Wheat outputs at soil scale (see :attr:`simulation.Simulation.SOILS_RUN_VARIABLES`)
    :param float delta_t: Delta t between 2 outputs (in seconds).

    :return: post-processing for each scale:
            * plant (see :attr:`PLANTS_RUN_POSTPROCESSING_VARIABLES`)
            * axis (see :attr:`AXES_RUN_POSTPROCESSING_VARIABLES`)
            * metamer (see :attr:`PHYTOMERS_RUN_POSTPROCESSING_VARIABLES`)
            * organ (see :attr:`ORGANS_RUN_POSTPROCESSING_VARIABLES`)
            * hidden zone (see :attr:`HIDDENZONE_RUN_POSTPROCESSING_VARIABLES`)
            * element (see :attr:`ELEMENTS_RUN_POSTPROCESSING_VARIABLES`)
            * and soil (see :attr:`SOILS_RUN_POSTPROCESSING_VARIABLES`)
        depending of the dataframes given as argument.
        For example, if user passes only dataframes `plants_df`, `axes_df` and `metamers_df`,
        then only post-processing dataframes of plants, axes and metamers are returned.
    :rtype: tuple [pandas.DataFrame]
    """

    returned_dataframes = []

    # plants
    if plants_df is not None:
        pp_plants_df = pd.concat([plants_df, pd.DataFrame(columns=PLANTS_POSTPROCESSING_VARIABLES)], sort=False)
        pp_plants_df = pp_plants_df.reindex(PLANTS_RUN_POSTPROCESSING_VARIABLES, axis=1, copy=False)
        pp_plants_df['plant'] = pp_plants_df['plant'].astype(int)
        returned_dataframes.append(pp_plants_df)
    else:
        returned_dataframes.append(pd.DataFrame({'A': []}))

    # metamers
    if metamers_df is not None:
        pp_metamers_df = pd.concat([metamers_df, pd.DataFrame(columns=PHYTOMERS_POSTPROCESSING_VARIABLES)], sort=False)
        pp_metamers_df = pp_metamers_df.reindex(PHYTOMERS_RUN_POSTPROCESSING_VARIABLES, axis=1, copy=False)
        pp_metamers_df[['plant', 'metamer']] = pp_metamers_df[['plant', 'metamer']].astype(int)
        returned_dataframes.append(pp_metamers_df)
    else:
        returned_dataframes.append(pd.DataFrame({'A': []}))

    # organs
    if organs_df is not None and axes_df is not None:
        axes_df = axes_df[axes_df['axis'] == 'MS'].copy()  # TODO : Temporary !

        pp_organs_df = pd.concat([organs_df, pd.DataFrame(columns=ORGANS_POSTPROCESSING_VARIABLES)], sort=False)

        organs_df['sum_dry_mass'] = (((organs_df.fillna(0)['structure'] + organs_df.fillna(0)[
            'starch']) * 1E-6 * cnwheat_model.EcophysiologicalConstants.C_MOLAR_MASS / cnwheat_model.EcophysiologicalConstants.HEXOSE_MOLAR_MASS_C_RATIO) +
                                     (organs_df.fillna(0)[
                                          'sucrose'] * 1E-6 * cnwheat_model.EcophysiologicalConstants.C_MOLAR_MASS) / cnwheat_model.EcophysiologicalConstants.HEXOSE_MOLAR_MASS_C_RATIO +
                                     (organs_df.fillna(0)['starch'] * 1E-6 * cnwheat_model.EcophysiologicalConstants.C_MOLAR_MASS) / cnwheat_model.EcophysiologicalConstants.HEXOSE_MOLAR_MASS_C_RATIO +
                                     (organs_df.fillna(0)[
                                          'nitrates'] * 1E-6 * cnwheat_model.EcophysiologicalConstants.N_MOLAR_MASS) / cnwheat_model.EcophysiologicalConstants.NITRATES_MOLAR_MASS_N_RATIO +
                                     (organs_df.fillna(0)[
                                          'amino_acids'] * 1E-6 * cnwheat_model.EcophysiologicalConstants.N_MOLAR_MASS) / cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                                     (organs_df.fillna(0)[
                                          'proteins'] * 1E-6 * cnwheat_model.EcophysiologicalConstants.N_MOLAR_MASS) / cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                                     organs_df.fillna(0)['mstruct'])
        organs_df['C_g'] = ((organs_df.fillna(0)['sucrose'] * 1E-6 * cnwheat_model.EcophysiologicalConstants.C_MOLAR_MASS) +
                            (organs_df.fillna(0)['starch'] * 1E-6 * cnwheat_model.EcophysiologicalConstants.C_MOLAR_MASS) +
                            (organs_df.fillna(0)[
                                 'amino_acids'] * 1E-6 * cnwheat_model.EcophysiologicalConstants.N_MOLAR_MASS) * cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_MOLAR_MASS_C_RATIO /
                            cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_MOLAR_MASS_N_RATIO + (organs_df.fillna(0)['proteins'] * 1E-6 * cnwheat_model.EcophysiologicalConstants.N_MOLAR_MASS) *
                            cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_MOLAR_MASS_C_RATIO / cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                            organs_df.fillna(0)['mstruct'] * cnwheat_model.EcophysiologicalConstants.RATIO_C_mstruct)
        organs_df['N_g'] = ((organs_df.fillna(0)['nitrates'] * 1E-6 * cnwheat_model.EcophysiologicalConstants.N_MOLAR_MASS) +
                            (organs_df.fillna(0)['amino_acids'] * 1E-6 * cnwheat_model.EcophysiologicalConstants.N_MOLAR_MASS) +
                            (organs_df.fillna(0)['proteins'] * 1E-6 * cnwheat_model.EcophysiologicalConstants.N_MOLAR_MASS) +
                            organs_df.fillna(0)['Nstruct'])

        pp_organs_df['N_tot'] = organs_df['N_g']

        # roots
        roots_df = organs_df.loc[organs_df.organ == 'roots']
        pp_organs_df.loc[pp_organs_df.organ == 'roots', 'WSC_g'] = Roots.calculate_WSC_g(roots_df['sucrose'])
        pp_organs_df.loc[pp_organs_df.organ == 'roots', 'Conc_Nitrates'] = Roots.calculate_Conc_Nitrates(roots_df['nitrates'], roots_df['mstruct'])
        pp_organs_df.loc[pp_organs_df.organ == 'roots', 'Conc_Amino_Acids'] = Roots.calculate_Conc_Amino_Acids(roots_df['amino_acids'], roots_df['mstruct'])
        pp_organs_df.loc[pp_organs_df.organ == 'roots', 'Conc_Sucrose'] = Roots.calculate_conc_sucrose(roots_df['sucrose'], roots_df['mstruct'])
        pp_organs_df.loc[pp_organs_df.organ == 'roots', 'Conc_cytokinins'] = Roots.calculate_conc_cytokinins(roots_df['cytokinins'], roots_df['mstruct'])
        pp_organs_df.loc[pp_organs_df.organ == 'roots', 'Conc_cytokinins'] = Roots.calculate_conc_cytokinins(roots_df['cytokinins'], roots_df['mstruct'])
        R_residual = np.array(
            list(map(respiwheat_model.RespirationModel.R_residual, roots_df['sucrose'], roots_df['mstruct'] * cnwheat_model.Roots.PARAMETERS.ALPHA, roots_df['Total_Organic_Nitrogen'],
                     soils_df['Tsoil'])))
        pp_organs_df.loc[pp_organs_df.organ == 'roots', 'R_residual'] = R_residual

        # phloem
        phloems_df = organs_df.loc[organs_df.organ == 'phloem']
        pp_organs_df.loc[pp_organs_df.organ == 'phloem', 'WSC_g'] = Phloem.calculate_WSC_g(phloems_df['sucrose'])
        # WARNING the computation of this concentration has been switched to a normalization by shoot dry mass only!
        pp_organs_df.loc[pp_organs_df.organ == 'phloem', 'Conc_Amino_Acids'] = Phloem.calculate_conc_amino_acids(phloems_df['amino_acids'], axes_df.set_index(phloems_df.index)['mstruct'] - roots_df.set_index(phloems_df.index)['mstruct'])
        pp_organs_df.loc[pp_organs_df.organ == 'phloem', 'Conc_Sucrose'] = Phloem.calculate_conc_sucrose(phloems_df['sucrose'], axes_df.set_index(phloems_df.index)['mstruct'] - roots_df.set_index(phloems_df.index)['mstruct'])

        # grains
        grains_df = organs_df.loc[organs_df.organ == 'grains']
        pp_organs_df.loc[pp_organs_df.organ == 'grains', 'Dry_Mass'] = Grains.calculate_dry_mass(grains_df['structure'], grains_df['starch'], grains_df['proteins'])
        pp_organs_df.loc[pp_organs_df.organ == 'grains', 'WSC_g'] = Grains.calculate_WSC_g(grains_df['sucrose'], grains_df['starch'])
        pp_organs_df.loc[pp_organs_df.organ == 'grains', 'Proteins_N_Mass'] = Grains.calculate_protein_N_mass(grains_df['proteins'])
        pp_organs_df = pp_organs_df.reindex(columns=ORGANS_RUN_POSTPROCESSING_VARIABLES, copy=False)
        pp_organs_df['plant'] = pp_organs_df['plant'].astype(int)
        returned_dataframes.append(pp_organs_df)
    else:
        returned_dataframes.append(pd.DataFrame({'A': []}))

    # elements
    if elements_df is not None:
        elements_df['sum_dry_mass'] = Element.calculate_dry_mass(elements_df.fillna(0)['triosesP'],
                                                                 elements_df.fillna(0)['sucrose'],
                                                                 elements_df.fillna(0)['starch'],
                                                                 elements_df.fillna(0)['fructan'],
                                                                 elements_df.fillna(0)['nitrates'],
                                                                 elements_df.fillna(0)['amino_acids'],
                                                                 elements_df.fillna(0)['proteins'],
                                                                 elements_df['mstruct'])
        elements_df['sum_dry_mass_total'] = Element.calculate_dry_mass(elements_df.fillna(0)['triosesP'],
                                                                       elements_df.fillna(0)['sucrose'],
                                                                       elements_df.fillna(0)['starch'],
                                                                       elements_df.fillna(0)['fructan'],
                                                                       elements_df.fillna(0)['nitrates'],
                                                                       elements_df.fillna(0)['amino_acids'],
                                                                       elements_df.fillna(0)['proteins'],
                                                                       elements_df['max_mstruct'] + elements_df['Nresidual'])
        elements_df['C_g'] = Element.calculate_C_g(elements_df.fillna(0)['triosesP'],
                                                   elements_df.fillna(0)['sucrose'],
                                                   elements_df.fillna(0)['starch'],
                                                   elements_df.fillna(0)['fructan'],
                                                   elements_df.fillna(0)['amino_acids'],
                                                   elements_df.fillna(0)['proteins'],
                                                   elements_df['mstruct'])
        elements_df['N_g'] = Element.calculate_N_g(elements_df.fillna(0)['nitrates'],
                                                   elements_df.fillna(0)['amino_acids'],
                                                   elements_df.fillna(0)['proteins'],
                                                   elements_df['Nstruct'])
        elements_df['N_g_total'] = Element.calculate_N_g(elements_df.fillna(0)['nitrates'],
                                                         elements_df.fillna(0)['amino_acids'],
                                                         elements_df.fillna(0)['proteins'],
                                                         elements_df['Nstruct'] + elements_df['Nresidual'])
        elements_df['WSC_g'] = Element.calculate_WSC_g(elements_df.fillna(0)['triosesP'],
                                                       elements_df.fillna(0)['sucrose'],
                                                       elements_df.fillna(0)['starch'],
                                                       elements_df.fillna(0)['fructan'], )
        pp_elements_df = pd.concat([elements_df, pd.DataFrame(columns=ELEMENTS_POSTPROCESSING_VARIABLES)], sort=False)
        pp_elements_df.loc[:, 'Conc_TriosesP'] = Element.calculate_conc_triosesP(elements_df['triosesP'], elements_df['mstruct'])
        pp_elements_df.loc[:, 'Conc_Starch'] = Element.calculate_conc_starch(elements_df['starch'], elements_df['mstruct'])
        pp_elements_df.loc[:, 'Conc_Sucrose'] = Element.calculate_conc_sucrose(elements_df['sucrose'], elements_df['mstruct'])
        pp_elements_df.loc[:, 'Conc_Fructan'] = Element.calculate_conc_fructan(elements_df['fructan'], elements_df['mstruct'])
        pp_elements_df.loc[:, 'Conc_Nitrates'] = Element.calculate_Conc_Nitrates(elements_df['nitrates'], elements_df['mstruct'])
        pp_elements_df.loc[:, 'Conc_Amino_Acids'] = Element.calculate_Conc_Amino_Acids(elements_df['amino_acids'], elements_df['mstruct'])
        pp_elements_df.loc[:, 'Conc_Proteins'] = Element.calculate_conc_proteins(elements_df['proteins'], elements_df['mstruct'])
        pp_elements_df.loc[:, 'Conc_cytokinins'] = Element.calculate_conc_cytokinins(elements_df['cytokinins'], elements_df['mstruct'])
        pp_elements_df.loc[:, 'Cont_Fructan_DM'] = Element.calculate_fructan_g(elements_df['fructan']) / elements_df['sum_dry_mass'] * 100
        pp_elements_df.loc[:, 'Cont_WSC_DM'] = elements_df['WSC_g'] / elements_df['sum_dry_mass'] * 100
        pp_elements_df.loc[:, 'Surfacic_NS'] = Element.calculate_surfacic_non_structural(elements_df['sum_dry_mass'],
                                                                                         elements_df['mstruct'], elements_df['green_area'])
        pp_elements_df.loc[:, 'NS'] = Element.calculate_ratio_non_structural(elements_df['sum_dry_mass'],
                                                                             elements_df['mstruct'])
        pp_elements_df.loc[:, 'N_content'] = elements_df['N_g'] / elements_df['sum_dry_mass'] * 100
        pp_elements_df.loc[:, 'N_content_total_DM'] = elements_df['N_g_total'] / elements_df['sum_dry_mass_total'] * 100
        pp_elements_df.loc[:, 'N_tot'] = elements_df['N_g_total']
        pp_elements_df.loc[:, 'SLN'] = Element.calculate_SLN(elements_df['nitrates'], elements_df['amino_acids'], elements_df['proteins'],
                                                             elements_df['Nstruct'], elements_df['green_area'])
        pp_elements_df.loc[pp_elements_df['is_growing'] == 1., 'SLN'] = np.nan
        pp_elements_df.loc[:, 'SLN_nonstruct'] = Element.calculate_SLN_nonstruct(elements_df['nitrates'], elements_df['amino_acids'], elements_df['proteins'],
                                                                                 elements_df['green_area'])
        pp_elements_df.loc[pp_elements_df['is_growing'] == 1., 'SLN_nonstruct'] = np.nan
        pp_elements_df.loc[:, 'SLA'] = Element.calculate_SLA(elements_df['sum_dry_mass'], elements_df['green_area'])
        pp_elements_df.loc[pp_elements_df['is_growing'] == 1., 'SLA'] = np.nan
        pp_elements_df.loc[:, 'Photosynthetic_efficiency'] = elements_df['Ag'] / elements_df['PARa']

        grouped = elements_df.groupby('organ')
        for organ_type, parameters_class in \
                (('ear', cnwheat_parameters.CHAFF_ELEMENT_PARAMETERS),
                 ('blade', cnwheat_parameters.LAMINA_ELEMENT_PARAMETERS),
                 ('internode', cnwheat_parameters.INTERNODE_ELEMENT_PARAMETERS),
                 ('peduncle', cnwheat_parameters.PEDUNCLE_ELEMENT_PARAMETERS),
                 ('sheath', cnwheat_parameters.SHEATH_ELEMENT_PARAMETERS)):
            if organ_type not in grouped.groups:
                continue
            group = grouped.get_group(organ_type)
            if len(group) == 0:  # TODO: faire mm tri que dans simulation de cnwheat (surface nulle)
                continue
            curr_organ_elements_df = elements_df.loc[group.index]
            pp_curr_organ_elements_df = pp_elements_df.loc[group.index]
            R_residual = np.vectorize(respiwheat_model.RespirationModel.R_residual)(curr_organ_elements_df['sucrose'], curr_organ_elements_df['mstruct'] * parameters_class.ALPHA,
                                                                                    curr_organ_elements_df['Total_Organic_Nitrogen'], curr_organ_elements_df['Ts'])
            pp_curr_organ_elements_df.loc[:, 'R_residual'] = R_residual
        pp_elements_df = pp_elements_df.reindex(columns=ELEMENTS_RUN_POSTPROCESSING_VARIABLES, copy=False)
        pp_elements_df[['plant', 'metamer']] = pp_elements_df[['plant', 'metamer']].astype(int)
        returned_dataframes.append(pp_elements_df)
    else:
        returned_dataframes.append(pd.DataFrame({'A': []}))

    # hidden zones
    if hiddenzones_df is not None:

        # hiddenzones_df = hiddenzones_df.merge(nb_replications_df, on='metamer')

        hiddenzones_df['sum_dry_mass'] = HiddenZone.calculate_dry_mass(hiddenzones_df.fillna(0)['sucrose'],
                                                                       0,  # hiddenzones_df.fillna(0)['starch'],
                                                                       hiddenzones_df.fillna(0)['fructan'],
                                                                       hiddenzones_df.fillna(0)['amino_acids'],
                                                                       hiddenzones_df.fillna(0)['proteins'],
                                                                       hiddenzones_df['mstruct'])
        hiddenzones_df['C_g'] = HiddenZone.calculate_C_g(hiddenzones_df.fillna(0)['sucrose'],
                                                         0,  # hiddenzones_df.fillna(0)['starch'],
                                                         hiddenzones_df.fillna(0)['fructan'],
                                                         hiddenzones_df.fillna(0)['amino_acids'],
                                                         hiddenzones_df.fillna(0)['proteins'],
                                                         hiddenzones_df['mstruct'])
        hiddenzones_df['N_g'] = HiddenZone.calculate_N_g(hiddenzones_df.fillna(0)['amino_acids'],
                                                         hiddenzones_df.fillna(0)['proteins'],
                                                         hiddenzones_df['leaf_enclosed_Nstruct'] + hiddenzones_df['internode_enclosed_Nstruct'])
        hiddenzones_df['WSC_g'] = HiddenZone.calculate_WSC_g(hiddenzones_df.fillna(0)['sucrose'],
                                                             hiddenzones_df.fillna(0)['fructan'])
        pp_hiddenzones_df = pd.concat([hiddenzones_df, pd.DataFrame(columns=HIDDENZONE_POSTPROCESSING_VARIABLES)], sort=False)
        pp_hiddenzones_df.loc[:, 'Conc_Amino_Acids'] = HiddenZone.calculate_Conc_Amino_Acids(hiddenzones_df['amino_acids'], hiddenzones_df['mstruct'])
        pp_hiddenzones_df.loc[:, 'Conc_Fructan'] = HiddenZone.calculate_conc_fructan(hiddenzones_df['fructan'], hiddenzones_df['mstruct'])
        pp_hiddenzones_df.loc[:, 'Conc_Proteins'] = HiddenZone.calculate_conc_protein(hiddenzones_df['proteins'], hiddenzones_df['mstruct'])
        pp_hiddenzones_df.loc[:, 'Conc_Sucrose'] = HiddenZone.calculate_conc_sucrose(hiddenzones_df['sucrose'], hiddenzones_df['mstruct'])
        pp_hiddenzones_df.loc[:, 'Cont_Fructan_DM'] = HiddenZone.calculate_fructan_g(hiddenzones_df['fructan']) / hiddenzones_df['sum_dry_mass'] * 100
        pp_hiddenzones_df.loc[:, 'Cont_Proteins_DM'] = HiddenZone.calculate_proteins_g(hiddenzones_df['proteins']) / hiddenzones_df['sum_dry_mass'] * 100
        pp_hiddenzones_df.loc[:, 'Cont_WSC_DM'] = hiddenzones_df['WSC_g'] / hiddenzones_df['sum_dry_mass'] * 100
        pp_hiddenzones_df.loc[:, 'N_content'] = pp_hiddenzones_df['N_g'] / pp_hiddenzones_df['sum_dry_mass'] * 100
        if set(hiddenzones_df.columns).issuperset(['delta_leaf_L', 'leaf_L']):
            # this is temporary: those post-processing should be done in model "elong-wheat"
            pp_hiddenzones_df.loc[:, 'RER'] = HiddenZone.calculate_RER(hiddenzones_df['delta_leaf_L'], hiddenzones_df['leaf_L'], delta_t)
        pp_hiddenzones_df = pp_hiddenzones_df.reindex(columns=HIDDENZONE_RUN_POSTPROCESSING_VARIABLES, copy=False)
        pp_hiddenzones_df[['plant', 'metamer']] = pp_hiddenzones_df[['plant', 'metamer']].astype(int)
        returned_dataframes.append(pp_hiddenzones_df)
    else:
        returned_dataframes.append(pd.DataFrame({'A': []}))

    # axes
    if axes_df is not None:
        axes_df = axes_df[axes_df['axis'] == 'MS'].copy()  # TODO : Temporary !
        pp_axes_df = pd.concat([axes_df, pd.DataFrame(columns=AXES_POSTPROCESSING_VARIABLES)], sort=False)

        # Integrated variables TODO : Homogeneiser la structure de ce bout de code
        if (hiddenzones_df is not None) and (organs_df is not None) and (elements_df is not None):
            # Add missing row if any
            axes_row_keys = organs_df[AXES_T_INDEXES].drop_duplicates()
            pp_axes_df = pp_axes_df.merge(axes_row_keys, how='outer', on=AXES_T_INDEXES)
            pp_axes_df.sort_values(AXES_T_INDEXES, inplace=True)  # Make sure axes_df is sorted
            pp_axes_df.reset_index(drop=True, inplace=True)

            # Roots
            dry_mass_roots = organs_df[(organs_df['organ'] == 'roots')].groupby(['t', 'plant', 'axis'])['sum_dry_mass'].agg('sum')

            # Total mstruct shoot and root
            hz_df_MS = hiddenzones_df[hiddenzones_df['axis'] == 'MS'].copy()
            elt_df_MS = elements_df[elements_df['axis'] == 'MS'].copy()
            hz_df_MS['mstruct_tillers'] = hz_df_MS['mstruct'] * hz_df_MS['nb_replications']
            elt_df_MS['mstruct_tillers'] = elt_df_MS['mstruct'] * elt_df_MS['nb_replications']
            sum_mstruct_shoot = hz_df_MS.groupby(['t', 'plant', 'axis'])['mstruct_tillers'].agg('sum') + elt_df_MS.groupby(['t', 'plant', 'axis'])['mstruct_tillers'].agg('sum')
            sum_mstruct_laminae = elt_df_MS[elt_df_MS.element == 'LeafElement1'].groupby(['t', 'plant', 'axis'])['mstruct_tillers'].agg('sum')
            sum_mstruct_stem = sum_mstruct_shoot - sum_mstruct_laminae
            sum_mstruct_roots = organs_df[(organs_df['organ'] == 'roots')].groupby(['t', 'plant', 'axis'])['mstruct'].agg('sum')

            shoot_roots_mstruct_ratio = sum_mstruct_shoot / sum_mstruct_roots

            # Phloem
            phloem_shoot_root = 1 / (1 + 1 / shoot_roots_mstruct_ratio)
            phloem_laminae = sum_mstruct_laminae / (sum_mstruct_shoot + sum_mstruct_roots)
            phloem_stem = sum_mstruct_stem / (sum_mstruct_shoot + sum_mstruct_roots)
            sum_dry_mass_phloem = organs_df[(organs_df['organ'] == 'phloem')].groupby(['t', 'plant', 'axis'])['sum_dry_mass'].agg('sum')
            sum_dry_mass_phloem_shoot = sum_dry_mass_phloem * phloem_shoot_root
            sum_dry_mass_phloem_laminae = sum_dry_mass_phloem * phloem_laminae
            sum_dry_mass_phloem_roots = sum_dry_mass_phloem * (1 - phloem_shoot_root)
            sum_N_g_phloem_shoot = organs_df[(organs_df['organ'] == 'phloem')].groupby(['t', 'plant', 'axis'])['N_g'].agg('sum') * phloem_shoot_root
            sum_C_g_phloem_shoot = organs_df[(organs_df['organ'] == 'phloem')].groupby(['t', 'plant', 'axis'])['C_g'].agg('sum') * phloem_shoot_root
            sum_WSC_g_phloem_shoot = pp_organs_df[(organs_df['organ'] == 'phloem')].groupby(['t', 'plant', 'axis'])['WSC_g'].agg('sum') * phloem_shoot_root
            sum_WSC_g_phloem_laminae = pp_organs_df[(organs_df['organ'] == 'phloem')].groupby(['t', 'plant', 'axis'])['WSC_g'].agg('sum') * phloem_laminae
            sum_WSC_g_phloem_stem = pp_organs_df[(organs_df['organ'] == 'phloem')].groupby(['t', 'plant', 'axis'])['WSC_g'].agg('sum') * phloem_stem
            sum_WSC_g_phloem_roots = pp_organs_df[(organs_df['organ'] == 'phloem')].groupby(['t', 'plant', 'axis'])['WSC_g'].agg('sum') * (1 - phloem_shoot_root)

            # Total shoot
            hz_df_MS['sum_dry_mass_tillers'] = hz_df_MS['sum_dry_mass'] * hz_df_MS['nb_replications']
            elt_df_MS['sum_dry_mass_tillers'] = elt_df_MS['sum_dry_mass'] * elt_df_MS['nb_replications']
            elt_df_MS['sum_dry_mass_total_tillers'] = elt_df_MS['sum_dry_mass_total'] * elt_df_MS['nb_replications']
            sum_dry_mass_shoot = sum_dry_mass_phloem_shoot + \
                                 hz_df_MS.groupby(['t', 'plant', 'axis'])['sum_dry_mass_tillers'].agg('sum') + \
                                 elt_df_MS.groupby(['t', 'plant', 'axis'])['sum_dry_mass_tillers'].agg('sum')
            sum_dry_mass_laminae = sum_dry_mass_phloem_laminae + \
                                 elt_df_MS[elt_df_MS.element == 'LeafElement1'].groupby(['t', 'plant', 'axis'])['sum_dry_mass_tillers'].agg('sum')
            sum_dry_mass_stem = sum_dry_mass_shoot - sum_dry_mass_laminae
            sum_dry_mass_total_shoot = sum_dry_mass_phloem_shoot + \
                                       hz_df_MS.groupby(['t', 'plant', 'axis'])['sum_dry_mass_tillers'].agg('sum') + \
                                       elt_df_MS.groupby(['t', 'plant', 'axis'])['sum_dry_mass_total_tillers'].agg('sum')
            # Total root
            sum_dry_mass_roots = sum_dry_mass_phloem_roots + dry_mass_roots

            # Total shoot + roots
            sum_dry_mass = sum_dry_mass_shoot + sum_dry_mass_roots
            sum_mstruct = sum_mstruct_roots + sum_mstruct_shoot

            # N content
            hz_df_MS['N_g_tillers'] = hz_df_MS['N_g'] * hz_df_MS['nb_replications']
            elt_df_MS['N_g_tillers'] = elt_df_MS['N_g'] * elt_df_MS['nb_replications']
            elt_df_MS['N_g_total_tillers'] = elt_df_MS['N_g_total'] * elt_df_MS['nb_replications']
            sum_N_g = (organs_df[organs_df['axis'] == 'MS'].groupby(['t', 'plant', 'axis'])['N_g'].agg('sum') +
                       hz_df_MS.groupby(['t', 'plant', 'axis'])['N_g_tillers'].agg('sum') +
                       elt_df_MS.groupby(['t', 'plant', 'axis'])['N_g_tillers'].agg('sum'))
            sum_N_g_total = (organs_df[organs_df['axis'] == 'MS'].groupby(['t', 'plant', 'axis'])['N_g'].agg('sum') +
                             hz_df_MS.groupby(['t', 'plant', 'axis'])['N_g_tillers'].agg('sum') +
                             elt_df_MS.groupby(['t', 'plant', 'axis'])['N_g_total_tillers'].agg('sum'))
            N_content = sum_N_g / sum_dry_mass * 100
            N_content_mstruct = sum_N_g / sum_mstruct * 100

            sum_N_g_shoot = (sum_N_g_phloem_shoot +
                             hz_df_MS.groupby(['t', 'plant', 'axis'])['N_g_tillers'].agg('sum') +
                             elt_df_MS.groupby(['t', 'plant', 'axis'])['N_g_tillers'].agg('sum'))
            sum_N_g_total_shoot = (sum_N_g_phloem_shoot +
                                   hz_df_MS.groupby(['t', 'plant', 'axis'])['N_g_tillers'].agg('sum') +
                                   elt_df_MS.groupby(['t', 'plant', 'axis'])['N_g_total_tillers'].agg('sum'))
            N_content_shoot = sum_N_g_shoot / sum_dry_mass_shoot * 100
            N_content_total_DM_shoot = sum_N_g_total_shoot / sum_dry_mass_total_shoot * 100
            N_content_mstruct_shoot = sum_N_g_shoot / sum_mstruct_shoot * 100

            N_content_roots = (N_content * sum_dry_mass - N_content_shoot * sum_dry_mass_shoot) / sum_dry_mass_roots
            N_content_mstruct_roots = (N_content_mstruct * sum_mstruct - N_content_mstruct_shoot * sum_mstruct_shoot) / sum_mstruct_roots

            # WSC
            hz_df_MS['WSC_g_tillers'] = hz_df_MS['WSC_g'] * hz_df_MS['nb_replications']
            elt_df_MS['WSC_g_tillers'] = elt_df_MS['WSC_g'] * elt_df_MS['nb_replications']

            WSC_g_plant = hz_df_MS.groupby(['t', 'plant', 'axis'])['WSC_g_tillers'].agg('sum') + elt_df_MS.groupby(['t', 'plant', 'axis'])['WSC_g_tillers'].agg('sum') + \
                          pp_organs_df.groupby(['t', 'plant', 'axis'])['WSC_g'].agg('sum')

            sum_WSC_g_shoot = (sum_WSC_g_phloem_shoot +
                             hz_df_MS.groupby(['t', 'plant', 'axis'])['WSC_g_tillers'].agg('sum') +
                             elt_df_MS.groupby(['t', 'plant', 'axis'])['WSC_g_tillers'].agg('sum'))
            sum_WSC_g_laminae = (sum_WSC_g_phloem_laminae +
                               elt_df_MS[elt_df_MS.element == 'LeafElement1' ].groupby(['t', 'plant', 'axis'])['WSC_g_tillers'].agg('sum'))
            sum_WSC_g_stem = (sum_WSC_g_phloem_stem +
                               hz_df_MS.groupby(['t', 'plant', 'axis'])['WSC_g_tillers'].agg('sum') +
                               elt_df_MS[~(elt_df_MS.element == 'LeafElement1')].groupby(['t', 'plant', 'axis'])['WSC_g_tillers'].agg('sum'))
            sum_WSC_g_roots = (sum_WSC_g_phloem_roots + pp_organs_df[(organs_df['organ'] == 'roots')].groupby(['t', 'plant', 'axis'])['WSC_g'].agg('sum'))


            # C/N ratio
            hz_df_MS['C_g_tillers'] = hz_df_MS['C_g'] * hz_df_MS['nb_replications']
            elt_df_MS['C_g_tillers'] = elt_df_MS['C_g'] * elt_df_MS['nb_replications']

            sum_C_g = (organs_df[organs_df['axis'] == 'MS'].groupby(['t', 'plant', 'axis'])['C_g'].agg('sum') +
                       hz_df_MS.groupby(['t', 'plant', 'axis'])['C_g_tillers'].agg('sum') +
                       elt_df_MS.groupby(['t', 'plant', 'axis'])['C_g_tillers'].agg('sum'))
            sum_C_g_shoot = (sum_C_g_phloem_shoot +
                             hz_df_MS.groupby(['t', 'plant', 'axis'])['C_g_tillers'].agg('sum') +
                             elt_df_MS.groupby(['t', 'plant', 'axis'])['C_g_tillers'].agg('sum'))

            C_N_ratio = sum_C_g / sum_N_g
            C_N_ratio_shoot = sum_C_g_shoot / sum_N_g_shoot

            # C non structural
            sum_NSC_g = sum_C_g.values - (axes_df['mstruct'] * cnwheat_model.EcophysiologicalConstants.RATIO_C_mstruct)

            # Photosyntheses
            elements_df['Tillers_Photosynthesis'] = elements_df['Photosynthesis'] * elements_df['nb_replications']
            elements_df['Tillers_Photosynthesis_An'] = elements_df['An'] * elements_df['green_area'] * 3600 * elements_df['nb_replications']
            tillers_photosynthesis = elements_df[elements_df['axis'] == 'MS'].groupby(['t', 'plant', 'axis'])['Tillers_Photosynthesis'].agg('sum')  # TEMPORARY : porter au niveau de la plante
            tillers_photosynthesis_An = elements_df[elements_df['axis'] == 'MS'].groupby(['t', 'plant', 'axis'])['Tillers_Photosynthesis_An'].agg('sum')
            tot_photosynthesis = elements_df[elements_df['axis'] == 'MS'].groupby(['t', 'plant', 'axis'])['Photosynthesis'].agg('sum')

            # INN
            DM_t_ha = sum_dry_mass_shoot * 250 * 10 ** -2  # convert from g.plant-1 to t.ha-1
            N_content_critical = np.where(DM_t_ha < 1.55, 4.4, 5.35 * DM_t_ha ** -0.442)  # from Justes 1994 : valid at field scale from Feekes 3 i.e. mid tillering
            NNI = N_content_shoot / N_content_critical

            # Ratio Non Structural Mass
            NS_shoot = (1 - sum_mstruct_shoot / sum_dry_mass_shoot) * 100
            NS_stem = (1 - sum_mstruct_stem / sum_dry_mass_stem) * 100
            NS_laminae = (1 - sum_mstruct_laminae / sum_dry_mass_laminae) * 100
            NS_roots = (1 - sum_mstruct_roots / sum_dry_mass_roots) * 100
            NS = (1 - sum_mstruct / sum_dry_mass) * 100

            # C_respired_shoot
            #hz_df_MS['Respi_growth_tillers'] = hz_df_MS['Respi_growth'].fillna(0) * hz_df_MS['nb_replications'].fillna(1)
            #C_respired_shoot = pp_axes_df.sum_respi_shoot.fillna(0.) + hz_df_MS.groupby(AXES_T_INDEXES, as_index=False)['Respi_growth_tillers'].sum().Respi_growth_tillers

            # C_respired_roots
            #C_respired_roots = pp_axes_df.sum_respi_roots.fillna(0.) + organs_df[(organs_df['organ'] == 'roots')]['Respi_growth'].reset_index(drop=True)

            # Add to axes df
            pp_axes_df.sort_values(AXES_T_INDEXES, inplace=True)  # Make sure axes_df is sorted
            pp_axes_df.loc[:, 'C_N_ratio'] = C_N_ratio.values
            pp_axes_df.loc[:, 'C_N_ratio_shoot'] = C_N_ratio_shoot.values
            pp_axes_df.loc[:, 'N_content'] = N_content.values
            pp_axes_df.loc[:, 'N_content_shoot'] = N_content_shoot.values
            pp_axes_df.loc[:, 'N_content_roots'] = N_content_roots.values
            pp_axes_df.loc[:, 'N_content_mstruct'] = N_content_mstruct.values
            pp_axes_df.loc[:, 'N_content_mstruct_shoot'] = N_content_mstruct_shoot.values
            pp_axes_df.loc[:, 'N_content_total_DM_shoot'] = N_content_total_DM_shoot.values
            pp_axes_df.loc[:, 'N_content_mstruct_roots'] = N_content_mstruct_roots.values
            pp_axes_df.loc[:, 'sum_N_g'] = sum_N_g.values
            pp_axes_df.loc[:, 'sum_N_g_shoot'] = sum_N_g_shoot.values
            pp_axes_df.loc[:, 'sum_dry_mass'] = sum_dry_mass.values
            pp_axes_df.loc[:, 'sum_dry_mass_shoot'] = sum_dry_mass_shoot.values
            pp_axes_df.loc[:, 'sum_dry_mass_laminae'] = sum_dry_mass_laminae.values
            pp_axes_df.loc[:, 'sum_dry_mass_stem'] = sum_dry_mass_stem.values
            pp_axes_df.loc[:, 'sum_dry_mass_roots'] = sum_dry_mass_roots.values
            pp_axes_df.loc[:, 'sum_C_g'] = sum_C_g.values
            pp_axes_df.loc[:, 'sum_C_g_shoot'] = sum_C_g_shoot.values
            pp_axes_df.loc[:, 'sum_NSC_g'] = sum_NSC_g.values
            pp_axes_df.loc[:, 'dry_mass_phloem'] = sum_dry_mass_phloem.values
            pp_axes_df.loc[:, 'shoot_roots_ratio'] = pp_axes_df['sum_dry_mass_shoot'] / pp_axes_df['sum_dry_mass_roots']
            pp_axes_df.loc[:, 'shoot_roots_mstruct_ratio'] = shoot_roots_mstruct_ratio.values
            pp_axes_df.loc[:, 'Total_Photosynthesis'] = tot_photosynthesis.values
            pp_axes_df.loc[:, 'Tillers_Photosynthesis'] = tillers_photosynthesis.values
            pp_axes_df.loc[:, 'Tillers_Photosynthesis_An'] = tillers_photosynthesis_An.values
            pp_axes_df.loc[:, 'NNI'] = NNI.values
            pp_axes_df.loc[:, 'NS_roots'] = NS_roots.values
            pp_axes_df.loc[:, 'NS_shoot'] = NS_shoot.values
            pp_axes_df.loc[:, 'NS_stem'] = NS_stem.values
            pp_axes_df.loc[:, 'NS_laminae'] = NS_laminae.values
            pp_axes_df.loc[:, 'NS'] = NS.values
            pp_axes_df.loc[:, 'mstruct_shoot'] = sum_mstruct_shoot.values
            pp_axes_df.loc[:, 'mstruct_laminae'] = sum_mstruct_laminae.values
            pp_axes_df.loc[:, 'mstruct_stem'] = sum_mstruct_stem.values
            # pp_axes_df.loc[:, 'C_respired_shoot'] = C_respired_shoot
            # pp_axes_df.loc[:, 'C_respired_roots'] = C_respired_roots
            pp_axes_df.loc[:, 'WSC_g'] = WSC_g_plant.values
            pp_axes_df.loc[:, 'Cont_WSC_DM'] = pp_axes_df['WSC_g'] / pp_axes_df['sum_dry_mass'] * 100
            pp_axes_df.loc[:, 'Cont_WSC_DM_shoot'] = sum_WSC_g_shoot.values / sum_dry_mass_shoot.values * 100
            pp_axes_df.loc[:, 'Cont_WSC_DM_laminae'] = sum_WSC_g_laminae.values / sum_dry_mass_laminae.values * 100
            pp_axes_df.loc[:, 'Cont_WSC_DM_stem'] = sum_WSC_g_stem.values / sum_dry_mass_stem.values * 100
            pp_axes_df.loc[:, 'Cont_WSC_DM_roots'] = sum_WSC_g_roots.values / sum_dry_mass_roots.values * 100

        pp_axes_df = pp_axes_df.reindex(AXES_RUN_POSTPROCESSING_VARIABLES, axis=1, copy=False)
        pp_axes_df['plant'] = pp_axes_df['plant'].astype(int)
        returned_dataframes.append(pp_axes_df)
    else:
        returned_dataframes.append(pd.DataFrame({'A': []}))

    # soils
    if soils_df is not None:
        pp_soils_df = pd.concat([soils_df, pd.DataFrame(columns=SOILS_POSTPROCESSING_VARIABLES)])
        pp_soils_df = pp_soils_df.reindex(columns=SOILS_RUN_POSTPROCESSING_VARIABLES, copy=False)
        pp_soils_df[['plant']] = pp_soils_df[['plant']].astype(int)
        returned_dataframes.append(pp_soils_df)
    else:
        returned_dataframes.append(pd.DataFrame({'A': []}))

    return tuple(returned_dataframes)


# -----------------------------------------------------------
# --------------- GRAPHS GENERATION FRONT-END ---------------
# -- PLEASE USE THIS FUNCTION FOR THE GENERATION OF GRAPHS --
# -----------------------------------------------------------

def generate_graphs(axes_df=None, hiddenzones_df=None, organs_df=None, elements_df=None, soils_df=None, graphs_dirpath='.'):
    """
    Generate graphs to validate the outputs of CN-Wheat, and save them in directory `graphs_dirpath`.

    :param pandas.DataFrame axes_df: CN-Wheat outputs and post-processing at axis scale (see :attr:`PLANTS_RUN_POSTPROCESSING_VARIABLES`)
    :param pandas.DataFrame hiddenzones_df: CN-Wheat outputs at hidden zone scale (see :attr:`HIDDENZONE_RUN_POSTPROCESSING_VARIABLES`)
    :param pandas.DataFrame organs_df: CN-Wheat outputs at organ scale (see :attr:`ORGANS_RUN_POSTPROCESSING_VARIABLES`)
    :param pandas.DataFrame elements_df: CN-Wheat outputs at element scale (see :attr:`ELEMENTS_RUN_POSTPROCESSING_VARIABLES`)
    :param pandas.DataFrame soils_df: CN-Wheat outputs at soil scale (see :attr:`SOILS_RUN_POSTPROCESSING_VARIABLES`)
    :param str graphs_dirpath: the path of the directory to save the generated graphs in
    """

    x_name = 't'
    x_label = 'Time (Hour)'

    colors = ['blue', 'darkorange', 'green', 'red', 'darkviolet', 'gold', 'magenta', 'brown', 'darkcyan', 'grey', 'lime']
    colors = colors + colors + colors + colors + colors

    # 1) Photosynthetic organs
    if elements_df is not None:
        elements_df = elements_df.loc[elements_df['mstruct'] != 0]
        graph_variables_ph_elements = {'Ag': u'Gross photosynthesis (�mol m$^{-2}$ s$^{-1}$)', 'Tr': u'Organ surfacic transpiration rate (mmol H$_{2}$0 m$^{-2}$ s$^{-1}$)',
                                       'Transpiration': u'Organ transpiration rate (mmol H$_{2}$0 s$^{-1}$)', 'Ts': u'Temperature surface (�C)', 'Conc_TriosesP': u'[TriosesP] (�mol g$^{-1}$ mstruct)',
                                       'Conc_Starch': u'[Starch] (�mol g$^{-1}$ mstruct)', 'Conc_Sucrose': u'[Sucrose] (�mol g$^{-1}$ mstruct)', 'Conc_Fructan': u'[Fructan] (�mol g$^{-1}$ mstruct)',
                                       'Conc_Nitrates': u'[Nitrates] (�mol g$^{-1}$ mstruct)', 'Conc_Amino_Acids': u'[Amino_Acids] (�mol g$^{-1}$ mstruct)',
                                       'Conc_Proteins': u'[Proteins] (g g$^{-1}$ mstruct)', 'Cont_Fructan_DM': u'Fructan content (% DM)', 'Cont_WSC_DM': u'WSC content (% DM)',
                                       'Nitrates_import': u'Total nitrates imported (�mol h$^{-1}$)',
                                       'Amino_Acids_import': u'Total amino acids imported (�mol N h$^{-1}$)', 'S_Amino_Acids': u'[Rate of amino acids synthesis] (�mol N g$^{-1}$ mstruct h$^{-1}$)',
                                       'S_Proteins': u'Rate of protein synthesis (�mol N g$^{-1}$ mstruct h$^{-1}$)', 'D_Proteins': u'Rate of protein degradation (�mol N g$^{-1}$ mstruct h$^{-1}$)',
                                       'Loading_Sucrose': u'Loading Sucrose (�mol C sucrose h$^{-1}$)', 'Loading_Amino_Acids': u'Loading Amino acids (�mol N amino acids h$^{-1}$)',
                                       'green_area': u'Green area (m$^{2}$)', 'R_phloem_loading': u'Respiration phloem loading (�mol C h$^{-1}$)',
                                       'R_Nnit_red': u'Respiration nitrate reduction (�mol C h$^{-1}$)', 'R_residual': u'Respiration residual (�mol C h$^{-1}$)', 'mstruct': u'Structural mass (g)',
                                       'Nstruct': u'Structural N mass (g)', 'Conc_cytokinins': u'[cytokinins] (UA g$^{-1}$ mstruct)', 'D_cytokinins': u'Cytokinin degradation (UA g$^{-1}$ mstruct)',
                                       'cytokinins_import': u'Cytokinin import (UA)', 'Surfacic_NS': u'Surfacic Non Structural mass (g m$^{-2}$)',
                                       'NS': u'Ratio of Non Structural mass', 'N_content': u'N content in the green tissues (% DM)',
                                       'N_content_total_DM': u'N content in the green + senesced tissues (% DM)', 'N_tot': u'N mass (g)',
                                       'SLA': u'Specific Leaf Area (m$^{2}$.kg$^{-1}$)',
                                       'SLN': u'Surfacic Leaf Nitrogen (g.m$^{-2}$)',
                                       'SLN_nonstruct': u'Surfacic Leaf Non-structural Nitrogen (g.m$^{-2}$)', 'length': u'Length (m)',
                                       'Photosynthetic_efficiency': u'Photosynthetic yield (�mol C/�mol PARa)'}

        for org_ph in (['blade'], ['sheath'], ['internode'], ['peduncle', 'ear']):
            for variable_name, variable_label in graph_variables_ph_elements.items():
                graph_name = variable_name + '_' + '_'.join(org_ph) + '.PNG'
                cnwheat_tools.plot_cnwheat_ouputs(elements_df,
                                                  x_name=x_name,
                                                  y_name=variable_name,
                                                  x_label=x_label,
                                                  y_label=variable_label,
                                                  colors=[colors[i - 1] for i in elements_df.metamer.unique().tolist()],
                                                  filters={'organ': org_ph},
                                                  plot_filepath=os.path.join(graphs_dirpath, graph_name),
                                                  explicit_label=False)

    # 2) Roots, grains and phloem
    if organs_df is not None:
        # 'R_growth': u'Growth respiration of roots (�mol C h$^{-1}$)',
        graph_variables_organs = {'Conc_Sucrose': u'[Sucrose] (�mol g$^{-1}$ mstruct)', 'Dry_Mass': 'Dry mass (g)', 'Conc_Nitrates': u'[Nitrates] (�mol g$^{-1}$ mstruct)',
                                  'Conc_Amino_Acids': u'[Amino Acids] (�mol g$^{-1}$ mstruct)', 'Proteins_N_Mass': u'[N Proteins] (g)', 'Uptake_Nitrates': u'Nitrates uptake (�mol h$^{-1}$)',
                                  'sucrose': u'Sucrose (�mol)', 'amino_acids': u'Amino Acids (�mol)', 'Unloading_Sucrose': u'Unloaded sucrose (�mol C g$^{-1}$ mstruct h$^{-1}$)',
                                  'Unloading_Amino_Acids': u'Unloaded Amino Acids (�mol N AA g$^{-1}$ mstruct h$^{-1}$)',
                                  'S_Amino_Acids': u'Rate of amino acids synthesis (�mol N g$^{-1}$ mstruct h$^{-1}$)', 'S_Proteins': u'Rate of protein synthesis (�mol N h$^{-1}$)',
                                  'Export_Nitrates': u'Total export of nitrates (�mol N h$^{-1}$)', 'Export_Amino_Acids': u'Total export of Amino acids (�mol N h$^{-1}$)',
                                  'R_Nnit_upt': u'Respiration nitrates uptake (�mol C h$^{-1}$)', 'R_Nnit_red': u'Respiration nitrate reduction (�mol C h$^{-1}$)',
                                  'R_residual': u'Respiration residual (�mol C h$^{-1}$)',
                                  'R_grain_growth_struct': u'Respiration grain structural growth (�mol C h$^{-1}$)', 'R_grain_growth_starch': u'Respiration grain starch growth (�mol C h$^{-1}$)',
                                  'mstruct': u'Structural mass (g)', 'C_exudation': u'Carbon lost by root exudation (�mol C g$^{-1}$ mstruct h$^{-1}$',
                                  'N_exudation': u'Nitrogen lost by root exudation (�mol N g$^{-1}$ mstruct h$^{-1}$', 'Conc_cytokinins': u'[cytokinins] (UA g$^{-1}$ mstruct)',
                                  'S_cytokinins': u'Rate of cytokinins synthesis (UA g$^{-1}$ mstruct)', 'Export_cytokinins': 'Export of cytokinins from roots (UA h$^{-1}$)',
                                  'HATS_LATS': u'Potential uptake (�mol h$^{-1}$)', 'regul_transpiration': u'Regulating transpiration function',
                                  'senesced_mstruct': u'Cumulated senesced mstruct (g)', 'synthetized_mstruct': u'Cumulated synthetized mstruct (g)'}

        for org in (['roots'], ['grains'], ['phloem']):
            for variable_name, variable_label in graph_variables_organs.items():
                graph_name = variable_name + '_' + '_'.join(org) + '.PNG'
                cnwheat_tools.plot_cnwheat_ouputs(organs_df,
                                                  x_name=x_name,
                                                  y_name=variable_name,
                                                  x_label=x_label,
                                                  y_label=variable_label,
                                                  colors=['blue'],
                                                  filters={'organ': org},
                                                  plot_filepath=os.path.join(graphs_dirpath, graph_name),
                                                  explicit_label=False)

    # 3) Soil
    if soils_df is not None:
        _, (ax1) = plt.subplots(1)
        Conc_Nitrates_soil = soils_df['Conc_Nitrates_Soil'] * 14E-6
        ax1.plot(soils_df['t'], Conc_Nitrates_soil)
        ax1.set_ylabel(u'[Nitrates] (g m$^{-3}$)')
        ax1.set_xlabel('Time from flowering (hour)')
        ax1.set_title = 'Conc Nitrates Soil'
        ax1.set_ylim(bottom=0)
        plt.savefig(os.path.join(graphs_dirpath, 'Conc_Nitrates_Soil.PNG'), format='PNG', bbox_inches='tight')
        plt.close()

    # 4) Hidden zones
    if hiddenzones_df is not None:
        graph_variables_hiddenzones = {'Conc_Sucrose': u'[Sucrose] (�mol g$^{-1}$ mstruct)', 'Conc_Amino_Acids': u'[Amino Acids] (�mol g$^{-1}$ mstruct)',
                                       'Conc_Proteins': u'[Proteins] (g g$^{-1}$ mstruct)', 'Conc_Fructan': u'[Fructan] (�mol g$^{-1}$ mstruct)', 'Cont_Fructan_DM': u'Fructan content (% DM)',
                                       'Cont_Proteins_DM': u'Protein content (% DM)', 'Cont_WSC_DM': u'WSC content (% DM)', 'Unloading_Sucrose': u'Rate of Sucrose unloading (�mol C h${-1}$)',
                                       'Unloading_Amino_Acids': u'Rate of Amino_acids unloading (�mol N h${-1}$)', 'mstruct': u'Structural mass (g)', 'Nstruct': u'Structural N mass (g)',
                                       'leaf_L': u'Leaf length in hz (m)', 'delta_leaf_L': u'delta of leaf length (m)', 'internode_L': u'Internode length in hz (m)',
                                       'leaf_pseudostem_length': u'leaf pseudostem length (m)'}

        for variable_name, variable_label in graph_variables_hiddenzones.items():
            graph_name = variable_name + '_hz' + '.PNG'
            cnwheat_tools.plot_cnwheat_ouputs(hiddenzones_df,
                                              x_name=x_name,
                                              y_name=variable_name,
                                              x_label=x_label,
                                              y_label=variable_label,
                                              colors=[colors[i - 1] for i in hiddenzones_df.metamer.unique().tolist()],
                                              filters={'plant': 1, 'axis': 'MS'},
                                              plot_filepath=os.path.join(graphs_dirpath, graph_name),
                                              explicit_label=False)

    # 4) Axes
    graph_variables_axes = {'mstruct': 'Axis mstruct (g)',
                            'C_N_ratio': u'C/N mass ratio', 'C_N_ratio_shoot': u'C/N mass ratio of the shoot',
                            'N_content': u'N content in the axis (% DM)', 'N_content_shoot': u'N content in the shoot (% DM)',
                            'N_content_total_DM_shoot': u'N content in the shoot acounting for partly senesced leaves (% DM)',
                            'N_content_roots': u'N content in the roots (% DM)',
                            'N_content_mstruct': u'N content in the axis (% mstruct)', 'N_content_mstruct_shoot': u'N content in the shoot (% mstruct)',
                            'N_content_mstruct_roots': u'N content in the roots (% mstruct)',
                            'sum_N_g': u'N mass (g)', 'sum_N_g_shoot': u'N mass in the shoot (g)',
                            'shoot_roots_ratio': u'Shoot/Roots dry mass ratio',
                            'shoot_roots_mstruct_ratio': u'Shoot/Roots mstruct ratio',
                            'sum_dry_mass': u'Total dry mass (g)', 'sum_dry_mass_shoot': u'Dry mass of the shoot (g)',
                            'sum_dry_mass_roots': u'Dry mass of the roots (g)', 'dry_mass_phloem': u'Dry mass of the phloem (g)',
                            'NNI': u'Nitrogen Nutrition Index',
                            'NS_shoot': u'Ratio of Non Structural Mass in the shoot', 'NS_roots': u'Ratio of Non Structural Mass in the roots',
                            'NS': u'Ratio of Non Structural Mass for the plant',
                            'mstruct_shoot': u'Structural Mass of the shoot (g)',
                            'Cont_WSC_DM': u'WSC content (% DM)'}

    for variable_name, variable_label in graph_variables_axes.items():
        graph_name = variable_name + '_axis' + '.PNG'
        cnwheat_tools.plot_cnwheat_ouputs(axes_df,
                                          x_name=x_name,
                                          y_name=variable_name,
                                          x_label=x_label,
                                          y_label=variable_label,
                                          colors=['blue'],
                                          filters={'plant': 1, 'axis': 'MS'},
                                          plot_filepath=os.path.join(graphs_dirpath, graph_name),
                                          explicit_label=False)
