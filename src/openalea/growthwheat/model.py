# -*- coding: latin-1 -*-

from __future__ import division  # use "//" to do integer division

from openalea.growthwheat import parameters as growthwheat_parameters

"""
    growthwheat.model
    ~~~~~~~~~~~~~

    The module :mod:`growthwheat.model` defines the equations of the kinetic of leaf growth (mass flows) according to leaf elongation. Also includes root growth.

"""

class GrowthWheatModel(object):
    def __init__(self):
        self.parameters = growthwheat_parameters.parameters_factory()
        self.organ_init = growthwheat_parameters.OrganInit()

        #: the inputs needed by GrowthWheat
        # -: hiddenzones
        self.hiddenzone_inputs = ['leaf_is_growing', 'internode_is_growing', 'leaf_pseudo_age', 'delta_leaf_pseudo_age',
                             'internode_pseudo_age', 'delta_internode_pseudo_age', 'leaf_L', 'delta_leaf_L',
                             'internode_L', 'delta_internode_L', 'leaf_pseudostem_length',
                             'delta_leaf_pseudostem_length', 'internode_distance_to_emerge',
                             'delta_internode_distance_to_emerge', 'SSLW', 'LSSW', 'LSIW', 'leaf_is_emerged',
                             'internode_is_visible',
                             'sucrose', 'amino_acids', 'fructan', 'proteins', 'leaf_enclosed_mstruct',
                             'leaf_enclosed_Nstruct', 'internode_enclosed_mstruct',
                             'internode_enclosed_Nstruct', 'mstruct', 'internode_Lmax', 'leaf_Lmax', 'sheath_Lmax',
                             'is_over', 'leaf_is_remobilizing', 'internode_is_remobilizing']

        self.element_inputs = ['is_growing', 'mstruct', 'senesced_mstruct', 'green_area', 'length', 'sucrose', 'amino_acids',
                          'fructan', 'proteins', 'cytokinins', 'Nstruct']
        self.root_inputs = ['sucrose', 'amino_acids', 'mstruct', 'Nstruct']
        self.axis_inputs = ['delta_teq', 'delta_teq_roots', 'nb_leaves']

        #: the outputs computed by GrowthWheat
        self.hiddenzone_outputs = ['sucrose', 'amino_acids', 'fructan', 'proteins', 'leaf_enclosed_mstruct',
                              'leaf_enclosed_Nstruct', 'internode_enclosed_mstruct', 'internode_enclosed_Nstruct',
                              'mstruct', 'Nstruct', 'Respi_growth', 'sucrose_consumption_mstruct', 'AA_consumption_mstruct', 'is_over',
                              'leaf_is_remobilizing', 'internode_is_remobilizing']
        self.element_outputs = ['sucrose', 'amino_acids', 'fructan', 'proteins', 'nitrates', 'mstruct', 'Nstruct',
                           'green_area', 'max_proteins', 'max_mstruct', 'Nresidual', 'senesced_length_element',
                           'senesced_mstruct']
        self.root_outputs = ['sucrose', 'amino_acids', 'mstruct', 'Nstruct', 'Respi_growth', 'R_min_upt', 'delta_mstruct_growth',
                        'sucrose_consumption_mstruct', 'AA_consumption_mstruct']
        self.axis_outputs = []

    @staticmethod
    def calculate_ratio_mstruct_DM(mstruct, sucrose, fructans, amino_acids, proteins):
        """
        Ratio mstruct/dry matter (dimensionless)

        :param float mstruct: Structural mass (g)
        :param float sucrose: Sucrose amount (µmol C)
        :param float fructans: Fructans amount (µmol C)
        :param float amino_acids: Amino acids amount (µmol N)
        :param float proteins: proteins amount (µmol N)

        :return: Ratio mstruct/dry matter (dimensionless)
        :rtype: float
        """
        C_MOLAR_MASS = growthwheat_parameters.C_MOLAR_MASS
        N_MOLAR_MASS = growthwheat_parameters.N_MOLAR_MASS
        dry_mass = ((sucrose * 1E-6 * C_MOLAR_MASS) / growthwheat_parameters.HEXOSE_MOLAR_MASS_C_RATIO +
                    (fructans * 1E-6 * C_MOLAR_MASS) / growthwheat_parameters.HEXOSE_MOLAR_MASS_C_RATIO +
                    (amino_acids * 1E-6 * N_MOLAR_MASS) / growthwheat_parameters.AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                    (proteins * 1E-6 * N_MOLAR_MASS) / growthwheat_parameters.AMINO_ACIDS_MOLAR_MASS_N_RATIO +
                    mstruct)

        return mstruct / dry_mass

    def calculate_delta_leaf_enclosed_mstruct(self, leaf_L, delta_leaf_L, ratio_mstruct_DM, init_leaf_L=None, leaf_pseudo_age=None):
        """ Relation between length and mstruct for the leaf segment located in the hidden zone during the exponential-like growth phase.
        Parameters alpha_mass_growth and beta_mass_growth estimated from Williams (1960) and expressed in g of dry mass
        The actual ratio_mstruct_DM is then used to convert in g of structural dry mass.

        :param float leaf_L: Total leaf length (m)
        :param float delta_leaf_L: delta of leaf length (m)
        :param float ratio_mstruct_DM: Ratio mstruct/dry matter (dimensionless)
        :param float init_leaf_L: Total leaf length before update in turgor-growth sub-model (m). Not used in this instance of the model
        :param float leaf_pseudo_age: Pseudo age of the leaf since beginning of automate elongation (s). Not used in this instance of the model

        :return: delta_leaf_enclosed_mstruct (g)
        :rtype: float
        """
        return self.parameters.ALPHA * self.parameters.BETA * leaf_L ** (self.parameters.BETA - 1) * delta_leaf_L * ratio_mstruct_DM


    def calculate_delta_leaf_enclosed_mstruct_postE(self, delta_leaf_pseudo_age, leaf_pseudo_age, leaf_pseudostem_L, enclosed_mstruct, LSSW, sucrose, mstruct):
        """ mstruct of the enclosed leaf from the emergence of the leaf to the end of elongation.
        Final mstruct of the enclosed leaf matches sheath mstruct calculation when it is mature.
        #TODO : Hiddenzone mstruct calculation is not correct for sheath shorten than previous one.
        :param float delta_leaf_pseudo_age: Delta of Pseudo age of the leaf since beginning of automate elongation (s)
        :param float leaf_pseudo_age: Pseudo age of the leaf since beginning of automate elongation (s)
        :param float leaf_pseudostem_L: Pseudostem length (m)
        :param float enclosed_mstruct: mstruct of the enclosed leaf (g)
        :param float LSSW: Lineic Structural Sheath Weight (g m-1)
        :param float sucrose: Sucrose amount (µmol C)
        :param float mstruct: Structural mass (g)

        :return: delta_leaf_enclosed_mstruct (g)
        :rtype: float
        """
        conc_sucrose_effective = max(0., sucrose / mstruct - self.parameters.conc_sucrose_offset)

        enclosed_mstruct_max = leaf_pseudostem_L * LSSW

        if leaf_pseudo_age < self.parameters.te and conc_sucrose_effective > 0:
            delta_enclosed_mstruct = (enclosed_mstruct_max - enclosed_mstruct) / (self.parameters.te - leaf_pseudo_age) * delta_leaf_pseudo_age
        else:
            delta_enclosed_mstruct = 0

        return max(0., delta_enclosed_mstruct)

    def calculate_delta_internode_enclosed_mstruct(self, internode_L, delta_internode_L, ratio_mstruct_DM):
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
        return self.parameters.RATIO_ENCLOSED_LEAF_INTERNODE * self.parameters.ALPHA * self.parameters.BETA * internode_L ** (self.parameters.BETA - 1) * delta_internode_L * ratio_mstruct_DM

    def calculate_delta_internode_enclosed_mstruct_postL(self, delta_internode_pseudo_age, internode_pseudo_age, internode_L, internode_pseudostem_L, internode_Lmax, LSIW, enclosed_mstruct):
        """ mstruct of the enclosed internode from the ligulation of the leaf to the end of elongation.
        Final mstruct of the enclosed internode matches internode mstruct calculation when it is mature.

        :param float delta_internode_pseudo_age: Delta of Pseudo age of the internode since beginning of automate elongation (s)
        :param float internode_pseudo_age: Pseudo age of the internode since beginning of automate elongation (s)
        :param float internode_L: Current length of the internode (m)
        :param float internode_pseudostem_L: Pseudostem length of the internode (m)
        :param float internode_Lmax: Final length of the internode (m)
        :param float LSIW: Lineic Structural Internode Weight (g m-1)
        :param float enclosed_mstruct: mstruct of the enclosed leaf (g)


        :return: delta_internode_enclosed_mstruct (g)
        :rtype: float
        """
        if not internode_Lmax:
            enclosed_mstruct_max = internode_L * LSIW
        else:
            enclosed_mstruct_max = min(internode_pseudostem_L, internode_Lmax) * LSIW

        if internode_pseudo_age < self.parameters.te_IN:
            delta_enclosed_mstruct = (enclosed_mstruct_max - enclosed_mstruct) / (self.parameters.te_IN - internode_pseudo_age) * delta_internode_pseudo_age
        else:
            delta_enclosed_mstruct = 0

        return max(0., delta_enclosed_mstruct)

    @staticmethod
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

    def calculate_delta_Nstruct(self, delta_mstruct):
        """ delta Nstruct of hidden zone and emerged tissue (lamina and sheath).

        :param float delta_mstruct: delta of mstruct (g)

        :return: delta Nstruct (g)
        :rtype: float
        """
        return delta_mstruct * self.parameters.RATIO_AMINO_ACIDS_MSTRUCT

    @staticmethod
    def calculate_export(delta_mstruct, metabolite, hiddenzone_mstruct):
        """Export of metabolite from the hidden zone towards the emerged part of the leaf integrated over delta_t.

        :param float delta_mstruct: Delta of structural dry mass of the emerged part of the leaf (g)
        :param float metabolite: Metabolite amount in the hidden zone (µmol C or N)
        :param float hiddenzone_mstruct: Structural mass of the hidden zone (g)

        :return: metabolite export (µmol N)
        :rtype: float
        """
        return delta_mstruct * max(0., (metabolite / hiddenzone_mstruct))

    def calculate_init_cytokinins_emerged_tissue(self, delta_mstruct):
        """Initial amount of cytokinins allocated in the mstruct of a newly emerged tissue.

        :param float delta_mstruct: Delta of structural dry mass of the emerged part of the leaf (g)

        :return: cytokinins addition (AU)
        :rtype: float
        """
        return delta_mstruct * self.parameters.INIT_CYTOKININS_EMERGED_TISSUE  # TODO: Set according to protein concentration ?

    @staticmethod
    def calculate_s_Nstruct_amino_acids(delta_hiddenzone_Nstruct, delta_lamina_Nstruct, delta_sheath_Nstruct, delta_internode_Nstruct):
        """Consumption of amino acids for the calculated mstruct growth (µmol N consumed by mstruct growth)

        :param float delta_hiddenzone_Nstruct: Nstruct growth of the hidden zone (g)
        :param float delta_lamina_Nstruct: Nstruct growth of the lamina (g)
        :param float delta_sheath_Nstruct: Nstruct growth of the sheath (g)
        :param float delta_internode_Nstruct: Nstruct growth of the internode (g)

        :return: Amino acid consumption (µmol N)
        :rtype: float
        """
        return (delta_hiddenzone_Nstruct + delta_lamina_Nstruct + delta_sheath_Nstruct + delta_internode_Nstruct) / growthwheat_parameters.N_MOLAR_MASS * 1E6

    def calculate_s_mstruct_sucrose(self, delta_hiddenzone_mstruct, delta_lamina_mstruct, delta_sheath_mstruct, s_Nstruct_amino_acids_N):
        """Consumption of sucrose for the calculated mstruct growth (µmol C consumed by mstruct growth)

        :param float delta_hiddenzone_mstruct: mstruct growth of the hidden zone (g)
        :param float delta_lamina_mstruct: mstruct growth of the lamina (g)
        :param float delta_sheath_mstruct: mstruct growth of the sheath (g)
        :param float s_Nstruct_amino_acids_N: Total amino acid consumption (µmol N) due to Nstruct (µmol N)

        :return: Sucrose consumption (µmol C)
        :rtype: float
        """
        s_Nstruct_amino_acids = s_Nstruct_amino_acids_N / self.parameters.AMINO_ACIDS_N_RATIO  #: µmol of AA
        s_mstruct_amino_acids_C = s_Nstruct_amino_acids * self.parameters.AMINO_ACIDS_C_RATIO  #: µmol of C coming from AA
        s_mstruct_C = (delta_hiddenzone_mstruct + delta_lamina_mstruct + delta_sheath_mstruct) * self.parameters.RATIO_SUCROSE_MSTRUCT / growthwheat_parameters.C_MOLAR_MASS * 1E6  #: Total C for mstruct growth (µmol C)
        s_mstruct_sucrose_C = s_mstruct_C - s_mstruct_amino_acids_C  #: µmol of C coming from sucrose

        return s_mstruct_sucrose_C

    @staticmethod
    def calculate_sheath_mstruct(sheath_L, LSSW):
        """ mstruct of the sheath.
          Final mstruct of the enclosed leaf matches sheath mstruct calculation when it is mature.

          :param float sheath_L: Sheath length (m)
          :param float LSSW: Lineic Structural Sheath Weight (g m-1).

          :return: Structural mass of the sheath (g)
          :rtype: float
          """
        return sheath_L * LSSW

    def calculate_roots_mstruct_growth(self, sucrose, amino_acids, mstruct, delta_teq, postflowering_stages, nb_leaves,
                                       xylem_water_potential=None):
        """Root structural dry mass growth integrated over delta_t

        :param float sucrose: Amount of sucrose in roots (µmol C)
        :param float amino_acids: Amount of amino acids in roots (µmol N)
        :param float mstruct: Root structural mass (g)
        :param float delta_teq: Time compensated for the effect of temperature - Time equivalent at Tref (s)
        :param bool postflowering_stages: Option : True to run a simulation with postflo parameter
        :param int nb_leaves: Current number of leaves on the axis
        :param float xylem_water_potential: Water potential of xylem (Mpa)

        :return: mstruct_C_growth (µmol C), mstruct_growth (g), Nstruct_growth (g), Nstruct_N_growth (µmol N)
        :rtype: (float, float, float, float)
        """
        conc_sucrose_effective = max(0., sucrose / mstruct - self.parameters.conc_sucrose_offset)

        if postflowering_stages:
            Vmax = self.parameters.VMAX_ROOTS_GROWTH_POSTFLO
            K_ROOTS_GROWTH = self.parameters.K_ROOTS_GROWTH_POSTFLO
        else:
            Vmax = self.parameters.Vmax_roots_growth(nb_leaves)
            K_ROOTS_GROWTH = self.parameters.K_ROOTS_GROWTH_PREFLO
        N = self.parameters.N_ROOTS_GROWTH

        if conc_sucrose_effective > 0.:
            mstruct_C_growth = max(0., ((conc_sucrose_effective ** N) * Vmax) / ((conc_sucrose_effective ** N) + (K_ROOTS_GROWTH ** N)) * delta_teq * mstruct)  #: root growth in C (µmol of C)
        else:
            mstruct_C_growth = 0.
        mstruct_growth = mstruct_C_growth * growthwheat_parameters.CONVERSION_MMOL_C_G_MSTRUCT_ROOTS  #: root growth (g of structural dry mass)

        Nstruct_growth = mstruct_growth * growthwheat_parameters.RATIO_N_MSTRUCT_ROOTS  #: root growth in N (g of structural dry mass)
        Nstruct_N_growth = min(amino_acids, (Nstruct_growth / growthwheat_parameters.N_MOLAR_MASS) * 1E6)  #: root growth in nitrogen (µmol N)

        return mstruct_C_growth, mstruct_growth, Nstruct_growth, Nstruct_N_growth

    def calculate_roots_s_mstruct_sucrose(self, delta_roots_mstruct, s_Nstruct_amino_acids_N):
        """Consumption of sucrose for the calculated mstruct growth (µmol C consumed by mstruct growth)

        :param float delta_roots_mstruct: mstruct growth of the roots (g)
        :param float s_Nstruct_amino_acids_N: Total amino acid consumption (µmol N) due to Nstruct (µmol N)

        :return: Sucrose consumption (µmol C)
        :rtype: float
        """
        s_Nstruct_amino_acids = s_Nstruct_amino_acids_N / self.parameters.AMINO_ACIDS_N_RATIO  #: µmol of AA
        s_mstruct_amino_acids_C = s_Nstruct_amino_acids * self.parameters.AMINO_ACIDS_C_RATIO  #: µmol of C coming from AA
        s_mstruct_C = delta_roots_mstruct * self.parameters.RATIO_SUCROSE_MSTRUCT / growthwheat_parameters.C_MOLAR_MASS * 1E6  #: Total C used for mstruct growth (µmol C)
        s_mstruct_sucrose_C = s_mstruct_C - s_mstruct_amino_acids_C  #: µmol of coming from sucrose

        return s_mstruct_sucrose_C

    def calculate_mineral_plant(self, mstruct, senesced_mstruct):
        """ Mineral mass.

        :param float mstruct: structural mass of the plant (g)
        :param float senesced_mstruct: senesced structural mass of the plant (g)

        :return: Mineral mass of the plant (g)
        :rtype: float
        """
        mineral_plant = (mstruct * self.parameters.MINERAL_LIVING_TISSUE) + (senesced_mstruct * self.parameters.MINERAL_SENESCED_TISSUE)
        return mineral_plant


class GrowthWheatModelHydraulics(GrowthWheatModel):
    def __init__(self):
        super().__init__()
        self.parameters = growthwheat_parameters.parameters_factory(hydraulics=True)

        #: the inputs needed by GrowthWheat
        # -: Hiddenzones
        self.hiddenzone_inputs = list(self.hiddenzone_inputs) + ['init_leaf_L']
        # -: Axes
        self.axis_inputs = list(self.axis_inputs) + ['xylem_water_potential']

    def calculate_delta_leaf_enclosed_mstruct(self, leaf_L, delta_leaf_L, ratio_mstruct_DM, init_leaf_L=0, leaf_pseudo_age=0):
        """ Relation between length and mstruct for the leaf segment located in the hidden zone during the exponential-like growth phase.
        Parameters alpha_mass_growth and beta_mass_growth estimated from Williams (1960) and expressed in g of dry mass
        The actual ratio_mstruct_DM is then used to convert in g of structural dry mass.

        :param float leaf_L: Total leaf length after update in turgor-growth sub-model (m)
        :param float delta_leaf_L: delta of leaf length (m)
        :param float ratio_mstruct_DM: Ratio mstruct/dry matter (dimensionless)
        :param float init_leaf_L: Total leaf length before update in turgor-growth sub-model (m)
        :param float leaf_pseudo_age: Pseudo age of the leaf since beginning of automate elongation (s)

        :return: delta_leaf_enclosed_mstruct (g)
        :rtype: float
        """
        if leaf_pseudo_age >= 0:
            delta_leaf_L_update = leaf_L - init_leaf_L
        else:
            delta_leaf_L_update = delta_leaf_L

        return self.parameters.ALPHA * self.parameters.BETA * leaf_L ** (
                    self.parameters.BETA - 1) * delta_leaf_L_update * ratio_mstruct_DM


    def calculate_roots_mstruct_growth(self, sucrose, amino_acids, mstruct, delta_teq, postflowering_stages,
                                       nb_leaves=None, xylem_water_potential=0):
        """Root structural dry mass growth integrated over delta_t

        :param float sucrose: Amount of sucrose in roots (µmol C)
        :param float amino_acids: Amount of amino acids in roots (µmol N)
        :param float mstruct: Root structural mass (g)
        :param float delta_teq: Time compensated for the effect of temperature - Time equivalent at Tref (s)
        :param bool postflowering_stages: Option : True to run a simulation with postflo parameter
        :param int nb_leaves: Current number of leaves on the axis
        :param float xylem_water_potential: Water potential of xylem (Mpa)

        :return: mstruct_C_growth (µmol C), mstruct_growth (g), Nstruct_growth (g), Nstruct_N_growth (µmol N)
        :rtype: (float, float, float, float)
        """
        conc_sucrose_effective = max(0., sucrose / mstruct - self.parameters.conc_sucrose_offset)

        if postflowering_stages:
            Vmax = self.parameters.VMAX_ROOTS_GROWTH_POSTFLO
        else:
            Vmax = self.parameters.VMAX_ROOTS_GROWTH_PREFLO
        N = self.parameters.N_ROOTS_GROWTH

        if conc_sucrose_effective > 0.:
            mstruct_C_growth = max(0., ((conc_sucrose_effective ** N) * Vmax) / ((conc_sucrose_effective ** N) + (
                        self.parameters.K_ROOTS_GROWTH ** N)) * delta_teq * mstruct)  #: root growth in C (µmol of C)
        else:
            mstruct_C_growth = 0.

        hydraulic_regulation = 1 / (1 + (
                    xylem_water_potential / self.parameters.water_potential_crit) ** self.parameters.n)  #: Regulation with plant water status

        mstruct_growth = mstruct_C_growth * growthwheat_parameters.CONVERSION_MMOL_C_G_MSTRUCT_ROOTS * hydraulic_regulation  #: root growth (g of structural dry mass)

        Nstruct_growth = mstruct_growth * growthwheat_parameters.RATIO_N_MSTRUCT_ROOTS  #: root growth in N (g of structural dry mass)
        Nstruct_N_growth = min(amino_acids,
                               (Nstruct_growth / growthwheat_parameters.N_MOLAR_MASS) * 1E6)  #: root growth in nitrogen (µmol N)

        return mstruct_C_growth, mstruct_growth, Nstruct_growth, Nstruct_N_growth