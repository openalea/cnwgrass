# -*- coding: latin-1 -*-

from __future__ import division  # use "//" to do integer division
import pandas as pd
from openalea.cnwgrass.morphogenesis import parameters as morphogenesis_parameters
from math import exp, log10

"""
    morphogenesis.model
    ~~~~~~~~~~~~~

    The module :mod:`morphogenesis.model` defines the equations of the kinetic of leaf and internode elongation according to C and N status.

"""


class MorphogenesisModel(object):
    def __init__(self):
        self.parameters = morphogenesis_parameters.parameters_factory()
        self.hiddenzone_init = morphogenesis_parameters.HiddenZoneInit()
        self.element_init = morphogenesis_parameters.ElementInit()

        #: the inputs needed by Morphogenesis
        # -: Axes
        self.axis_inputs = ['SAM_temperature', 'delta_teq', 'teq_since_primordium', 'status', 'nb_leaves', 'GA',
                            'SAM_height', 'cohort', 'sum_TT']
        #-: Hiddenzones
        self.hiddenzone_inputs = ['leaf_is_growing', 'internode_is_growing', 'leaf_pseudo_age', 'internode_pseudo_age',
                             'leaf_pseudostem_length', 'internode_distance_to_emerge', 'leaf_L', 'internode_L',
                             'hiddenzone_age', 'leaf_Lmax', 'leaf_Lmax_em', 'lamina_Lmax', 'sheath_Lmax', 'leaf_Wmax',
                             'SSLW', 'LSSW', 'leaf_is_emerged', 'internode_Lmax', 'internode_Lmax_lig', 'LSIW',
                             'internode_is_visible', 'sucrose', 'amino_acids', 'fructan', 'proteins',
                             'leaf_enclosed_mstruct', 'leaf_enclosed_Nstruct', 'internode_enclosed_mstruct',
                             'internode_enclosed_Nstruct', 'mstruct', 'mean_conc_sucrose', 'leaf_is_remobilizing', 'internode_is_remobilizing']
        # -: Elements
        self.element_inputs = ['length', 'Wmax', 'is_growing', 'age', 'age_teq']


        #: the outputs computed by Morphogenesis
        # TODO : add be default all the attributes of the class HiddenZoneInit and ElementInit, and define which attribute is set by growth.parameters or morphogenesis.parameters
        # -: Axes
        self.axis_outputs = ['SAM_temperature', 'delta_teq', 'delta_teq_roots', 'teq_since_primordium', 'status',
                             'nb_leaves', 'GA', 'SAM_height', 'cohort', 'sum_TT']

        # -: Hiddenzones
        self.hiddenzone_outputs = ['leaf_is_growing', 'internode_is_growing', 'leaf_pseudo_age', 'delta_leaf_pseudo_age', 'internode_pseudo_age', 'delta_internode_pseudo_age', 'leaf_pseudostem_length',
                      'hiddenzone_age', 'delta_leaf_pseudostem_length', 'internode_distance_to_emerge',
                      'delta_internode_distance_to_emerge', 'leaf_L', 'delta_leaf_L', 'internode_L', 'delta_internode_L', 'leaf_Lmax', 'leaf_Lmax_em', 'lamina_Lmax', 'sheath_Lmax', 'leaf_Wmax',
                      'SSLW', 'LSSW', 'leaf_is_emerged', 'internode_Lmax', 'internode_Lmax_lig', 'LSIW', 'internode_is_visible', 'leaf_enclosed_mstruct',
                      'ratio_DZ', 'mean_conc_sucrose', 'leaf_is_remobilizing', 'internode_is_remobilizing']
        # -: Elements
        self.element_outputs = ['age', 'age_teq']


        #: topology columns for ligule height dataframe
        self.ligule_topology_columns = ['axis_id', 'phytomer', 'ligule height']



    # -------------------------------------------------------------------------------------------------------------------
    # --- SAM
    # -------------------------------------------------------------------------------------------------------------------
    @staticmethod
    def calculate_growing_temperature(Tair, Tsol, SAM_height, sowing_depth):
        """ Return temperature to be used for growth zone

        :param float Tair: Air temperature at t (degree Celsius)
        :param float Tsol: Soil temperature at t (degree Celsius)
        :param float SAM_height: Height of SAM, calculated from internode length (m)
        :param float sowing_depth: Sowing depth (m)

        :return: Return temperature to be used for growth zone at t (degree Celsius)
        :rtype: float
        """

        if SAM_height > sowing_depth:
            growth_temperature = Tair
        else:
            growth_temperature = Tsol
        return growth_temperature

    def modified_Arrhenius_equation(self, temperature):
        """ Return value of equation from Johnson and Lewin (1946) for temperature. The equation is modified to return zero below zero degree.

        :param float temperature: organ temperature (degree Celsius)

        :return: Return value of Eyring equation from Johnson and Lewin (1946) for temperature. The equation is modified to return zero below zero degree.
        :rtype: float
        """

        def Arrhenius_equation(T):
            """
            :param float T: temperature (Kelvin)
            """
            return T * exp(-self.parameters.Temp_Ea_R / T) / (
                        1 + exp(self.parameters.Temp_DS_R - self.parameters.Temp_DH_R / T))

        temperature_K = temperature + 273.15

        if temperature < 0:
            res = 0
        elif temperature < self.parameters.Temp_Ttransition:
            res = temperature * Arrhenius_equation(
                self.parameters.Temp_Ttransition + 273.15) / self.parameters.Temp_Ttransition
        else:
            res = Arrhenius_equation(temperature_K)

        return res

    def calculate_time_equivalent_Tref(self, temperature, time):
        """ Return the time equivalent to a reference temperature i.e. temperature-compensated time (Parent, 2010).

        :param float temperature: temperature (degree Celsius)
        :param float time: time duration (s)

        :return: temperature-compensated time (s)
        :rtype: float
        """
        return time * self.modified_Arrhenius_equation(temperature) / self.modified_Arrhenius_equation(
            self.parameters.Temp_Tref)

    def calculate_cumulated_thermal_time(self, sum_TT, temperature, delta_teq):
        """ Return cumulated thermal time (used by Adel-Wheat model to calculate leaf geometry).

        :param float sum_TT: cumulated thermal time (degree-days)
        :param float temperature: temperature (degree Celsius)
        :param float delta_teq: time duration equivalent at Tref(s)

        :return: cumulated thermal time (s)
        :rtype: float
        """
        if temperature > 0:
            return sum_TT + delta_teq * self.parameters.Temp_Tref / 24.0 / 3600  # Conversion to days at 12°C
        else:
            return sum_TT

    @staticmethod
    def calculate_element_age_teq(age_teq, delta_teq):
        """ Age time equivalent of the element (s)

        :param float age_teq: Previous age_teq of the element since emergence (s)
        :param float delta_teq: Temperature-compensated time = time duration at a reference temperature (s)

        :return: Updated age_teq (s)
        :rtype: float
        """
        return age_teq + delta_teq

    def calculate_SAM_primodia(self, status, teq_since_primordium, delta_teq, nb_leaves, cohort_id):
        """ Update SAM status, leaf number

        :param str status: SAM status ('vegetative', if emitting leaf primordia or 'reproductive')
        :param float teq_since_primordium: Time since last primordium initiation (in time equivalent to a reference temperature) (s)
        :param float delta_teq: time increment (in time equivalent to a reference temperature) (s)
        :param int nb_leaves: Number of leaves already emitted by the SAM.
        :param int cohort_id: Corresponding leaf on the Main Stem for the first leaf of a tiller

        :return: Number of leaf to be initiated (should be 0 or 1), updated leaf number on the SAM, status, time since last primordium initiation (in time equivalent to a reference temperature, s)
        :rtype: (int, int, str, float)
        """
        init_leaf = 0
        teq_since_primordium += delta_teq

        while teq_since_primordium > self.parameters.PLASTOCHRONE and status == 'vegetative':
            init_leaf += 1
            teq_since_primordium -= self.parameters.PLASTOCHRONE
            if init_leaf > 1:
                raise ValueError(
                    'Error : {} leaf primordia have been created in one time step. Consider shorter timestep'.format(
                        init_leaf))
            if (nb_leaves + init_leaf) >= (self.parameters.max_nb_leaves - cohort_id + 1):
                status = 'reproductive'

        nb_leaves += init_leaf

        return init_leaf, nb_leaves, status, teq_since_primordium

    def calculate_SAM_GA(self, status, teq_since_primordium):
        """ Synthesis of GA by the SAM according to its stage

        :param str status: SAM status ('vegetative', if emitting leaf primordia or 'reproductive')
        :param float teq_since_primordium: Time since last primordium initiation (in time equivalent to a reference temperature) (s)

        :return: whether GA production or not
        :rtype: bool
        """
        is_producing_GA = status == 'reproductive' and teq_since_primordium > self.parameters.delta_TT_GA
        return is_producing_GA

    # -------------------------------------------------------------------------------------------------------------------
    # --- Leaves
    # -------------------------------------------------------------------------------------------------------------------
    @staticmethod
    def calculate_ligule_height(sheath_internode_length, all_element_inputs, SAM_id, all_ligule_height_df):
        """ Calculate ligule heights below each phytomer of an axis from lengths of internodes and sheaths.

        The result is formated into a dataframe which is concatenated with a shared dataframe storing all ligule heights

        :param dict sheath_internode_length: Dictionary with sheath and internode lengths for a given axis, (m)
        :param dict all_element_inputs: Dictionary of all element inputs
        :param tuple SAM_id: Id of the SAM (plant_id, axis_id)
        :param pandas.DataFrame all_ligule_height_df: Shared dataframe used to store all ligule heights

        :return: A dataframe with ligule height of a given axis
        :rtype: pandas.DataFrame
        """

        for phytomer_id, lengths in sheath_internode_length.items():
            lamina_id = SAM_id + (phytomer_id, 'blade', 'LeafElement1')

            if lamina_id in all_element_inputs.keys() and not all_element_inputs[lamina_id]['is_growing']:
                ligule_height = sum(lengths['sheath'] + lengths['cumulated_internode'])
                ligule_height_df = pd.DataFrame([(SAM_id, phytomer_id, ligule_height)],
                                                columns=list(all_ligule_height_df))
                all_ligule_height_df = pd.concat((all_ligule_height_df, ligule_height_df))

            elif phytomer_id == 0 and len(sheath_internode_length[phytomer_id]['sheath']) != 0:  # Coleoptile
                coleoptile_length_df = pd.DataFrame(
                    [(SAM_id, phytomer_id, sheath_internode_length[phytomer_id]['sheath'][0])],
                    columns=list(all_ligule_height_df))
                all_ligule_height_df = pd.concat((all_ligule_height_df, coleoptile_length_df))

        return all_ligule_height_df

    @staticmethod
    def calculate_leaf_pseudostem_length(ligule_heights, bottom_hiddenzone_height, phytomer_id):
        """ Distance between the bottom of the hidden zone and the highest previous ligule

        :param pandas.Series ligule_heights: height of all ligules of a given axis (m)
        :param float bottom_hiddenzone_height: height of the bottom of the hidden zone given by the cumulated lengths of below internodes (m)
        :param int phytomer_id: phytomer number

        :return: Distance for the leaf to emerge out of the pseudostem (m)
        :rtype: float
        """
        top_ligule_height = float(
            max(ligule_heights[ligule_heights['phytomer'] < phytomer_id]['ligule height']))  # highest previous ligule
        leaf_pseudostem_length = top_ligule_height - bottom_hiddenzone_height

        return max(leaf_pseudostem_length, 0)

    def calculate_deltaL_preE(self, sucrose, leaf_L, amino_acids, mstruct, delta_teq, leaf_rank, optimal_growth_option,
                              xylem_water_potential=None):
        """ Delta of leaf length over delta_t as a function of sucrose and amino acids, from initiation to the emergence of the previous leaf.

        :param float sucrose: Amount of sucrose (µmol C)
        :param float leaf_L: Total leaf length (m)
        :param float amino_acids: Amount of amino acids (µmol N)
        :param float mstruct: Structural mass (g)
        :param float delta_teq: Temperature-compensated time = time duration at a reference temperature (s)
        :param int leaf_rank: leaf phytomer number
        :param bool optimal_growth_option: if True the function will calculate leaf elongation assuming optimal growth conditions (except if sucrose and amino acids are zero)
        :param None xylem_water_potential: not used in this instance of the model

        :return: delta delta_leaf_L (m)
        :rtype: float
        """

        conc_sucrose_effective = max(0., sucrose / mstruct - self.parameters.conc_sucrose_offset)

        if optimal_growth_option:
            RER_max = self.parameters.RERmax_Ljutovac_fit.get(leaf_rank, self.parameters.RERmax_Ljutovac_fit[
                max(self.parameters.RERmax_Ljutovac_fit.keys())])
            delta_leaf_L = leaf_L * RER_max * delta_teq
        elif conc_sucrose_effective > 0 and amino_acids > 0:
            RER_max = self.parameters.RERmax.get(leaf_rank, self.parameters.RERmax[max(self.parameters.RERmax.keys())])
            # Enzymatic rate for bi-substrates with random fixation
            conc_amino_acids = (amino_acids / mstruct)
            delta_leaf_L = delta_teq * leaf_L * RER_max / ((1 + self.parameters.RER_Kc / conc_sucrose_effective) * (
                        1 + self.parameters.RER_Kn / conc_amino_acids))
        else:
            delta_leaf_L = 0

        return delta_leaf_L

    def calculate_leaf_pseudo_age(self, leaf_pseudo_age, delta_teq, turgor_water_potential=None):
        """ Pseudo age of the leaf since beginning of automate elongation (s)

        :param float leaf_pseudo_age: Previous pseudo age of the leaf since beginning of automate elongation (s)
        :param float delta_teq: Temperature-compensated time = time duration at a reference temperature (s)
        :param float turgor_water_potential: Not used in this instance of the model

        :return: Updated leaf_pseudo_age (s)
        :rtype: float
        """
        return leaf_pseudo_age + delta_teq

    def Beta_function(self, leaf_pseudo_age, leaf_rank):
        """ Normalized leaf length from the emergence of the previous leaf to the end of elongation (automate function depending on leaf pseudo age).

        :param float leaf_pseudo_age: Pseudo age of the leaf since beginning of automate elongation (s)
        :param int leaf_rank: rank of the leaf

        :return: Normalized leaf length (m)
        :rtype: float
        """
        if leaf_rank in (1, 2):
            leaf_L = abs((1 + (max(0., (self.parameters.te_init - leaf_pseudo_age)) / (
                        self.parameters.te_init - self.parameters.tm_init))) *
                         (min(1.0, float(leaf_pseudo_age - self.parameters.tb_init) / float(
                             self.parameters.te_init - self.parameters.tb_init)) **
                          ((self.parameters.te_init - self.parameters.tb_init) / (
                                      self.parameters.te_init - self.parameters.tm_init))))
        else:
            leaf_L = abs(
                (1 + (max(0., (self.parameters.te - leaf_pseudo_age)) / (self.parameters.te - self.parameters.tm))) *
                (min(1.0,
                     float(leaf_pseudo_age - self.parameters.tb) / float(self.parameters.te - self.parameters.tb)) **
                 ((self.parameters.te - self.parameters.tb) / (self.parameters.te - self.parameters.tm))))
        return leaf_L

    def calculate_deltaL_postE(self, prev_leaf_pseudo_age, leaf_pseudo_age, prev_leaf_L, leaf_Lmax, sucrose,
                               amino_acids, mstruct, leaf_rank, optimal_growth_option=False):
        """ Leaf length from the emergence of the previous leaf to the end of elongation (automate function depending on leaf pseudo age and final length).

        :param float prev_leaf_pseudo_age: Pseudo age of the leaf since beginning of automate elongation at previous time step (s)
        :param float leaf_pseudo_age: Pseudo age of the leaf since beginning of automate elongation (s)
        :param float prev_leaf_L: Leaf length at previous time step (m)
        :param float leaf_Lmax: Final leaf length (m)
        :param float sucrose: Amount of sucrose (µmol C)
        :param float amino_acids: Amount of amino acids (µmol N)
        :param float mstruct: Structural mass (g)
        :param int leaf_rank: rank of the leaf
        :param bool optimal_growth_option: if True the function will calculate leaf elongation assuming optimal growth conditions (except if sucrose and amino acids are zero)

        :return: delta_leaf_L (m)
        :rtype: float
        """

        conc_sucrose_effective = max(0., sucrose / mstruct - self.parameters.conc_sucrose_offset)

        if leaf_pseudo_age <= self.parameters.tb and conc_sucrose_effective > 0 and amino_acids > 0:
            delta_leaf_L = prev_leaf_L - self.Beta_function(0., leaf_rank) * leaf_Lmax
        elif leaf_pseudo_age < self.parameters.te:
            # Beta function
            delta_leaf_L_Beta_0 = min(leaf_Lmax, leaf_Lmax * (
                        self.Beta_function(leaf_pseudo_age, leaf_rank) - self.Beta_function(prev_leaf_pseudo_age,
                                                                                            leaf_rank)))

            if optimal_growth_option:
                # Current leaf length
                delta_leaf_L = delta_leaf_L_Beta_0
            elif conc_sucrose_effective > 0 and amino_acids > 0:
                # Regulation by C and N
                conc_amino_acids = amino_acids / mstruct
                regul = self.parameters.leaf_pseudo_age_Vmax / (
                            1 + self.parameters.leaf_pseudo_age_Kc / conc_sucrose_effective) / (
                                    1 + self.parameters.leaf_pseudo_age_Kn / conc_amino_acids)
                # Actual leaf elongation
                delta_leaf_L = delta_leaf_L_Beta_0 * regul
            else:  # empty sucrose or amino acids
                delta_leaf_L = 0
        else:  # end of leaf elongation
            delta_leaf_L = 0

        return max(0., delta_leaf_L)

    def calculate_update_leaf_Lmax(self, leaf_Lmax_em, leaf_L, leaf_pseudo_age, leaf_rank):
        """ Update leaf_Lmax following a reduction of delta_leaf_L due to C and N regulation.
        Updated final length is calculated as the sum of the theoretical remaining length to elongate (leaf_Lmax_em * (1 - Beta_function(leaf_pseudo_age)))
        and the actual elongation at the end of the time step. This could lead to shorter or longer leaves, but the duration of elongation is not modified.

        :param float leaf_Lmax_em: Estimate of final leaf length at previous leaf emergence (m)
        :param float leaf_L: actual leaf length at the end of the time step, calculated according to CN concentration (m)
        :param float leaf_pseudo_age: Pseudo age of the leaf since beginning of automate elongation at the end of the time step (s)
        :param int leaf_rank: rank of the leaf

        :return: updated leaf_Lmax (m)
        :rtype: float
        """
        return leaf_L + leaf_Lmax_em * (1 - self.Beta_function(leaf_pseudo_age, leaf_rank))

    def calculate_ratio_DZ_postE(self, leaf_L, leaf_Lmax, leaf_pseudostem_length):
        """ Ratio of the hiddenzone length which is made of division zone.
        Calculated from an inverse beta function representing the ratio of the leaf composed by the division zone according to its relative length in log.
        The model was fitted on literature data on wheat.

        :param float leaf_L: Leaf length (m)
        :param float leaf_Lmax: Final leaf length (m)
        :param float leaf_pseudostem_length: Length of the pseudostem (m)

        :return: ratio_DZ (dimensionless)
        :rtype: float
        """
        log_L_init = log10(self.parameters.ratio_DZ_l_init)
        log_L_mid = log10(self.parameters.ratio_DZ_l_mid)
        log_Lend = log10(self.parameters.ratio_DZ_l_end)

        leaf_L_normalised = leaf_L / leaf_Lmax
        log_leaf_L_normalised = log10(leaf_L_normalised)

        if log_leaf_L_normalised >= log_Lend:
            return 0
        elif log_leaf_L_normalised <= log_L_init:
            return 1
        else:
            return min((1 - (1 + (log_Lend - log_leaf_L_normalised) / (log_Lend - log_L_mid)) *
                        (((log_leaf_L_normalised - log_L_init) / (log_Lend - log_L_init)) ** (
                                    (log_Lend - log_L_init) / (log_Lend - log_L_mid)))) * leaf_L / min(leaf_L,
                                                                                                       leaf_pseudostem_length),
                       1.)

    @staticmethod
    def calculate_leaf_emergence(leaf_L, leaf_pseudostem_length):
        """Calculate if a given leaf has emerged from its pseudostem

        :param float leaf_L: Leaf length (m)
        :param float leaf_pseudostem_length: Length of the pseudostem (m)

        :return: Specifies if the leaf has emerged (True) or not (False)
        :rtype: bool
        """
        return leaf_L > leaf_pseudostem_length

    @staticmethod
    def calculate_lamina_L(leaf_L, leaf_pseudostem_length, lamina_Lmax):
        """ Emerged lamina length given by the difference between leaf length and leaf pseudostem length.

        :param float leaf_L: Total leaf length (m)
        :param float leaf_pseudostem_length: Length of the pseudostem (m)
        :param float lamina_Lmax: Final lamina length (m)

        :return: lamina length (m)
        :rtype: float
        """
        lamina_L = leaf_L - leaf_pseudostem_length

        return max(10 ** -5,
                   min(lamina_L,
                       lamina_Lmax))  # Minimum length set to 10^-5 m to make sure growth can run even if the lamina turns back hidden (case when an older sheath elongates faster)

    def calculate_leaf_Lmax(self, leaf_Lem_prev, leaf_rank):
        """ Final leaf length.

        :param float leaf_Lem_prev: Leaf length at the emergence of the previous leaf (m)
        :param int leaf_rank: rank of the leaf

        :return: Final leaf length (m)
        :rtype: float
        """
        return min(leaf_Lem_prev / self.Beta_function(0., leaf_rank), self.parameters.leaf_Lmax_MAX)

    def calculate_SL_ratio(self, leaf_rank):
        """ Sheath:Lamina final length ratio according to the rank. Parameters from Dornbush (2011).

        :param int leaf_rank: leaf phytomer number

        :return: Sheath:Lamina ratio (dimensionless)
        :rtype: float
        """
        return self.parameters.SL_ratio_a * leaf_rank ** 3 + self.parameters.SL_ratio_b * leaf_rank ** 2 + self.parameters.SL_ratio_c * leaf_rank + self.parameters.SL_ratio_d

    @staticmethod
    def calculate_lamina_Lmax(leaf_Lmax, sheath_lamina_ratio):
        """ Final lamina length.

        :param float leaf_Lmax: Final leaf length (m)
        :param float sheath_lamina_ratio: Sheath:Lamina ratio (dimensionless)

        :return: Final lamina length (m)
        :rtype: float
        """
        return leaf_Lmax / (1 + sheath_lamina_ratio)

    @staticmethod
    def calculate_sheath_Lmax(leaf_Lmax, lamina_Lmax):
        """ Final sheath length.

        :param float leaf_Lmax: Final leaf length (m)
        :param float lamina_Lmax: Final lamina length (m)

        :return: Final sheath length (m)
        :rtype: float
        """
        return leaf_Lmax - lamina_Lmax

    @staticmethod
    def calculate_mean_conc_sucrose(prev_mean_conc_sucrose, time_prev_leaf2_emergence, delta_teq, sucrose, mstruct):
        """ Update the mean sucrose concentration of the hiddenzone since leaf n-2 emergence.
        Calculation starts at leaf n-2 emergence, the updated mean accounts for the [sucrose] of current time step weighted by a function of temperature.

        :param float prev_mean_conc_sucrose: Mean sucrose concentration of the hiddenzone at the end of previous simulation time step (µmol C g-1).
        :param float time_prev_leaf2_emergence: Time elapsed since leaf n-2 emergence (s at Tref).
        :param float delta_teq: Duration of the current simulation time step (s at Tref).
        :param float sucrose: Sucrose in the hidden zone (µmol C).
        :param float mstruct: Mstruct of the hidden zone (g).

        :return: Updated mean sucrose concentration of the hiddenzone since leaf n-2 emergence (µmol C g-1)
        :rtype: float
        """
        conc_sucrose = sucrose / mstruct
        #: Just after leaf initiation
        if time_prev_leaf2_emergence == 0 or prev_mean_conc_sucrose is None:
            new_integral_conc_sucrose = conc_sucrose
        #: Else
        else:
            new_integral_conc_sucrose = (
                                                    prev_mean_conc_sucrose * time_prev_leaf2_emergence + conc_sucrose * delta_teq) / (
                                                    time_prev_leaf2_emergence + delta_teq)
        return new_integral_conc_sucrose

    def calculate_leaf_Wmax(self, lamina_Lmax, leaf_rank, integral_conc_sucr, optimal_growth_option=False):
        """ Maximal lamina width.

        :param float lamina_Lmax: Maximal lamina length (m)
        :param int leaf_rank: leaf phytomer number
        :param float integral_conc_sucr: mean sucrose concentration of the hiddenzone since leaf n-2 emergence (µmol C g-1)
        :param bool optimal_growth_option: if True the function will calculate leaf Wmax assuming optimal growth conditions

        :return: Maximal leaf width (m)
        :rtype: float
        """
        if optimal_growth_option or leaf_rank in (1, 2):
            Wmax = self.parameters.leaf_Wmax_dict[leaf_rank]

        else:
            #: Width:length ratio
            W_L_ratio = max(self.parameters.leaf_W_L_MIN,
                            self.parameters.leaf_W_L_a - (self.parameters.leaf_W_L_b / self.parameters.leaf_W_L_c) * (
                                        1 - exp(-self.parameters.leaf_W_L_c * integral_conc_sucr)))

            #: Maximal width (m)
            Wmax = lamina_Lmax * W_L_ratio
        return Wmax

    def calculate_SSLW(self, leaf_rank, integral_conc_sucr, optimal_growth_option=False):
        """ Structural Specific Lamina Weight.

        :param int leaf_rank: leaf phytomer number
        :param float integral_conc_sucr: mean sucrose concentration of the hiddenzone since leaf n-2 emergence (µmol C g-1)
        :param bool optimal_growth_option: if True the function will calculate SSLW assuming optimal growth conditions


        :return: Structural Specific Leaf Weight (g m-2)
        :rtype: float
        """
        SSLW_min = self.parameters.leaf_SSLW_MIN
        SSLW_max = self.parameters.leaf_SSLW_MAX

        if optimal_growth_option or leaf_rank in (1, 2):
            SSLW = self.parameters.leaf_SSLW[leaf_rank]
        else:
            SSLW = (self.parameters.leaf_SSLW_a * integral_conc_sucr) / (
                        self.parameters.leaf_SSLW_b + integral_conc_sucr)
        return max(min(SSLW, SSLW_max), SSLW_min)

    def calculate_LSSW(self, leaf_rank, integral_conc_sucr, optimal_growth_option=False):
        """ Lineic Structural Sheath Weight.

        :param int leaf_rank: leaf phytomer number
        :param float integral_conc_sucr: mean sucrose concentration of the hiddenzone since leaf n-2 emergence (µmol C g-1)
        :param bool optimal_growth_option: if True the function will calculate LSLW assuming optimal growth conditions

        :return: Lineic Structural Sheath Weight (g m-1)
        :rtype: float
        """
        if optimal_growth_option or leaf_rank in (1, 2):
            LSSW = self.parameters.leaf_LSSW_dict[leaf_rank]
        else:
            LSSW_nominal = self.parameters.leaf_LSSW_nominal_A * leaf_rank + self.parameters.leaf_LSSW_nominal_B
            LSSW = LSSW_nominal + self.parameters.leaf_LSSW_a * (
                        integral_conc_sucr - self.parameters.leaf_LSSW_integral_min)

        return max(min(LSSW, self.parameters.leaf_LSSW_MAX), self.parameters.leaf_LSSW_MIN)

    @staticmethod
    def calculate_emerged_sheath_L(leaf_L, leaf_pseudostem_length, lamina_L, sheath_Lmax):
        """ Emerged sheath length. Assumes that leaf_L = leaf_pseudostem_length + emerged_sheath_L + lamina_L

        :param float leaf_L: Total leaf length (m)
        :param float leaf_pseudostem_length: Length of the pseudostem (m)
        :param float lamina_L: Lamina length (m)
        :param float sheath_Lmax: Final sheath length (m)

        :return: Sheath length (m)
        :rtype: float
        """
        return max(min(leaf_L - leaf_pseudostem_length - lamina_L, sheath_Lmax), 0.)

    @staticmethod
    def calculate_hidden_lamina_L(lamina_L, lamina_Lmax):
        """ Hidden lamina length at the end of lamina growth.

        :param float lamina_L: Length of the emerged lamina at the end of its growth (m)
        :param float lamina_Lmax: Final lamina length (m)

        :return: Hidden lamina length (m)
        :rtype: float
        """
        return max(lamina_Lmax - lamina_L, 0.)

    # -------------------------------------------------------------------------------------------------------------------
    # --- Internodes
    # -------------------------------------------------------------------------------------------------------------------
    @staticmethod
    def calculate_cumulated_internode_length(internode_lengths):
        """ Calculate cumulative internode length

        :param list internode_lengths: list of internode lengths below a given phytomer index (m)

        :return: Cumulated internode length below a given phytomer index (m)
        :rtype: float
        """

        return sum(internode_lengths)

    @staticmethod
    def calculate_internode_distance_to_emerge(ligule_heights, bottom_hiddenzone_height, phytomer_rank,
                                               curr_internode_L):
        """ Distance for the internode to be visible given by the pseudostem length.

        :param pandas.Series ligule_heights: height of all ligules of a given axis (m)
        :param float bottom_hiddenzone_height: height of the bottom of the hidden zone given by the cumulated lengths of below internodes (m)
        :param int phytomer_rank: phytomer rank
        :param float curr_internode_L: Internode length of the current phytomer (m).

        :return: Distance for the internode to be visible (m)
        :rtype: float
        """

        top_ligule_height = float(
            max(ligule_heights[ligule_heights['phytomer'] < phytomer_rank]['ligule height']))  # highest previous ligule
        leaf_pseudostem_length = top_ligule_height - bottom_hiddenzone_height

        internode_distance_to_emerge = leaf_pseudostem_length + curr_internode_L

        return max(0, internode_distance_to_emerge)

    def calculate_internode_Lmax(self, internode_L_lig):
        """ Final internode length.

        :param float internode_L_lig: Internode length at the ligulation of the previous leaf (m)

        :return: Final internode length (m)
        :rtype: float
        """

        internode_Lmax = internode_L_lig / self.Beta_function_internode(0.)

        return internode_Lmax

    def calculate_LSIW(self, LSSW, phytomer_rank, optimal_growth_option=False):
        """ Lineic Structural Internode Weight.

        :param float LSSW: Lineic Structural Sheath Weight (g m-1).
        :param int phytomer_rank: phytomer rank
        :param bool optimal_growth_option: if True the function will calculate LSIW assuming optimal growth conditions


        :return: Lineic Structural Internode Weight (g m-1)
        :rtype: float
        """
        if optimal_growth_option:
            LSIW = self.parameters.internode_LSIW_dict.get(phytomer_rank, self.parameters.internode_LSIW_dict[
                max(self.parameters.internode_LSIW_dict.keys())])
        else:
            LSIW = LSSW * self.parameters.ratio_LSIW_LSSW  # TODO : change this calculation as the ratio is not really stable according to the phytomer
        return LSIW

    def calculate_init_internode_elongation(self, hiddenzone_age, leaf_rank):
        """Initialize internode elongation.

        :param float hiddenzone_age: Time since primordium initiation (in time equivalent to a reference temperature) (s)
        :param int leaf_rank: rank of the leaf

        :return: Specifies if the internode has started the elongation (True) or not (False), and initialize internode_L
        :rtype: (bool, float)
        """
        if leaf_rank in (1, 2, 3):
            if hiddenzone_age > 5 * self.parameters.PLASTOCHRONE:  #6
                is_growing = True
                internode_L = self.parameters.internode_L_init
            else:
                is_growing = False
                internode_L = 0
        else:
            if hiddenzone_age > self.parameters.nb_PLASTO_internode_init * self.parameters.PLASTOCHRONE:
                is_growing = True
                internode_L = self.parameters.internode_L_init
            else:
                is_growing = False
                internode_L = 0

        return is_growing, internode_L

    def calculate_delta_internode_L_preL(self, phytomer_rank, sucrose, internode_L, amino_acids, mstruct, delta_teq,
                                         optimal_growth_option=False):
        """ delta of internode length over delta_t as a function of sucrose and amino acids, from initiation to the ligulation of the previous leaf.

        :param int phytomer_rank: phytomer rank
        :param float sucrose: Amount of sucrose (µmol C)
        :param float internode_L: Total internode length (m)
        :param float amino_acids: Amount of amino acids (µmol N)
        :param float mstruct: Structural mass of the hidden zone(g)
        :param float delta_teq: Temperature-compensated time = time duration at a reference temperature (s)
        :param bool optimal_growth_option: if True the function will calculate delta of internode length assuming optimal growth conditions

        :return: delta delta_internode_L (m)
        :rtype: float
        """

        conc_sucrose_effective = max(0., sucrose / mstruct - self.parameters.conc_sucrose_offset)

        if optimal_growth_option:
            RER_max = self.parameters.RERmax_dict_IN.get(phytomer_rank, self.parameters.RERmax_dict_IN[
                max(self.parameters.RERmax_dict_IN.keys())])
            delta_internode_L = internode_L * RER_max * delta_teq

        else:
            if conc_sucrose_effective > 0 and amino_acids > 0:
                RER_max = self.parameters.RERmax_dict_IN.get(phytomer_rank, self.parameters.RERmax_dict_IN[
                    max(self.parameters.RERmax_dict_IN.keys())])
                # Enzymatic rate for bi-substrates with random fixation
                conc_amino_acids = (amino_acids / mstruct)
                delta_internode_L = internode_L * RER_max * delta_teq / (
                            1 + self.parameters.RER_Kc / conc_sucrose_effective) / (
                                                1 + self.parameters.RER_Kn / conc_amino_acids)
            else:
                delta_internode_L = 0

        return delta_internode_L

    @staticmethod
    def calculate_internode_pseudo_age(internode_pseudo_age, delta_teq):
        """ Pseudo age of the internode since beginning of automate elongation (s)

        :param float internode_pseudo_age: Pseudo age of the leaf since beginning of automate elongation (s)
        :param float delta_teq: Temperature-compensated time = time duration at a reference temperature (s)

        :return: internode_pseudo_age (s)
        :rtype: float
        """
        return internode_pseudo_age + delta_teq

    def calculate_short_internode_Lmax(self, internode_L_lig, internode_pseudo_age):
        """ Final internode length.

        :param float internode_L_lig: Internode length at the ligulation of the previous leaf (m)
        :param float internode_pseudo_age: Internode pseudo age since previous leaf ligulation (s)

        :return: Final internode length (m)
        :rtype: float
        """

        L0 = 1 / self.Beta_function_internode(
            internode_pseudo_age)  #: Initial relative length of the short internode according to its pseudo age
        internode_Lmax = internode_L_lig * L0

        return internode_Lmax

    def Beta_function_internode(self, internode_pseudo_age):
        """ Normalized internode length from the emergence of the previous leaf to the end of elongation (automate function depending on internode pseudo age).

        :param float internode_pseudo_age: Pseudo age of the leaf since beginning of automate elongation (s)

        :return: normalized internode_L (m)
        :rtype: float
        """
        return (abs((1 + (max(0, (self.parameters.te_IN - internode_pseudo_age)) / (
                    self.parameters.te_IN - self.parameters.tm_IN))) *
                    (min(1.0, float(internode_pseudo_age - self.parameters.tb_IN) /
                         float(self.parameters.te_IN - self.parameters.tb_IN)) ** (
                                 (self.parameters.te_IN - self.parameters.tb_IN) /
                                 (self.parameters.te_IN - self.parameters.tm_IN)))))

    def calculate_delta_internode_L_postL(self, prev_internode_pseudo_age, internode_pseudo_age, prev_internode_L,
                                          internode_Lmax_lig, sucrose, amino_acids, mstruct,
                                          optimal_growth_option=False):
        """ Internode length, from the ligulation of the previous leaf to the end of elongation (automate function depending on leaf pseudo age and final length).

        :param float prev_internode_pseudo_age: Pseudo age of the internode since beginning of automate elongation at previous time step (s)
        :param float internode_pseudo_age: Pseudo age of the internode since beginning of automate elongation (s)
        :param float prev_internode_L: Internode length before elongation (m)
        :param float internode_Lmax_lig: Estimate of final internode length at previous leaf ligulation (m)
        :param float sucrose: Amount of sucrose (µmol C)
        :param float amino_acids: Amount of amino acids (µmol N)
        :param float mstruct: Structural mass (µmol N)
        :param bool optimal_growth_option: if True the function will calculate delta of internode length assuming optimal growth conditions

        :return: internode_L (m)
        :rtype: float
        """
        conc_sucrose_effective = max(0., sucrose / mstruct - self.parameters.conc_sucrose_offset)

        if internode_pseudo_age <= self.parameters.tb_IN and conc_sucrose_effective > 0 and amino_acids > 0:
            delta_internode_L = prev_internode_L - self.Beta_function_internode(0.) * internode_Lmax_lig
        elif internode_pseudo_age < self.parameters.te_IN:
            # Beta function
            delta_internode_L_Beta_0 = min(internode_Lmax_lig, internode_Lmax_lig * (
                        self.Beta_function_internode(internode_pseudo_age) - self.Beta_function_internode(
                    prev_internode_pseudo_age)))

            if optimal_growth_option:
                # Current internode length
                delta_internode_L = delta_internode_L_Beta_0
            elif conc_sucrose_effective > 0 and amino_acids > 0:  # TODO: not tested yet
                # Regulation by C and N
                conc_amino_acids = amino_acids / mstruct
                regul = self.parameters.leaf_pseudo_age_Vmax / (
                            1 + self.parameters.leaf_pseudo_age_Kc / conc_sucrose_effective) / (
                                    1 + self.parameters.leaf_pseudo_age_Kn / conc_amino_acids)
                # Current internode length
                delta_internode_L = regul * delta_internode_L_Beta_0
            else:
                delta_internode_L = 0
        else:
            delta_internode_L = 0

        return delta_internode_L

    def calculate_update_internode_Lmax(self, internode_Lmax_lig, internode_L, internode_pseudo_age):
        """ Update internode_Lmax following a reduction of delta_leaf_L due to C and N regulation

        :param float internode_Lmax_lig: Estimate of final internode length at previous leaf ligulation (m)
        :param float internode_L: actual internode length at the end of the time step, calculated according to CN concentration  (m)
        :param float internode_pseudo_age: Pseudo age of the internode since beginning of automate elongation at the end of the time step (s)

        :return: Updated internode Lmax (m)
        :rtype: float
        """
        return internode_L + internode_Lmax_lig * (1 - self.Beta_function_internode(internode_pseudo_age))

    @staticmethod
    def calculate_internode_visibility(internode_L, internode_distance_to_emerge):
        """Calculate if a given internode is visible

        :param float internode_L: Total internode length (m)
        :param float internode_distance_to_emerge: Length of the pseudostem + internode length (m)

        :return: Specifies if the internode is visible (True) or not (False)
        :rtype: float
        """
        return internode_L > internode_distance_to_emerge

    @staticmethod
    def calculate_emerged_internode_L(internode_L, internode_distance_to_emerge):
        """ Emerged internode length.

        :param float internode_L: Total internode length (m)
        :param float internode_distance_to_emerge: Length of the pseudostem + internode length (m)

        :return: Emerged internode length (m)
        :rtype: float
        """
        return internode_L - internode_distance_to_emerge

    def calculate_end_internode_elongation(self, internode_L, internode_Lmax, internode_pseudo_age):
        """Calculate if a given internode has finished elongating

        :param float internode_L: Total internode length (m)
        :param float internode_Lmax: Maximum internode length (m)
        :param float internode_pseudo_age: Internode pseudo age (s)

        :return: Specifies if the internode has completed elongation (True) or not (False)
        :rtype: float
        """
        condition_A = (
                    internode_pseudo_age >= self.parameters.te_IN)  # test if pseudo age is higher than the maximal duration of IN elongation
        condition_B = False  # test if IN length is longer than maximal length
        # condition_C = (internode_L == self.parameters.internode_L_init)  # test if IN length is still equal to the initial length (ie IN has not grown)
        if internode_Lmax:
            condition_B = (internode_L >= internode_Lmax)
        return condition_A or condition_B  #or condition_C

    # -------------------
    # --- COLEOPTILES ---
    # -------------------
    def calculate_coleo_Lmax(self, sowing_depth):
        """Calculate final length of coleoptile from sowing depth

        :param float sowing_depth: Sowing depth (m)

        :return: final length of coleoptile
        :rtype: float
        """
        return sowing_depth + self.parameters.delta_L_emergence

    @staticmethod
    def calculate_coleoptile_emergence(coleoptile_L, sowing_depth):
        """Calculate if a given coleoptile has emerged above the soil

        :param float coleoptile_L: Coleoptile length (m)
        :param float sowing_depth: Sowing depth (m)

        :return: Specifies if the coleoptile has emerged (True) or not (False)
        :rtype: bool
        """
        return coleoptile_L > sowing_depth

    @staticmethod
    def calculate_emerged_coleo_L(coleoptile_L, sowing_depth):
        """Calculate coleoptile length emerged above the soil

        :param float coleoptile_L: Coleoptile length (m)
        :param float sowing_depth: Sowing depth (m)

        :return: Visible coleoptile length (m)
        :rtype: float
        """
        return coleoptile_L - sowing_depth


class MorphogenesisModelHydraulics(MorphogenesisModel):
    def __init__(self):
        super().__init__()
        self.parameters = morphogenesis_parameters.parameters_factory(hydraulics=True)
        self.hiddenzone_init = morphogenesis_parameters.HiddenZoneInitHydraulics()

        variables_to_delete = {'leaf_Lmax', 'leaf_Lmax_em', 'leaf_Wmax', 'Wmax'}  # Not needed in the hydraulics version

        #: the inputs needed by Morphogenesis
        # -: Axes
        self.axis_inputs = list(self.axis_inputs) + ['xylem_water_potential']
        #-: Hiddenzones
        self.hiddenzone_inputs = [var for var in self.hiddenzone_inputs if var not in variables_to_delete]
        self.hiddenzone_inputs = list(self.hiddenzone_inputs) + ['width', 'thickness',
                                                                 'turgor_water_potential']
        # -: Elements
        self.element_inputs = [var for var in self.element_inputs if var not in variables_to_delete]


        #: the outputs computed by Morphogenesis
        # -: Hiddenzones
        self.hiddenzone_outputs = [var for var in self.hiddenzone_outputs if var not in variables_to_delete]
        self.hiddenzone_outputs = list(self.hiddenzone_outputs) + ['init_leaf_L']
        # -: Elements
        self.element_outputs = [var for var in self.element_outputs if var not in variables_to_delete]
        self.element_outputs = list(self.element_outputs) + ['width', 'thickness'] #todo voir avec VA


    def calculate_deltaL_preE(self, sucrose, leaf_L, amino_acids, mstruct, delta_teq, leaf_rank, optimal_growth_option,
                              xylem_water_potential=0.):
        """ Delta of leaf length over delta_t as a function of sucrose and amino acids, from initiation to the emergence of the previous leaf.

        :param float sucrose: Amount of sucrose (µmol C)
        :param float leaf_L: Total leaf length (m)
        :param float amino_acids: Amount of amino acids (µmol N)
        :param float mstruct: Structural mass (g)
        :param float delta_teq: Temperature-compensated time = time duration at a reference temperature (s)
        :param int leaf_rank: leaf phytomer number
        :param bool optimal_growth_option: if True the function will calculate leaf elongation assuming optimal growth conditions (except if sucrose and amino acids are zero)
        :param float xylem_water_potential: xylem water potential (Mpa)

        :return: delta delta_leaf_L (m)
        :rtype: float
        """

        conc_sucrose_effective = max(0., sucrose / mstruct - self.parameters.conc_sucrose_offset)

        if optimal_growth_option:
            RER_max = self.parameters.RERmax_Ljutovac_fit.get(leaf_rank, self.parameters.RERmax_Ljutovac_fit[
                max(self.parameters.RERmax_Ljutovac_fit.keys())])
            delta_leaf_L = leaf_L * RER_max * delta_teq
        elif conc_sucrose_effective > 0 and amino_acids > 0:
            RER_max = self.parameters.RERmax.get(leaf_rank, self.parameters.RERmax[max(self.parameters.RERmax.keys())])
            # Enzymatic rate for bi-substrates with random fixation
            conc_amino_acids = (amino_acids / mstruct)
            regul_water = 1 / (1 + (
                        xylem_water_potential / self.parameters.psi_ref) ** self.parameters.n)  # Sigmoidal function of xylem water potential control of RER
            delta_leaf_L = delta_teq * leaf_L * RER_max / ((1 + self.parameters.RER_Kc / conc_sucrose_effective) * (
                        1 + self.parameters.RER_Kn / conc_amino_acids)) * regul_water
        else:
            delta_leaf_L = 0

        return delta_leaf_L

    def calculate_leaf_pseudo_age(self, leaf_pseudo_age, delta_teq, turgor_water_potential=0):
        """ Pseudo age of the leaf since beginning of automate elongation (s)

        :param float leaf_pseudo_age: Previous pseudo age of the leaf since beginning of automate elongation (s)
        :param float delta_teq: Temperature-compensated time = time duration at a reference temperature (s)
        :param float turgor_water_potential: turgor potential of the hz (MPa)

        :return: Updated leaf_pseudo_age (s)
        :rtype: float
        """
        if turgor_water_potential <= self.parameters.GAMMA:
            return leaf_pseudo_age
        else:
            return leaf_pseudo_age + delta_teq

    def calculate_ratio_DZ_postE(self, leaf_L, leaf_pseudo_age, leaf_pseudostem_length):

        """ Ratio of the hiddenzone length which is made of division zone.
        Calculated from an inverse beta function representing the ratio of the leaf composed by the division zone according to its relative length in log.
        The model was fitted on literature data on wheat.

        :param float leaf_L: Leaf length (m)
        :param float leaf_pseudostem_length: Length of the pseudostem (m)
        :param float leaf_pseudo_age: Pseudo age of the leaf since beginning of automate elongation (s)

        :return: ratio_DZ (dimensionless)
        :rtype: float
        """
        log_L_init = log10(self.parameters.ratio_DZ_l_init)
        log_L_mid = log10(self.parameters.ratio_DZ_l_mid)
        log_Lend = log10(self.parameters.ratio_DZ_l_end)

        age_normalised = leaf_pseudo_age / self.parameters.te
        log_age_normalised = log10(age_normalised)

        if log_age_normalised >= log_Lend:
            return 0
        elif log_age_normalised <= log_L_init:
            return 1
        else:
            return min((1 - (1 + (log_Lend - log_age_normalised) / (log_Lend - log_L_mid)) *
                        (((log_age_normalised - log_L_init) / (log_Lend - log_L_init)) ** (
                                (log_Lend - log_L_init) / (log_Lend - log_L_mid)))) * leaf_L / min(leaf_L,
                                                                                                   leaf_pseudostem_length),
                       1.)

    def calculate_lamina_Lmax(self, leaf_rank):

        """ Final lamina length.

        :param int leaf_rank: rank of the leaf

        :return: "Dummy" final lamina length (m)
        :rtype: float
        """
        lamina_Lmax = self.parameters.lamina_Lmax_dict.get(leaf_rank, self.parameters.lamina_Lmax_dict[max(self.parameters.lamina_Lmax_dict.keys())])
        if leaf_rank > 11:
            lamina_Lmax = 0.45
        return lamina_Lmax