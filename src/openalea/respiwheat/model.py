# -*- coding: latin-1 -*-

from __future__ import division  # use '//' to do integer division

"""
    respiwheat.model
    ~~~~~~~~~~~~~~~~~~~

    Model of respiration based on Thornley and Cannell, 2000.
    The model computes the respiration associated with the main biological processes.

    R_total = sum(R_growth) + R_phloem + R_Namm_upt + R_Nnit_upt + R_Nnit_red(shoot + root) + R_N2fix + R_min_upt + sum(R_residual)

    :copyright: Copyright 2014-2017 INRA-ECOSYS, see AUTHORS.
    :license: CeCILL-C, see LICENSE for details.
    
    **Acknowledgments**: The research leading these results has received funding through the 
    Investment for the Future programme managed by the Research National Agency 
    (BreedWheat project ANR-10-BTBR-03).
    
    .. seealso:: Barillot et al. 2016.
"""


class RespirationModel(object):
    SECOND_TO_HOUR_RATE_CONVERSION = 3600

    # R_growth#
    YG = 0.75  # 0.8               # Growth yield (units of C appearing in new biomass per unit of C substrate utilized for growth). Range 0±75 to 0±85 in Cannell and Thornley, 2000.

    YG_GRAINS = 0.71  # Growth yield (units of C appearing in new biomass per unit of C substrate utilized for growth)

    # R_phloem#
    CPHLOEM = 0.006  # Units C respired per unit C substrate loaded into the phloem

    # R_Namm_upt#
    C_AMM_UPT = 0.198  # µmol of C substrate respired per µmol of N ammonium taken up

    # R_Nnit_upt#
    C_NIT_UPT = 0.397  # µmol of C substrate respired per µmol of N nitrates taken up

    # R_Nnit_red#
    F_NIT_RED_SH_CS = 0.5  # fraction of nitrate reduced in shoot using C substrate rather than using excess ATP and reducing power obtained directly from photosynthesis
    C_NIT_RED = 1.98  # µmol of C substrate per µmol of N nitrates reduced

    # R_N2fix #
    C_NFIX = 6  # kg substrate C respired (kg N fixed)-1 (in the range  5 to 12)

    # R_min_upt #
    CMIN_UPT = 5000  # µmol of C substrate respired per g of minerals taken up

    # R_residual #
    KM_MAX = 4.1E-6  # 8E-6 # 4.1E-6       # Maximum value of the maintenance constant when C is much greater than KM (µmol of C substrate respired per µmol N s-1)
    KM = 1.67E3  # The Michaelis-Menten constant affinity i.e. the C substrate concentration at half the value of KM_MAX (µmol of C substrate per g of structural mass)

    @classmethod
    def R_growth(cls, mstruct_growth):
        """ Local growth respiration

        :param float mstruct_growth: gross growth of mstruct (µmol C added in mstruct)

        :return: R_growth (µmol C respired)
        :rtype: float
        """
        R_growth = ((1 - cls.YG) / cls.YG) * mstruct_growth
        return R_growth

    @classmethod
    def R_grain_growth(cls, mstruct_growth, starch_filling, mstruct):
        """ Grain growth respiration

        :param float mstruct_growth: gross growth of grain structure (µmol C added in grain structure)
        :param float starch_filling: gross growth of grain starch (µmol C added in grain starch g-1 mstruct)
        :param float mstruct: structural dry mass of organ (g)

        :return: R_grain_growth_struct, R_grain_growth_starch (µmol C respired)
        :rtype: (float, float)
        """
        R_grain_growth_struct = ((1 - cls.YG_GRAINS) / cls.YG_GRAINS) * mstruct_growth
        R_grain_growth_starch = ((1 - cls.YG_GRAINS) / cls.YG_GRAINS) * (starch_filling * mstruct)
        return R_grain_growth_struct, R_grain_growth_starch

    @classmethod
    def R_phloem(cls, sucrose_loading, mstruct):
        """ Phloem loading respiration

        :param float sucrose_loading: Loading flux from the C substrate pool to phloem (µmol C g-1 mstruct)
        :param float mstruct: structural dry mass of organ (g)

        :return: R_phloem, sucrose_loading (µmol C respired, µmol C)
        :rtype: (float, float)
        """
        R_phloem = max(0., cls.CPHLOEM * sucrose_loading * mstruct)  #: Do not count a respiratory cost for negative loading i.e. unloading (assumed to be passive)
        return R_phloem, sucrose_loading

    @classmethod
    def R_Namm_upt(cls, U_Namm):
        """ Ammonium uptake respiration

        :param float U_Namm: uptake of N ammonium (µmol N)

        :return: R_Namm (µmol C respired)
        :rtype: float
        """
        R_Namm_upt = cls.C_AMM_UPT * U_Namm
        return R_Namm_upt

    @classmethod
    def R_Nnit_upt(cls, U_Nnit, sucrose):
        """ Nitrate uptake respiration

        :param float U_Nnit: uptake of N nitrates (µmol N)
        :param float sucrose: amount of C sucrose in organ (µmol C)

        :return: R_Nnit_upt (µmol C respired)
        :rtype: float
        """
        if sucrose > 0:
            R_Nnit_upt = cls.C_NIT_UPT * U_Nnit
        else:
            R_Nnit_upt = 0
        return R_Nnit_upt

    @classmethod
    def R_Nnit_red(cls, s_amino_acids, sucrose, mstruct, root=False):
        """ Nitrate reduction-linked respiration
        Distinction is made between nitrate realised in roots or in shoots where a part of the energy required is derived from ATP
        and reducing power obtained directly from photosynthesis (rather than C substrate)

        :param float s_amino_acids: consumption of N for the synthesis of amino acids (µmol N g-1 mstruct)
        (in the present version, this is used to approximate nitrate reduction needed in the original model of Thornley and Cannell, 2000)
        :param float sucrose: amount of C sucrose in organ (µmol C)
        :param float mstruct: structural dry mass of organ (g)
        :param bool root: specifies if the nitrate reduction-linked respiration is computed for shoot (False) or root (True) tissues.

        :return: R_Nnit_upt, s_amino_acids (µmol C respired, µmol N g-1 mstruct)
        :rtype: (float, float)
        """
        if not root:
            R_Nnit_red = cls.F_NIT_RED_SH_CS * cls.C_NIT_RED * s_amino_acids * mstruct  # Respiration in shoot tissues
        else:
            R_Nnit_red = cls.C_NIT_RED * s_amino_acids * mstruct  # Respiration in root tissues
            if sucrose < R_Nnit_red:
                R_Nnit_red = 0
                s_amino_acids = 0

        return R_Nnit_red, s_amino_acids

    @classmethod
    def R_N2fix(cls, I_Nfix):
        """ N2-fixation respiration

        :param float I_Nfix: flux of fixed N into the root substrate N pool (kg fixed N)

        :return: R_N2fix (µmol C respired)
        :rtype: float`
        """
        if I_Nfix <= 0:
            R_N2fix = 0
        else:
            R_N2fix = cls.C_NFIX * I_Nfix
        return R_N2fix

    @classmethod
    def R_min_upt(cls, delta_mineral_plant):
        """ Mineral ion (other than N) uptake-linked respiration

        :param float delta_mineral_plant: Uptake of mineral by the plant (g)

        :return: R_min_upt (µmol C respired)
        :rtype float
        """
        # Uptake of minerals (g)
        Umin = delta_mineral_plant
        # Respiratory cost
        R_min_upt = cls.CMIN_UPT * Umin
        return max(0., R_min_upt)

    @classmethod
    def R_residual(cls, sucrose, mstruct, Ntot, Ts):
        """ Residual respiration (cost from protein turn-over, cell ion gradients, futile cycles...)

        :param float sucrose: amount of C sucrose (µmol C)
        :param float mstruct: structural dry mass of organ (g)
        :param float Ntot: total N in organ (µmol N)
        :param float Ts : organ temperature (°C)

        :return: R_residual (µmol C respired h-1)
        :rtype: float
        """

        Q10 = 2.
        T_ref = 20.

        if sucrose <= 0. or mstruct <= 0.:
            R_residual = 0.
        else:
            conc_sucrose = sucrose / mstruct
            R_residual = ((cls.KM_MAX * conc_sucrose) / (cls.KM + conc_sucrose)) * Ntot * Q10 ** ((Ts - T_ref) / 10) * cls.SECOND_TO_HOUR_RATE_CONVERSION

        return R_residual
