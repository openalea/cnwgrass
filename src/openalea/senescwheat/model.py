# -*- coding: latin-1 -*-

from __future__ import division  # use '//' to do integer division
from senescwheat import parameters

"""
    senescwheat.model
    ~~~~~~~~~~~~~~~~~~~

    Model of senescence.

    :copyright: Copyright 2014-2015 INRA-ECOSYS, see AUTHORS.
    :license: see LICENSE for details.

"""


class SenescenceModel(object):

    @classmethod
    def calculate_N_content_total(cls, proteins, amino_acids, nitrates, Nstruct, max_mstruct, Nresidual):
        """ N content in the whole element (both green and senesced tissues).

        :param float proteins: protein concentration (�mol N proteins g-1 mstruct)
        :param float amino_acids: amino acids concentration (�mol N amino acids g-1 mstruct)
        :param float nitrates: nitrates concentration (�mol N nitrates g-1 mstruct)
        :param float Nstruct: structural N mass (g). Should be constant during leaf life.
        :param float max_mstruct: structural mass maximal of the element i.e. structural mass of the whole element before senescence (g)
        :param float Nresidual: residual mass of N in the senescent tissu (g)

        :return: N_content_total (between 0 and 1)
        :rtype: float
        """
        return ((proteins + amino_acids + nitrates) * 1E-6 * parameters.N_MOLAR_MASS + Nresidual + Nstruct) / max_mstruct

    @classmethod
    def calculate_forced_relative_delta_green_area(cls, green_area_df, group_id, prev_green_area):
        """relative green_area variation due to senescence

        :param pandas.DataFrame green_area_df: a pandas DataFrame containing the green area values for each photosynthetic element at each time
        :param tuple group_id: the group id to be used to select data in the DataFrame
        :param float prev_green_area: previous value of an organ green area (m-2)

        :return: new_green_area (m-2), relative_delta_green_area (dimensionless)
        :rtype: tuple [float, float]
        """
        new_green_area = green_area_df.get_group(group_id).green_area.values[0]
        relative_delta_green_area = (prev_green_area - new_green_area) / prev_green_area
        return new_green_area, relative_delta_green_area

    @classmethod
    def calculate_relative_delta_green_area(cls, organ_name, prev_green_area, proteins, max_proteins, delta_t, update_max_protein):
        """relative green_area variation due to senescence

        :param str organ_name: name of the organ to which belongs the element (used to distinguish lamina from stem organs)
        :param float prev_green_area: previous value of an organ green area (m-2)
        :param float proteins: protein concentration (�mol N proteins g-1 mstruct)
        :param float max_proteins: maximal protein concentrations experienced by the organ (�mol N proteins g-1 mstruct)
        :param float delta_t: value of the timestep (s)
        :param bool update_max_protein: whether to update the max proteins or not.

        :return: new_green_area (m-2), relative_delta_green_area (dimensionless)
        :rtype: tuple [float, float]

        .. todo:: remove update_max_protein
        """

        if organ_name == 'blade':
            fraction_N_max = parameters.FRACTION_N_MAX['blade']
        else:
            fraction_N_max = parameters.FRACTION_N_MAX['stem']

        # Overwrite max proteins
        if max_proteins < proteins and update_max_protein:
            max_proteins = proteins
            new_green_area = prev_green_area
            relative_delta_green_area = 0
        # Senescence if (actual proteins/max_proteins) < fraction_N_max
        elif max_proteins == 0 or (proteins / max_proteins) < fraction_N_max:
            senesced_area = min(prev_green_area, parameters.SENESCENCE_MAX_RATE * delta_t)
            new_green_area = max(0., prev_green_area - senesced_area)
            relative_delta_green_area = senesced_area / prev_green_area
        else:
            new_green_area = prev_green_area
            relative_delta_green_area = 0
        return new_green_area, relative_delta_green_area, max_proteins

    # Temporaire
    @classmethod
    def calculate_relative_delta_senesced_length(cls, organ_name, prev_senesced_length, length, proteins, max_proteins, delta_t, update_max_protein):
        """relative senesced length variation

        :param str organ_name: name of the organ to which belongs the element (used to distinguish lamina from stem organs)
        :param float prev_senesced_length: previous senesced length of an organ (m-2)
        :param float length: organ length (m)
        :param float proteins: protein concentration (�mol N proteins g-1 mstruct)
        :param float max_proteins: maximal protein concentrations experienced by the organ (�mol N proteins g-1 mstruct)
        :param float delta_t: value of the timestep (s)
        :param bool update_max_protein: whether to update the max proteins or not.

        :return: new_senesced_length (m), relative_delta_senesced_length (dimensionless), max_proteins (�mol N proteins g-1 mstruct)
        :rtype: tuple [float, float, float]
        
        .. todo:: remove update_max_protein
        """

        if organ_name == 'blade':
            fraction_N_max = parameters.FRACTION_N_MAX['blade']
        else:
            fraction_N_max = parameters.FRACTION_N_MAX['stem']

        # Overwrite max proteins
        if max_proteins < proteins and update_max_protein:
            max_proteins = proteins
            new_senesced_length = prev_senesced_length
            relative_delta_senesced_length = 0
        # Senescence if (actual proteins/max_proteins) < fraction_N_max
        elif max_proteins == 0 or (proteins / max_proteins) < fraction_N_max:
            senesced_length = parameters.SENESCENCE_LENGTH_MAX_RATE * delta_t
            new_senesced_length = min(length, prev_senesced_length + senesced_length)
            if length == new_senesced_length:
                relative_delta_senesced_length = 1
            else:
                relative_delta_senesced_length = 1 - (length - new_senesced_length) / (length - prev_senesced_length)
        else:
            new_senesced_length = prev_senesced_length
            relative_delta_senesced_length = 0
        return new_senesced_length, relative_delta_senesced_length, max_proteins

    @classmethod
    def calculate_delta_mstruct_shoot(cls, relative_delta_green_area, prev_mstruct, prev_Nstruct):
        """delta of structural mass due to senescence of photosynthetic elements

        :param float relative_delta_green_area: relative variation of a photosynthetic element green area (dimensionless)
        :param float prev_mstruct: previous value of an organ structural mass (g)
        :param float prev_Nstruct: previous value of an organ structural N (g)

        :return: delta_mstruct (g), delta_Nstruct (g)
        :rtype: tuple [float, float]
        """
        delta_mstruct = prev_mstruct * relative_delta_green_area
        delta_Nstruct = prev_Nstruct * relative_delta_green_area
        return delta_mstruct, delta_Nstruct

    @classmethod
    def calculate_remobilisation(cls, metabolite, relative_delta_structure):
        """Metabolite remobilisation due to senescence over DELTA_T (�mol).
        
        :param float metabolite: amount of any metabolite to be remobilised (�mol) 
        :param float relative_delta_structure: could be relative variation of a photosynthetic element green area or relative variation of mstruct
        
        :return: metabolite remobilisation (�mol)
        :rtype: float
        """
        return metabolite * relative_delta_structure

    @classmethod
    def calculate_if_element_is_over(cls, green_area, is_growing, mstruct):
        """Define is an element is fully senescent

        :param float green_area: Green area of the element (m2)
        :param bool is_growing: flag is the element is still growing
        :param float mstruct: Strucural mass of the element (g)

        :return: is_over which indicates if the element is fully senescent
        :rtype: bool
        """
        is_over = False
        if (green_area < parameters.MIN_GREEN_AREA or mstruct == 0) and not is_growing:
            is_over = True
        return is_over

    @classmethod
    def calculate_remobilisation_proteins(cls, organ, element_index, proteins, relative_delta_green_area, ratio_N_mstruct_max, full_remob):
        """Protein remobilisation due to senescence over DELTA_T. Part is remobilized as amino_acids (�mol N), the rest is increasing Nresidual (g).
        
        :param str organ: name of the organ
        :param int element_index: phytomer rank
        :param float proteins: amount of proteins (�mol N)
        :param float relative_delta_green_area: relative variation of a photosynthetic element green area
        :param float ratio_N_mstruct_max: N content in the whole element (both green and senesced tissues).
        :param bool full_remob: whether all proteins should be remobilised
        
        :return: Quantity of proteins remobilised either in amino acids, either in residual N (�mol),
                 Quantity of proteins converted into amino_acids (�mol N), 
                 Increment of Nresidual (g)
        :rtype: tuple [float, float, float]
        """

        if full_remob or organ != 'blade':
            remob_proteins = delta_amino_acids = proteins * relative_delta_green_area
            delta_Nresidual = 0
        else:
            if ratio_N_mstruct_max <= parameters.RATIO_N_MSTRUCT.get(element_index, parameters.DEFAULT_RATIO_N_MSTRUCT):  # all the proteins are converted into Nresidual
                remob_proteins = proteins
                delta_Nresidual = remob_proteins * 1E-6 * parameters.N_MOLAR_MASS
                delta_amino_acids = 0
            else:  # part of the proteins are converted into amino_acids
                remob_proteins = proteins * relative_delta_green_area
                delta_amino_acids = remob_proteins
                delta_Nresidual = 0
        return remob_proteins, delta_amino_acids, delta_Nresidual

    @classmethod
    def calculate_roots_senescence(cls, mstruct, Nstruct, postflowering_stages):
        """Root senescence
        :param float mstruct: structural mass (g)
        :param float Nstruct: structural N (g)
        :return: Rate of mstruct loss by root senescence (g mstruct s-1), rate of Nstruct loss by root senescence (g Nstruct s-1)
        :rtype: tuple [float, float]
        """
        if postflowering_stages:
            rate_senescence = parameters.SENESCENCE_ROOTS_POSTFLOWERING
        else:
            rate_senescence = parameters.SENESCENCE_ROOTS_PREFLOWERING
        return mstruct * rate_senescence, Nstruct * rate_senescence

    @classmethod
    def calculate_relative_delta_mstruct_roots(cls, rate_mstruct_death, root_mstruct, delta_teq):
        """Relative delta of root structural dry matter (g) over delta_t

        :param float rate_mstruct_death: Rate of mstruct loss by root senescence (g mstruct s-1)
        :param float root_mstruct: actual mstruct of roots (g)
        :param float delta_teq: Temperature-consensated time = time duration at a reference temperature (s)

        :return: relative_delta_mstruct (dimensionless)
        :rtype: float
        """
        return (rate_mstruct_death * delta_teq) / root_mstruct

    @classmethod
    def calculate_delta_mstruct_root(cls, rate_mstruct_death, rate_Nstruct_death, delta_teq):
        """delta of structural mass due to senescence of roots

        :param float rate_mstruct_death: relative delta of root structural mass over delta_t (g s-1)
        :param float rate_Nstruct_death: relative delta of root N structural mass over delta_t (g s-1)
        :param float delta_teq: Temperature-consensated time = time duration at a reference temperature (s)

        :return: delta_mstruct (g), delta_Nstruct (g)
        :rtype: tuple [float, float]
        """
        delta_mstruct = rate_mstruct_death * delta_teq
        delta_Nstruct = rate_Nstruct_death * delta_teq

        return delta_mstruct, delta_Nstruct
