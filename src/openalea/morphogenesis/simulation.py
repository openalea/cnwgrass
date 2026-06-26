# -*- coding: latin-1 -*-

from __future__ import division  # use "//" to do integer division

import warnings
import copy

import pandas as pd
import numpy as np

from openalea.morphogenesis import model

"""
    morphogenesis.simulation
    ~~~~~~~~~~~~~~~~~~

    The module :mod:`morphogenesis.simulation` is the front-end to run the Morphogenesis model.

"""


class SimulationError(Exception):
    pass


class SimulationRunError(SimulationError):
    pass


class Simulation(object):
    """The Simulation class allows to initialize and run a simulation.
    """

    def __init__(self, delta_t=1, hydraulics=False, optimal_growth_option=False, update_parameters=None):
        """
        :param int delta_t: the delta t of the simulation (in seconds)
        :param bool hydraulics: if True the model will assume the coupling to the turgor-driven growth model
        :param bool optimal_growth_option: if True the model will assume optimal growth conditions
        :param None or dict update_parameters: if a dict is provided, the specified parameters in keys will be updated
        """

        #: `inputs` is a dictionary of dictionaries:
        #:     {'hiddenzone': {(plant_index, axis_label, metamer_index): {hiddenzone_input_name: hiddenzone_input_value, ...}, ...},
        #:      'elements': {(plant_index, axis_label, metamer_index, organ_label, element): {element_input_name: element_input_value, ...}, ...},
        #:      'axes': {(plant_index, axis_label): {axis_input_name: axis_input_value, ...}, ...},
        #:      'sheath_internode_lengths': {(plant_index, axis_label, metamer_index): {'sheath': [list of sheath length belonging to the phytomer],
        #:                                                                             'cumulated_internode': [list of internode lengths cumulated from phytomer 1 to n]}, ...}}
        self.inputs = {}

        #: `outputs` is a dictionary of dictionaries:
        #:     {'hiddenzone': {(plant_index, axis_label, metamer_index): {hiddenzone_input_name: hiddenzone_input_value, ...}, ...},
        #:      'elements': {(plant_index, axis_label, metamer_index, organ_label, element): {element_output_name: element_output_value, ...}, ...},
        #:      'axes': {(plant_index, axis_label): {axis_output_name: axis_output_value, ...}, ...}}
        self.outputs = {}

        #: the delta t of the simulation (in seconds)
        self.delta_t = delta_t

        #: Checks whether the Hydraulic version should be used
        self.hydraulics = hydraulics
        if not self.hydraulics:
            self.model = model.MorphogenesisModel()
        else:
            self.model = model.MorphogenesisModelHydraulics()
        self.hiddenzone_init = self.model.hiddenzone_init
        self.element_init = self.model.element_init

        #: The inputs and outputs of Morphogenesis at each scale.
        self.axis_inputs = self.model.axis_inputs
        self.axis_outputs = self.model.axis_outputs
        self.axis_inputs_outputs = sorted(set(self.axis_inputs + self.axis_outputs))
        self.hiddenzone_inputs = self.model.hiddenzone_inputs
        self.hiddenzone_outputs = self.model.hiddenzone_outputs
        self.hiddenzone_inputs_outputs = sorted(set(self.hiddenzone_inputs + self.hiddenzone_outputs))
        self.element_inputs = self.model.element_inputs
        self.element_outputs = self.model.element_outputs
        self.element_inputs_outputs = sorted(set(self.element_inputs + self.element_outputs))

        #: Optimal growth option
        self.optimal_growth_option = optimal_growth_option

        #: Update parameters if specified
        if update_parameters:
            for key, value in update_parameters.items():
                if hasattr(self.model.parameters, key):
                    setattr(self.model.parameters, key, value)
                else:
                    warnings.warn(f"Parameter '{key}' is not defined in class self.model.parameters.")

    def initialize(self, inputs):
        """
        Initialize :attr:`inputs` from `inputs`.

        :param dict inputs: must be a dictionary with the same structure as :attr:`inputs`.
        """
        self.inputs.clear()
        self.inputs.update(inputs)

    def run(self, Tair, Tsoil, Zsowing=0.025):
        """
        Run the simulation.

        :param float Tair: Air temperature at t (degree Celsius)
        :param float Tsoil: Soil temperature at t (degree Celsius)
        :param float Zsowing: Sowing depth (m)
        """

        # Copy the inputs into the output dict
        self.outputs.update({inputs_type: copy.deepcopy(all_inputs) for inputs_type, all_inputs in self.inputs.items()
                             if inputs_type in {'hiddenzone', 'elements', 'axes', 'sheath_internode_lengths'}})

        # Hidden zones
        all_hiddenzone_inputs = self.inputs['hiddenzone']
        all_hiddenzone_outputs = self.outputs['hiddenzone']

        # Elements
        all_element_inputs = self.inputs['elements']
        all_element_outputs = self.outputs['elements']

        # Axes
        all_axes_inputs = self.inputs['axes']
        all_axes_outputs = self.outputs['axes']

        # sheath and internode lengths
        all_sheath_internode_lengths = self.inputs['sheath_internode_lengths']

        # Ligule heights
        all_ligule_height_df = pd.DataFrame(columns=self.model.ligule_topology_columns)

        # --  Beginning of computations

        # -----------------------------
        # ---------- SAM --------------
        # -----------------------------

        for axis_id, axis_inputs in sorted(all_axes_inputs.items()):
            curr_axis_outputs = all_axes_outputs[axis_id]
            nb_leaves = curr_axis_outputs['nb_leaves']

            # height of the SAM
            below_internode_lengths = all_sheath_internode_lengths[axis_id][nb_leaves]['cumulated_internode']
            SAM_height = self.model.calculate_cumulated_internode_length(below_internode_lengths)
            curr_axis_outputs['SAM_height'] = SAM_height

            # SAM temperature
            growth_temperature = self.model.calculate_growing_temperature(Tair, Tsoil, SAM_height, Zsowing)
            curr_axis_outputs['SAM_temperature'] = growth_temperature

            # temperature-compensated time
            curr_axis_outputs['delta_teq'] = self.model.calculate_time_equivalent_Tref(growth_temperature, self.delta_t)
            curr_axis_outputs['delta_teq_roots'] = self.model.calculate_time_equivalent_Tref(Tsoil, self.delta_t)

            # cumulated thermal time
            curr_axis_outputs['sum_TT'] = self.model.calculate_cumulated_thermal_time(curr_axis_outputs['sum_TT'], growth_temperature, curr_axis_outputs['delta_teq'])

            # update SAM status, leaf number and
            init_leaf, curr_axis_outputs['nb_leaves'], curr_axis_outputs['status'], curr_axis_outputs['teq_since_primordium'] = self.model.calculate_SAM_primodia(axis_inputs['status'],
                                                                                                                                                             curr_axis_outputs['teq_since_primordium'],
                                                                                                                                                             curr_axis_outputs['delta_teq'], nb_leaves,
                                                                                                                                                             curr_axis_outputs['cohort'])

            # GA production
            curr_axis_outputs['GA'] = self.model.calculate_SAM_GA(curr_axis_outputs['status'], curr_axis_outputs['teq_since_primordium'])

            # hiddenzone initiation
            for i in range(0, init_leaf):
                # Initialise hiddenzone
                hiddenzone_id = axis_id + tuple([1 + i + curr_axis_outputs['nb_leaves'] - init_leaf])  # TODO: peut etre simplifié tant que 'calculate_SAM_status' renvoie 1 erreur si init_leaf>1
                new_hiddenzone = self.hiddenzone_init.__dict__.copy()
                self.outputs['hiddenzone'][hiddenzone_id] = new_hiddenzone

            # Ligule height
            all_ligule_height_df = self.model.calculate_ligule_height(all_sheath_internode_lengths[axis_id], all_element_inputs, axis_id, all_ligule_height_df)

            self.outputs['axes'][axis_id] = curr_axis_outputs

        # -----------------------------
        # ---------- Elements ---------
        # -----------------------------

        for element_id, element_inputs in sorted(all_element_inputs.items()):
            # Update element's age, only used by ADEL to adapt element's geometry (so far lamina curvature, could be used to adapt stem geometry too)
            # TODO : the calculation of element's age must be extracted from morphogenesis as it is run even for mature leaves
            curr_age = all_element_inputs[element_id]['age']
            curr_age_teq = all_element_inputs[element_id]['age_teq']
            axis_id = element_id[:2]
            curr_axis_outputs = all_axes_outputs[axis_id]
            self.outputs['elements'][element_id]['age'] = self.model.calculate_cumulated_thermal_time(curr_age, curr_axis_outputs['SAM_temperature'], curr_axis_outputs['delta_teq'])
            self.outputs['elements'][element_id]['age_teq'] = self.model.calculate_element_age_teq(curr_age_teq, curr_axis_outputs['delta_teq'])

        # -----------------------------
        # ---------- Hiddenzones ------
        # -----------------------------
        for hiddenzone_id, hiddenzone_inputs in sorted(all_hiddenzone_inputs.items()):

            axe_label = hiddenzone_id[1]
            axis_id = hiddenzone_id[:2]
            phytomer_id = hiddenzone_id[2]

            if phytomer_id == 0:  #: coleoptile
                continue

            # Get xylem water potential if hydraulics
            xylem_water_potential = all_axes_inputs[axis_id].get('xylem_water_potential')  # Set at None if hydraulics is False
            hz_turgor_water_potential = hiddenzone_inputs.get('turgor_water_potential')    # Set at None if hydraulics is False

            #: Tillers: in this version tillers functioning is replicated from corresponding elements of MS
            if axe_label != 'MS':
                tiller_to_MS_phytomer_id = tuple([axis_id[0], 'MS', all_axes_outputs[axis_id]['cohort'] + phytomer_id - 1])

                if tiller_to_MS_phytomer_id in all_hiddenzone_outputs.keys():
                    self.outputs['hiddenzone'][hiddenzone_id] = all_hiddenzone_outputs[tiller_to_MS_phytomer_id]

                    if all_hiddenzone_outputs[tiller_to_MS_phytomer_id]['leaf_is_emerged']:
                        # Lamina
                        tiller_to_MS_lamina_id = tiller_to_MS_phytomer_id + tuple(['blade', 'LeafElement1'])
                        tiller_lamina_id = hiddenzone_id + tuple(['blade', 'LeafElement1'])
                        if tiller_to_MS_lamina_id in all_element_outputs.keys():
                            self.outputs['elements'][tiller_lamina_id] = all_element_outputs[tiller_to_MS_lamina_id]
                        # else:
                        #     warnings.warn('No leaf found on main stem for tiller {}.'.format(tiller_to_MS_lamina_id))

                        # Emerged Sheath
                        tiller_to_MS_emerged_sheath_id = tiller_to_MS_phytomer_id + tuple(['sheath', 'StemElement'])
                        tiller_emerged_sheath_id = hiddenzone_id + tuple(['sheath', 'StemElement'])
                        if tiller_to_MS_emerged_sheath_id in all_element_outputs.keys():
                            self.outputs['elements'][tiller_emerged_sheath_id] = all_element_outputs[tiller_to_MS_emerged_sheath_id]
                        # else:
                        #     warnings.warn('No emerged sheath found on main stem for tiller {}.'.format(tiller_to_MS_emerged_sheath_id))

                        # Enclosed Sheath
                        tiller_to_MS_enclosed_sheath_id = tiller_to_MS_phytomer_id + tuple(['sheath', 'HiddenElement'])
                        tiller_enclosed_sheath_id = hiddenzone_id + tuple(['sheath', 'HiddenElement'])
                        if tiller_to_MS_enclosed_sheath_id in all_element_outputs.keys():
                            self.outputs['elements'][tiller_enclosed_sheath_id] = all_element_outputs[tiller_to_MS_enclosed_sheath_id]
                        # else:
                        #     warnings.warn('No enclosed sheath found on main stem for tiller {}.'.format(tiller_to_MS_enclosed_sheath_id))

                        # Emerged internode
                        tiller_to_MS_emerged_internode_id = tiller_to_MS_phytomer_id + tuple(['internode', 'StemElement'])
                        tiller_emerged_internode_id = hiddenzone_id + tuple(['internode', 'StemElement'])
                        if tiller_to_MS_emerged_internode_id in all_element_outputs.keys():
                            self.outputs['elements'][tiller_emerged_internode_id] = all_element_outputs[tiller_to_MS_emerged_internode_id]
                        # else:
                        #     warnings.warn('No emerged internode found on main stem for tiller {}.'.format(tiller_to_MS_emerged_internode_id))

                        # Enclosed internode
                        tiller_to_MS_enclosed_internode_id = tiller_to_MS_phytomer_id + tuple(['internode', 'HiddenElement'])
                        tiller_enclosed_internode_id = hiddenzone_id + tuple(['internode', 'HiddenElement'])
                        if tiller_to_MS_enclosed_internode_id in all_element_outputs.keys():
                            self.outputs['elements'][tiller_enclosed_internode_id] = all_element_outputs[tiller_to_MS_enclosed_internode_id]
                        # else:
                        #     warnings.warn('No enclosed internode found on main stem for tiller {}.'.format(tiller_to_MS_enclosed_internode_id))
                # else:
                #     warnings.warn('No main stem found for tiller {}.'.format(tiller_to_MS_phytomer_id))

            #: Main Stem
            else:
                curr_hiddenzone_outputs = all_hiddenzone_outputs[hiddenzone_id]
                curr_axis_outputs = all_axes_outputs[axis_id]

                curr_hiddenzone_outputs['hiddenzone_age'] += curr_axis_outputs['delta_teq']
                if self.hydraulics:
                    curr_hiddenzone_outputs['init_leaf_L'] = hiddenzone_inputs['leaf_L']    # leaf length before Hydraulics sub-model

                hidden_sheath_id = hiddenzone_id + tuple(['sheath', 'HiddenElement'])
                visible_sheath_id = hiddenzone_id + tuple(['sheath', 'StemElement'])
                hidden_lamina_id = hiddenzone_id + tuple(['blade', 'HiddenElement'])
                hidden_internode_id = hiddenzone_id + tuple(['internode', 'HiddenElement'])
                visible_internode_id = hiddenzone_id + tuple(['internode', 'StemElement'])
                next_hiddenzone_id = tuple(list(hiddenzone_id[:2]) + [hiddenzone_id[2] + 1])

                # Found previous hidden zone
                prev_hiddenzone_id = tuple(list(axis_id) + [phytomer_id - 1])
                if prev_hiddenzone_id in all_hiddenzone_inputs:
                    prev_leaf_emerged = all_hiddenzone_inputs[prev_hiddenzone_id]['leaf_is_emerged']
                    prev_hiddenzone_inputs = all_hiddenzone_inputs[prev_hiddenzone_id]
                else:
                    prev_hiddenzone_inputs = None
                    prev_leaf_emerged = True

                prev_leaf2_hiddenzone_id = tuple(list(axis_id) + [phytomer_id - 2])
                if prev_leaf2_hiddenzone_id in all_hiddenzone_inputs:
                    prev_leaf2_emerged = all_hiddenzone_inputs[prev_leaf2_hiddenzone_id]['leaf_is_emerged']
                else:
                    prev_leaf2_emerged = True

                # Cumulated length of internodes up to the hidden zone
                below_internode_lengths = all_sheath_internode_lengths[axis_id][phytomer_id]['cumulated_internode']
                bottom_hiddenzone_height = self.model.calculate_cumulated_internode_length(below_internode_lengths)

                # Distance between the bottom of the hiddenzone and the highest previous ligule
                leaf_pseudostem_length = self.model.calculate_leaf_pseudostem_length(all_ligule_height_df[all_ligule_height_df['axis_id'] == axis_id], bottom_hiddenzone_height, phytomer_id)
                curr_hiddenzone_outputs['leaf_pseudostem_length'] = leaf_pseudostem_length
                curr_hiddenzone_outputs['delta_leaf_pseudostem_length'] = leaf_pseudostem_length - hiddenzone_inputs['leaf_pseudostem_length']  # Variable used in growth

                # Calculate the internode pseudostem length
                curr_internode_L = hiddenzone_inputs['internode_L']
                internode_distance_to_emerge = self.model.calculate_internode_distance_to_emerge(all_ligule_height_df[all_ligule_height_df['axis_id'] == axis_id], bottom_hiddenzone_height, phytomer_id,
                                                                                            curr_internode_L)
                curr_hiddenzone_outputs['internode_distance_to_emerge'] = internode_distance_to_emerge
                curr_hiddenzone_outputs['delta_internode_distance_to_emerge'] = internode_distance_to_emerge - hiddenzone_inputs['internode_distance_to_emerge']  # Variable used in growth

                # In case leaf is already mature but internode is growing, we update sheath visible and hidden lengths.
                if not curr_hiddenzone_outputs['leaf_is_growing'] and curr_hiddenzone_outputs['leaf_is_emerged']:
                    if hidden_sheath_id in self.inputs['elements'].keys():
                        sheath_hidden_length = self.inputs['elements'][hidden_sheath_id]['length']
                    else:
                        sheath_hidden_length = 0.
                        new_sheath = self.element_init.__dict__.copy()
                        self.outputs['elements'][hidden_sheath_id] = new_sheath
                    if visible_sheath_id not in self.outputs['elements'].keys():
                        new_sheath = self.element_init.__dict__.copy()
                        self.outputs['elements'][visible_sheath_id] = new_sheath
                    total_sheath_L = sheath_hidden_length + self.outputs['elements'][visible_sheath_id]['length']
                    updated_sheath_hidden_length = min(total_sheath_L, leaf_pseudostem_length)
                    updated_sheath_visible_length = max(0, total_sheath_L - updated_sheath_hidden_length)
                    self.outputs['elements'][hidden_sheath_id]['length'] = updated_sheath_hidden_length
                    self.outputs['elements'][visible_sheath_id]['length'] = updated_sheath_visible_length

                ###################
                # Leaf elongation #
                ###################
                if curr_hiddenzone_outputs['leaf_is_growing']:

                    if leaf_pseudostem_length < 0:
                        warnings.warn('Pseudostem length of {} decreased while leaf growing.'.format(hiddenzone_id))

                    if prev_leaf2_emerged and not curr_hiddenzone_outputs['leaf_is_emerged'] and phytomer_id != 1:
                        prev_leaf2_lamina_id = tuple(list(axis_id) + [phytomer_id - 2] + ['blade'] + ['LeafElement1'])
                        # Two alternatives to find the time of leaf n-2 emergence: either from the age (time equivalent) of lamina n-2 or from the pseudo_age of hz n-1
                        if prev_leaf2_lamina_id in all_element_inputs:
                            time_prev_leaf2_emergence = all_element_inputs[prev_leaf2_lamina_id]['age_teq']
                        elif prev_hiddenzone_inputs is not None:
                            time_prev_leaf2_emergence = prev_hiddenzone_inputs['leaf_pseudo_age']
                        else:
                            raise Exception("Can't find previous hz for current hz {}".format(hiddenzone_id))
                        curr_hiddenzone_outputs['mean_conc_sucrose'] = self.model.calculate_mean_conc_sucrose(hiddenzone_inputs['mean_conc_sucrose'],
                                                                                                         time_prev_leaf2_emergence,
                                                                                                         curr_axis_outputs['delta_teq'],
                                                                                                         hiddenzone_inputs['sucrose'],
                                                                                                         hiddenzone_inputs['mstruct'])

                    if not prev_leaf_emerged and phytomer_id != 1:  #: Before the emergence of the previous leaf. Exponential-like elongation. Leaf1 starts directly in phase 2
                        # delta leaf length
                        delta_leaf_L = self.model.calculate_deltaL_preE(hiddenzone_inputs['sucrose'],
                                                                   hiddenzone_inputs['leaf_L'],
                                                                   hiddenzone_inputs['amino_acids'],
                                                                   hiddenzone_inputs['mstruct'],
                                                                   curr_axis_outputs['delta_teq'], phytomer_id,
                                                                   self.optimal_growth_option,
                                                                   xylem_water_potential)

                        leaf_L = hiddenzone_inputs['leaf_L'] + delta_leaf_L
                        curr_hiddenzone_outputs['ratio_DZ'] = 1

                    else:  #: After the emergence of the previous leaf.

                        # Leaf length
                        leaf_pseudo_age = self.model.calculate_leaf_pseudo_age(hiddenzone_inputs['leaf_pseudo_age'], curr_axis_outputs['delta_teq'], hz_turgor_water_potential)
                        curr_hiddenzone_outputs['leaf_pseudo_age'] = leaf_pseudo_age
                        curr_hiddenzone_outputs['delta_leaf_pseudo_age'] = leaf_pseudo_age - hiddenzone_inputs['leaf_pseudo_age']

                        if not self.hydraulics:
                            delta_leaf_L = self.model.calculate_deltaL_postE(hiddenzone_inputs['leaf_pseudo_age'], leaf_pseudo_age, hiddenzone_inputs['leaf_L'], hiddenzone_inputs['leaf_Lmax_em'],
                                                                    hiddenzone_inputs['sucrose'], hiddenzone_inputs['amino_acids'], hiddenzone_inputs['mstruct'], phytomer_id, self.optimal_growth_option)
                            leaf_L = hiddenzone_inputs['leaf_L'] + delta_leaf_L

                            # Update leaf_Lmax. Subsequently, lamina_Lmax and sheath_Lmax will be updated depending on each element status (growing or mature)
                            curr_hiddenzone_outputs['leaf_Lmax'] = self.model.calculate_update_leaf_Lmax(hiddenzone_inputs['leaf_Lmax_em'], leaf_L, leaf_pseudo_age, phytomer_id)

                            # Ratio (mass) of Division Zone in the hiddenzone
                            curr_hiddenzone_outputs['ratio_DZ'] = self.model.calculate_ratio_DZ_postE(leaf_L, curr_hiddenzone_outputs['leaf_Lmax'], leaf_pseudostem_length)
                        else:
                            delta_leaf_L = 0
                            leaf_L = hiddenzone_inputs['leaf_L'] + delta_leaf_L
                            # Ratio (mass) of Division Zone in the hiddenzone
                            curr_hiddenzone_outputs['ratio_DZ'] = self.model.calculate_ratio_DZ_postE(hiddenzone_inputs['leaf_L'], leaf_pseudo_age, leaf_pseudostem_length)

                        #: Lamina : sheath distinction
                        lamina_id = hiddenzone_id + tuple(['blade', 'LeafElement1'])

                        ##: Lamina has not emerged
                        if not curr_hiddenzone_outputs['leaf_is_emerged']:

                            #: Test of leaf emergence against distance to leaf emergence. Assumes that a leaf cannot emerge before the previous one
                            #  TODO: besoin correction pour savoir a quel pas de temps exact??
                            curr_hiddenzone_outputs['leaf_is_emerged'] = self.model.calculate_leaf_emergence(hiddenzone_inputs['leaf_L'], leaf_pseudostem_length)
                            if curr_hiddenzone_outputs['leaf_is_emerged']:  # Initialise lamina outputs
                                print('leaf emergence', hiddenzone_id)
                                new_lamina = self.element_init.__dict__.copy()
                                self.outputs['elements'][lamina_id] = new_lamina

                                # Length of the HZ at the emergence
                                curr_hiddenzone_outputs['length_hz_En'] = curr_hiddenzone_outputs['leaf_L']

                                curr_lamina_outputs = all_element_outputs[lamina_id]
                                # Length of emerged lamina
                                lamina_L = self.model.calculate_lamina_L(leaf_L, leaf_pseudostem_length, curr_hiddenzone_outputs['lamina_Lmax'])
                                curr_lamina_outputs['length'] = lamina_L

                                # Update of lamina outputs
                                self.outputs['elements'][lamina_id] = curr_lamina_outputs

                                # Initialise variables for the next hidden zone as its previous leaf has now emerged
                                if next_hiddenzone_id in all_hiddenzone_inputs:
                                    next_hiddenzone_inputs = all_hiddenzone_inputs[next_hiddenzone_id]
                                    next_hiddenzone_outputs = all_hiddenzone_outputs[next_hiddenzone_id]

                                    if not self.hydraulics:
                                        next_hiddenzone_outputs['leaf_Lmax'] = self.model.calculate_leaf_Lmax(next_hiddenzone_inputs['leaf_L'], phytomer_id+1)  #: Final leaf length
                                        next_hiddenzone_outputs['leaf_Lmax_em'] = next_hiddenzone_outputs['leaf_Lmax']  #: Final leaf length at Ln-1 em
                                        sheath_lamina_ratio = self.model.calculate_SL_ratio(next_hiddenzone_id[2])  #: Sheath:Lamina final length ratio
                                        next_hiddenzone_outputs['lamina_Lmax'] = self.model.calculate_lamina_Lmax(next_hiddenzone_outputs['leaf_Lmax'], sheath_lamina_ratio)  #: Final lamina length
                                        next_hiddenzone_outputs['sheath_Lmax'] = self.model.calculate_sheath_Lmax(next_hiddenzone_outputs['leaf_Lmax'],
                                                                                                             next_hiddenzone_outputs['lamina_Lmax'])  #: Final sheath length
                                    else:
                                        next_hiddenzone_outputs['lamina_Lmax'] = self.model.calculate_lamina_Lmax(hiddenzone_id[2])  #: Maximal lamina length

                                    next_hiddenzone_outputs['leaf_pseudo_age'] = 0  #: Pseudo age of the leaf since beginning of automate growth (s)
                                    self.outputs['hiddenzone'][next_hiddenzone_id] = next_hiddenzone_outputs
                                else:
                                    warnings.warn('No next hidden zone found for hiddenzone {}.'.format(hiddenzone_id))

                                # Define lamina_Wmax and structural weight of the current sheath and lamina
                                if not self.hydraulics:
                                    curr_hiddenzone_outputs['leaf_Wmax'] = self.outputs['elements'][lamina_id]['Wmax'] = self.model.calculate_leaf_Wmax(curr_hiddenzone_outputs['lamina_Lmax'], hiddenzone_id[2],
                                                                                                                                                   curr_hiddenzone_outputs['mean_conc_sucrose'],
                                                                                                                                                   self.optimal_growth_option)
                                curr_hiddenzone_outputs['SSLW'] = self.model.calculate_SSLW(hiddenzone_id[2], curr_hiddenzone_outputs['mean_conc_sucrose'], self.optimal_growth_option)
                                curr_hiddenzone_outputs['LSSW'] = self.model.calculate_LSSW(hiddenzone_id[2], curr_hiddenzone_outputs['mean_conc_sucrose'], self.optimal_growth_option)

                            #: Test end of elongation when a leaf stops elongation inside the pseudostem (extreme stress)
                            if not self.hydraulics and leaf_L >= curr_hiddenzone_outputs['leaf_Lmax']:
                                # Update lamina_Lmax and sheath_Lmax based on updates of leaf_Lmax
                                sheath_lamina_ratio = self.model.calculate_SL_ratio(hiddenzone_id[2])
                                # Initialise hidden lamina
                                new_hidden_lamina = self.element_init.__dict__.copy()
                                self.outputs['elements'][hidden_lamina_id] = new_hidden_lamina
                                self.outputs['elements'][hidden_lamina_id]['length'] = curr_hiddenzone_outputs['lamina_Lmax'] = self.model.calculate_lamina_Lmax(curr_hiddenzone_outputs['leaf_Lmax'],
                                                                                                                                                            sheath_lamina_ratio)
                                self.outputs['elements'][hidden_lamina_id]['Wmax'] = curr_hiddenzone_outputs['leaf_Wmax']
                                self.outputs['elements'][hidden_lamina_id]['is_growing'] = False
                                # Initialise hidden sheath outputs
                                new_sheath = self.element_init.__dict__.copy()
                                self.outputs['elements'][hidden_sheath_id] = new_sheath
                                self.outputs['elements'][hidden_sheath_id]['length'] = curr_hiddenzone_outputs['sheath_Lmax'] = self.model.calculate_sheath_Lmax(curr_hiddenzone_outputs['leaf_Lmax'],
                                                                                                                                                            curr_hiddenzone_outputs['lamina_Lmax'])
                                self.outputs['elements'][hidden_sheath_id]['is_growing'] = False
                                # End of leaf growth
                                curr_hiddenzone_outputs['leaf_is_growing'] = False
                                curr_hiddenzone_outputs['leaf_is_remobilizing'] = True

                        else:
                            transition_lamina_to_sheath = False
                            end_leaf_elongation = False
                            if self.hydraulics:
                                # hydraulics -> the whole leaf is considered
                                if curr_hiddenzone_outputs['leaf_pseudo_age'] >= self.model.parameters.te:
                                    end_leaf_elongation = True

                            else:
                                #: Lamina has emerged and is growing
                                if all_element_inputs[lamina_id]['is_growing']:
                                    curr_lamina_outputs = all_element_outputs[lamina_id]
                                    # Update lamina_Lmax and sheath_Lmax based on updates of leaf_Lmax
                                    sheath_lamina_ratio = self.model.calculate_SL_ratio(hiddenzone_id[2])
                                    curr_hiddenzone_outputs['lamina_Lmax'] = self.model.calculate_lamina_Lmax(curr_hiddenzone_outputs['leaf_Lmax'], sheath_lamina_ratio)
                                    curr_hiddenzone_outputs['sheath_Lmax'] = self.model.calculate_sheath_Lmax(curr_hiddenzone_outputs['leaf_Lmax'], curr_hiddenzone_outputs['lamina_Lmax'])

                                    # Length of emerged lamina
                                    lamina_L = self.model.calculate_lamina_L(leaf_L, leaf_pseudostem_length, curr_hiddenzone_outputs['lamina_Lmax'])
                                    curr_lamina_outputs['length'] = lamina_L
                                    # Update of lamina outputs
                                    self.outputs['elements'][lamina_id] = curr_lamina_outputs

                                    # Test end of elongation
                                    if (lamina_L >= curr_hiddenzone_outputs['lamina_Lmax']) or (leaf_L >= curr_hiddenzone_outputs['leaf_Lmax']):
                                        transition_lamina_to_sheath = True

                                # Mature lamina, growing sheath
                                else:
                                    curr_visible_sheath_outputs = all_element_outputs[visible_sheath_id]
                                    curr_hidden_sheath_outputs = all_element_outputs[hidden_sheath_id]

                                    # Update only sheath_Lmax based on updates of leaf_Lmax
                                    curr_hiddenzone_outputs['sheath_Lmax'] = curr_hiddenzone_outputs['leaf_Lmax'] - \
                                                                             curr_hiddenzone_outputs['lamina_Lmax']

                                    # Length of emerged sheath
                                    lamina_L = \
                                        self.outputs['elements'][hiddenzone_id + tuple(['blade', 'LeafElement1'])][
                                            'length']
                                    emerged_sheath_L = self.model.calculate_emerged_sheath_L(leaf_L,
                                                                                             leaf_pseudostem_length,
                                                                                             lamina_L,
                                                                                             curr_hiddenzone_outputs[
                                                                                                 'sheath_Lmax'])
                                    curr_visible_sheath_outputs['length'] = emerged_sheath_L

                                    # Length of hidden sheath
                                    curr_hidden_sheath_outputs['length'] = leaf_pseudostem_length

                                    #: Test end of elongation
                                    if leaf_L >= curr_hiddenzone_outputs['leaf_Lmax']:
                                        end_leaf_elongation = True

                            # Deals with the initiation of the sheath element and the end of lamina elongation (only in the non-hydraulics version)
                            if transition_lamina_to_sheath:
                                curr_lamina_outputs['is_growing'] = False
                                curr_lamina_outputs['length'] = min(curr_hiddenzone_outputs['lamina_Lmax'], lamina_L)

                                # Initialise hidden sheath outputs
                                self.outputs['elements'][hidden_sheath_id] = self.element_init.__dict__.copy()
                                self.outputs['elements'][hidden_sheath_id]['length'] = min(curr_hiddenzone_outputs['sheath_Lmax'], leaf_pseudostem_length)
                                # Initialise visible sheath outputs
                                self.outputs['elements'][visible_sheath_id] = self.element_init.__dict__.copy()
                                self.outputs['elements'][visible_sheath_id]['length'] = self.model.calculate_emerged_sheath_L(leaf_L, leaf_pseudostem_length,
                                                                                      lamina_L, curr_hiddenzone_outputs['sheath_Lmax'])

                                # Case a mature sheath is shorter than the previous one.
                                hidden_lamina_L = self.model.calculate_hidden_lamina_L(lamina_L, curr_hiddenzone_outputs['lamina_Lmax'])
                                if hidden_lamina_L > 0:
                                    # Initialise hidden lamina if any
                                    self.outputs['elements'][hidden_lamina_id] = self.element_init.__dict__.copy()
                                    self.outputs['elements'][hidden_lamina_id]['length'] = hidden_lamina_L
                                    self.outputs['elements'][hidden_lamina_id]['Wmax'] = curr_hiddenzone_outputs['leaf_Wmax']
                                    self.outputs['elements'][hidden_lamina_id]['is_growing'] = False
                                    # End of sheath elongation
                                    self.outputs['elements'][visible_sheath_id]['is_growing'] = False
                                    self.outputs['elements'][hidden_sheath_id]['is_growing'] = False
                                    # End of leaf growth
                                    end_leaf_elongation = True

                                # Initialise variables for the next internode
                                if next_hiddenzone_id in all_hiddenzone_inputs:
                                    next_hiddenzone_outputs = all_hiddenzone_outputs[next_hiddenzone_id]
                                    if curr_axis_outputs['GA']:
                                        next_hiddenzone_outputs['internode_Lmax'] = self.model.calculate_internode_Lmax(next_hiddenzone_outputs['internode_L'])  #: Estimate of final internode length
                                        next_hiddenzone_outputs['internode_Lmax_lig'] = next_hiddenzone_outputs['internode_Lmax']  #: Estimate of final internode length at previous leaf ligulation
                                    next_hiddenzone_outputs['LSIW'] = self.model.calculate_LSIW(next_hiddenzone_outputs['LSSW'], next_hiddenzone_id[2],
                                                                                           optimal_growth_option=True)  #: Lineic Structural Internode Weight
                                    next_hiddenzone_outputs['internode_pseudo_age'] = 0  #: Pseudo age of the internode since beginning of automate growth (s)
                                    self.outputs['hiddenzone'][next_hiddenzone_id] = next_hiddenzone_outputs
                                else:
                                    warnings.warn('No next hidden zone found for hiddenzone {}.'.format(hiddenzone_id))

                            # End of leaf growth
                            if end_leaf_elongation:
                                curr_hiddenzone_outputs['leaf_is_growing'] = False
                                curr_hiddenzone_outputs['leaf_is_remobilizing'] = True

                                if self.hydraulics:
                                    all_element_outputs[lamina_id]['is_growing'] = False
                                    curr_hiddenzone_outputs['sheath_is_emerged'] = True

                                    # Initialise hidden and visible sheath outputs
                                    self.outputs['elements'][hidden_sheath_id] = self.element_init.__dict__.copy()
                                    self.outputs['elements'][visible_sheath_id] = self.element_init.__dict__.copy()
                                    #- Hidden sheath
                                    self.outputs['elements'][hidden_sheath_id]['length'] = leaf_pseudostem_length
                                    self.outputs['elements'][hidden_sheath_id]['width'] = all_hiddenzone_inputs[hiddenzone_id]['width']
                                    self.outputs['elements'][hidden_sheath_id]['thickness'] = all_hiddenzone_inputs[hiddenzone_id]['thickness']
                                    self.outputs['elements'][hidden_sheath_id]['is_growing'] = False

                                    # - Visible sheath
                                    SL_ratio = self.model.calculate_SL_ratio(hiddenzone_id[2])
                                    all_element_outputs[lamina_id]['length'] = curr_hiddenzone_outputs['leaf_L'] / (1 + SL_ratio)
                                    self.outputs['elements'][visible_sheath_id]['length'] = max(1E-03,
                                                                                                curr_hiddenzone_outputs[
                                                                                                    'leaf_L'] -
                                                                                                all_element_outputs[
                                                                                                    lamina_id][
                                                                                                    'length'] - leaf_pseudostem_length)

                                    self.outputs['elements'][visible_sheath_id]['is_growing'] = False
                                    self.outputs['elements'][visible_sheath_id]['width'] = all_hiddenzone_inputs[hiddenzone_id]['width']
                                    self.outputs['elements'][visible_sheath_id]['thickness'] = all_hiddenzone_inputs[hiddenzone_id]['thickness']

                                    # Hiddenzone
                                    curr_hiddenzone_outputs['lamina_Lmax'] = all_element_outputs[lamina_id]['length']
                                    curr_hiddenzone_outputs['sheath_Lmax'] = self.outputs['elements'][visible_sheath_id]['length'] + leaf_pseudostem_length
                                    curr_hiddenzone_outputs['SL_ratio'] = SL_ratio

                                    # Initialise variables for the next internode todo: temporary version for hydraulics : new internodes should initiate at ligulation
                                    if next_hiddenzone_id in all_hiddenzone_inputs:
                                        next_hiddenzone_outputs = all_hiddenzone_outputs[next_hiddenzone_id]
                                        if curr_axis_outputs['GA']:
                                            next_hiddenzone_outputs['internode_Lmax'] = self.model.calculate_internode_Lmax(next_hiddenzone_outputs['internode_L'])  #: Estimate of final internode length
                                            next_hiddenzone_outputs['internode_Lmax_lig'] = next_hiddenzone_outputs['internode_Lmax']  #: Estimate of final internode length at previous leaf ligulation
                                        next_hiddenzone_outputs['LSIW'] = self.model.calculate_LSIW(next_hiddenzone_outputs['LSSW'], next_hiddenzone_id[2],
                                                                                               optimal_growth_option=True)  #: Lineic Structural Internode Weight
                                        next_hiddenzone_outputs['internode_pseudo_age'] = 0  #: Pseudo age of the internode since beginning of automate growth (s)
                                        self.outputs['hiddenzone'][next_hiddenzone_id] = next_hiddenzone_outputs
                                    else:
                                        warnings.warn('No next hidden zone found for hiddenzone {}.'.format(hiddenzone_id))

                                else: # Without hydraulics
                                    self.outputs['elements'][visible_sheath_id]['is_growing'] = False
                                    self.outputs['elements'][hidden_sheath_id]['is_growing'] = False



                #: Leaf is mature but internode may be growing
                else:
                    leaf_L = hiddenzone_inputs['leaf_L']
                    delta_leaf_L = 0

                # Update of leaf outputs
                curr_hiddenzone_outputs['leaf_L'] = leaf_L
                curr_hiddenzone_outputs['delta_leaf_L'] = delta_leaf_L

                ########################
                # Internode elongation #
                ########################

                #: Initialisation of internode elongation
                if (not curr_hiddenzone_outputs['internode_is_growing']) and (curr_hiddenzone_outputs['internode_L'] == 0):
                    #: As for leaf primordia, we neglect CN growth due to IN length initialisation
                    curr_hiddenzone_outputs['internode_is_growing'], curr_hiddenzone_outputs['internode_L'] = self.model.calculate_init_internode_elongation(curr_hiddenzone_outputs['hiddenzone_age'],
                                                                                                                                                        phytomer_id)
                    if curr_hiddenzone_outputs['internode_is_growing']:
                        new_internode = self.element_init.__dict__.copy()
                        self.outputs['elements'][hidden_internode_id] = new_internode
                        self.outputs['elements'][hidden_internode_id]['length'] = curr_hiddenzone_outputs['internode_L']

                if curr_hiddenzone_outputs['internode_is_growing']:

                    #: Found previous lamina to know if the previous leaf is ligulated.
                    prev_lamina_id = tuple(list(hiddenzone_id[:2]) + [hiddenzone_id[2] - 1]) + tuple(['blade', 'LeafElement1'])
                    if phytomer_id == 1 and tuple(list(hiddenzone_id[:2]) + [hiddenzone_id[2] - 1]) + tuple(['sheath', 'StemElement']) in all_element_inputs:
                        prev_leaf_ligulated = not all_element_inputs[tuple(list(hiddenzone_id[:2]) + [hiddenzone_id[2] - 1]) + tuple(['sheath', 'StemElement'])]['is_growing']
                    elif prev_lamina_id in all_element_inputs:
                        prev_leaf_ligulated = not all_element_inputs[prev_lamina_id]['is_growing']
                    else:
                        prev_leaf_ligulated = False

                    #: Before ligulation of the leaf on the previous phytomer. Exponential-like elong. cf. Kirby 1988, Malvoisin 1984 II
                    if not prev_leaf_ligulated:
                        delta_internode_L = self.model.calculate_delta_internode_L_preL(phytomer_id, curr_hiddenzone_outputs['sucrose'], curr_hiddenzone_outputs['internode_L'],
                                                                                   curr_hiddenzone_outputs['amino_acids'], curr_hiddenzone_outputs['mstruct'],
                                                                                   curr_axis_outputs['delta_teq'])
                        internode_L = curr_hiddenzone_outputs['internode_L'] + delta_internode_L  # TODO: Check internode_L is not too large (in the case of long delta_t)

                        curr_hiddenzone_outputs['internode_L'] = internode_L
                        curr_hiddenzone_outputs['delta_internode_L'] = delta_internode_L
                        # Hidden internode
                        if hidden_internode_id not in self.outputs['elements'].keys():
                            new_internode = self.element_init.__dict__.copy()
                            self.outputs['elements'][hidden_internode_id] = new_internode
                        self.outputs['elements'][hidden_internode_id]['length'] = internode_L

                    #: After ligulation of the leaf on the previous phytomer.
                    else:
                        prev_internode_pseudo_age = curr_hiddenzone_outputs['internode_pseudo_age']
                        internode_pseudo_age = self.model.calculate_internode_pseudo_age(prev_internode_pseudo_age, curr_axis_outputs['delta_teq'])

                        curr_hiddenzone_outputs['internode_pseudo_age'] = internode_pseudo_age
                        curr_hiddenzone_outputs['delta_internode_pseudo_age'] = internode_pseudo_age - curr_hiddenzone_outputs['internode_pseudo_age']

                        #: Elongation only if Gibberelin production by SAM
                        if curr_axis_outputs['GA']:
                            #: Case of internodes that will not fully elongate, GA synthesis started after their previous leaf ligulation (i.e. no Lmax defined)
                            if pd.isnull(curr_hiddenzone_outputs['internode_Lmax']) or pd.isnull(curr_hiddenzone_outputs['internode_Lmax_lig']):
                                curr_hiddenzone_outputs['internode_Lmax'] = curr_hiddenzone_outputs['internode_Lmax_lig'] = self.model.calculate_short_internode_Lmax(curr_hiddenzone_outputs['internode_L'],
                                                                                                                                                                 curr_hiddenzone_outputs[
                                                                                                                                                                     'internode_pseudo_age'])
                            delta_internode_L = self.model.calculate_delta_internode_L_postL(prev_internode_pseudo_age, curr_hiddenzone_outputs['internode_pseudo_age'],
                                                                                        hiddenzone_inputs['internode_L'], curr_hiddenzone_outputs['internode_Lmax_lig'],
                                                                                        hiddenzone_inputs['sucrose'], hiddenzone_inputs['amino_acids'],
                                                                                        hiddenzone_inputs['mstruct'], optimal_growth_option=True)
                            internode_L = hiddenzone_inputs['internode_L'] + delta_internode_L

                            # Update internode_Lmax
                            if hiddenzone_inputs['internode_Lmax']:
                                curr_hiddenzone_outputs['internode_Lmax'] = self.model.calculate_update_internode_Lmax(curr_hiddenzone_outputs['internode_Lmax_lig'], internode_L, internode_pseudo_age)

                            #: Internode is not visible
                            if not curr_hiddenzone_outputs['internode_is_visible']:
                                #: Test of internode emergence.
                                curr_hiddenzone_outputs['internode_is_visible'] = self.model.calculate_internode_visibility(curr_hiddenzone_outputs['internode_L'], internode_distance_to_emerge)
                                if curr_hiddenzone_outputs['internode_is_visible']:  #: Initialise internode outputs
                                    new_internode_outputs = self.element_init.__dict__.copy()
                                    self.outputs['elements'][visible_internode_id] = new_internode_outputs
                                    self.outputs['elements'][visible_internode_id]['length'] = min(curr_hiddenzone_outputs['internode_Lmax'],
                                                                                                   self.model.calculate_emerged_internode_L(internode_L, internode_distance_to_emerge))
                                    self.outputs['elements'][hidden_internode_id]['length'] = internode_distance_to_emerge
                                else:
                                    self.outputs['elements'][hidden_internode_id]['length'] = internode_L

                            #: Internode is visible
                            else:
                                self.outputs['elements'][visible_internode_id]['length'] = min(curr_hiddenzone_outputs['internode_Lmax'],
                                                                                               self.model.calculate_emerged_internode_L(internode_L, internode_distance_to_emerge))
                                self.outputs['elements'][hidden_internode_id]['length'] = internode_distance_to_emerge

                            #: Test end of elongation
                            if self.model.calculate_end_internode_elongation(internode_L, curr_hiddenzone_outputs['internode_Lmax'], curr_hiddenzone_outputs['internode_pseudo_age']):
                                curr_hiddenzone_outputs['internode_is_growing'] = False
                                curr_hiddenzone_outputs['internode_is_remobilizing'] = True
                                # Visible internode
                                if curr_hiddenzone_outputs['internode_is_visible']:
                                    self.outputs['elements'][visible_internode_id]['is_growing'] = False
                                # Hidden internode
                                self.outputs['elements'][hidden_internode_id]['is_growing'] = False

                        #: Internode not elongating because no GA
                        else:
                            internode_L = curr_hiddenzone_outputs['internode_L']
                            delta_internode_L = 0
                            # Test end of elongation
                            if self.model.calculate_end_internode_elongation(internode_L, curr_hiddenzone_outputs['internode_Lmax'], curr_hiddenzone_outputs['internode_pseudo_age']):
                                curr_hiddenzone_outputs['internode_is_growing'] = False
                                curr_hiddenzone_outputs['internode_is_remobilizing'] = True
                                curr_hiddenzone_outputs['internode_Lmax'] = curr_hiddenzone_outputs['internode_L']
                                # Visible internode
                                if curr_hiddenzone_outputs['internode_is_visible']:
                                    self.outputs['elements'][visible_internode_id]['is_growing'] = False
                                # Hidden internode
                                if hidden_internode_id not in self.outputs['elements'].keys():
                                    new_internode = self.element_init.__dict__.copy()
                                    self.outputs['elements'][hidden_internode_id] = new_internode
                                self.outputs['elements'][hidden_internode_id]['is_growing'] = False
                                self.outputs['elements'][hidden_internode_id]['length'] = min(internode_L, internode_distance_to_emerge)
                                print('internode {} is over'.format(hidden_internode_id))

                #: Internode not elongating (not yet or already mature)
                else:
                    internode_L = curr_hiddenzone_outputs['internode_L']
                    delta_internode_L = 0

                # Update internodes outputs
                curr_hiddenzone_outputs['internode_L'] = internode_L
                curr_hiddenzone_outputs['delta_internode_L'] = delta_internode_L

                # Only growing hiddenzones are sent
                # after end of elongation (leaf and/or internode), it should :
                #       - pass by Growth for remobilisation
                #       - pass another time by morphogenesis for update of curr_element_outputs['final_hidden_length']
                # the hiddenzone will then be deleted since both growing flags are False and both delta_L are zeros.
                if hiddenzone_inputs['internode_is_growing'] or curr_hiddenzone_outputs['leaf_is_growing'] or \
                        curr_hiddenzone_outputs['leaf_is_remobilizing'] or curr_hiddenzone_outputs['internode_is_remobilizing']:
                    self.outputs['hiddenzone'][hiddenzone_id] = curr_hiddenzone_outputs
                else:  # End of internode elongation
                    del self.outputs['hiddenzone'][hiddenzone_id]
                    print('hgz {} is over'.format(hiddenzone_id))

        # -----------------------------
        # -------- COLEOPTILE ---------
        # -----------------------------
        # todo : do not use with hydraulics version
        # todo: add coleoptile internode (seems to elongate only if sowing depth is too high ; nodes/internodes of leaf 1-2 can also elongate in such cases)
        coleoptiles = {coleo_id: coleo_inputs for coleo_id, coleo_inputs in all_hiddenzone_inputs.items() if coleo_id[1:] == ('MS', 0)}
        for coleo_id, coleo_inputs in coleoptiles.items():
            curr_axis_outputs = all_axes_outputs[coleo_id[:2]]
            curr_coleoptile_outputs = all_hiddenzone_outputs[coleo_id]

            leaf1_id = coleo_id[:2] + (1,)
            if leaf1_id in all_hiddenzone_outputs.keys():
                leaf1_outputs = all_hiddenzone_outputs[leaf1_id]
                if not leaf1_outputs['leaf_is_growing']:
                    continue
            else:
                continue
            #: Maximal coleoptile length provided in inputs or calculated from sowing depth
            if np.isnan(curr_coleoptile_outputs['leaf_Lmax']):
                curr_coleoptile_outputs['leaf_Lmax'] = self.model.calculate_coleo_Lmax(Zsowing)

            #: Coleoptile elongation equals that of leaf 1
            curr_coleoptile_outputs['delta_leaf_L'] = min(leaf1_outputs['delta_leaf_L'], curr_coleoptile_outputs['leaf_Lmax'] - curr_coleoptile_outputs['leaf_L'])
            curr_coleoptile_outputs['leaf_L'] += curr_coleoptile_outputs['delta_leaf_L']
            curr_coleoptile_outputs['leaf_pseudostem_length'] = curr_coleoptile_outputs['leaf_L']

            curr_coleoptile_outputs['leaf_pseudo_age'] += curr_axis_outputs['delta_teq']
            curr_coleoptile_outputs['delta_leaf_pseudo_age'] = leaf1_outputs['delta_leaf_pseudo_age']

            #: Coleoptile emergence
            if not curr_coleoptile_outputs['leaf_is_emerged']:
                curr_coleoptile_outputs['leaf_is_emerged'] = self.model.calculate_coleoptile_emergence(curr_coleoptile_outputs['leaf_L'], Zsowing)
                if curr_coleoptile_outputs['leaf_is_emerged']:
                    print('coleoptile has emerged', curr_coleoptile_outputs['leaf_L'])
                    # Initialise visible coleoptile outputs
                    visible_coleo_id = coleo_id + tuple(['sheath', 'StemElement'])
                    visible_coleoptile = self.element_init.__dict__.copy()
                    self.outputs['elements'][visible_coleo_id] = visible_coleoptile
                    visible_coleoptile_outputs = all_element_outputs[visible_coleo_id]
                    emerged_coleoptile_L = self.model.calculate_emerged_coleo_L(curr_coleoptile_outputs['leaf_L'], Zsowing)
                    visible_coleoptile_outputs['length'] = emerged_coleoptile_L
                    self.outputs['elements'][visible_coleo_id] = visible_coleoptile_outputs  # Update of coleoptile outputs

                    # Initialise hidden coleoptile outputs
                    hidden_coleo_id = coleo_id + tuple(['sheath', 'HiddenElement'])
                    hidden_coleoptile = self.element_init.__dict__.copy()
                    self.outputs['elements'][hidden_coleo_id] = hidden_coleoptile
                    self.outputs['elements'][hidden_coleo_id]['length'] = Zsowing  # Length of hidden coleoptile

            else:  # Coleoptile has emerged
                visible_coleo_id = coleo_id + tuple(['sheath', 'StemElement'])
                visible_coleoptile_outputs = all_element_outputs[visible_coleo_id]
                emerged_coleoptile_L = self.model.calculate_emerged_coleo_L(curr_coleoptile_outputs['leaf_L'], Zsowing)
                visible_coleoptile_outputs['length'] = emerged_coleoptile_L
                self.outputs['elements'][visible_coleo_id] = visible_coleoptile_outputs  # Update of coleoptile outputs

            # Update of coleoptile outputs
            self.outputs['hiddenzone'][coleo_id] = curr_coleoptile_outputs

            #: Test end of elongation (independently of coleoptile emergence)
            if curr_coleoptile_outputs['leaf_L'] >= curr_coleoptile_outputs['leaf_Lmax'] and curr_coleoptile_outputs['leaf_is_growing']:
                print('end of coleoptile', curr_coleoptile_outputs['leaf_Lmax'])
                curr_coleoptile_outputs['leaf_is_growing'] = False
                curr_coleoptile_outputs['leaf_is_remobilizing'] = True
                # Update of coleoptile outputs
                self.outputs['hiddenzone'][coleo_id] = curr_coleoptile_outputs

                # Update of visible coleoptile outputs
                visible_coleo_id = coleo_id + tuple(['sheath', 'StemElement'])
                if visible_coleo_id in all_element_outputs.keys():
                    visible_coleoptile_outputs = all_element_outputs[visible_coleo_id]
                    visible_coleoptile_outputs['is_growing'] = False
                    self.outputs['elements'][visible_coleo_id] = visible_coleoptile_outputs
                # Update of hidden coleoptile outputs
                hidden_coleo_id = coleo_id + tuple(['sheath', 'HiddenElement'])
                hidden_coleoptile_outputs = all_element_outputs[hidden_coleo_id]
                hidden_coleoptile_outputs['is_growing'] = False
                self.outputs['elements'][hidden_coleo_id] = hidden_coleoptile_outputs

            if not curr_coleoptile_outputs['leaf_is_growing'] and not curr_coleoptile_outputs['leaf_is_remobilizing']:  # to be adapted if accounting for coleoptile's internode elongation
                del self.outputs['hiddenzone'][coleo_id]
