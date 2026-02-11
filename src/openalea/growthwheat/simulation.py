# -*- coding: latin-1 -*-

from __future__ import division  # use "//" to do integer division

import copy

from growthwheat import model
from growthwheat import parameters

from respiwheat.model import RespirationModel

"""
    growthwheat.simulation
    ~~~~~~~~~~~~~~~~~~

    The module :mod:`growthwheat.simulation`.

    :copyright: Copyright 2014-2015 INRA-ECOSYS, see AUTHORS.
    :license: LICENSE for details.

"""

#: the inputs needed by GrowthWheat
HIDDENZONE_INPUTS = ['leaf_is_growing', 'internode_is_growing', 'leaf_pseudo_age', 'delta_leaf_pseudo_age', 'internode_pseudo_age', 'delta_internode_pseudo_age', 'leaf_L', 'delta_leaf_L',
                     'internode_L', 'delta_internode_L', 'leaf_pseudostem_length',
                     'delta_leaf_pseudostem_length', 'internode_distance_to_emerge', 'delta_internode_distance_to_emerge', 'SSLW', 'LSSW', 'LSIW', 'leaf_is_emerged', 'internode_is_visible',
                     'sucrose', 'amino_acids', 'fructan', 'proteins', 'leaf_enclosed_mstruct', 'leaf_enclosed_Nstruct', 'internode_enclosed_mstruct',
                     'internode_enclosed_Nstruct', 'mstruct', 'internode_Lmax', 'leaf_Lmax', 'sheath_Lmax', 'is_over', 'leaf_is_remobilizing', 'internode_is_remobilizing']
ELEMENT_INPUTS = ['is_growing', 'mstruct', 'senesced_mstruct', 'green_area', 'length', 'sucrose', 'amino_acids', 'fructan', 'proteins', 'cytokinins', 'Nstruct']
ROOT_INPUTS = ['sucrose', 'amino_acids', 'mstruct', 'Nstruct']
AXIS_INPUTS = ['delta_teq', 'delta_teq_roots']

#: the outputs computed by GrowthWheat
HIDDENZONE_OUTPUTS = ['sucrose', 'amino_acids', 'fructan', 'proteins', 'leaf_enclosed_mstruct', 'leaf_enclosed_Nstruct', 'internode_enclosed_mstruct', 'internode_enclosed_Nstruct', 'mstruct',
                      'Nstruct', 'Respi_growth', 'sucrose_consumption_mstruct', 'AA_consumption_mstruct', 'is_over', 'leaf_is_remobilizing', 'internode_is_remobilizing']
ELEMENT_OUTPUTS = ['sucrose', 'amino_acids', 'fructan', 'proteins', 'nitrates', 'mstruct', 'Nstruct',
                   'green_area', 'max_proteins', 'max_mstruct', 'Nresidual', 'senesced_length_element', 'senesced_mstruct']
ROOT_OUTPUTS = ['sucrose', 'amino_acids', 'mstruct', 'Nstruct', 'Respi_growth', 'R_min_upt', 'delta_mstruct_growth', 'sucrose_consumption_mstruct', 'AA_consumption_mstruct']
AXIS_OUTPUTS = []

#: the inputs and outputs of GrowthWheat.
HIDDENZONE_INPUTS_OUTPUTS = sorted(set(HIDDENZONE_INPUTS + HIDDENZONE_OUTPUTS))
ELEMENT_INPUTS_OUTPUTS = sorted(set(ELEMENT_INPUTS + ELEMENT_OUTPUTS))
ROOT_INPUTS_OUTPUTS = sorted(set(ROOT_INPUTS + ROOT_OUTPUTS))
AXIS_INPUTS_OUTPUTS = sorted(set(AXIS_INPUTS + AXIS_OUTPUTS))


class SimulationError(Exception):
    pass


class SimulationRunError(SimulationError):
    pass


class Simulation(object):
    """The Simulation class permits to initialize and run a simulation.
    """

    def __init__(self, delta_t=1, update_parameters=None, cnwheat_roots=True):

        #: The inputs of growth-Wheat.
        #:
        #: `inputs` is a dictionary of dictionaries:
        #:     {'hiddenzone': {(plant_index, axis_label, metamer_index): {hiddenzone_input_name: hiddenzone_input_value, ...}, ...},
        #:      'elements': {(plant_index, axis_label, metamer_index, organ_label, element_label): {organ_input_name: organ_input_value, ...}, ...},
        #:      'roots': {(plant_index, axis_label): {root_input_name: root_input_value, ...}, ...}}
        #: See :TODO?
        #: for more information about the inputs.
        self.inputs = {}

        #: The outputs of growth-Wheat.
        #:
        #: `outputs` is a dictionary of dictionaries:
        #:     {'hiddenzone': {(plant_index, axis_label, metamer_index): {hiddenzone_input_name: hiddenzone_input_value, ...}, ...},
        #:      'elements': {(plant_index, axis_label, metamer_index, organ_label, element_label): {organ_input_name: organ_input_value, ...}, ...}
        #:      'roots': {(plant_index, axis_label): {root_input_name: root_input_value, ...}, ...}}
        #: See :TODO?
        #: for more information about the inputs.
        self.outputs = {}

        #: the delta t of the simulation (in seconds)
        self.delta_t = delta_t

        #: Update parameters if specified
        if update_parameters:
            parameters.__dict__.update(update_parameters)

        self.cnwheat_roots = cnwheat_roots

    def initialize(self, inputs):
        """
        Initialize :attr:`inputs` from `inputs`.

        :param dict inputs: must be a dictionary with the same structure as :attr:`inputs`.
        """
        self.inputs.clear()
        self.inputs.update(inputs)

    def run(self, postflowering_stages=False):
        """
        Run the simulation.

        :param bool postflowering_stages: if True the model will calculate root growth with the parameters calibrated for post flowering stages
        """
        # Copy the inputs into the output dict
        self.outputs.update({inputs_type: copy.deepcopy(all_inputs) for inputs_type, all_inputs in self.inputs.items() if inputs_type in {'hiddenzone', 'elements', 'roots', 'axes'}})

        # Hidden growing zones
        all_hiddenzone_inputs = self.inputs['hiddenzone']
        all_hiddenzone_outputs = self.outputs['hiddenzone']

        # elements
        all_elements_inputs = self.inputs['elements']
        all_elements_outputs = self.outputs['elements']

        # roots
        all_roots_inputs = self.inputs['roots']
        all_roots_outputs = self.outputs['roots']

        # axes
        all_axes_inputs = self.inputs['axes']
        all_axes_outputs = self.outputs['axes']

        # ----------------------------------------------
        # ----------- Hiddenzones and elements ---------
        # ----------------------------------------------
        for hiddenzone_id, hiddenzone_inputs in sorted(all_hiddenzone_inputs.items()):

            curr_hiddenzone_outputs = all_hiddenzone_outputs[hiddenzone_id]

            axe_label = hiddenzone_id[1]
            #: Tillers (we copy corresponding elements of MS)
            if axe_label != 'MS':  # TODO: temporary or should be an option at least
                pass

            #: Main stem
            else:
                # Initialisation of the exports towards the growing lamina or sheath
                delta_leaf_enclosed_mstruct = delta_leaf_enclosed_Nstruct = delta_lamina_mstruct = delta_sheath_mstruct = delta_lamina_Nstruct = delta_sheath_Nstruct = leaf_export_sucrose = \
                    delta_internode_mstruct = delta_internode_Nstruct = leaf_export_amino_acids = leaf_remob_fructan = leaf_export_proteins = internode_export_sucrose = \
                    internode_export_amino_acids = internode_remob_fructan = internode_export_proteins = 0.

                # -- Delta Growth internode

                if hiddenzone_inputs['internode_pseudo_age'] < parameters.internode_rapid_growth_t:  #: Internode is not yet in rapide growth stage TODO : tester sur une variable "is_ligulated"
                    # delta mstruct of the internode
                    ratio_mstruct_DM = model.calculate_ratio_mstruct_DM(hiddenzone_inputs['mstruct'], hiddenzone_inputs['sucrose'], hiddenzone_inputs['fructan'],
                                                                        hiddenzone_inputs['amino_acids'], hiddenzone_inputs['proteins'])
                    delta_internode_enclosed_mstruct = model.calculate_delta_internode_enclosed_mstruct(hiddenzone_inputs['internode_L'], hiddenzone_inputs['delta_internode_L'], ratio_mstruct_DM)
                    # delta Nstruct of the internode
                    delta_internode_enclosed_Nstruct = model.calculate_delta_Nstruct(delta_internode_enclosed_mstruct)
                else:
                    # delta mstruct of the enclosed internode
                    delta_internode_enclosed_mstruct = model.calculate_delta_internode_enclosed_mstruct_postL(hiddenzone_inputs['delta_internode_pseudo_age'],
                                                                                                              hiddenzone_inputs['internode_pseudo_age'],
                                                                                                              hiddenzone_inputs['internode_L'],
                                                                                                              hiddenzone_inputs['internode_distance_to_emerge'],
                                                                                                              hiddenzone_inputs['internode_Lmax'],
                                                                                                              hiddenzone_inputs['LSIW'],
                                                                                                              hiddenzone_inputs['internode_enclosed_mstruct'])
                    # delta Nstruct of the enclosed internode
                    delta_internode_enclosed_Nstruct = model.calculate_delta_Nstruct(delta_internode_enclosed_mstruct)

                if hiddenzone_inputs['internode_is_visible']:  #: Internode is visible
                    visible_internode_id = hiddenzone_id + tuple(['internode', 'StemElement'])
                    curr_visible_internode_inputs = all_elements_inputs[visible_internode_id]
                    curr_visible_internode_outputs = all_elements_outputs[visible_internode_id]
                    # Delta mstruct of the emerged internode
                    delta_internode_mstruct = model.calculate_delta_emerged_tissue_mstruct(hiddenzone_inputs['LSIW'], curr_visible_internode_inputs['mstruct'], curr_visible_internode_inputs['length'])
                    # Delta Nstruct of the emerged internode
                    delta_internode_Nstruct = model.calculate_delta_Nstruct(delta_internode_mstruct)
                    # Export of sucrose from hiddenzone towards emerged internode
                    internode_export_sucrose = model.calculate_export(delta_internode_mstruct, hiddenzone_inputs['sucrose'], hiddenzone_inputs['mstruct'])
                    # Export of amino acids from hiddenzone towards emerged internode
                    internode_export_amino_acids = model.calculate_export(delta_internode_mstruct, hiddenzone_inputs['amino_acids'], hiddenzone_inputs['mstruct'])
                    internode_remob_fructan = model.calculate_export(delta_internode_mstruct, hiddenzone_inputs['fructan'], hiddenzone_inputs['mstruct'])
                    internode_export_proteins = model.calculate_export(delta_internode_mstruct, hiddenzone_inputs['proteins'], hiddenzone_inputs['mstruct'])

                    # Update of internode outputs
                    curr_visible_internode_outputs['mstruct'] += delta_internode_mstruct
                    curr_visible_internode_outputs['max_mstruct'] = curr_visible_internode_outputs['mstruct']
                    curr_visible_internode_outputs['Nstruct'] += delta_internode_Nstruct
                    curr_visible_internode_outputs['sucrose'] += internode_export_sucrose + internode_remob_fructan
                    curr_visible_internode_outputs['amino_acids'] += internode_export_amino_acids
                    curr_visible_internode_outputs['proteins'] += internode_export_proteins
                    self.outputs['elements'][visible_internode_id] = curr_visible_internode_outputs

                # -- Delta Growth leaf
                if not hiddenzone_inputs['leaf_is_emerged']:  #: Leaf is not emerged
                    # delta mstruct of the hidden leaf
                    ratio_mstruct_DM = model.calculate_ratio_mstruct_DM(hiddenzone_inputs['mstruct'], hiddenzone_inputs['sucrose'], hiddenzone_inputs['fructan'],
                                                                        hiddenzone_inputs['amino_acids'], hiddenzone_inputs['proteins'])
                    delta_leaf_enclosed_mstruct = model.calculate_delta_leaf_enclosed_mstruct(hiddenzone_inputs['leaf_L'], hiddenzone_inputs['delta_leaf_L'], ratio_mstruct_DM)
                    # delta Nstruct of the hidden leaf
                    delta_leaf_enclosed_Nstruct = model.calculate_delta_Nstruct(delta_leaf_enclosed_mstruct)
                elif hiddenzone_inputs['leaf_is_growing']:  #: Leaf has emerged and growing
                    # delta mstruct of the enclosed leaf (which length is assumed to equal the length of the pseudostem)
                    delta_leaf_enclosed_mstruct = model.calculate_delta_leaf_enclosed_mstruct_postE(hiddenzone_inputs['delta_leaf_pseudo_age'],
                                                                                                    hiddenzone_inputs['leaf_pseudo_age'],
                                                                                                    hiddenzone_inputs['leaf_pseudostem_length'],
                                                                                                    hiddenzone_inputs['leaf_enclosed_mstruct'],
                                                                                                    hiddenzone_inputs['LSSW'])
                    # delta Nstruct of the enclosed en leaf
                    delta_leaf_enclosed_Nstruct = model.calculate_delta_Nstruct(delta_leaf_enclosed_mstruct)

                    # leaf has emerged and still growing
                    visible_lamina_id = hiddenzone_id + tuple(['blade', 'LeafElement1'])
                    #: Lamina is growing
                    if all_elements_inputs[visible_lamina_id]['is_growing']:
                        curr_visible_lamina_inputs = all_elements_inputs[visible_lamina_id]
                        curr_visible_lamina_outputs = all_elements_outputs[visible_lamina_id]
                        # Delta mstruct of the emerged lamina
                        delta_lamina_mstruct = model.calculate_delta_emerged_tissue_mstruct(hiddenzone_inputs['SSLW'], curr_visible_lamina_inputs['mstruct'], curr_visible_lamina_inputs['green_area'])
                        # Delta Nstruct of the emerged lamina
                        delta_lamina_Nstruct = model.calculate_delta_Nstruct(delta_lamina_mstruct)
                        # Export of metabolite from hiddenzone towards emerged lamina
                        leaf_export_sucrose = model.calculate_export(delta_lamina_mstruct, hiddenzone_inputs['sucrose'], hiddenzone_inputs['mstruct'])
                        leaf_export_amino_acids = model.calculate_export(delta_lamina_mstruct, hiddenzone_inputs['amino_acids'], hiddenzone_inputs['mstruct'])
                        leaf_remob_fructan = model.calculate_export(delta_lamina_mstruct, hiddenzone_inputs['fructan'], hiddenzone_inputs['mstruct'])
                        leaf_export_proteins = model.calculate_export(delta_lamina_mstruct, hiddenzone_inputs['proteins'], hiddenzone_inputs['mstruct'])
                        # Cytokinins in the newly visible mstruct
                        addition_cytokinins = model.calculate_init_cytokinins_emerged_tissue(delta_lamina_mstruct)

                        # Update of lamina outputs
                        curr_visible_lamina_outputs['mstruct'] += delta_lamina_mstruct
                        curr_visible_lamina_outputs['max_mstruct'] = curr_visible_lamina_outputs['mstruct']
                        curr_visible_lamina_outputs['Nstruct'] += delta_lamina_Nstruct
                        curr_visible_lamina_outputs['sucrose'] += leaf_export_sucrose + leaf_remob_fructan
                        curr_visible_lamina_outputs['amino_acids'] += leaf_export_amino_acids
                        curr_visible_lamina_outputs['proteins'] += leaf_export_proteins
                        curr_visible_lamina_outputs['cytokinins'] += addition_cytokinins

                        self.outputs['elements'][visible_lamina_id] = curr_visible_lamina_outputs

                    else:  #: Mature lamina, growing sheath
                        # The hidden part of the sheath is only updated once, at the end of leaf elongation, by remobilisation from the hiddenzone
                        visible_sheath_id = hiddenzone_id + tuple(['sheath', 'StemElement'])
                        curr_visible_sheath_inputs = all_elements_inputs[visible_sheath_id]
                        curr_visible_sheath_outputs = all_elements_outputs[visible_sheath_id]
                        # Delta mstruct of the emerged sheath
                        delta_sheath_mstruct = model.calculate_delta_emerged_tissue_mstruct(hiddenzone_inputs['LSSW'], curr_visible_sheath_inputs['mstruct'], curr_visible_sheath_inputs['length'])
                        # Delta Nstruct of the emerged sheath
                        delta_sheath_Nstruct = model.calculate_delta_Nstruct(delta_sheath_mstruct)
                        # Export of metabolite from hiddenzone towards emerged sheath
                        leaf_export_sucrose = model.calculate_export(delta_sheath_mstruct, hiddenzone_inputs['sucrose'], hiddenzone_inputs['mstruct'])
                        leaf_export_amino_acids = model.calculate_export(delta_sheath_mstruct, hiddenzone_inputs['amino_acids'], hiddenzone_inputs['mstruct'])
                        leaf_remob_fructan = model.calculate_export(delta_sheath_mstruct, hiddenzone_inputs['fructan'], hiddenzone_inputs['mstruct'])
                        leaf_export_proteins = model.calculate_export(delta_sheath_mstruct, hiddenzone_inputs['proteins'], hiddenzone_inputs['mstruct'])
                        addition_cytokinins = model.calculate_init_cytokinins_emerged_tissue(delta_sheath_mstruct)

                        # Update of sheath outputs
                        curr_visible_sheath_outputs['mstruct'] += delta_sheath_mstruct
                        curr_visible_sheath_outputs['max_mstruct'] = curr_visible_sheath_outputs['mstruct']
                        curr_visible_sheath_outputs['Nstruct'] += delta_sheath_Nstruct
                        curr_visible_sheath_outputs['sucrose'] += leaf_export_sucrose + leaf_remob_fructan
                        curr_visible_sheath_outputs['amino_acids'] += leaf_export_amino_acids
                        curr_visible_sheath_outputs['proteins'] += leaf_export_proteins
                        curr_visible_sheath_outputs['cytokinins'] += addition_cytokinins
                        self.outputs['elements'][visible_sheath_id] = curr_visible_sheath_outputs

                # -- CN consumption due to mstruct/Nstruct growth of the enclosed leaf and of the internode
                curr_hiddenzone_outputs['AA_consumption_mstruct'] = model.calculate_s_Nstruct_amino_acids((delta_leaf_enclosed_Nstruct + delta_internode_enclosed_Nstruct),
                                                                                                          delta_lamina_Nstruct,
                                                                                                          delta_sheath_Nstruct,
                                                                                                          delta_internode_Nstruct)  #: Consumption of amino acids due to mstruct growth (µmol N)
                curr_hiddenzone_outputs['sucrose_consumption_mstruct'] = model.calculate_s_mstruct_sucrose((delta_leaf_enclosed_mstruct + delta_internode_enclosed_mstruct),
                                                                                                           delta_lamina_mstruct,
                                                                                                           delta_sheath_mstruct,
                                                                                                           curr_hiddenzone_outputs[
                                                                                                               'AA_consumption_mstruct'])  #: Consumption of sucrose due to mstruct growth (µmol C)
                curr_hiddenzone_outputs['Respi_growth'] = RespirationModel.R_growth(curr_hiddenzone_outputs['sucrose_consumption_mstruct'])  #: Respiration growth (µµmol C)

                # -- Update of hiddenzone outputs
                curr_hiddenzone_outputs['leaf_enclosed_mstruct'] += delta_leaf_enclosed_mstruct
                curr_hiddenzone_outputs['leaf_enclosed_Nstruct'] += delta_leaf_enclosed_Nstruct
                curr_hiddenzone_outputs['internode_enclosed_mstruct'] += delta_internode_enclosed_mstruct
                curr_hiddenzone_outputs['internode_enclosed_Nstruct'] += delta_internode_enclosed_Nstruct
                curr_hiddenzone_outputs['mstruct'] = curr_hiddenzone_outputs['leaf_enclosed_mstruct'] + curr_hiddenzone_outputs['internode_enclosed_mstruct']
                curr_hiddenzone_outputs['Nstruct'] = curr_hiddenzone_outputs['leaf_enclosed_Nstruct'] + curr_hiddenzone_outputs['internode_enclosed_Nstruct']
                curr_hiddenzone_outputs['sucrose'] -= (
                            curr_hiddenzone_outputs['sucrose_consumption_mstruct'] + curr_hiddenzone_outputs['Respi_growth'] + leaf_export_sucrose + internode_export_sucrose)
                curr_hiddenzone_outputs['fructan'] -= (leaf_remob_fructan + internode_remob_fructan)
                curr_hiddenzone_outputs['amino_acids'] -= (curr_hiddenzone_outputs['AA_consumption_mstruct'] + leaf_export_amino_acids + internode_export_amino_acids)
                curr_hiddenzone_outputs['proteins'] -= (leaf_export_proteins + internode_export_proteins)
                self.outputs['hiddenzone'][hiddenzone_id] = curr_hiddenzone_outputs

                # -- Remobilisation at the end of leaf elongation
                if hiddenzone_inputs['leaf_is_remobilizing']:
                    share_leaf = curr_hiddenzone_outputs['leaf_enclosed_mstruct'] / curr_hiddenzone_outputs['mstruct']

                    # Case when the hiddenzone contains a hidden part of lamina
                    if hiddenzone_inputs['leaf_pseudostem_length'] > hiddenzone_inputs['sheath_Lmax']:
                        hidden_sheath_mstruct = model.calculate_sheath_mstruct(hiddenzone_inputs['sheath_Lmax'], hiddenzone_inputs['LSSW'])
                        share_hidden_sheath = hidden_sheath_mstruct / curr_hiddenzone_outputs['leaf_enclosed_mstruct']
                    else:
                        share_hidden_sheath = 1

                    # Add to hidden part of the sheath
                    hidden_sheath_id = hiddenzone_id + tuple(['sheath', 'HiddenElement'])
                    if hidden_sheath_id not in self.outputs['elements'].keys():
                        new_sheath_outputs = parameters.OrganInit().__dict__
                        self.outputs['elements'][hidden_sheath_id] = new_sheath_outputs
                    curr_hidden_sheath_outputs = self.outputs['elements'][hidden_sheath_id]
                    curr_hidden_sheath_outputs['mstruct'] = curr_hiddenzone_outputs['leaf_enclosed_mstruct'] * share_hidden_sheath
                    curr_hidden_sheath_outputs['max_mstruct'] = curr_hiddenzone_outputs['leaf_enclosed_mstruct'] * share_hidden_sheath
                    curr_hidden_sheath_outputs['Nstruct'] = curr_hiddenzone_outputs['leaf_enclosed_Nstruct'] * share_hidden_sheath
                    curr_hidden_sheath_outputs['sucrose'] = curr_hiddenzone_outputs['sucrose'] * share_leaf * share_hidden_sheath
                    curr_hidden_sheath_outputs['amino_acids'] = curr_hiddenzone_outputs['amino_acids'] * share_leaf * share_hidden_sheath
                    curr_hidden_sheath_outputs['fructan'] = curr_hiddenzone_outputs['fructan'] * share_leaf * share_hidden_sheath
                    curr_hidden_sheath_outputs['proteins'] = curr_hiddenzone_outputs['proteins'] * share_leaf * share_hidden_sheath
                    self.outputs['elements'][hidden_sheath_id] = curr_hidden_sheath_outputs

                    # Add to hidden part of the lamina, if any
                    if share_hidden_sheath < 1:
                        hidden_lamina_id = hiddenzone_id + tuple(['blade', 'HiddenElement'])
                        curr_hidden_lamina_outputs = self.outputs['elements'][hidden_lamina_id]
                        curr_hidden_lamina_outputs['mstruct'] = curr_hiddenzone_outputs['leaf_enclosed_mstruct'] * (1 - share_hidden_sheath)
                        curr_hidden_lamina_outputs['max_mstruct'] = curr_hiddenzone_outputs['leaf_enclosed_mstruct'] * (1 - share_hidden_sheath)
                        curr_hidden_lamina_outputs['Nstruct'] = curr_hiddenzone_outputs['leaf_enclosed_Nstruct'] * (1 - share_hidden_sheath)
                        curr_hidden_lamina_outputs['sucrose'] = curr_hiddenzone_outputs['sucrose'] * share_leaf * (1 - share_hidden_sheath)
                        curr_hidden_lamina_outputs['amino_acids'] = curr_hiddenzone_outputs['amino_acids'] * share_leaf * (1 - share_hidden_sheath)
                        curr_hidden_lamina_outputs['fructan'] = curr_hiddenzone_outputs['fructan'] * share_leaf * (1 - share_hidden_sheath)
                        curr_hidden_lamina_outputs['proteins'] = curr_hiddenzone_outputs['proteins'] * share_leaf * (1 - share_hidden_sheath)
                        self.outputs['elements'][hidden_lamina_id] = curr_hidden_lamina_outputs

                    # Remove in hiddenzone
                    curr_hiddenzone_outputs = self.outputs['hiddenzone'][hiddenzone_id]
                    curr_hiddenzone_outputs['leaf_enclosed_mstruct'] = 0
                    curr_hiddenzone_outputs['leaf_enclosed_Nstruct'] = 0
                    curr_hiddenzone_outputs['mstruct'] = curr_hiddenzone_outputs['internode_enclosed_mstruct']
                    curr_hiddenzone_outputs['Nstruct'] = curr_hiddenzone_outputs['internode_enclosed_Nstruct']
                    curr_hiddenzone_outputs['sucrose'] -= curr_hiddenzone_outputs['sucrose'] * share_leaf
                    curr_hiddenzone_outputs['amino_acids'] -= curr_hiddenzone_outputs['amino_acids'] * share_leaf
                    curr_hiddenzone_outputs['fructan'] -= curr_hiddenzone_outputs['fructan'] * share_leaf
                    curr_hiddenzone_outputs['proteins'] -= curr_hiddenzone_outputs['proteins'] * share_leaf
                    self.outputs['hiddenzone'][hiddenzone_id] = curr_hiddenzone_outputs

                    # Turn remobilizing flag to False
                    self.outputs['hiddenzone'][hiddenzone_id]['leaf_is_remobilizing'] = False

                # -- Remobilisation at the end of internode elongation
                # Internodes stop to elongate after leaves. We cannot test delta_internode_L > 0 for the cases of short internodes which are mature before GA production.
                if hiddenzone_inputs['internode_is_remobilizing']:

                    # Add to hidden part of the internode
                    hidden_internode_id = hiddenzone_id + tuple(['internode', 'HiddenElement'])
                    if hidden_internode_id not in self.outputs['elements'].keys():
                        new_internode_outputs = parameters.OrganInit().__dict__
                        self.outputs['elements'][hidden_internode_id] = new_internode_outputs
                    curr_hidden_internode_outputs = self.outputs['elements'][hidden_internode_id]
                    curr_hidden_internode_outputs['mstruct'] += curr_hiddenzone_outputs['internode_enclosed_mstruct']
                    curr_hidden_internode_outputs['max_mstruct'] = curr_hidden_internode_outputs['mstruct']
                    curr_hidden_internode_outputs['Nstruct'] += curr_hiddenzone_outputs['internode_enclosed_Nstruct']
                    curr_hidden_internode_outputs['sucrose'] += curr_hiddenzone_outputs['sucrose']
                    curr_hidden_internode_outputs['amino_acids'] += curr_hiddenzone_outputs['amino_acids']
                    curr_hidden_internode_outputs['fructan'] += curr_hiddenzone_outputs['fructan']
                    curr_hidden_internode_outputs['proteins'] += curr_hiddenzone_outputs['proteins']
                    curr_hidden_internode_outputs['is_growing'] = False
                    self.outputs['elements'][hidden_internode_id] = curr_hidden_internode_outputs

                    # Turn remobilizing flag to False
                    self.outputs['hiddenzone'][hiddenzone_id]['internode_is_remobilizing'] = False

                    #: Turn the flag to true after remobilisation in order to Delete Hiddenzone in both MTG and shared_outputs
                    self.outputs['hiddenzone'][hiddenzone_id]['is_over'] = True

        # --------------------------------
        # -------------- Roots -----------
        # --------------------------------
        if self.cnwheat_roots:
            for root_id, root_inputs in all_roots_inputs.items():
                curr_root_outputs = all_roots_outputs[root_id]

                axe_id = root_id[:2]
                curr_axis_outputs = all_axes_outputs[axe_id]

                # Temperature-compensated time (delta_teq)
                delta_teq = all_axes_inputs[axe_id]['delta_teq_roots']

                # Growth
                mstruct_C_growth, mstruct_growth, Nstruct_growth, Nstruct_N_growth = model.calculate_roots_mstruct_growth(root_inputs['sucrose'], root_inputs['amino_acids'],
                                                                                                                          root_inputs['mstruct'], delta_teq, postflowering_stages)
                # Respiration growth
                curr_root_outputs['Respi_growth'] = RespirationModel.R_growth(mstruct_C_growth)

                # Update of root outputs
                curr_root_outputs['mstruct'] += mstruct_growth
                curr_root_outputs['AA_consumption_mstruct'] = Nstruct_N_growth
                curr_root_outputs['sucrose_consumption_mstruct'] = model.calculate_roots_s_mstruct_sucrose(mstruct_growth, Nstruct_N_growth)
                curr_root_outputs['sucrose'] -= (curr_root_outputs['sucrose_consumption_mstruct'] + curr_root_outputs['Respi_growth'])
                curr_root_outputs['Nstruct'] += Nstruct_growth
                curr_root_outputs['amino_acids'] -= curr_root_outputs['AA_consumption_mstruct']
                curr_root_outputs['delta_mstruct_growth'] = mstruct_growth
                self.outputs['roots'][root_id] = curr_root_outputs

                # Update of axis outputs
                self.outputs['axes'][axe_id] = curr_axis_outputs
