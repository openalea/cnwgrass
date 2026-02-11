# -*- coding: latin-1 -*-

from __future__ import division  # use "//" to do integer division

from senescwheat import model
from senescwheat import parameters

"""
    senescwheat.simulation
    ~~~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`senescwheat.simulation` is the front-end to run the Senesc-Wheat :mod:`model <senescwheat.model>`.

    :copyright: Copyright 2014-2015 INRA-ECOSYS, see AUTHORS.
    :license: see LICENSE for details.

"""


class Simulation(object):
    """The Simulation class permits to initialize and run a simulation.
    """

    def __init__(self, delta_t=1, update_parameters=None, cnwheat_roots=True):

        #: The inputs of Senesc-Wheat.
        #:
        #: `inputs` is a dictionary of dictionaries:
        #:     {'roots': {(plant_index, axis_label): {roots_input_name: roots_input_value, ...}, ...},
        #:      'elements': {(plant_index, axis_label, metamer_index, organ_label, element_label): {element_input_name: element_input_value, ...}, ...}}
        self.inputs = {}

        #: The outputs of Senesc-Wheat.
        #:
        #: `outputs` is a dictionary of dictionaries:
        #:     {'roots': {(plant_index, axis_label): {roots_output_name: roots_output_value, ...}, ...},
        #:      'elements': {(plant_index, axis_label, metamer_index, organ_label, element_label): {element_output_name: element_output_value, ...}, ...}}
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

        :param dict inputs: The inputs by roots and element. `inputs` must be a dictionary with the same structure as :attr:`inputs`.
        """
        self.inputs.clear()
        self.inputs.update(inputs)

    def run(self, forced_max_protein_elements=None, opt_full_remob=False, postflowering_stages=False):
        """
        Compute Senesc-Wheat outputs from :attr:`inputs`, and update :attr:`outputs`.

        :param set forced_max_protein_elements: The elements ids with fixed max proteins.
        :param bool postflowering_stages: True to run a simulation with postflo parameter
        :param bool opt_full_remob: whether all proteins should be remobilised

        .. todo:: remove forced_max_protein_elements

        """

        if postflowering_stages:
            opt_full_remob = True

        self.outputs.update({inputs_type: {} for inputs_type in self.inputs.keys()})

        # axes
        all_axes_inputs = self.inputs['axes']

        # Roots
        if self.cnwheat_roots:
            all_roots_inputs = self.inputs['roots']
            all_roots_outputs = self.outputs['roots']
            for roots_inputs_id, roots_inputs_dict in all_roots_inputs.items():
                # Temperature-compensated time (delta_teq)
                delta_teq = all_axes_inputs[roots_inputs_id]['delta_teq_roots']

                # loss of mstruct and Nstruct
                rate_mstruct_death, rate_Nstruct_death = model.SenescenceModel.calculate_roots_senescence(roots_inputs_dict['mstruct'], roots_inputs_dict['Nstruct'], postflowering_stages)
                relative_delta_mstruct = model.SenescenceModel.calculate_relative_delta_mstruct_roots(rate_mstruct_death, roots_inputs_dict['mstruct'], delta_teq)
                delta_mstruct, delta_Nstruct = model.SenescenceModel.calculate_delta_mstruct_root(rate_mstruct_death, rate_Nstruct_death, delta_teq)
                # loss of cytokinins (losses of nitrates, amino acids and sucrose are neglected)
                loss_cytokinins = model.SenescenceModel.calculate_remobilisation(roots_inputs_dict['cytokinins'], relative_delta_mstruct)
                # Update of root outputs
                all_roots_outputs[roots_inputs_id] = {'mstruct': roots_inputs_dict['mstruct'] - delta_mstruct,
                                                      'senesced_mstruct': roots_inputs_dict['senesced_mstruct'] + delta_mstruct,
                                                      'rate_mstruct_death': rate_mstruct_death,
                                                      'Nstruct': roots_inputs_dict['Nstruct'] - delta_Nstruct,
                                                      'cytokinins': roots_inputs_dict['cytokinins'] - loss_cytokinins}

        # Elements
        all_elements_inputs = self.inputs['elements']
        all_elements_outputs = self.outputs['elements']
        for element_inputs_id, element_inputs_dict in all_elements_inputs.items():

            axe_label = element_inputs_id[1]
            if axe_label != 'MS':  # TODO: Calculation only for the main stem
                continue

            # Temperature-compensated time (delta_teq)
            axe_id = element_inputs_id[:2]
            delta_teq = all_axes_inputs[axe_id]['delta_teq']

            # Senescence
            element_outputs_dict = element_inputs_dict.copy()

            if model.SenescenceModel.calculate_if_element_is_over(element_inputs_dict['green_area'], element_inputs_dict['is_growing'], element_inputs_dict['mstruct']):
                element_outputs_dict['green_area'] = 0.0
                element_outputs_dict['senesced_length_element'] = element_inputs_dict['length']
                element_outputs_dict['mstruct'] = 0
                element_outputs_dict['senesced_mstruct'] += element_inputs_dict['mstruct']
                element_outputs_dict['is_over'] = True
            elif not element_inputs_dict['is_growing']:
                update_max_protein = forced_max_protein_elements is None or element_inputs_id not in forced_max_protein_elements

                if postflowering_stages:
                    new_green_area, relative_delta_green_area, max_proteins = model.SenescenceModel.calculate_relative_delta_green_area(element_inputs_id[3], element_inputs_dict['green_area'],
                                                                                                                                        element_inputs_dict['proteins'] / element_inputs_dict[
                                                                                                                                            'mstruct'],
                                                                                                                                        element_inputs_dict['max_proteins'], delta_teq,
                                                                                                                                        update_max_protein)

                    # Temporaire
                    new_senesced_length = relative_delta_green_area * (element_inputs_dict['length'] - element_inputs_dict.get('senesced_length_element', 0))

                else:
                    # Temporaire
                    new_senesced_length, relative_delta_senesced_length, max_proteins = model.SenescenceModel.calculate_relative_delta_senesced_length(element_inputs_id[3],
                                                                                                                                                       element_inputs_dict['senesced_length_element'],
                                                                                                                                                       element_inputs_dict['length'],
                                                                                                                                                       element_inputs_dict['proteins'] /
                                                                                                                                                       element_inputs_dict['mstruct'],
                                                                                                                                                       element_inputs_dict['max_proteins'], delta_teq,
                                                                                                                                                       update_max_protein)
                    # Senescence with element age
                    if element_inputs_id[3] != 'internode' and relative_delta_senesced_length == 0 and element_inputs_dict['age'] > parameters.AGE_EFFECT_SENESCENCE:
                        new_senesced_length, relative_delta_senesced_length, max_proteins = model.SenescenceModel.calculate_relative_delta_senesced_length(element_inputs_id[3],
                                                                                                                                                           element_inputs_dict['senesced_length_element'],
                                                                                                                                                           element_inputs_dict['length'],
                                                                                                                                                           0,
                                                                                                                                                           max_proteins, delta_teq,
                                                                                                                                                           update_max_protein)
                    # Temporaire :
                    relative_delta_green_area = relative_delta_senesced_length
                    new_green_area = element_inputs_dict['green_area'] * (1 - relative_delta_green_area)

                # Remobilisation
                N_content_total = model.SenescenceModel.calculate_N_content_total(element_inputs_dict['proteins'], element_inputs_dict['amino_acids'], element_inputs_dict['nitrates'],
                                                                                  element_inputs_dict['Nstruct'], element_inputs_dict['max_mstruct'], element_inputs_dict['Nresidual'])

                remob_starch = model.SenescenceModel.calculate_remobilisation(element_inputs_dict['starch'], relative_delta_green_area)
                remob_fructan = model.SenescenceModel.calculate_remobilisation(element_inputs_dict['fructan'], relative_delta_green_area)
                remob_proteins, delta_aa, delta_Nresidual = model.SenescenceModel.calculate_remobilisation_proteins(element_inputs_id[3], element_inputs_id[2], element_inputs_dict['proteins'],
                                                                                                                    relative_delta_green_area, N_content_total, opt_full_remob)
                loss_cytokinins = model.SenescenceModel.calculate_remobilisation(element_inputs_dict['cytokinins'], relative_delta_green_area)
                loss_nitrates = model.SenescenceModel.calculate_remobilisation(element_inputs_dict['nitrates'], relative_delta_green_area)

                # Loss of mstruct and Nstruct
                delta_mstruct, delta_Nstruct = model.SenescenceModel.calculate_delta_mstruct_shoot(relative_delta_green_area, element_inputs_dict['mstruct'], element_inputs_dict['Nstruct'])
                new_mstruct = element_inputs_dict['mstruct'] - delta_mstruct
                new_Nstruct = element_inputs_dict['Nstruct'] - delta_Nstruct

                delta_Nresidual += element_inputs_dict['Nstruct'] - new_Nstruct

                if new_mstruct == 0:
                    is_over = True
                else:
                    is_over = False

                # Turn 'is_over' to True when the element is fully senescent (to delete the element in the shared elements inputs/outputs)
                element_outputs_dict = {'green_area': new_green_area,
                                        'senesced_length_element': new_senesced_length,
                                        'mstruct': new_mstruct,
                                        'senesced_mstruct': element_inputs_dict['senesced_mstruct'] + delta_mstruct,
                                        'Nstruct': new_Nstruct,
                                        'starch': element_inputs_dict['starch'] - remob_starch,
                                        'sucrose': element_inputs_dict['sucrose'] + remob_starch + remob_fructan,
                                        'fructan': element_inputs_dict['fructan'] - remob_fructan,
                                        'proteins': element_inputs_dict['proteins'] - remob_proteins,
                                        'amino_acids': element_inputs_dict['amino_acids'] + delta_aa,
                                        'cytokinins': element_inputs_dict['cytokinins'] - loss_cytokinins,
                                        'nitrates': element_inputs_dict['nitrates'] - loss_nitrates,
                                        'max_proteins': max_proteins,
                                        'Nresidual': element_inputs_dict['Nresidual'] + delta_Nresidual,
                                        'N_content_total': N_content_total,
                                        'is_over': is_over}

            all_elements_outputs[element_inputs_id] = element_outputs_dict
