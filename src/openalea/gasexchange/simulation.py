# -*- coding: latin-1 -*-

from __future__ import division  # use "//" to do integer division

import warnings
import inspect

import numpy as np

from openalea.gasexchange import model
from openalea.gasexchange import parameters

"""
    gasexchange.simulation
    ~~~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`gasexchange.simulation` is the front-end to run the Gas-Exchange :mod:`model <gasexchange.model>`.

"""


class SimulationError(Exception):
    pass


class SimulationInputsError(SimulationError):
    pass


class Simulation:
    """The Simulation class permits to initialize and run a simulation.
    """

    def __init__(self, update_parameters=None, stomatal_model_name='BWB', hydraulics=False):
        """
        :param None or dict update_parameters: if a dict is provided, the specified parameters in keys will be updated.
        :param str stomatal_model_name: the model of stomatal conductance. Should be one of 'BWB', 'Leuning', 'Tuzet' or 'hydraulics'.
        :param bool hydraulics: if True the model will assume the coupling to the turgor-driven growth model.
        """

        #: `inputs` is a dictionary of dictionaries:
        #:     {(plant_index, axis_label, metamer_index, organ_label, element_label): {element_input_name: element_input_value, ...}, ...}
        self.inputs = {}
        #: the inputs needed by Gas-Exchange at element scale
        self.elements_inputs = ['width', 'height', 'PARa', 'nitrates', 'amino_acids', 'proteins', 'Nstruct',
                                         'green_area', 'sucrose', 'starch', 'fructan', 'PARa_prim', 'area_prim', 'Ci']
        if hydraulics:
            self.elements_inputs.append('water_potential')
        #: the inputs needed by Gas-Exchange at axis scale
        self.axes_inputs = ['SAM_temperature', 'height_canopy']

        #: `outputs` is a dictionary of dictionaries:
        #:     {(plant_index, axis_label, metamer_index, organ_label, element_label): {element_output_name: element_output_value, ...}, ...}
        self.outputs = {}
        #: the outputs computed by Gas-Exchange
        self.elements_outputs = ['Ag', 'An', 'Rd', 'Tr', 'Ts', 'gs', 'Ci', 'width', 'height']

        self.elements_inputs_outputs = set(self.elements_inputs + self.elements_outputs)

        #: Update parameters if specified
        if update_parameters:
            for key, value in update_parameters.items():
                if hasattr(parameters, key):
                    setattr(parameters, key, value)
                else:
                    warnings.warn(f"Parameter '{key}' is not defined in class self.model.parameters.")

        #: Stomatal conductance model
        if stomatal_model_name not in model.STOMATAL_MODELS_MAPPING:
            raise ValueError(
                f"Unknown stomatal conductance model."
                f"Must be one of : {list(model.STOMATAL_MODELS_MAPPING.keys())}"
            )
        HYDRAULIC_DEPENDENT_MODELS = ['Tuzet', 'hydraulics']
        if stomatal_model_name in HYDRAULIC_DEPENDENT_MODELS and not hydraulics:
            raise ValueError(
                f"Configuration error: the stomatal conductance model '{stomatal_model_name}' "
                f"requires the hydraulics option to be True in order to calculate the water potentials."
            )

        # Gets the selected model
        self.stomatal_model = model.STOMATAL_MODELS_MAPPING[stomatal_model_name]
        self.gs_args_keys = inspect.signature(self.stomatal_model).parameters
        self.hydraulics = hydraulics

    def initialize(self, inputs):
        """
        Initialize :attr:`inputs` from `inputs`.

        :param dict inputs: Dictionary of two dictionaries :
                    - `elements` : The inputs by element.
                    - `axes` : The inputs by axis.
              `inputs` must be a dictionary with the same structure as :attr:`inputs`.

            See :meth:`Model.run <gasexchange.model.run>`
               for more information about the inputs.
        """
        self.inputs.clear()
        self.inputs.update(inputs)

    def run(self, Ta, ambient_CO2, RH, Ur):
        """
        Compute Farquhar variables for each element in :attr:`inputs` and put
        the results in :attr:`outputs`.

        :param float Ta: air temperature at t (degree Celsius)
        :param float ambient_CO2: air CO2 at t (µmol mol-1)
        :param float RH: relative humidity at t (decimal fraction)
        :param float Ur: wind speed at the top of the canopy at t (m s-1)
        """

        self.outputs.update({inputs_type: {} for inputs_type in self.inputs['elements'].keys()})

        for (element_id, element_inputs) in self.inputs['elements'].items():

            axis_id = element_id[:2]
            organ_label = element_id[3]

            axe_label = axis_id[1]
            if axe_label != 'MS':  # Calculation only for the main stem
                continue
            # In case it is an HiddenElement, we need temperature calculation.
            # Cases of Visible Element without geometry property (because too small) don't have photosynthesis calculation neither.
            if element_inputs['height'] is None or np.isnan(element_inputs['height']):
                Ag, An, Rd, Tr, gsw, Ci = 0., 0., 0., 0., 0., 0.

                Ts = self.inputs['axes'][axis_id]['SAM_temperature']
            else:
                Ts = Ta # Initial value of Ts (°C)
                Ci = element_inputs['Ci']      #: previous organ internal CO2 concentration (µmol mol-1) todo Ci = parameters.Ci_init_ratio * ambient_CO2 see with Victoria if we keep this
                height_canopy = self.inputs['axes'][axis_id]['height_canopy']
                water_potential = element_inputs.get('water_potential', None)

                if parameters.SurfacicProteins:
                    surfacic_photosynthetic_proteins = model.calculate_surfacic_photosynthetic_proteins(element_inputs['proteins'],
                                                                                                        element_inputs['green_area'])

                    surfacic_nitrogen = model.calculate_surfacic_nonstructural_nitrogen_Farquhar(surfacic_photosynthetic_proteins)

                else:
                    surfacic_nitrogen = model.calculate_surfacic_nitrogen(element_inputs['nitrates'],
                                                                          element_inputs['amino_acids'],
                                                                          element_inputs['proteins'],
                                                                          element_inputs['Nstruct'],
                                                                          element_inputs['green_area'])

                surfacic_NSC = model.calculate_surfacic_WSC(element_inputs['sucrose'], element_inputs['starch'], element_inputs['fructan'], element_inputs['green_area'])

                # Checks if the calculation are made at the whole element scale or at the primitive scale
                if parameters.prim_scale:
                    PARa_list = element_inputs['PARa_prim']
                    areas_list = element_inputs['area_prim']
                else:
                    PARa_list = [element_inputs['PARa']]  # Unique item in this case
                    areas_list = [1.0]  # Set to 1 for mean PARa calculation below

                # Run the model
                Ag_prim_list = []
                for PARa in PARa_list:
                    count = 0
                    # Iteration until convergence
                    while True:
                        prec_Ci, prec_Ts = Ci, Ts

                        # Farquhar model calculations
                        Ag, An, Rd = model.calculate_photosynthesis(PARa,
                                                                    surfacic_nitrogen, parameters.NSC_Retroinhibition, surfacic_NSC,
                                                                    Ts, Ci)

                        # Stomatal conductance to water
                        potentials_gs_args = {'Ag': Ag, 'An': An,
                                              'surfacic_nitrogen': surfacic_nitrogen,
                                              'ambient_CO2': ambient_CO2, 'RH': RH,
                                              'Ta': Ta, 'water_potential': water_potential}
                        filtered_gs_args = {k: v for k, v in potentials_gs_args.items() if k in self.gs_args_keys}
                        gsw = self.stomatal_model(**filtered_gs_args)

                        # Calculation of Ci
                        Ci = model.calculate_Ci(ambient_CO2, An, gsw)

                        # Calculation of Ts and Tr
                        Ts, Tr = model.organ_temperature(element_inputs['width'], element_inputs['height'], height_canopy,
                                                          Ur, PARa, gsw, Ta, Ts, RH, organ_label)
                        count += 1

                        if count >= 30:
                            if abs((Ci - prec_Ci) / prec_Ci) >= parameters.DELTA_CONVERGENCE:
                                print('{}, Ci cannot converge, prec_Ci= {}, Ci= {}'.format(organ_label, prec_Ci, Ci))
                            if prec_Ts != 0 and abs((Ts - prec_Ts) / prec_Ts) >= parameters.DELTA_CONVERGENCE:
                                print('{}, Ts cannot converge, prec_Ts= {}, Ts= {}'.format(organ_label, prec_Ts, Ts))
                            break
                        if abs((Ci - prec_Ci) / prec_Ci) < parameters.DELTA_CONVERGENCE and (
                                (prec_Ts == 0 and (Ts - prec_Ts) == 0) or abs(
                                (Ts - prec_Ts) / prec_Ts) < parameters.DELTA_CONVERGENCE):
                            break

                    #: Conversion of Tr from mm s-1 to mmol m-2 s-1 (more suitable for further use of Tr)
                    Tr = (Tr * 1E6) / parameters.MM_WATER  # Using 1 mm = 1kg m-2
                    #: Decrease efficiency of non-lamina organs
                    if organ_label != 'blade':
                        Ag = Ag * parameters.EFFICENCY_STEM
                    Ag_prim_list.append(Ag)

                    # Aggregation of photosynthesis at element scale
                    if not Ag_prim_list:
                        Ag = 0
                    else:
                        Ag = sum([Ag * area for Ag, area in zip(Ag_prim_list, areas_list)]) / sum(areas_list)

            element_outputs = {'Ag': Ag, 'An': An, 'Rd': Rd,
                               'Tr': Tr, 'Ts': Ts, 'gs': gsw, 'Ci': Ci,
                               'width': element_inputs['width'], 'height': element_inputs['height']}

            self.outputs[element_id] = element_outputs
