# -*- coding: latin-1 -*-

"""
    senescwheat.parameters
    ~~~~~~~~~~~~~~~~~~~

    Parameters used in the model of senescence.

    :copyright: Copyright 2014-2015 INRA-ECOSYS, see AUTHORS.
    :license: see LICENSE for details.

"""

CONVERSION_FACTOR_20_TO_12 = 0.45  # modified_Arrhenius_equation(12)/modified_Arrhenius_equation(20)

N_MOLAR_MASS = 14  #: Molar mass of nitrogen (g mol-1)

SENESCENCE_ROOTS_POSTFLOWERING = 3.5E-7 * CONVERSION_FACTOR_20_TO_12  #: Rate of root turnover at 12°C (s-1). Value at 20°C coming from Johnson and Thornley (1985), see also Asseng et al. (1997).
SENESCENCE_ROOTS_PREFLOWERING = 0  # TODO: should be ontogenic for vegetative stages, 0 in Asseng 1997, but not null in Johnson and Thornley

FRACTION_N_MAX = {'blade': 0.5, 'stem': 0.425}  # Threshold of ([proteins]/[proteins]max) below which tissue death is triggered

SENESCENCE_MAX_RATE = 0.2E-8 * CONVERSION_FACTOR_20_TO_12  # maximal senescence m² s-1 at 12°C (Tref)
SENESCENCE_LENGTH_MAX_RATE = SENESCENCE_MAX_RATE / 3.5e-3  # maximal senescence m s-1 at 12°C (Tref)

RATIO_N_MSTRUCT = {1: 0.02, 2: 0.02, 3: 0.02, 4: 0.02, 5: 0.0175, 6: 0.015, 7: 0.01, 8: 0.005, 9: 0.005, 10: 0.005, 11: 0.005}  # Residual Mass of N in 1 g of mstruct at full senescence of the blade (from experiment NEMA)
DEFAULT_RATIO_N_MSTRUCT = 0.005  #: default N content in total organ mass (senesced + green) if phytomer rank not found above

AGE_EFFECT_SENESCENCE = 400  #: Age-induced senescence (degree-day since leaf emergence calculated from elong-wheat as equivalent at 12°C)

MIN_GREEN_AREA = 0.5E-8  #: Minimal green area of an element (m2). Below this area, set green_area to 0.0.
