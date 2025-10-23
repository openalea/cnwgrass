# -*- coding: latin-1 -*-
import logging

"""
    cnwheat
    ~~~~~~~

    The model CN-Wheat.
    
    CN-Wheat computes the CN exchanges in a wheat architecture. See:
    
        * :mod:`cnwheat.simulation`: the simulator (front-end) to run the model,
        * :mod:`cnwheat.model`: the state and the equations of the model,
        * :mod:`cnwheat.parameters`: the parameters of the model,
        * :mod:`cnwheat.postprocessing`: the post-processing and graph functions,
        * :mod:`cnwheat.tools`: tools to help for the validation of the outputs,
        * and :mod:`cnwheat.converter`: functions to convert CN-Wheat inputs/outputs to/from Pandas dataframes.

"""

__version__ = '3.0'

# Add a do-nothing handler to prevent an error message being output to sys.stderr in the absence of logging configuration
logging.getLogger(__name__).addHandler(logging.NullHandler())
