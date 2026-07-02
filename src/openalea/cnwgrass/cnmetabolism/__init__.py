# -*- coding: latin-1 -*-
import logging

"""
    cnmetabolism
    ~~~~~~~

    The model CN-Metabolism.
    
    CN-Metabolism computes the CN exchanges in a wheat architecture. See:
    
        * :mod:`cnmetabolism.simulation`: the simulator (front-end) to run the model,
        * :mod:`cnmetabolism.model`: the state and the equations of the model,
        * :mod:`cnmetabolism.parameters`: the parameters of the model,
        * :mod:`cnmetabolism.postprocessing`: the post-processing and graph functions,
        * :mod:`cnmetabolism.tools`: tools to help for the validation of the outputs,
        * and :mod:`cnmetabolism.converter`: functions to convert CN-Metabolism inputs/outputs to/from Pandas dataframes.

"""

# Add a do-nothing handler to prevent an error message being output to sys.stderr in the absence of logging configuration
logging.getLogger(__name__).addHandler(logging.NullHandler())
