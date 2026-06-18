"""
    turgorgrowth
    ~~~~~~~

    Turgor-Growth is a turgor-driven model of leaf growth. See:

        * :mod:`turgorgrowth.simulation`: the simulator (front-end) to run the model,
        * :mod:`turgorgrowth.model`: the state and the equations of the model,
        * :mod:`turgorgrowth.parameters`: the parameters of the model,
        * :mod:`turgorgrowth.postprocessing`: the post-processing and graph functions,
        * :mod:`turgorgrowth.tools`: tools to help for the validation of the outputs,
        * and :mod:`turgorgrowth.converter`: functions to convert Turgor-Growth inputs/outputs to/from Pandas dataframes.

"""

# Add a do-nothing handler to prevent an error message being output to sys.stderr in the absence of logging configuration
import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())
