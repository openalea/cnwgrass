"""
    hydraulics
    ~~~~~~~

    Hydraulics is a turgor-driven model of leaf growth. See:

        * :mod:`hydraulics.simulation`: the simulator (front-end) to run the model,
        * :mod:`hydraulics.model`: the state and the equations of the model,
        * :mod:`hydraulics.parameters`: the parameters of the model,
        * :mod:`hydraulics.postprocessing`: the post-processing and graph functions,
        * :mod:`hydraulics.tools`: tools to help for the validation of the outputs,
        * and :mod:`hydraulics.converter`: functions to convert Hydraulics inputs/outputs to/from Pandas dataframes.

"""

# Add a do-nothing handler to prevent an error message being output to sys.stderr in the absence of logging configuration
import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())
