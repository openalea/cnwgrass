# -*- coding: latin-1 -*-

import os
import logging
import logging.config
import json


"""
    cnmetabolism.tools
    ~~~~~~~~~~~~~

    This module provides tools to help for the validation of the outputs: 
    
        * plot of multiple variables on the same graph, 
        * set up of loggers,
        * quantitative comparison test,
        * and progress-bar to follow the evolution of long simulations.  

"""

OUTPUTS_INDEXES = ['t', 'plant', 'axis', 'metamer', 'organ', 'element']  #: All the possible indexes of CN-Metabolism outputs


def setup_logging(config_filepath='logging.json', level=logging.INFO,
                  log_model=False, log_compartments=False, log_derivatives=False, 
                  remove_old_logs=False):
    """Setup logging configuration.

    :param str config_filepath: The file path of the logging configuration.
    :param int level: The global level of the logging. Use either
          `logging.DEBUG`, `logging.INFO`, `logging.WARNING`, `logging.ERROR` or
          `logging.CRITICAL`.
    :param bool log_model: if `True`, log the messages from :mod:`cnmetabolism.model`. `False` otherwise.
    :param bool log_compartments: if `True`, log the values of the compartments. `False` otherwise.
    :param bool log_derivatives: if `True`, log the values of the derivatives. `False` otherwise.
    :param bool remove_old_logs: if `True`, remove all files in the logs directory documented in `config_filepath`.
    """
    if os.path.exists(config_filepath):
        with open(config_filepath, 'r') as f:
            config = json.load(f)
        if remove_old_logs:
            logs_dir = os.path.dirname(os.path.abspath(config['handlers']['file_info']['filename']))
            for logs_file in os.listdir(logs_dir):
                os.remove(os.path.join(logs_dir, logs_file))
        logging.config.dictConfig(config)
    else:
        logging.basicConfig()
    root_logger = logging.getLogger()
    root_logger.setLevel(level)

    cnmetabolism_model_logger = logging.getLogger('cnmetabolism.model')
    cnmetabolism_model_logger.disabled = not log_model  # set to False to log messages from openalea.cnwgrass.cnmetabolism.model
    logging.getLogger('cnmetabolism.compartments').disabled = not log_compartments  # set to False to log the compartments
    logging.getLogger('cnmetabolism.derivatives').disabled = not log_derivatives  # set to False to log the derivatives
