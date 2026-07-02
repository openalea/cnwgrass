# -*- coding: latin-1 -*-

"""
    hydraulics.tools
    ~~~~~~~~~~~~~

    This module provides tools to help for the validation of the outputs:

        * set up of loggers,

"""

import os
import logging
import logging.config
import json


OUTPUTS_INDEXES = ['t', 'plant', 'axis', 'metamer', 'organ', 'element']  #: All the possible indexes of Turgor_Growth outputs


def setup_logging(config_filepath='logging.json', level=logging.INFO,
                  log_model=False, log_compartments=False, log_derivatives=False, 
                  remove_old_logs=False):
    """Setup logging configuration.

    :param str config_filepath: the file path of the logging configuration.
    :param int level: the global level of the logging. Use either `logging.DEBUG`, `logging.INFO`, `logging.WARNING`, `logging.ERROR` or `logging.CRITICAL`.
    :param bool log_model: if `True`, log the messages from :mod:`hydraulics.model`. `False` otherwise.
    :param bool log_compartments: if `True`, log the values of the compartments. `False` otherwise.
    :param bool log_derivatives: if `True`, log the values of the derivatives. `False` otherwise.
    :param bool remove_old_logs: if `True` remove old logs. `False` otherwise.
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

    hydraulics_model_logger = logging.getLogger('hydraulics.model')
    hydraulics_model_logger.disabled = not log_model  # set to False to log messages from hydraulics.model
    logging.getLogger('hydraulics.compartments').disabled = not log_compartments  # set to False to log the compartments
    logging.getLogger('hydraulics.derivatives').disabled = not log_derivatives  # set to False to log the derivatives
