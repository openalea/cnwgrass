# -*- coding: latin-1 -*-
import logging
"""
    integration
    ~~~~~~~

The integration and coupling of the different sub-models.

"""

# Add a do-nothing handler to prevent an error message being output to sys.stderr in the absence of logging configuration
logging.getLogger(__name__).addHandler(logging.NullHandler())
