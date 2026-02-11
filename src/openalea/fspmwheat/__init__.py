# -*- coding: latin-1 -*-
import logging
"""
    fspmwheat
    ~~~~~~~

    The model FSPM-Wheat.

    A Functional Structural Plant Model of Wheat.

    :copyright: Copyright 2014-2015 INRA-ECOSYS, see AUTHORS.
    :license: see LICENSE for details.

"""

__version__ = '2.0'

# Add a do-nothing handler to prevent an error message being output to sys.stderr in the absence of logging configuration
logging.getLogger(__name__).addHandler(logging.NullHandler())
