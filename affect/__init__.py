"""
Tools for processing results of a finite element analysis.
"""

import sys

__author__ = 'kevincopps'

if sys.version_info[0] == 3 and sys.version_info[1] >= 6:
    pass
else:
    raise ImportError('Python version 3.6 or above is required for the affect module.')

from . import arithmetic
from . import connect
from . import display
from . import dynamics
from . import exodus
from . import util

del sys
