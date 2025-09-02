"""
collisional_ot
==============
Main entry point of the library. Provides hierarchical access to different backends.
"""

from . import numpy
from . import pytorch
from . import cython
from . import utility

__all__ = ["numpy", "pytorch", "cython", "utility"]

