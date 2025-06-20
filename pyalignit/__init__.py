# Copyright (C) 2021 OliverBScott
"""
pyalignit
---------
Python wrappers for the Align-it™ tool from Silicos-it.

"""
from .cpyalignit import *

from .core import *
from .draw import *
from .pymol import *
from .exclusion import *
from .pharfile import *

# Get version information from cpp
__version__ = GetVersion()
