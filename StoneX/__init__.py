#!/opt/local/bin/env python
# -*- coding: utf-8 -*-
"""
    StoneX program.
    Numerical calculation of temperature dependent reversal cycles.

    Copyright (C) 2016  Jérôme Richy

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

__version__ = "1.4.3"
__date__ = "2016-07-28"
__author__ = "Jérôme RICHY"
__copyright__ = "Copyright 2016"
__license__ = "GPL"
__maintainer__ = "Jérôme RICHY"
__email__ = "jerome.richy@mines-nancy.org"
__status__ = "Development"


# General modules
import shutil
import time

# Modules for physics
import numpy as np
from matplotlib import pylab as pl

# StoneX modules import
from StoneX.Logging import *            # Log functions & logger creation
from StoneX.Physics import *            # Physical constants and unit functions
from StoneX.VSM import *                # VSM class
from StoneX.Cycles import *
from StoneX.Sample import *             # Sample classes
from StoneX.Models import *             # Model classes


disclaimer = """
    StoneX  Copyright (C) 2016  Jérôme Richy
    This program comes with ABSOLUTELY NO WARRANTY.
    This is free software, and you are welcome to redistribute it
    under certain conditions.
    See the LICENSE file for details.
"""

# Activate logging (cleaning the log on the first logger)
logger = first_init_log(
    __name__,
    console_level='debug',
    file_level='debug',
    mode='a'
)
logger.info(
    "Program version {}, {}\n{}".format(
        __version__,
        __date__,
        disclaimer
    )
)
