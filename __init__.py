#!/opt/local/bin/python-3.4.3
# -*- coding: utf-8 -*-

# Version of the program
__version__ = '1.0, 2015-06-08'

## General modules
import shutil
import time
# Can't work without
import numpy as np
from matplotlib import pylab as pl

## StoneX modules
from StoneX.Logging import *            # Log functions & logger creation
from StoneX.Physics import *            # Physical constants and unit functions
from StoneX.VSM import *                # VSM class
from StoneX.Cycles import *
from StoneX.Sample import *             # Sample classes
from StoneX.Models import *             # Model classes


## Activate logging (cleaning the log on the first logger)
logger = first_init_log(__name__, console_level='debug', file_level='debug', mode='w')
logger.info("Program version {}".format(__version__))
