#!/opt/local/bin/python-3.4.3
# -*- coding: utf-8 -*-

# Version of the program
__version__ = '0.4, 2015-02-27'
version='ani'

## General modules
# Can't work without
import numpy as np

## StoneX modules
from StoneX.Logging import *            # Log functions & logger creation
from StoneX.Physics import *            # Physical constants and unit functions
from StoneX.VSM import *                # VSM class
from StoneX.Cycles import *
from StoneX.Sample import *             # Sample classes
from StoneX.Models import *             # Model classes
