#!/opt/local/bin/ipython
# -*- coding: utf-8 -*-

"""
    Main program.
    StoneX V0.4 compatible
"""
# Path to the module
module_path = '/Users/zorg/These/Documents/Programmes/Python_Modules/'

# Define module path
import sys
sys.path.append(module_path)

# Importing StoneX module
from StoneX.New_classes import *
from StoneX import logging_init

# Activate logging (cleaning the log on the first logger)
logger = logging_init.init_log(__name__, console_level='debug', file_level='info', mode='w')
#logger.info("Version du programme {}".format(StoneX.__version__))
logger.info("Nouvelle version")


# First, lets create the VSM
vsm = VSM()

# Set the vsm parameters
vsm.H = (5, 0.1, 'cgs')
vsm.phi = (10, 50, 50, 'deg')
print(vsm)

# Then the sample
#sample = create_sample(Stoner_Wohlfarth)
#sample = create_sample(Meiklejohn_Bean)
#sample = create_sample(Garcia_Otero)
sample = create_sample(Franco_Conde)
#sample = create_sample(Rotatable_AF)

# Set the sample parameters
sample.theta = (1, 'deg')
sample.T = 300

print(sample)

# Loading the sample into the VSM
vsm.load(sample)

vsm.measure()
