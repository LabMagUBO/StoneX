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


# Fist, lets create the VSM
vsm = VSM()

# VSM settings
vsm.H_field = 10

# Then the sample
#sample = create_sample(Stoner_Wohlfarth)
#sample = create_sample(Meiklejohn_Bean)
sample = create_sample(Rotatable_AF)

# Loading the sample into the VSM
vsm.load(sample)




print(vsm.measure())
