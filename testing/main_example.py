#!/opt/local/bin/ipython
# -*- coding: utf-8 -*-

"""
    Main program.
    StoneX V1.0 compatible
"""

################################################################################
# MODULES
################################################################################
## Path to the module
module_path = '/Users/zorg/These/Documents/Programmes/Python_Modules/'

## Define module path
import sys
sys.path.append(module_path)

# Importing StoneX module
from StoneX import *

## Activate logging (cleaning the log on the first logger)
logger = init_log(__name__, console_level='debug', file_level='info', mode='w')
#logger.info("Version du programme {}".format(StoneX.__version__))
logger.info("Nouvelle version")

################################################################################
# VSM & SAMPLE CREATION
################################################################################
# First, lets create the VSM
vsm = VSM()

# Set the vsm parameters
vsm.H = (20, 0.1, 'cgs')
vsm.phi = (1, 360, 5, 'deg')
print(vsm)

# Then the sample
sample = create_sample(Stoner_Wohlfarth)
#sample = create_sample(Meiklejohn_Bean)
#sample = create_sample(Garcia_Otero)
#sample = create_sample(Franco_Conde)
#sample = create_sample(Rotatable_AF, 'sample')

# Set the sample parameters
sample.theta = (0.1, 'deg')
sample.alpha = (0.1, 'deg')

# Anisotropy
sample.K_f = sample.K_f * 5
#sample.V_f = sample.V_f / 10
#sample.K_af = sample.K_af * 6
#sample.J_ex = 0
#print(sample.J_ex)
sample.T = 300

print(sample)



################################################################################
# MEASUREMENTS
################################################################################
# Loading the sample into the VSM
vsm.load(sample)
# Measuring
vsm.measure()

print(type(sample.rotation))


################################################################################
# PLOTTING
################################################################################
