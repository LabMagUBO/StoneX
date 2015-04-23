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
logger = init_log(__name__, console_level='debug', file_level='debug', mode='w')
#logger.info("Version du programme {}".format(StoneX.__version__))
logger.info("Nouvelle version")

################################################################################
# VSM & SAMPLE CREATION
################################################################################
# First, lets create the VSM
vsm = VSM()

# Set the vsm parameters
vsm.H = (20, 0.5, 'cgs')
vsm.phi = (0, 91, 10, 'deg')
vsm.T = (300, 301, 50, 'K')

# Plotting
vsm.plot_cycles = True
vsm.plot_azimuthal = True
vsm.plot_energyPath = True
vsm.plot_energyLandscape = True
vsm.plot_T = True

# Displaying parameters
logger.info(vsm)



# Then the sample
#sample = create_sample(Stoner_Wohlfarth, 'sample')
#sample = create_sample(Meiklejohn_Bean, 'sample')
#sample = create_sample(Garcia_Otero, 'sample')
sample = create_sample(Franco_Conde, 'sample')
#sample = create_sample(Rotatable_AF, 'sample')
#sample = create_sample(Double_Macrospin, 'sample')

# Set the sample parameters
sample.theta = (1, 'deg')
sample.alpha = (1, 'deg')

# Anisotropy
#sample.K_f = 0#sample.K_f * 5

#sample.V_f = sample.V_f / 10
#sample.K_af = sample.K_af * 4
#sample.J_ex *= 2
#sample.K_af /=
#print(sample.J_ex)
sample.T = 100

print(sample)



################################################################################
# MEASUREMENTS
################################################################################
# Loading the sample into the VSM
vsm.load(sample)
# Measuring
vsm.measure()



################################################################################
# PLOTTING
################################################################################
