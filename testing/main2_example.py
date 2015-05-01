#!/opt/local/bin/python
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
import StoneX

## Activate logging (cleaning the log on the first logger)
logger = StoneX.init_log(__name__, console_level='debug', file_level='debug', mode='w')
logger.info("Program version {}".format(StoneX.__version__))


################################################################################
# VSM CREATION
################################################################################
# First, lets create the VSM
vsm = StoneX.VSM()

# Set the vsm parameters
vsm.H = (20, 0.1, 'cgs')
vsm.phi = (0, 91, 100, 'deg')
vsm.T = (300, 1001, 1000, 'K')

# Plotting
vsm.plot_cycles = True
vsm.plot_azimuthal = False
vsm.plot_energyPath = False
vsm.plot_energyLandscape = False    #Takes a lot of time
vsm.plot_T = False

# Export
vsm.export_data = True

# Displaying parameters
logger.info(vsm)


################################################################################
# SAMPLE CREATION
################################################################################
# First, we create a specific domain
domain = StoneX.create_domain(StoneX.Meiklejohn_Bean, 'sample', mkdir=True)
domain.T = 100
#domain.K_f = 0
domain.J_ex *= 1

# Then we create a sample based on the domain
sample = StoneX.create_sample(domain, 10)

sample.domains[0].J_ex = 0
sample.domains[1].J_ex *= 1
sample.domains[2].J_ex *= 2
sample.domains[3].J_ex *= 3

################################################################################
# MEASUREMENTS
################################################################################
#sample.measure(vsm)
vsm.load(sample)
vsm.measure()


##
sample.sum_domains()

sample.rotations[0].cycles[0].plot('sample')

# END OF PROGRAM
