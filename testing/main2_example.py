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
import numpy as np

## Activate logging (cleaning the log on the first logger)
logger = StoneX.init_log(__name__, console_level='debug', file_level='debug', mode='w')
logger.info("Program version {}".format(StoneX.__version__))


################################################################################
# VSM CREATION
################################################################################
# First, lets create the VSM
vsm = StoneX.VSM()

# Set the vsm parameters
vsm.H = (30, 0.5, 'cgs')
vsm.phi = (5, 91, 100, 'deg')
vsm.T = (10, 1001, 100, 'K')

# Plotting
vsm.plot_cycles = True
vsm.plot_azimuthal = False
vsm.plot_energyPath = False
vsm.plot_energyLandscape = False    #Takes a lot of time
vsm.plot_T = True

# Export
vsm.export_data = True

# Displaying parameters
logger.info(vsm)


################################################################################
# SAMPLE CREATION
################################################################################
# First, we create a specific domain, changing if necessary the parameters
# Models available : Stoner_Wohlfarth, Meiklejohn_Bean, Garcia_Otero, Franco_Conde, Rotatable_AF, Double_MacroSpin
domain = StoneX.create_domain(StoneX.Rotatable_AF, 'sample')
domain.T = 100
domain.theta = (1, 'deg')

if domain.model == 'Rotatable_AF':
    domain.alpha = (2, 'deg')
    domain.S = (200e-9)**2
    domain.V_f = 10e-9 * domain.S
    domain.V_af = 80e-9 * domain.S
    domain.K_f = 400
    domain.K_af = 110
    domain.J_ex = 1e-6

    domain.Ms = 400 * 1e-6 * 1e-3 / (1e-4 * 10 * 1e-9)
    #sample.M_af = 1/1000 * sample.Ms * sample.V_af / sample.V_f

print(domain.V_af)

# Then we create a sample based on the domain
n = 50
sample = StoneX.create_sample(domain, n)

for i, t in enumerate(np.random.normal(50, 20, n)):
    logger.info("N: {}, thickness t = {}nm".format(i, t))
    sample.domains[i].V_af = domain.S * t *1e-9



################################################################################
# MEASUREMENTS
################################################################################
# We can measure the sample or one domain only
vsm.load(sample)
vsm.measure()




# END OF PROGRAM
