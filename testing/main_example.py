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
vsm.H = (40, 1, 'cgs')
vsm.phi = (0, 91, 40, 'deg')
vsm.T = (100, 800, 5000, 'K')

# Plotting
vsm.plot_cycles = True
vsm.plot_azimuthal = True
vsm.plot_energyPath = True
vsm.plot_energyLandscape = False    #Takes a lot of time
vsm.plot_T = True

# Export
vsm.export_data = True

# Displaying parameters
logger.info(vsm)



# Then the sample
#sample = create_sample(Stoner_Wohlfarth, 'sample')
#sample = create_sample(Meiklejohn_Bean, 'sample')
#sample = create_sample(Garcia_Otero, 'sample')
#sample = create_sample(Franco_Conde, 'sample')
#sample = create_sample(Rotatable_AF, 'sample')
sample = create_sample(Double_MacroSpin, 'sample2')

# Set the sample parameters
sample.theta = (1, 'deg')
sample.alpha = (2, 'deg')

# Anisotropy
#sample.K_f = sample.K_f * 2

#sample.V_f = sample.V_f / 10
#sample.K_af = sample.K_af * 2
#sample.J_ex *= 1
#sample.K_af /=
#print(sample.J_ex)
#sample.T = 100


# For thermal models
if True:
    sample.S = (200e-9)**2
    sample.V_f = 10e-9 * sample.S
    sample.V_af = 50e-9 * sample.S
    sample.K_f = 400
    sample.K_af = 110
    sample.J_ex = 1e-6

    sample.Ms = 400 * 1e-6 * 1e-3 / (1e-4 * 10 * 1e-9)
    sample.M_af = 1/1000 * sample.Ms * sample.V_af / sample.V_f


print(sample)
print("V_f", sample.V_f)
print("V_af", sample.V_af)
print("J_ex =", sample.J_ex)
print("K_f=", sample.K_f )
print("K_af =", sample.K_af)

print("M_af", sample.M_af)
print("sample.Ms", sample.Ms)

print("mu_af", sample.M_af * sample.V_af)
print("mu_f", sample.Ms * sample.V_f)

print("Energy")
print("25 k_B T(300K)", 300*k_B * np.log(tau_mes * f0))
print("mu0 H Ms V_f", mu_0 * vsm.H[0] * sample.Ms * sample.V_f)
print("K_f V_f", sample.K_f * sample.V_f)
print("K_af V_af", sample.K_af * sample.V_af)
print("J_ex S", sample.J_ex * sample.S)

################################################################################
# MEASUREMENTS
################################################################################
# Loading the sample into the VSM
vsm.load(sample)
# Measuring
vsm.measure()


# END OF PROGRAM
