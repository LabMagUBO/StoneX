#!/opt/local/bin/python
# -*- coding: utf-8 -*-

"""
    Main program file.
    Example of Stoner Wohlfarth.

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
from matplotlib import pylab as pl

## Activate logging (cleaning the log on the first logger)
logger = StoneX.init_log(__name__, console_level='debug', file_level='debug', mode='w')
logger.info("Program version {}".format(StoneX.__version__))


################################################################################
# VSM CREATION
################################################################################
# First, lets create the VSM
vsm = StoneX.VSM()

# Set the vsm parameters
vsm.H = (0.5/StoneX.mu_0, 0.05/StoneX.mu_0, 0.15/StoneX.mu_0, 0.005/StoneX.mu_0, 'si')
vsm.phi = (0, 91, 15, 'deg')
vsm.T = (300, 1001, 1000, 'K')

# Plotting
vsm.plot_cycles = True
vsm.plot_azimuthal = True
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
# First, we create a specific domain, changing if necessary the parameters
# Models available : Stoner_Wohlfarth, Meiklejohn_Bean, Garcia_Otero, Franco_Conde, Rotatable_AF, Double_MacroSpin
domain = StoneX.create_domain(StoneX.Stoner_Wohlfarth, 'Xmpl_SW')

# Setting the temperature
domain.T = 0

# Theta step
domain.theta = (0.1, 'deg')

# Setting the physical parameters
domain.Ms = 1.4e6       # Magnetization in A/m
r = 2e-9               # radius of the particule, in meter
domain.V_f = 4/3*np.pi * r**3
domain.K_f =  1e5     # Uniaxial anisotropy, in J/m**3
domain.gamma_f = 0        # Uniaxial anisotropy direction, in radians
domain.K_bq = 0
domain.gamma_bq = 0
domain.K_iso = 0

print(domain)

## Then we create a sample based on the domain
# Distribution parameters
N = 1000
xMin = np.log(0.03)/np.log(10)
xMax = np.log(5) / np.log(10)
nx = 10
mu = 2
sigma = 10

# LogNormal distribution
R = np.logspace(xMin, xMax, nx)
X, m, s = StoneX.lognormale(R, mu, sigma)

logger.info("""Distribution de probabilité
        Normale : mu = {}, sigma = {}
        Log-normale : m = {}, s = {}""".format(mu, sigma, m, s))

# Creating the sample
Density = N * X
sample = StoneX.create_sample(domain, Density)

if True:
    pl.plot(R, np.around(X * N), '-ro')
    pl.savefig('{}/distrib.pdf'.format(domain.name), dpi=100)


for i, radius in enumerate(R):
    sample.domains[i].V_f = 4/3 * np.pi * (radius * 1e-9)**3
    #print(sample.domains[i].V_f)

################################################################################
# MEASUREMENTS
################################################################################
# We can measure the sample or one domain only
vsm.load(domain)
vsm.measure()




# END OF PROGRAM
