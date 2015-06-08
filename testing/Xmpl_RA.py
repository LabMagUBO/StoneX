#!/opt/local/bin/python
# -*- coding: utf-8 -*-

"""
    Main program file.
    Example of Franco-Conde Model.

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
vsm.H = (250, 5, 'cgs')
vsm.phi = (0, 361, 1000, 'deg')
vsm.T = (10, 400, 50, 'K')

# Plotting
vsm.plot_cycles = True
vsm.plot_azimuthal = False
vsm.plot_energyPath = True
vsm.plot_energyLandscape = False    #Takes a lot of time and memory
vsm.plot_T = True

# Export
vsm.export_data = False

# Displaying parameters
logger.info(vsm)


################################################################################
# SAMPLE CREATION
################################################################################
# First, we create a specific domain, changing if necessary the parameters
# Models available : Stoner_Wohlfarth, Meiklejohn_Bean, Garcia_Otero, Franco_Conde, Rotatable_AF, Double_MacroSpin
domain = StoneX.create_domain(StoneX.Rotatable_AF, 'sampleRA')

# Setting the temperature
domain.T = 1

# Theta step
domain.theta = (1, 'deg')

# Alpha step
domain.alpha = (2, 'deg')

# Setting the physical parameters
domain.Ms = 1.1e6*5#1.4e6       # Magnetization in A/m
r_f = 4e-9               # radius of the particule, in meter
domain.V_f = 4/3*np.pi * r_f**3
domain.K_f =  3e4 /5    # Uniaxial anisotropy, in J/m**3
domain.gamma_f = 0        # Uniaxial anisotropy direction, in radians
domain.K_bq = 0
domain.gamma_bq = 0
domain.K_iso = 0

# AF parameters
d_af = 1                                        # Thickness of the AF shell, in nm
r_af = d_af * 1e-9 + r_f                        # radius of the AF shell, in m
domain.V_af = 4/3*np.pi * (r_af**3 - r_f**3)
domain.J_ex = 3e-7#3e-6#11e-3 * 1e-7 * 1e4
domain.S = 4 * np.pi * r_f**2
domain.K_af = 2e4

logger.info(domain)

## Then we create a sample based on the domain
# Distribution parameters
N = 10000                           # number of particules (average)
xMin = np.log(1.5)/np.log(10)
xMax = np.log(8) / np.log(10)
nx = 5                               # number of domain
mu = 4                              # mean particule radius
sigma = 2                          # standard deviation

# LogNormal distribution
R = np.logspace(xMin, xMax, nx)
X, m, s = StoneX.lognormale(R, mu, sigma)

logger.info("""Probability distribution : X
        Normale : mu = {}, sigma = {}
        Log-normale : m = {}, s = {}
        Mean : mean(X) = {}
        Number of particules : N = {}
        """.format(mu,
                    sigma,
                    m,
                    s,
                    np.sum(R * X)/np.sum(X),
                    np.sum(X * N)
                )
        )

# Creating the sample
Density = N * X
sample = StoneX.create_sample(domain, Density)

if True:
    pl.plot(R, X * N, '-ro')
    pl.savefig('{}/distrib.pdf'.format(domain.name), dpi=100)


for i, radius in enumerate(R):
    #print(i, radius)
    sample.domains[i].V_f = 4/3 * np.pi * (radius * 1e-9)**3
    sample.domains[i].S = 4 * np.pi * (radius * 1e-9)**2
    #sample.domains[i].V_af = 4/3 * np.pi * ( ((radius+d_af) * 1e-9)**3 - (radius * 1e-9)**3   )
    sample.domains[i].V_af = 4/3 * np.pi * ( (15 * 1e-9)**3 - (radius * 1e-9)**3   )
    logger.info(domain)

################################################################################
# DOMAIN'S PROPERTIES
################################################################################
n = sample.domains.size
#volume
#for i, dom in enumerate(sample.domains):
#    domains.



################################################################################
# MEASUREMENTS
################################################################################
# We can measure the sample or one domain only
vsm.load(sample)
vsm.measure()




# END OF PROGRAM
