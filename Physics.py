#!/opt/local/bin/python-3.4.3
# -*- coding: utf-8 -*-

import numpy as np
#from functions import *
from StoneX.Logging import *

################################################
# Constants and variables definitions
################################################
# Recording the file name of the main script
#import __main__ as main
#main_file = main.__file__


# Physical constants
mu_0 = 4*np.pi*1e-7 #H/m
k_B = 1.3806488 * 1e-23 #m2 kg s-2 K-1

## Néel relaxation parameters
f0 = 1e9          #attempt frequency, in Hz
tau_mes = 100       #measurement mean time, in s


## Distribution function
def normale(x, m, s):
    """
        Normal distribution function.
        x is the variable, m the mean, s the standard deviation
    """

    return 1/np.sqrt(2*np.pi*s**2) * np.exp( -(x - m)**2 / 2 / s**2 )


def lognormale(x, m, s):
    """
        Lognormal distribution function.
        x is the variable, m and s respectively the normal mean and standard deviation
        mu the mean, sigma the lognormal standard deviation

        Return the log expected value and standard deviation as well.
    """
    if np.any(x == 0):
        logger = init_log(__name__, console_level='debug', file_level='info')
        logger.error("lognormale distribution: unable to divide by 0.")

    mu = np.log(m) - np.log(1 + s / m**2) / 2
    sigma = np.sqrt(np.log(1 + s / m**2))

    return 1/(x * sigma * np.sqrt(2 * np.pi)) * np.exp( -(np.log(x) - mu)**2 / 2 / sigma**2 ), mu, sigma

## Conversion unit functions
def convert_to(value, factor, system, logger):
    if system == 'si':
        return value * factor
    elif system == 'cgs':
        return value / factor
    else:
        logger.warn("this system does not exists")


def convert_field(value, system):
    """
        Convert «value» to the «system»
    """
    logger = init_log(__name__, console_level='debug', file_level='info')
    factor = 1e3 / 4 / np.pi
    return convert_to(value, factor, system, logger)



def convert_moment(value, system):
    """
        Convert «value» to the «system»
    """
    logger = init_log(__name__, console_level='debug', file_level='info')
    factor = 1e-3       #1 emu = 1e-3 A m**2
    return convert_to(value, factor, system, logger)
