#!/opt/local/bin/env python
# -*- coding: utf-8 -*-

"""
    Physics tools.

    Copyright (C) 2016  Jérôme Richy
"""

import numpy as np
import time
# from functions import *
from StoneX.Logging import *

################################################
# Constants and variables definitions
################################################
# Recording the file name of the main script
# import __main__ as main
# main_file = main.__file__


# Physical constants
mu_0 = 4 * np.pi * 1e-7  # H/m
k_B = 1.3806488 * 1e-23  # m2 kg s-2 K-1

# Néel relaxation parameters
f0 = 1e9                 # attempt frequency, in Hz
tau_mes = 100            # measurement mean time, in s

# SI unit prefixes
prefixes = {
    'atto': 1e-18,
    'femto': 1e-15,
    'pico': 1e-12,
    'nano': 1e-9,
    'micro': 1e-6,
    'milli': 1e-3,
    'centi': 1e-2,
    'deci': 1e-1,
    '': 1,
    'deca': 1e1,
    'hecto': 1e2,
    'kilo': 1e3,
    'mega': 1e6,
    'giga': 1e9,
    'tera': 1e12
}


# Distribution function
def normale(x, m, s):
    """
        Normal distribution function.
        x is the variable, m the mean, s the standard deviation
    """

    return 1 / np.sqrt(2 * np.pi * s**2) * np.exp(-(x - m)**2 / 2 / s**2)


def lognormale(x, m, s):
    """
        Lognormal distribution function.
        Parameters:
            – x: variable
            — m: mean of x distribution
            — s: standard deviation of x distribution
            — mu: mean of log(x) distribution
            — sigma: standard deviation of log(x) distribution

        Returns:
            lognormal distribution x, mean of log(x),
            standard deviation of log(x)
    """
    # Logging info
    if np.any(x == 0):
        logger = init_log(__name__, console_level='debug', file_level='info')
        logger.error("lognormale distribution: unable to divide by 0.")

    # Parameters for log(x)
    mu = np.log(m) - np.log(1 + s**2 / m**2) / 2
    sigma = np.sqrt(np.log(1 + s**2 / m**2))

    # Results
    return 1/(
        x * sigma * np.sqrt(2 * np.pi)
        ) * np.exp(-(np.log(x) - mu)**2 / 2 / sigma**2), mu, sigma


# Conversion unit functions
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
    factor = 1e-3       # 1 emu = 1e-3 A m**2
    return convert_to(value, factor, system, logger)


def timer(startstop):
    """
        Timer to record or display the time, depending on the boolean startstop
        startstop=True : start the time
        startstop = False : stop the timer and return the time
    """
    global start_time
    if startstop:
        start_time = time.time()
    else:
        stop_time = time.time()
        dt = stop_time - start_time
        s = dt % 60
        m = np.int(((dt - s) / 60) % 60)
        h = np.int(((dt - 60 * m - s) / 3600) % 60)
        return " Total time : {}hour, {} minutes, {:.3} seconds".format(
            h, m, s
        )
