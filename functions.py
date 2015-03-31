#!/opt/local/bin/ipython-2.7
# -*- coding: utf-8 -*-

import os
import numpy as np
import pylab as pl
from StoneX.logging_init import *
#from scipy import optimize, integrate

#from Stoner.constants import *
#from objects import *

###########################
# Folder management
#############################
def create_folders(folders):
    """ Creation of the folders
    """
    for dir in folders:
        if not os.path.exists(dir):
            os.makedirs(dir)
def clean_folders(folders):
    """ Clean all the given folders"""

    for dir in folders:
        for file in os.listdir(dir):
            file_path = os.path.join(dir, file)
            try:
                os.remove(file_path)
            except Exception as e:
                print(e)



def latex_fonts(test):
    """ Fonction to activate latex fonts.
        The generation of the graphics are much longer.

        Usage : latex_fonts(True) to activate.
    """
    if test:
        print("Loading latex font...")
        from matplotlib import rc
        pl.rc('text', usetex=True)
        pl.rc('font', family='serif', serif='cm10', style='normal', weight='medium')
        print("done.")



def search_eqProb(sample):
    """ Function searching the equilibrium probability and returning the mean magnetization depending of the sample's temperature.
        If T=0K, search the equilibrium using search_eq(sample).
        Usage : search_eqProb(sample)
        Output : M_l, M_t, theta_eq
    """

    if sample.T == 0:
        th_eq = search_eq(sample)
        return sample.M_f * np.cos(th_eq), sample.M_f * np.sin(th_eq), th_eq
    else:
        ener_red = lambda theta: np.exp(-sample.energy(theta) / k_B / sample.T)      # Reduced energy E/(k_b*T)
        Z, error_Z = integrate.quad(ener_red, 0, 2*np.pi)      # Partition function
        print("Z", Z)
        f = lambda theta: ener_red(theta) / Z
        prob_Ml = lambda theta: sample.M_f * f(theta) * np.cos(theta)
        prob_Mt = lambda theta: sample.M_f * f(theta) * np.sin(theta)
        prob_theta_eq = lambda theta: f(theta) * theta
        M_l, error_Ml = integrate.quad(prob_Ml, 0, 2*np.pi)
        M_t, error_Mt = integrate.quad(prob_Mt, 0, 2*np.pi)
        theta_eq, error_theta_eq = integrate.quad(prob_theta_eq, 0, 2*np.pi)

        return M_l, M_t, np.degrees(theta_eq % (2*np.pi)), f

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
