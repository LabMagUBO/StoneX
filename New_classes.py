#!/opt/local/bin/ipython-2.7
# -*- coding: utf-8 -*-

################################
# Importing modules
###############################
# General modules
import numpy as np
import sympy as sp
from sympy.abc import x

# My own modules
#from constants import *
#from functions import *
from StoneX.logging_init import *
#from StoneX.models import *
from StoneX.constants import *




# Domain class
class VSM(object):
    def __init__(self):
        self.logger = init_log(__name__, console_level='debug', file_level='info', mode='w')
        self.logger.debug("Init. Class VSM")

        self.H_field = 0

        self.phi = 1

        self.logger.info("VSM loaded")

    @property
    def H(self, unit='si'):
        return self._H

    @H.setter
    def H(self, H):
        self._H = H

    # to be used?
    def load_sample(sample):
        self.sample = sample


class Domain(object):
    def __init__(self):
        # Creating logger for all children classes
        self.logger = init_log(__name__, console_level='debug', file_level='info', mode='w')

        # Initiating class
        self.logger.debug("Init. Class Domain")

        # Volume
        self.V = 10e-9 * (100e-9)**2


class Ferro(Domain):
    """
        Ferromagnetic domain.
        Different energies are coded :
        — Zeeman
        — Uniaxial anisotropy
        — Biquadratic anisotropy
        — ...
    """
    def __init__(self):
        super().__init__()

        self.logger.debug("Init. Class Ferro")
        # Magnetization
        self.Ms = 100

        # M possible directions
        self.nb_theta = 2**3
        self.theta = np.linspace(0, 2*np.pi, self.nb_theta)

    def energy(self, vsm):
        zeeman = lambda x: x
        uniaxial = lambda x: x
        biquadratic = lambda x: x
        return lambda x: zeeman(x) + uniaxial(x) + biquadratic(x)


class AntiFerro(Domain):
    def __init__(self):
        super().__init__()
        self.logger.debug("Init. Class AntiFerro")

        self.J_ex = 1

    # No energy function.

class AntiFerro_rotatable(AntiFerro):
    def __init__(self):
        super().__init__()

        self.logger.debug("Init. Class AntiFerro_rotatable")

    def energy(self):
        exchange = lambda x, y: x * 0.5 + y**2
        uniaxial = lambda x, y: 0.1 * x + y*10
        return lambda x, y: exchange(x, y) + uniaxial(x, y)


class Stoner_Wohlfarth(Ferro):
    def __init__(self):
        super().__init__()

        self.logger.debug("Init. Class Stoner_Wohlfarth")
        # Naming the model
        self.model = "Stoner_Wohlfarth"

        # Intiiating subclasses

        self.logger.info("Model loaded : {}".format(self.model))

    # No need to define energy

    def relax(self, vsm):
        pass



class Meiklejohn_Bean(AntiFerro, Ferro):
    def __init__(self):
        super().__init__()
        self.model = "Meiklejohn_Bean"
        self.logger.debug("Création modèle {}".format(self.model))

        self.J_ex = 1

    def energy(self, vsm):
        exchange = lambda x: self.J_ex * x
        return lambda x: Ferro.energy(self, vsm)(x) + AntiFerro.energy(self)(x)

class Rotatable_AF(AntiFerro_rotatable, Ferro):
    def __init__(self):
        # Initiate all parents classes
        super().__init__()

        # Naming the model
        self.model = "Rotatable_AF"

        # Self parameters
        self.T = 300

    def energy(self, vsm):
        return lambda x, y: Ferro.energy(self, vsm)(x) + AntiFerro_rotatable.energy(self)(x, y)


def create_sample(model):
    """
        List of available models
        — Stoner_Wohlfarth
        — Meiklejohn_Bean
        — Garcia_Otero
        — Franco_Conde
        — Rotatable_AF
    """

    class Sample(model):
        def __init__(self):
            # Initialization of the subclass
            super().__init__()

        def __str__(self):
            return "Modèle de l'échantillon : {}".format(self.model)


    return Sample()
