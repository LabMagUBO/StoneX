#!/opt/local/bin/python-3.4.3
# -*- coding: utf-8 -*-

"""
    Sample Class.
    Used to create the sample

    Attributes :


    Properties :

    Methods :

"""
import os
import shutil
import numpy as np
from StoneX.Physics import *
from StoneX.Logging import *
from StoneX.Cycles import *



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
        self.Ms = 400 * 1e-6 * 1e-3 / (1e-4 * 10 * 1e-9)      #magnetization of Py in A/m, 400 μemu / 1cm2 / 10nm
        self.V_f = (500 * 1e-9)**2 * 10e-9
        self.K_f = convert_field(2, 'si') * mu_0 * self.Ms / 2   #2 Oe uniaxial anisotropy
        self.gamma_f = 0
        self.K_bq = 0
        self.gamma_bq = 0

        # Theta
        self.theta = (1, 'deg')

    @property
    def theta(self):
        return self._theta

    @theta.setter
    def theta(self, val):
        """
            theta setter. Need a tuple as argument with three variables : step, unit
        """
        try:
            step, unit = val
        except ValueError:
            raise ValueError("vsm.theta => need the following tuple : (step, 'deg' or 'rad')")
        else:
            # If no exception
            if unit == 'deg':
                step = np.radians(step)
            elif unit != 'rad':
                self.logger.warn("theta.setter : unknown «{0}» system, using 'rad' by default.".format(unit))

            self._theta = np.arange(0, 2*np.pi, step)

    def energy(self):
        zeeman = lambda phi, H, th: - mu_0 * H * self.V_f * self.Ms * np.cos(th - phi)
        uniaxial = lambda th: self.K_f * self.V_f * np.sin(th - self.gamma_f)**2
        biquadratic = lambda th: self.K_bq * self.V_f * np.sin(th - self.gamma_bq)**2 * np.cos(th - self.gamma_bq)**2
        return lambda phi, H, th: zeeman(phi, H, th) + uniaxial(th) + biquadratic(th)


class AntiFerro(Domain):
    def __init__(self):
        super().__init__()
        self.logger.debug("Init. Class AntiFerro")

        self.gamma_af = 0
        self.V_af = self.V * 0.9

    # No energy function.

class AntiFerro_rotatable(AntiFerro):
    def __init__(self):
        super().__init__()

        self.logger.debug("Init. Class AntiFerro_rotatable")

        print("setting alpha array")
        self.alpha = (1, 'deg')
        self.K_af = convert_field(2, 'si') * mu_0 * 400 * 1e-6 * 1e-3 / (1e-4 * 10 * 1e-9) / 2 * 100 # K_f * 10

    @property
    def alpha(self):
        return self._alpha

    @alpha.setter
    def alpha(self, val):
        """
            alpha setter. Need a tuple as argument with three variables : step, unit
        """
        try:
            step, unit = val
        except ValueError:
            raise ValueError("vsm.alpha => need the following tuple : (step, 'deg' or 'rad')")
        else:
            # If no exception
            if unit == 'deg':
                step = np.radians(step)
            elif unit != 'rad':
                self.logger.warn("alpha.setter : unknown «{0}» system, using 'rad' by default.".format(unit))

            self._alpha = np.arange(0, 2*np.pi, step)

    def energy(self):
        uniaxial = lambda alph: self.K_af * self.V_af * np.sin(alph - self.gamma_af)**2
        return uniaxial


def create_sample(model, name='test'):
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

            self.name = name

            self.create_folder(name)

        def __str__(self):
            return "Sample's model: {}".format(self.model)

        def create_folder(self, name):
            """
                Create the folder for data exportation.
            """
            if not os.path.exists(name):
                os.makedirs(name)
                self.logger.info("Creating folder for plotting '{}'.".format(name))
            else:
                self.logger.info("Folder '{}' exists.".format(name))
                # On demande quoi faire
                yes_choice = np.array(['O', 'o', '', 'oui'])
                no_choice = np.array(['n', 'N', 'non'])
                answer = input("Voulez-vous l'effacer? [Y/n]\n?")

                if answer in yes_choice:
                    self.logger.info("Okayy... To the trash.")
                    shutil.rmtree(name)
                    self.logger.info("Creating folder '{}'.".format(name))
                    os.makedirs(name)
                elif answer in no_choice:
                    self.logger.info("Promise, I keep it.")
                else:
                    self.logger.error("Not the right answer.")
                    self.logger.warn("I quit.")
                    sys.exit(0)


    return Sample()
