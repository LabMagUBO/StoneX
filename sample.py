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
import copy
import numpy as np
from StoneX.Physics import *
from StoneX.Logging import *
from StoneX.Cycles import *



class Domain_base():
    def __init__(self):
        # Creating logger for all children classes
        self.logger = init_log(__name__)

        # Initiating class
        self.logger.debug("Init. Class Domain_base")

        # Volume
        self.V = 10e-9 * (100e-9)**2
        self.T = 0

    def __str__(self):
        txt = "Domain status :\n"
        txt += """
        Base
            V = {:.2} m**3 = {:.2} nm**3
            T = {} K
            25 k_B T = {} J
            """.format(self.V, self.V * 1e27, self.T, 25*k_B*self.T)
        return txt



class Ferro(Domain_base):
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
        self.K_iso = 0      #isotropy energy

        # Theta
        self.theta = (1, 'deg')

    def __str__(self):
        txt = super().__str__()
        txt += """
Ferro :
    theta : {} deg step
    Ms = {} A/m
    V_f = {:.3} m**3 = {:.3} nm**3
    K_f = {} J/m**3
    gamma_f = {} deg
    K_bq = {} J/m**3
    gamma_bq = {} deg
    K_f V_f = {} J
    K_f V_f / (25k_B) = T_B = {} K
                """.format(
                            np.degrees(self.theta[1] - self.theta[0]),
                            self.Ms, self.V_f,
                            self.V_f * 1e27,
                            self.K_f,
                            np.degrees(self.gamma_f),
                            self.K_bq,
                            np.degrees(self.gamma_bq),
                            (self.K_f * self.V_f / 2),
                            self.K_f * self.V_f / ( 25 * k_B)
                        )
        return txt


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
        isotropic = lambda phi, th: self.K_iso * self.V_f * np.sin(th - phi)**2
        zeeman = lambda phi, H, th: - mu_0 * H * self.V_f * self.Ms * np.cos(th - phi)
        uniaxial = lambda th: self.K_f * self.V_f * np.sin(th - self.gamma_f)**2
        biquadratic = lambda th: self.K_bq * self.V_f * np.sin(th - self.gamma_bq)**2 * np.cos(th - self.gamma_bq)**2
        return lambda phi, H, th: zeeman(phi, H, th) + uniaxial(th) + biquadratic(th) + isotropic(phi, th)


class AntiFerro(Domain_base):
    def __init__(self):
        super().__init__()
        self.logger.debug("Init. Class AntiFerro")

        self.gamma_af = 0
        self.V_af = self.V * 0.9

    def __str__(self):
        txt = super().__str__()
        txt += """
Anti-Ferro :
    V_af = {:.3} m**3
    gamma_af = {} deg
                """.format(self.V_af, np.degrees(self.gamma_af))
        return txt

    # No energy function.

class AntiFerro_Rotatable(AntiFerro):
    def __init__(self):
        super().__init__()

        self.logger.debug("Init. Class AntiFerro_rotatable")

        self.alpha = (1, 'deg')
        self.K_af = convert_field(2, 'si') * mu_0 * 400 * 1e-6 * 1e-3 / (1e-4 * 10 * 1e-9) / 2 * 100 # K_f * 10
        self.Bq_af = 0      # Biquadratic anisotropy

    def __str__(self):
        txt = super().__str__()
        txt +="""
    Rotatable:
        alpha = {} deg step
        K_af = {} J/m**3

        V_af K_af = {} J
        """.format(np.degrees(self.alpha[1] - self.alpha[0]), self.K_af, self.K_af * self.V_af)
        return txt

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
        biquadratic = lambda alph: self.Bq_af * self.V_af * np.sin(2*(alph - self.gamma_af))**2
        return lambda alph: uniaxial(alph) + biquadratic(alph)


class AntiFerro_Spin(AntiFerro_Rotatable):
    def __init__(self):
        super().__init__()

        self.logger.debug("Init. Class AntiFerro_rotatable")

        self.M_af = self.Ms * 50/100        # 1% of the F magnetization

    def __str__(self):
        txt = super().__str__()
        txt +="""
    AF magnetization:
        M_af = {} A/m
        """.format(self.M_af)
        return txt

    # Same property as AntiFerro_rotatable
    def energy(self):
        zeeman = lambda phi, H, alph: - mu_0 * H * self.V_af * self.M_af * np.cos(alph - phi)
        return lambda phi, H, alph: zeeman(phi, H, alph) + AntiFerro_Rotatable.energy(self)(alph)


def create_domain(model, name='test', mkdir=True):
    """
        Domain creation.
        List of available models
        — Stoner_Wohlfarth
        — Meiklejohn_Bean
        — Garcia_Otero
        — Franco_Conde
        — Rotatable_AF
        — ...
    """

    class Domain(model):
        def __init__(self):
            # Initialization of the subclass
            super().__init__()

            self.name = name

            if mkdir: self.create_folder(name)

        def __str__(self):
            txt = super().__str__()
            txt += """
Model : {}
                """.format(self.model)
            return txt

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


    return Domain()

def create_sample(domain, d):
    """
        Sample creation
        d : density of domain
    """
    class Sample(object):
        def __init__(self, density):
            # Creating logger
            self.logger = init_log(__name__)
            self.logger.info("Init. Sample class")

            self.model = domain.model
            self.name = domain.name
            self.domains = np.zeros(d.size, dtype=object)

            # Density: number for each domain
            self.density = density

            self.logger.info("Creating {} domain folders, under '{}/'".format(self.domains.size, self.name))
            for i in np.arange(self.domains.size):
                folder = "{}/domain{}".format(self.name, i)
                self.create_folder(folder, verbose=False)
                self.domains[i] = copy.copy(domain)
                self.domains[i].name = folder

        def __str__(self):
            txt = """
    Sample :
        Model : {}
        Number of domains : {}
        Folder : {}
                """.format(self.model, self.domains.size, self.name)
            return txt

        def measure(self, vsm):
            for i, domain in enumerate(self.domains):
                # Loading the domain into the VSM
                vsm.load(domain)

                # Measuring
                vsm.measure()

        def create_folder(self, name, verbose=True):
            """
                Create the folder for data exportation.
            """
            if not os.path.exists(name):
                os.makedirs(name)
                if verbose:
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

        def sum_domains(self):
            """
                Sum all the domains cycles.
            """
            # Creating a new rotations array for the sample from the first domain rotations (deepcopy needed for recurcive copy)
            self.rotations = copy.copy(self.domains[0].rotations)

            # Nulling the array by self-substracting
            #self.rotations.sum(self.domains[0].rotations, -1)

            # Adding each domains rotation cycles
            for i, domain in enumerate(self.domains):
                for j, rot in enumerate(self.rotations):
                    if i == 0:
                        # First nulling the array
                        rot.sum(self.domains[0].rotations[j], self.density[0], 0)
                    else:
                        rot.sum(self.domains[i].rotations[j], 1, self.density[i])

            # Then, recalculating each rotation parameters
            for i, rot in enumerate(self.rotations):
                rot.process()

    return Sample(d)
