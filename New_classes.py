#!/opt/local/bin/ipython-2.7
# -*- coding: utf-8 -*-

################################
# Importing modules
###############################
# General modules
import sys
import numpy as np
import sympy as sp
import scipy.ndimage as nd
from sympy.abc import x
from matplotlib import pylab as pl

# My own modules
#from constants import *
from StoneX.functions import *
from StoneX.logging_init import *
#from StoneX.models import *
from StoneX.constants import *




# Domain class
class VSM(object):
    def __init__(self):
        self.logger = init_log(__name__, console_level='debug', file_level='info', mode='w')
        self.logger.debug("Init. Class VSM")

        ## Cycle's parameters
        # H and phi store ararys for cycle and rotations
        # Field values, default
        self.H = (100 * 1e3 / 4 / np.pi, 2**4, 'si')

        # Orientation
        self.phi = (0, 95, 10, 'deg')

        #self.logger.info("""VSM loaded
        """Parameters :
            H = {0} A/m
            phi = {1} deg
            H_max = {2} A/m
            n_step = {3}

            No sample loaded. Use «vsm.load(sample)».
        """#.format(self.H_field, self.phi, self.H_field_max, self.n_H_field))

    def __str__(self):
        """
            Return the VSM status when using print(vsm)
        """

        str = """
        Field : max = {} Oe, nb step = {}, step = {} Oe
        Phi : start = {}, stop = {}deg, step = {}deg
        """.format( convert_field(self._H[0], 'cgs'),
                    self._H.size,
                    convert_field(self._H[0] - self._H[1], 'cgs'),
                    np.degrees(self._phi[0]),
                    np.degrees(self._phi[-1]),
                    np.degrees((self._phi[-1]-self._phi[0])/(self._phi.size-1)) if (self._phi.size > 1) else "???"
        )

        return str

    @property
    def H(self, unit='si'):
        return self._H

    @H.setter
    def H(self, val):
        """
            H setter. Need a tuple as argument with three variables : Hmax, Hstep, 'units'
        """
        try:
            Hmax, Hstep, unit = val
        except ValueError:
            raise ValueError("vsm.H => need the following tuple : (Hmax, Hstep, 'si' or 'cgs')")
        else:
            # If no exception
            if unit == 'cgs':
                Hmax, Hstep = convert_field(np.array([Hmax, Hstep]), 'si')
            elif unit != 'si':
                self.logger.warn("H.setter : unknown «{0}» system, using 'si' by default.".format(unit))

            if Hstep == 0:
                self.logger.warn("H.setter : need at least one value for H. Changing Hstep to 1 Oe.".format(unit))
                Hstep = convert_field(1, 'si')

            half = np.arange(- Hmax, Hmax + Hstep, Hstep)
            self._H = np.append(-half, half[1:])

    @property
    def phi(self, unit='deg'):
        return self._phi

    @phi.setter
    def phi(self, val):
        """
            Phi setter. Need a tuple as argument with three variables : start, stop, step, unit
        """
        try:
            start, stop, step, unit = val
        except ValueError:
            raise ValueError("vsm.phi => need the following tuple : (start, stop, step, 'deg' or 'rad')")
        else:
            # If no exception
            if unit == 'deg':
                start, stop, step = np.radians((start, stop, step))
            elif unit != 'rad':
                self.logger.warn("phi.setter : unknown «{0}» system, using 'rad' by default.".format(unit))

            if step == 0:
                self.logger.warn("phi.setter : step cannot be null. Setting step = 1 deg")
                step = np.radians(1)
            if stop == start:
                self.logger.warn("vsm.phi cannot be empty. Use start < step strictly. Using stop=step")
                stop = start + step
            self._phi = np.arange(start, stop, step)

    def load(self, sample):
        """
            Used to load the sample.
            «sample» is now callable inside the VSM class.
        """
        self.sample = sample

    def measure(self):
        self.sample.calculate_energy(self)
        self.sample.analyse_energy(self)
        self.process_cycles()

    def process_cycles(self):
        """
            Analyse (ie calculate the coercive fields and others) the sample's cycles and plot the graphs.
        """
        pass


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

        self.J_ex = 1
        self.alpha = 0

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

        # Initiating subclasses

        self.logger.info("Model loaded : {}".format(self.model))

    # No need to define energy

    def calculate_energy(self, vsm):
        """
            Calculate the energy for each parameters' values given by the vsm.
            The parameters are, in order : phi, H, theta
        """

        # sparse=True for saving memory
        # indexing : place the element in the index order
        Ph, Hf, Th = np.meshgrid(vsm.phi, vsm.H, self.theta, sparse=True, indexing='ij')

        # Creating a masked array. All values unmasked.
        self.E = np.ma.array(self.energy()(Ph, Hf, Th))

    def search_eq(self, E, eq):
        """
            Search the equilibrium state «eq» depending on the initial state «ini» and the energy landscape «E».
            Return the index of equilibrium
        """
        # stopping criteria
        found = False

        # Overflow security
        k = 0

        while not found:
            k += 1

            barrier = np.array([E[ (eq - 1) % E.size] - E[eq], E[ (eq + 1) % E.size] - E[eq]])

            if (barrier >= 0).all():
                # local minimum
                found = True
            elif barrier[0] == barrier[1]:
                self.logger.warn("Symmetrical unstable equilibrium. Choosing next largeur index (theta[idx+1]).")
                eq = (eq + 1) % E.size
            else:
                idx = np.argmin(barrier)
                eq = (eq + idx * 2 - 1) % E.size

            if k == 1000:
                self.logger.error("search_eq : Not converging.")
                break

        return eq

    def analyse_energy(self, vsm):
        """
            After calculating the energy depending on the parameters, search in which state is the magnetization.
            Returns the magnetization path : tab[phi, [H, theta, Mt, Ml]]
        """

        # Array containing all the data. Indexes : (phi, H, [H, theta_eq, Mt, Ml])
        self.data = np.zeros((vsm.phi.size, vsm.H.size, 4))

        # Index to access E(H) (Hmax -> Hmin -> Hmax)
        #Hindex = np.append(np.arange(vsm.H.size)[::-1], np.arange(vsm.H.size)[1:])

        # Index of the magnetization equilibrium, to start. Correspond to the global energy minimum
        eq = np.argmin(self.E[0, 0])

        # Loop over phi
        for i, phi in enumerate(vsm.phi):
            self.logger.debug("i= {}, Phi = {}deg".format(i, np.degrees(phi)))

            # Loop over H (max -> min -> max)
            for j, H in enumerate(vsm.H):
                self.logger.debug("j = {}, H = {}Oe".format(j, convert_field(H, 'cgs')))

                eq = self.search_eq(self.E[i, j, :], eq)

                self.data[i, j] = np.array([H, self.theta[eq], self.Ms * np.sin(self.theta[eq] - phi), self.Ms * np.cos(self.theta[eq] - phi)])

            pl.plot(convert_field(self.data[i, :, 0], 'cgs'), self.data[i, :, 3], '-ro')
            pl.plot(convert_field(self.data[i, :, 0], 'cgs'), self.data[i, :, 2], '-go')
            pl.show()


class Meiklejohn_Bean(Stoner_Wohlfarth, AntiFerro, Ferro):
    def __init__(self):
        super().__init__()
        self.model = "Meiklejohn_Bean"
        self.logger.debug("Creating model {}".format(self.model))

        self.S = (0.2*1e-9)**2
        self.J_ex = 0#11e-3 * 1e-7 * 1e4             #J/m**2

    def energy(self):
        exchange = lambda th: - self.J_ex * self.S * np.cos(th - self.alpha)
        return lambda phi, H, th: Ferro.energy(self)(phi, H, th) + exchange(th)


class Garcia_Otero(Meiklejohn_Bean):
    def __init__(self):
        # Initiating subclasses
        super().__init__()

        self.logger.debug("Init. Class Garcia_Otero")
        # Naming the model
        self.model = "Garcia_Otero"

        ## Additionnal parameters
        # Temperature
        self.T = 300     #K

        self.logger.info("Model loaded : {}".format(self.model))

    # No need to define energy, same as Meiklejohn_Bean

    def masking_energy(self, phiIdx, HIdx, eqIdx):
        """
            Method for masking all non-reachable states above the equilibrium.
            After searching all the available states, return the new equilibrium position.
            The non-reachable states are masqued in self.E
        """
        ok = False      # Test result
        while not ok:
            # Masking all the values above energy state equilibrium
            mask = self.E[phiIdx, HIdx, :] > self.E[phiIdx, HIdx, eqIdx] + 25 * k_B * self.T

            mask_lab, mask_num = nd.measurements.label(np.logical_not(mask))

            left = mask_lab[0]
            right = mask_lab[-1]
            # if the two extrem are reachable
            if left and right and (left != right):
                #replacement mask
                mask_replace = mask_lab == left
                mask_lab[mask_replace] = right

                # relabelling the mask
                labels = np.unique(mask_lab)
                mask_lab = np.searchsorted(labels, mask_lab)

            #Zone where the equilibrium is
            zone_eq = mask_lab[eqIdx]
            # Creating a new mask
            new_mask = mask_lab != zone_eq

            #Applying the mask
            self.E[phiIdx, HIdx, new_mask] = np.ma.masked

            # Searching the minimum
            mini = np.ma.min(self.E[phiIdx, HIdx, :])
            # If the minimum over the zone is the same that the eq energy
            if mini == self.E[phiIdx, HIdx, eqIdx]:
                # finish the loop
                ok = True
            else:
                # new equilibrium
                self.logger.warn("old eq {}".format(eqIdx))
                eqIdx = self.E[phiIdx, HIdx, :].argmin()
                self.logger.warn("new eq {}".format(eqIdx))



        if False:
            fig = pl.figure()
            ax = fig.add_subplot(111, aspect='equal')

            cax = ax.imshow(self.E[0], interpolation = 'nearest', origin='upper')
            cbar = fig.colorbar(cax)

            pl.show()

        return eqIdx

    def calculate_magnetization(self, phi, phiIdx, HIdx, eqIdx):
        """
            Calculate the magnetization, only using the local minimum state.
            Will be redefined in Franco_Conde to integrate all the available states.
            Arguments : self, phi, phiIdx, HIdx, eqIdx.
            Return Mt, Ml.
        """
        return self.Ms * np.sin(self.theta[eqIdx] - phi), self.Ms * np.cos(self.theta[eqIdx] - phi)

    # Redefining analyse_energy from Stoner_Wohlfarth
    def analyse_energy(self, vsm):
        """
            After calculating the energy depending on the parameters, search in which state is the magnetization.
            Returns the magnetization path : tab[phi, [H, theta, Mt, Ml]]
        """
        # Array containing all the data. Indexes : (phi, H, [H, theta_eq, Mt, Ml])
        self.data = np.zeros((vsm.phi.size, vsm.H.size, 4))

        # Index of the magnetization equilibrium, to start. Correspond to the global energy minimum
        eq = np.argmin(self.E[0, 0])

        # Loop over phi
        for i, phi in enumerate(vsm.phi):
            self.logger.debug("i= {}, Phi = {}deg".format(i, np.degrees(phi)))

            # Loop over H (max -> min -> max)
            for j, H in enumerate(vsm.H):
                # A little verbose
                self.logger.debug("j = {}, H = {}Oe".format(j, convert_field(H, 'cgs')))

                # Adjusting equilibrium
                eq = self.search_eq(self.E[i, j], eq)

                # Defining all the accessible states, depending on the temperature (changing equilibrium if necessary)
                eq = self.masking_energy(i, j, eq)

                # Calculating the magnetization
                Mt, Ml = self.calculate_magnetization(phi, i, j, eq)

                # Storing the results
                self.data[i, j] = np.array([H, self.theta[eq], Mt, Ml])


            print(self.E.mask)
            if True:
                pl.plot(convert_field(self.data[i, :, 0], 'cgs'), self.data[i, :, 3], '-ro')
                pl.plot(convert_field(self.data[i, :, 0], 'cgs'), self.data[i, :, 2], '-go')
                pl.grid()
                pl.show()
            if True:
                fig = pl.figure()
                ax = fig.add_subplot(111, aspect='equal')

                cax = ax.imshow(self.E[0], interpolation = 'nearest', origin='upper')
                cbar = fig.colorbar(cax)

                pl.show()


class Franco_Conde(Garcia_Otero):
    def __init__(self):
        # Initiating subclasses
        super().__init__()

        self.logger.debug("Init. Class Franco_Conde")
        # Naming the model
        self.model = "Franco_Conde"


        self.logger.info("Model loaded : {}".format(self.model))

    # Same energy function as Garcia_Otero

    # Redifining only calculate_magnetization
    def calculate_magnetization(self, phi, phiIdx, HIdx, eqIdx):
        """
            Calculate the magnetization, by integrating all the available states.
            Arguments : self, phi, phiIdx, HIdx, eqIdx.
            Return Mt, Ml.
        """
        # Energy array
        E = self.E[phiIdx, HIdx]

        # Probability array
        P = np.exp(- E / k_B / self.T)

        # Calculating magnetization
        Ml = np.ma.sum(P * np.cos(self.theta - phi)) / np.ma.sum(P)
        Mt = np.ma.sum(P * np.sin(self.theta - phi)) / np.ma.sum(P)

        return Mt, Ml


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
            return "Sample's model: {}".format(self.model)


    return Sample()
