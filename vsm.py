#!/opt/local/bin/python-3.4.3
# -*- coding: utf-8 -*-

"""
    VSM Class.

    Attributes :
        — H : returns the cycle field values
        — H(Hmax, Hstep, 'si' or 'cgs') : create H array [Hmax,...,-Hmax,..., Hmax] using Hstep in 'si' or 'cgs' system

    Properties :
        — phi : vsm direction of field, array type
        — phi(start, stop, step, 'deg' or 'rad') : create phi array [start,...,stop[ using step in 'deg' or 'rad' unit

    Methods :
        — self.load(sample) : load the sample into the VSM
        — self.measure() : measure all energy states and determine the magnetization cycle
        — self.process_cycles() : plot all the graphs
"""
import sys
import numpy as np
from StoneX.Physics import *
from StoneX.Logging import *
from StoneX.Cycles import *


class VSM(object):
    def __init__(self):
        self.logger = init_log(__name__)
        self.logger.debug("Init. Class VSM")

        ## Cycle's parameters
        # H and phi store ararys for cycle and rotations
        # Field values, default
        self.H = (100 * 1e3 / 4 / np.pi, 2**4, 'si')

        # Orientation
        self.phi = (0, 95, 10, 'deg')

        # Temperature
        self.T = (300, 301, 1, 'K')

        # Plotting parameters
        self.plot_cycles = True
        self.plot_azimuthal = True
        self.plot_energyLandscape = True
        self.plot_energyPath = True
        self.plot_T = True

        # Export parameters
        self.export_data = True
        #self.export_azimuthal = True
        #self.export_T = True

    def __str__(self):
        """
            Return the VSM status when using print(vsm)
        """

        str = """
        Field : max = {} Oe, nb step = {}, step = {} Oe
        Phi : start = {}, stop = {}deg, step = {}deg
        T : start = {}K, stop = {}K, step = {}K
        """.format( convert_field(self._H[0], 'cgs'),
                    self._H.size,
                    convert_field(self._H[0] - self._H[1], 'cgs'),
                    np.degrees(self._phi[0]),
                    np.degrees(self._phi[-1]),
                    np.degrees((self._phi[-1]-self._phi[0])/(self._phi.size-1)) if (self._phi.size > 1) else "???",
                    self._T[0],
                    self._T[-1],
                    (self._T[-1]-self._T[0])/(self._T.size-1) if (self._T.size > 1) else "???"
        )

        str += "\nPlotting :\n\t cycles = {} \n\t azimuthal = {} \n\t energy Path = {} \n\t energy Landscape = {} \n\t T evol = {}".format(self.plot_cycles, self.plot_azimuthal, self.plot_energyPath, self.plot_energyLandscape, self.plot_T)

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
                unit = 'si'

            if Hstep == 0:
                self.logger.warn("H.setter : need at least one value for H. Changing Hstep to 1 Oe.".format(unit))
                Hstep = convert_field(1, 'si')

            half = np.arange(- Hmax, Hmax + Hstep, Hstep)
            self._H = np.append(-half, half[1:])

    @property
    def phi(self, unit='rad'):
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
                unit = 'rad'

            if step == 0:
                self.logger.warn("phi.setter : step cannot be null. Setting step = 1 deg")
                step = np.radians(1)
            if stop == start:
                self.logger.warn("vsm.phi cannot be empty. Use start < step strictly. Using stop=step")
                stop = start + step
            self._phi = np.arange(start, stop, step)

    @property
    def T(self, unit='K'):
        return self._T

    @T.setter
    def T(self, val):
        """
            T setter. Need a tuple as argument with four variables (1 optionnal): start, stop, step, unit
        """
        try:
            start, stop, step, unit = val
        except ValueError:
            raise ValueError("vsm.T => need the following tuple : (start, stop, step, 'K' or 'degC')")
        else:
            # If no exception
            if unit == 'degC':
                start, stop, step = np.array([start + 273, stop + 273, step])
            elif unit != 'K':
                self.logger.warn("T.setter : unknown «{0}» system, using 'K' by default.".format(unit))
                unit = 'K'

            if step == 0:
                self.logger.warn("T.setter : step cannot be null. Setting step = 1 K")
                step = 1
            if stop == start:
                self.logger.warn("vsm.T cannot be empty. Use start < step strictly. Using stop=step")
                stop = start + step
            self._T = np.arange(start, stop, step)

    def load(self, sample):
        """
            Used to load the sample.
            «sample» is now callable inside the VSM class.
        """
        self.sample = sample

    def measure(self):
        #
        self.logger.info("Calculating energy...")
        #self.logger.info("Number of state to calculate : {}".format(self.H.size*self.phi.size*self.sample.theta.size*self.sample.alpha.size))
        self.sample.calculate_energy(self)

        self.logger.debug("Memory usage : sample {}Mib".format(sys.getsizeof(self.sample.E) * self.sample.E.size / 1024**2))
        self.logger.warn(type(self.sample.E[0, 0, 0, 0]))
        self.logger.warn(self.sample.E.shape)
        self.logger.warn(self.sample.E.size)
        #
        self.logger.info("Analysing_energy...")
        self.sample.analyse_energy(self)
        #
        self.process_cycles()

    @profile
    def process_cycles(self):
        """
            Analyse (ie calculate the coercive fields and others) the sample's cycles.
            All the calculated properties are stored in sample.attrib
            Calculated properties :
                — Hc1, Hc2
                — Mt1, Mt2
                – Mr1, Mr2
        """
        self.logger.info("Processing cycles...")
        for k, rot in enumerate(self.sample.rotations):
            rot.process()
            rot.plot(self.sample.name, plot_cycles=self.plot_cycles, plot_azimuthal=self.plot_azimuthal, plot_energyPath=self.plot_energyPath, plot_energyLandscape=self.plot_energyLandscape)

            if self.export_data:
                rot.export(self.sample.name)

        if self.plot_T :
            # Creating a sample's attribute to contain the (Hc, He, Mr, Mt…)
            self.sample.evolT = Tevol(self)

            self.sample.evolT.extract_data(self.sample.rotations)

            self.sample.evolT.plot_H(self.sample.name)

            if self.export_data:
                self.sample.evolT.export(self.sample.name)
