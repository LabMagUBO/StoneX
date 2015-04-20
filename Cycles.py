#!/opt/local/bin/python-3.4.3
# -*- coding: utf-8 -*-
"""
    Cycle class.
"""
## Global module
import numpy as np
import matplotlib.pylab as pl

## StoneX module
from StoneX.Logging import *


class Cycle(object):
    def __init__(self, Htab, cols):
        """
            Initiate the cycle.
            Arguments :
                — Htab : field array
                — cols : number of columns for the data (minimum 3 : H, Mt, Ml)

            self.data : H, Mt, Ml,
        """
        # Creating logger for all children classes
        self.logger = init_log(__name__, console_level='debug', file_level='info', mode='w')

        ## Check the number of cols
        if cols < 3:
            self.logger.error("Unsufficient number of columns.")
            self.logger.warn("Setting 3 columns instead.")
            cols = 3

        self.data = np.zeros((Htab.size, cols))

        # Storing the fields
        self.data[:, 0] = Htab

    def info(self, sample, phi):
        """
            Save cycles conditions : T, phi, model
        """
        self.T = sample.T
        self.phi = phi
        self.model = sample.model
        self.Ms = sample.Ms


    def calc_properties(self):
        """
            Calculate, depending of self.data :
                — Hc1, Hc2 : coercive fields
                — Mr1, Mr2 : remanent magnetizations
                — max(|Mt|) (decreasing and increasing field)
            H for the field, M the magnetic moment
            Zeros are calculated using linear regression.

            Store :
                — H_coer : two coercive fields (left, right)
                — Mr     : remanent magnetization
                — Mt_max : max of transverse
        """
        # Properties array
        self.H_coer = np.zeros(2)
        self.Mr = np.zeros(2)
        self.Mt_max = np.zeros(2)

        # Making views of the data
        H = self.data[:, 0]
        Mt = self.data[:, 1]
        Ml = self.data[:, 2]

        # Initial sign of magnetization
        sign_Ml = (Ml[0] > 0)
        sign_H = (H[0] > 0)

        # Loop over the field
        for i in np.arange(1, H.size):
            # if change of sign for Ml -> coercive field
            if (Ml[i] > 0) != sign_Ml:
                sign_Ml = not sign_Ml
                self.H_coer[sign_Ml] = H[i-1] - (H[i] - H[i-1])/(Ml[i] - Ml[i-1]) * Ml[i-1]

            # if change of sign for H -> remanent magnetization moment
            if (H[i] > 0) != sign_H:
                sign_H = not sign_H
                self.Mr[sign_H] = Ml[i-1] - (Ml[i] - Ml[i-1])/(H[i] - H[i-1]) * H[i-1]

        # Calcul du maximum, aller retour (doit être symétrique)
        half = H.size / 2
        self.Mt_max[0] = np.amax(np.absolute(Mt[:half]))
        self.Mt_max[1] = np.amax(np.absolute(Mt[half:]))


    def plot(self, path):
        """
            Trace et export le cycle Mt(H) et Ml(H)
        """
        # Setting the file name.
        file = "{0}/cycle_{1}_T{2}_phi{3}.pdf".format(path, self.model, round(self.T, 0), round(np.degrees(self.phi), 1) )

        Hc = round((self.H_coer[1] - self.H_coer[0]) / 2, 2)
        He = round((self.H_coer[1] + self.H_coer[0]) / 2, 2)

        #Création figure
        fig = pl.figure()
        fig.set_size_inches(18.5,10.5)
        #fig.suptitle("{0}".format(file))

        ax = fig.add_subplot(111)
        ax.grid(True)
        # Renseignements
        ax.plot(self.data[:, 0], self.data[:, 2], 'ro-', label='Ml')
        ax.plot(self.data[:, 0], self.data[:, 1], 'go-', label='Mt')
        ax.legend()

        y_lim = ax.get_ylim()
        x_lim = ax.get_xlim()
        x_text = (x_lim[1] - x_lim[0]) * 0.15 + x_lim[0]
        y_text = (y_lim[1] - y_lim[0]) * 0.8 + y_lim[0]
        ax.text(x_text, y_text, "Hc = {0}Oe\nHe = {1}Oe".format(Hc, He), style='italic', bbox={'facecolor':'white', 'alpha':1, 'pad':10})

        #On trace en exportant
        #self.logger.info("é : {}".format(file))
        pl.savefig(file, dpi=100)

        # Not forgetting to close figure (saves memory)
        pl.close(fig)


    def export(self, path):
        """
            Export the data.
        """
        pass

class Rotation(object):
    def __init__(self, phiTab):
        """
            Initiate the rotation.
            Contains all the cycles.
            Arguments :
                — phiTab : phi array
        """
        # Creating logger for all children classes
        self.logger = init_log(__name__, console_level='debug', file_level='info', mode='w')

        # Data rotation (phi, Hc, He, Mr1, Mr2, Mt_max, Mt_min)
        self.data = np.zeros((phiTab.size, 7))
        # Storing phiTab
        self.data[:, 0] = phiTab
        self.cycles = np.empty(phiTab.size, dtype='object')


    def info(self, sample):
        """
            Save cycles conditions : T, model, ...
        """
        self.T = sample.T
        self.model = sample.model
        self.Ms = sample.Ms


    def plot(self, path, plot_cycles=True):
        """
            Plot the azimutal properties.
        """
        self.logger.info("Plotting rotationnal data in {} folder".format(path))

        ## Plotting cycles
        if plot_cycles:
            for i, cycle in enumerate(self.cycles):

                cycle.plot(path)

        # Plotting azimutal data
        file = "{0}/azimuthal_{1}_T{2}.pdf".format(path, self.model, round(self.T, 0))
        data = self.data

        fig = pl.figure()
        fig.set_size_inches(18.5,18.5)
        coer = fig.add_subplot(221, polar=True)
        coer.grid(True)
        coer.plot(data[:, 0], np.abs(data[:, 1]), 'ro-', label='Hc (Oe)')
        coer.legend()

        ex = fig.add_subplot(222, polar=True)
        ex.grid(True)
        ex.plot(data[:, 0], np.abs(data[:, 2]), 'bo-', label='He (Oe)')
        ex.legend()

        rem = fig.add_subplot(224, polar = True)
        rem.grid(True)
        rem.plot(data[:,0], np.abs(data[:, 3] / self.Ms), 'mo-', label='Mr_1 / Ms')
        rem.plot(data[:,0], np.abs(data[:, 4] / self.Ms), 'co-', label='Mr_2 / Ms')
        rem.legend()

        trans = fig.add_subplot(223, polar = True)
        trans.grid(True)
        trans.plot(data[:, 0], np.abs(data[:, 5] / self.Ms), 'go-', label='max(Mt)1 (A m**2)')
        trans.plot(data[:, 0], np.abs(data[:, 6] / self.Ms), 'yo-', label='max(Mt)2 (A m**2)')
        trans.legend()

        #On trace en exportant
        self.logger.info("Exporting azimuthal plot : {}".format(file))
        pl.savefig(file, dpi=100)

        pl.close(fig)

    def process(self):
        """

        """
        for i, cycle in enumerate(self.cycles):
            cycle.calc_properties()

            # Storing cycles properties in rotation object
            Hc = (cycle.H_coer[1] - cycle.H_coer[0])/2
            He = (cycle.H_coer[1] + cycle.H_coer[0])/2
            Mr1 = cycle.Mr[0]
            Mr2 = cycle.Mr[1]
            Mt1 = cycle.Mt_max[0]
            Mt2 = cycle.Mt_max[1]
            self.data[i] = np.array([cycle.phi, Hc, He, Mr1, Mr2, Mt1, Mt2])

    def export(self, path):
        """
            Export the data.
        """
        pass
