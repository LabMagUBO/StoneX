#!/opt/local/bin/python-3.4.3
# -*- coding: utf-8 -*-
"""
    Cycle class.
"""
## Global module
import numpy as np
import matplotlib.pylab as pl
from matplotlib.backends.backend_pdf import PdfPages       # Multiple page pdf

## StoneX module
from StoneX.Logging import *
from StoneX.Physics import *


class Cycle(object):
    """
        Class containing the cycle's data.

        self.model : model used to calculate the data
        self.T
        self.phi
        self.Ms
        self.data : array [H, Mt, Ml, theta_eq, alpha_eq (eventually)]
    """
    def __init__(self, Htab, cols):
        """
            Initiate the cycle.
            Arguments :
                — Htab : field array
                — cols : number of columns for the data (minimum 3 : H, Mt, Ml)

            self.data : H, Mt, Ml,
        """
        # Creating logger for all children classes
        self.logger = init_log(__name__)

        ## Check the number of cols
        if cols < 3:
            self.logger.error("Unsufficient number of columns.")
            self.logger.warn("Setting 3 columns instead.")
            cols = 3

        self.data = np.zeros((Htab.size, cols))

        # Storing the fields
        self.data[:, 0] = Htab

        # Energy
        self.energy = np.zeros(Htab.size, dtype=object)

    def info(self, sample, phi):
        """
            Save cycles conditions : T, phi, model
        """
        self.T = sample.T
        self.phi = phi
        self.model = sample.model
        self.Ms = sample.Ms

    def sum(self, cycl, selfdensity, density):
        """
            Method for summing two Cycle objects.
            Only sums Mt & Ml.
            The density is the weight factor for the second cycle.
        """
        try:
            if not (self.data[:, 0] == cycl.data[:, 0]).all():
                self.logger.error("The cycles does not have the same field values.")
                self.logger.warn("Summing anyway.")

        except AttributeError as err:
            self.logger.error("Error when summing the cycles. The data does not have the same shape.\n{}".format(err))
            self.logger.critical("Stopping the program.")
            sys.exit(0)

        # Summing Mt & Ml
        self.data[:, 1:3] = self.data[:, 1:3] * selfdensity + cycl.data[:, 1:3] * density

        # Erasing the rest
        self.data[:, 3:] = None

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
        ## Cycle figure
        # Creating figure
        fig = pl.figure()
        fig.set_size_inches(18.5,10.5)
        fig.suptitle("Model : {}, T = {}K, phi= {}deg".format(
            self.model, self.T, np.degrees(self.phi)
        ))

        # Initiating axis
        ax = fig.add_subplot(111)
        ax.grid(True)
        # ax.set_title("Model : {}, T = {}K, phi= {}deg".format(self.model, self.T, np.degrees(self.phi)))
        ax.set_ylabel("Mag. Moment (A.m²)")

        # First axis
        ax.plot(self.data[:, 0], self.data[:, 2], 'ro-', label='Ml')
        ax.plot(self.data[:, 0], self.data[:, 1], 'go-', label='Mt')
        ax.set_xlabel('Mag. Field H (A/m)')
        ax.legend()

        # Second axis (magnetic induction B)
        x1, x2 = ax.get_xlim()
        ax2 = ax.twiny()
        ax2.set_xlim(mu_0 * x1 * 1e3, mu_0 * x2 * 1e3)
        ax2.set_xlabel('Mag. Induction B (mT)')

        # New y axis
        ax3 = ax.twinx()
        y1, y2 = ax.get_ylim()
        ax3.set_ylim(y1 * 1e3 * 1e9, y2 * 1e3 * 1e9)
        ax3.set_ylabel('Mag. Moment (micro emu)')

        # Displaying mu_0 Hc1 and mu_0 Hc2 in millitestla
        Bc = np.abs((self.H_coer[1] - self.H_coer[0]) / 2) * mu_0
        Be = (self.H_coer[1] + self.H_coer[0]) / 2 * mu_0
        y_lim = ax.get_ylim()
        x_lim = ax.get_xlim()
        x_text = (x_lim[1] - x_lim[0]) * 0.15 + x_lim[0]
        y_text = (y_lim[1] - y_lim[0]) * 0.8 + y_lim[0]
        ax.text(
            x_text, y_text,
            "Hc = {:.3}mT\nHe = {:.3}mT".format(Bc * 1e3, Be * 1e3),
            style='italic',
            bbox={'facecolor': 'white', 'alpha': 1, 'pad': 10}
        )

        # Exporting graph as pdf
        file = "{0}/cycle_{1}_T{2}_phi{3}.pdf".format(
            path,
            self.model,
            round(self.T, 2),
            round(np.degrees(self.phi), 2)
        )
        pl.savefig(file, dpi=100)

        # Not forgetting to close figure (saves memory)
        pl.close(fig)

    def plot_energyPath(self, path):
        """
            Plotting the path taken by the magnetic moments depending on H.
        """
        # Create the PdfPages object to which we will save the pages:
        # The with statement makes sure that the PdfPages object is closed properly at
        # the end of the block, even if an Exception occurs.
        pdfFile = "{}/path_T{}_phi{}.pdf".format(path, self.T, np.round(np.degrees(self.phi), 2))
        with PdfPages(pdfFile) as pdf:

            # Determine if the energy landscape is 2D or 1d
            if len(self.energy[0].shape) == 2:  # 2D
                ## Plotting (theta, alpha) path
                Alpha = np.degrees(self.data[:, 4])
                Theta = np.degrees(self.data[:, 3])
                Num = Alpha.size/2

                # Creating figure, with title
                fig = pl.figure()
                fig.set_size_inches(18.5,10.5)
                fig.suptitle("Model : {}, T = {}K, phi= {}deg".format(self.model, np.round(self.T, 2), np.degrees(self.phi)))

                # Axis
                ax = fig.add_subplot(121, aspect='equal')

                ax.plot(Theta[:Num], Alpha[:Num], '-ro', label="Decreasing H")
                ax.plot(Theta[Num:], Alpha[Num:], '-bo', label="Increasing H")

                ax.set_xlabel("Theta M_f(deg)")
                ax.set_ylabel("Alpha M_af(deg)")

                ax.set_xlim(0, 360)
                ax.set_ylim(0, 360)

                ax.legend()
                ax.grid()

                # Independant angles
                ax2 = fig.add_subplot(122)
                ax2.plot(Theta, '-ro', label='θ Mf')
                ax2.plot(Alpha, '-bo', label='α Maf')

                ax2.set_xlabel("Field index")
                ax2.set_ylabel("Angle (deg)")
                ax2.grid()
                ax2.legend()

                # Saving the figure
                pdf.savefig()
                # Closing
                pl.close()

            elif len(self.energy[0].shape) == 1:
                # Resizing the energy table (in order to have a 2D-array)
                # Transform a 1Darray of 1Darrays in a single 2D array
                E = np.ma.zeros((self.energy.size, self.energy[0].size))
                for i, val in enumerate(self.energy):
                    E[i, :] = val

                # Field
                Hmax = self.data[0, 0]

                # Number of step
                Num = self.data[:, 0].size

                # Creating figure, with title
                fig = pl.figure()
                fig.set_size_inches(25, 10.5)
                fig.suptitle("Model : {}, T = {}K, phi= {}deg".format(self.model, np.round(self.T, 2) , np.round(np.degrees(self.phi), 2)))

                ## Plotting energy : increasing field
                ax = fig.add_subplot(121, aspect='equal')
                # All energy states with transparency
                cax_all = ax.imshow(E[:Num/2, :].data, alpha=0.5, label="Reachable states", interpolation = 'nearest', origin='upper', extent=(0, 360, -Hmax, Hmax), aspect='auto')
                # Reachable states
                cax_reach = ax.imshow(E[:Num/2, :], label="Reachable states", interpolation = 'nearest', origin='upper', extent=(0, 360, -Hmax, Hmax), aspect='auto')
                cbar = fig.colorbar(cax_reach)

                C = ax.contour(E[:Num/2, :].data, 10, colors='black', linewidth=.5, extent=(0, 360, Hmax, -Hmax))
                #ax.clabel(C, inline=1, fontsize=10)

                ax.plot(np.degrees(self.data[:Num/2, 3]), self.data[:Num/2, 0], 'ro', label='Eq.')

                #ax.set_title("id:{}".format(HIdx))
                ax.set_xlabel("Theta M_f(deg)")
                ax.set_ylabel("H")

                ## Plotting energy : decreasing field
                ax = fig.add_subplot(122, aspect='equal')
                # All energy states with transparency
                cax_all = ax.imshow(E[Num/2:, :].data, alpha=0.5, label="Reachable states", interpolation = 'nearest', origin='upper', extent=(0, 360, Hmax, -Hmax), aspect='auto')
                # Reachable states
                cax_reach = ax.imshow(E[Num/2:, :], label="Reachable states", interpolation = 'nearest', origin='upper', extent=(0, 360, Hmax, -Hmax), aspect='auto')
                cbar = fig.colorbar(cax_reach)

                C = ax.contour(E[Num/2:, :].data, 10, colors='black', linewidth=.5, extent=(0, 360, -Hmax, Hmax))
                #ax.clabel(C, inline=1, fontsize=10)

                #ax.imshow(E[Num/2:, :].mask, label="Reachable states", interpolation = 'nearest', origin='upper', extent=(0, 360, Hmax, -Hmax), aspect='auto')

                ax.plot(np.degrees(self.data[Num/2:, 3]), self.data[Num/2:, 0], 'ro', label='Eq.')

                #ax.set_title("id:{}".format(HIdx))
                ax.set_xlabel("Theta M_f(deg)")
                ax.set_ylabel("H")

                # Saving the figure
                pdf.savefig()
                # Closing
                pl.close()

            else:
                self.logger.warn("Cannot plot energy landscape.")

    def plot_energyLandscape(self, path):
        """
            Plotting the energy landscape of the cycles.
            If the landscape have more than 2 variables (H, theta), plotting the landscape for each H value.
        """

        # Create the PdfPages object to which we will save the pages:
        # The with statement makes sure that the PdfPages object is closed properly at
        # the end of the block, even if an Exception occurs.
        pdfFile = "{}/landscape_T{}_phi{}.pdf".format(path, self.T, np.round(np.degrees(self.phi), 2))
        with PdfPages(pdfFile) as pdf:

            # Determine if the energy landscape is 2D or 1d
            if len(self.energy[0].shape) == 2:  # 2D
                for i, line in enumerate(self.data):
                    H = line[0] #field
                    E = self.energy[i]      #2D array
                    theta = np.degrees(line[3])
                    alpha = np.degrees(line[4])

                    # Creating figure, with title
                    fig = pl.figure()
                    fig.set_size_inches(18.5,10.5)
                    fig.suptitle("Model : {}, T = {}K, phi= {}deg, H = {} Oe".format(self.model, self.T, np.degrees(self.phi), np.round(convert_field(H, 'cgs'), 2)))

                    # Axis
                    ax = fig.add_subplot(111, aspect='equal')

                    cax1 = ax.imshow(E.data, label="Energy landscape", interpolation = 'nearest', origin='upper', alpha=0.6, extent=(0, 360, 360, 0))
                    cax2 = ax.imshow(E, label="Reachable states", interpolation = 'nearest', origin='upper', extent=(0, 360, 360, 0))
                    cbar1 = fig.colorbar(cax1)
                    cbar2 = fig.colorbar(cax2)

                    C = ax.contour(E.data, 10, colors='black', linewidth=.5, extent=(0, 360, 0, 360))
                    #ax.clabel(C, inline=1, fontsize=10)

                    ax.plot(theta, alpha, 'ro', label='Eq.')

                    #ax.set_title("id:{}".format(HIdx))
                    ax.set_xlabel("Theta M_f(deg)")
                    ax.set_ylabel("Alpha M_af(deg)")

                    # Saving the figure
                    pdf.savefig()
                    # Closing
                    pl.close()


            elif len(self.energy[0].shape) == 1:
                Num = self.energy[0].size
                step = 360 / Num    # Step in degrees
                Theta = np.arange(Num) * step
                for i, line in enumerate(self.data):
                    H = line[0] #field

                    E = self.energy[i]     #1D array
                    theta = np.degrees(line[3])

                    # Creating figure, with title
                    fig = pl.figure()
                    fig.set_size_inches(18.5,10.5)
                    fig.suptitle("Model : {}, T = {}K, phi= {}deg, H = {} Oe".format(self.model, self.T, np.degrees(self.phi), np.round(convert_field(H, 'cgs'), 2)))

                    # Axis
                    ax = fig.add_subplot(111)

                    # Energy landscape
                    ax.plot(Theta, E.data, 'bo', label="Energy")
                    # Accessible states
                    ax.plot(Theta, E, 'go', label="Accessible")
                    # Position of the local minimum
                    ax.plot(theta, E[theta / 360 * Num], 'ro', label='Eq.')
                    # Thermal energy from the equilibrium position
                    ax.plot(Theta, np.zeros(Num) + E[theta / 360 * Num] + np.log(f0 * tau_mes) * k_B * self.T, '-y', label='25k_B T')

                    # Legend / Axis
                    ax.set_xlabel("Theta M_f(deg)")
                    ax.set_ylabel("E (J)")
                    ax.legend()
                    pl.grid()

                    # Saving and closing
                    pdf.savefig()
                    pl.close()
            else:
                self.logger.warn("Cannot plot energy landscape.")

    def export(self, path):
        """
            Export the data cycle to the path.
        """
        # Define the filename
        file = "{0}/cycle_{1}_T{2}_phi{3}.dat".format(path, self.model, round(self.T, 3), round(np.degrees(self.phi), 3) )

        # Info
        self.logger.info("Exporting cycle data : {}".format(file))

        # Exporting
        header = "Model = {0} \nT = {1}K \nphi = {2} deg \nMs = {3} A/m\n\n".format(self.model, self.T, np.degrees(self.phi), self.Ms)
        header +="H (A/m) \t\tMt (Am**2) \t\tMl (Am**2) \t\t(theta, alpha, ...) (rad)"
        np.savetxt(file, self.data, delimiter='\t', header=header, comments='# ')

    def import_data(self, path):
        """
            Import the data cycle from the given path.
        """
        # Define the filename
        file = "{0}/cycle_{1}_T{2}_phi{3}.dat".format(
            path,
            self.model,
            round(self.T, 3),
            round(np.degrees(self.phi), 3)
        )

        # Info
        self.logger.info("Importing cycle data : {}".format(file))

        # Exporting
        self.data = np.loadtxt(file)


class Rotation(object):
    def __init__(self, phiTab):
        """
            Initiate the rotation.
            Contains all the cycles.
            Arguments :
                — phiTab : phi array

            self.data : [cycle.phi, Hc, He, Mr1, Mr2, Mt1, Mt2]
        """
        # Creating logger for all children classes
        self.logger = init_log(__name__)

        # Storing phiTab
        self.phi = phiTab

        # Creating the cycles array
        self.cycles = np.empty(phiTab.size, dtype='object')

    def info(self, sample):
        """
            Save cycles conditions : T, model, ...
        """
        self.T = sample.T
        self.model = sample.model
        self.Ms = sample.Ms

    def sum(self, rot, selfdensity, density):
        """
            Method for summing two Rotation objects.
        """
        print("rot sum", density)
        # Checking if the Rotations are at the same temperature
        if (self.T != rot.T):
            self.logger.warn("Sum two rotations with different temperature.")

        try:
            if not (self.phi == rot.phi).all():
                self.logger.error("The rotations does not have the same phi values.")
                self.logger.warn("Summing anyway.")

        except AttributeError as err:
            self.logger.error("The rotations does not seem tohave the same phi array.\n{}".format(err))

        #if self is rot:
            #print("same rot")

        # Summing
        for i, cycle in enumerate(self.cycles):
            cycle.sum(rot.cycles[i], selfdensity, density)

    def plot(
        self,
        path, plot_azimuthal=True,
        plot_cycles=True,
        plot_energyPath=True,
        plot_energyLandscape=True
            ):
        """
            Plot the azimutal properties.
        """
        self.logger.info("Plotting rotationnal data in {} folder".format(path))

        # Plotting cycles graph
        for i, cycle in enumerate(self.cycles):
            # Cycle
            if plot_cycles:
                self.logger.info("Plotting cycles...")
                cycle.plot(path)
            # Plotting path taken by the magnetization
            if plot_energyPath:
                self.logger.info("Plotting energy path...")
                cycle.plot_energyPath(path)
            # Energy landscape
            if plot_energyLandscape:
                self.logger.info("Plotting energy landscape...")
                cycle.plot_energyLandscape(path)

            # Freeing memory
            if hasattr(cycle, 'energy') and not plot_energyPath:
                self.logger.debug("Freeing memory.")
                del(cycle.energy)

        # Plotting azimuthal
        if plot_azimuthal:
            # Plotting azimutal data
            file = "{0}/azimuthal_{1}_T{2}.pdf".format(
                path, self.model, round(self.T, 0)
            )
            data = self.data

            fig = pl.figure()
            fig.set_size_inches(18.5, 18.5)
            coer = fig.add_subplot(221, polar=True)
            coer.grid(True)
            coer.plot(
                data[:, 0], np.abs(data[:, 1]) * mu_0 * 1e3,
                'ro-', label='Bc (mT)'
            )
            coer.legend()

            ex = fig.add_subplot(222, polar=True)
            ex.grid(True)
            ex.plot(
                data[:, 0], np.abs(data[:, 2]) * mu_0 * 1e3,
                'bo-', label='Be (mT)'
            )
            ex.legend()

            rem = fig.add_subplot(224, polar=True)
            rem.grid(True)
            rem.plot(
                data[:, 0], np.abs(data[:, 3] / self.Ms),
                'mo-', label='Mr_1 / Ms'
            )
            rem.plot(
                data[:, 0], np.abs(data[:, 4] / self.Ms),
                'co-', label='Mr_2 / Ms'
            )
            rem.legend()

            trans = fig.add_subplot(223, polar=True)
            trans.grid(True)
            trans.plot(
                data[:, 0], np.abs(data[:, 5] / self.Ms),
                'go-', label='max(Mt)1 (A m**2)'
            )
            trans.plot(
                data[:, 0], np.abs(data[:, 6] / self.Ms),
                'yo-', label='max(Mt)2 (A m**2)'
            )
            trans.legend()

            #On trace en exportant
            self.logger.info("Plotting azimuthal graph : {}".format(file))
            pl.savefig(file, dpi=100)

            pl.close(fig)

    def process(self):
        """
            Calculate the properties for each cycle of the rotation.
            The results are stored in self.data
        """
        # Data rotation (phi, Hc, He, Mr1, Mr2, Mt_max, Mt_min)
        self.data = np.zeros((self.phi.size, 7))
        self.data[:, 0] = self.phi

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
        # Exporting the cycles
        for i, cycle in enumerate(self.cycles):
            cycle.export(path)

        # Exporting the azimuthal data
        file = "{0}/azimuthal_{1}_T{2}.dat".format(
            path,
            self.model,
            round(self.T, 0)
        )

        # Verbose
        self.logger.info("Exporting azimuthal data : {}".format(file))

        # Export
        header = "Model = {0} \nT = {1}K \nMs = {2} A/m".format(self.model, self.T, self.Ms)
        header += "phi(rad) \t\tHc (A/m) \t\tHe (A/m) \t\tMr1 (Am**2) \t\tMr2(Am**2) \t\tMt1 (Am**2) \t\tMt2 (Am**2)\n"
        np.savetxt(file, self.data, delimiter='\t', header=header, comments='# ')

    def import_data(self, path):
        """
            Import a data file
        """
        # Importing the cycles
        for i, cycle in enumerate(self.cycles):
            cycle.import_data(path)

        # Importing the azimuthal data
        file = "{0}/azimuthal_{1}_T{2}.dat".format(
            path,
            self.model,
            round(self.T, 0)
        )

        # Verbose
        self.logger.info("Importing azimuthal data : {}".format(file))

        # Import
        self.data = np.loadtxt(file)


class Tevol(object):
    """
        Class containing the thermal evolution of the different properties.
        This class is associated with a sample.

        data : (T, phi, [cycle.phi, Hc, He, Mr1, Mr2, Mt1, Mt2])
    """
    def __init__(self, vsm):
        """
            Initiating the data
        """
        # Creating logger for all children classes
        self.logger = init_log(__name__)

        # Attributes
        self.T = vsm.T
        self.data = np.zeros((vsm.T.size, vsm.phi.size, 7))
        self.model = vsm.sample.model

    def extract_data(self, rotations):
        # Loop over T
        for k, rot in enumerate(rotations):
            self.data[k, :, :] = rot.data

    def plot_H(self, path):
        """
            Plotting the T evolution of Hc and He for each field's direction
        """

        # Create the PdfPages object to which we will save the pages:
        # The with statement makes sure that the PdfPages object is closed properly at
        # the end of the block, even if an Exception occurs.
        pdfFile = "{}/Tevol_HcHe.pdf".format(path)
        with PdfPages(pdfFile) as pdf:
            # Loop over phi
            for k, phi in enumerate(self.data[0, :, 0]):
                # Creating figure
                fig = pl.figure()
                fig.set_size_inches(25,10.5)
                fig.suptitle("Model : {}, phi= {}deg".format(self.model, np.round(np.degrees(phi), 2) ))

                # Hc
                ax = fig.add_subplot(121)
                #ax.plot(self.T, convert_field(self.data[:, k, 1], 'cgs'), 'ro', label="Hc")
                ax.plot(self.T, self.data[:, k, 1] * mu_0 * 1e3, 'ro', label="Hc")
                ax.set_xlabel("T (K)")
                ax.set_ylabel("Hc (mT)")
                ax.legend()
                ax.grid()

                # He
                ax2 = fig.add_subplot(122)
                #ax2.plot(self.T, np.round(convert_field(self.data[:, k, 2], 'cgs'), 2), 'go', label="He")
                ax2.plot(self.T, self.data[:, k, 2] * mu_0 * 1e3, 'go', label="He")
                ax2.set_xlabel("T (K)")
                ax2.set_ylabel("He (mT)")
                ax2.legend()
                ax2.grid()

                # Saving and closing
                pdf.savefig()
                pl.close()

    def plot_M(self, path):
        """
            Plotting the T evolution of Mr and Mt for each field's direction
        """

        # Create the PdfPages object to which we will save the pages:
        # The with statement makes sure that the PdfPages object is closed properly at
        # the end of the block, even if an Exception occurs.
        pdfFile = "{}/Tevol_MrMt.pdf".format(path)
        with PdfPages(pdfFile) as pdf:
            # Loop over phi
            for k, phi in enumerate(self.data[0, :, 0]):
                # Creating figure
                fig = pl.figure()
                fig.set_size_inches(25,10.5)
                fig.suptitle("Model : {}, phi= {}deg".format(self.model, np.round(np.degrees(phi), 2) ))

                # Hc
                ax = fig.add_subplot(121)
                ax.set_title("Mr")
                ax.plot(self.T, np.abs(self.data[:, k, 3]), 'ro', label="|Mr1|")
                ax.plot(self.T, np.abs(self.data[:, k, 4]), 'go', label="|Mr2|")
                ax.set_xlabel("T (K)")
                ax.set_ylabel("Mr (A.m²)")
                ax.legend()
                ax.grid()

                # He
                ax2 = fig.add_subplot(122)
                ax2.set_title("Mt")
                ax2.plot(self.T, np.abs(self.data[:, k, 5]), 'ro', label="Mt1")
                ax2.plot(self.T, np.abs(self.data[:, k, 6]), 'go', label="Mt2")
                ax2.set_xlabel("T (K)")
                ax2.set_ylabel("Mt (A.m²)")
                ax2.legend()
                ax2.grid()

                # Saving and closing
                pdf.savefig()
                pl.close()

    def export(self, path):
        """
            Exporting the temperature evolution for each field direction
        """
        self.logger.info("Exporting temperature evolution.")
        # Loop over phi
        for k, phi in enumerate(self.data[0, :, 0]):
            # Verbose
            self.logger.info("\t phi = {} deg".format(np.degrees(phi)))

            # Filename
            file = "{0}/Tevol_{1}_phi{2}.dat".format(path, self.model, round(np.degrees(phi), 1))

            # Joining data
            T_vert = self.T.reshape((self.T.size, 1))
            tab = np.hstack((T_vert, self.data[:, k, 1:]))

            # Exporting
            header = "Model = {0} \nphi = {1} deg\n".format(self.model, np.degrees(phi))
            header +="Hc (A/m) \t\tHe (A/m) \t\tMr1 (Am**2) \t\tMr2 (Am**2) \t\tMt1 (Am**2) \t\tMt2 (Am**2)\n"
            np.savetxt(file, tab, delimiter='\t', header=header, comments='# ')
