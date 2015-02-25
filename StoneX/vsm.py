#!/opt/local/bin/ipython-2.7
# -*- coding: utf-8 -*-

################################
# Importing modules
###############################
## General modules
# Folder management
import os
# Arrays
import numpy as np
#import sympy as sp
#from sympy.abc import x

# My own modules
from StoneX.logging_init import *
#from constants import *
from StoneX.functions import *
from StoneX.models import *


##############################
# Defining classes
##############################


class VSM(object):
    """
        VSM apparatus. All the magnetometer parameters are defined here.
        The angle phi is stored in radians, but is given in degrees.
        All the others units are form the Internatial System.
    """
    def __init__(self):
        """
            Initialisation function, with default values.
            H_field = 300 Oe        Actual field value
            phi = 0                 Sample angle in the VSM
            H_field_max = 100 Oe    Max field for hysteresis cycles
            n_H_field = 2**6        Cycles' number of points
        """
        # Setting logger
        self.logger = init_log(__name__, console_level='debug', file_level='info')

        # Main attributes
        self.H_field = 300 * 1e3 / 4 / np.pi
        self.phi = 0
        self.H_field_max = 100 * 1e3 / 4 / np.pi
        self.n_H_field = 2**6

        # Output folder (empty by default)
        self.dos_out = ''

        #A little bit chatty
        self.logger.info("VSM created.")
        self.logger.warning("Don't forget to define the output folder")

    def __str__(self):
        """
            Return the VSM status
        """
        texte = "\n#############\n# VSM\n#############\n"
        texte += "VSM status\n  Field : {0} A/m ou {1} Oe\n  Angle phi : {2} deg\n".format(self.H_field, convert_field(self.H_field, 'cgs'), np.degrees(self.phi))
        texte += "VSM Parameters\n  H_field_max = {0} A/m ou {1} Oe\n  n_H_field = {2} \n".format(self.H_field_max, convert_field(self.H_field_max, 'cgs'), self.n_H_field)
        return texte

    def set_T(self):
        return "Ne fonctionne pas"

    def set_field(self, H):
        """ Set the H_field parameters (need to be in Oe)"""
        self.H_field = convert_field(H, 'si')

    def set_angle(self, phi):
        """ Set the Phi angle (convert it directly in degrees) """
        self.phi = np.radians(phi)

    def set_output(self, folder='output'):
        """
            Define the output folder for data and pdf files.
            Need to be called after the VSM.
        """
        # Defining folders varibles
        # Need to be
        self.dos_out = folder
        self.dos_cycles = folder + "/cycles"
        self.dos_data = folder + "/dat"
        self.dos_pdf = folder + "/pdf"
        self.dos_energy = folder + "/energy"

        # Logging
        self.logger.debug("Creating/Cleaning folders. Output folder: {0}".format(folder))

        # Creating/cleaning
        self.create_folders([self.dos_cycles, self.dos_data, self.dos_pdf, self.dos_energy])
        self.clean_folders([self.dos_cycles, self.dos_data, self.dos_pdf, self.dos_energy])

    def create_folders(self, folders):
        """ Creation of the folders
        """
        for dir in folders:
            if not os.path.exists(dir):
                os.makedirs(dir)

    def clean_folders(self, folders):
        """ Clean all the given folders"""

        for dir in folders:
            for file in os.listdir(dir):
                file_path = os.path.join(dir, file)
                try:
                    os.remove(file_path)
                except Exception as e:
                    print(e)



    # Measure the sample
    def get_magnetization(self, sample):
        """
            Return the magnetic moment m_l and m_t, according to the VSM orientation, and depending of the probability distribution.
            Variables inside:
                - L : longitudinal moment
                - T : transverse moment
        """

        L = np.sum( np.cos(self.phi + sample.theta) * sample.probM ) * sample.M_f * sample.V_f
        T = np.sum( np.sin(self.phi + sample.theta) * sample.probM ) * sample.M_f * sample.V_f

        return L, T

    # Do measurements
    def do_cycle(self, sample, idx=0, display=False, export=True, display_energy=False, export_energy=False, verbose=True):
        """
            Do a hysteresis cycle on the sample, and return the table (H, ml, mt), with the coercitive fields (Hc_l, Hc_r),
            the transversal moment extrema (min(mt), max(mt)), the remanent moments (ml_rem_left, ml_rem_right).
            All are in SI units (A*m**2 for magnetic moments).
        """
        self.logger.info("Beginning cycle. Phi=%s deg", np.degrees(self.phi))
        self.logger.debug("Cycle parameters\n\t display={0}\n\t export={1}\n\t display_energy={2}\n\t export_energy={3}".format(display, export, display_energy, export_energy))
        self.logger.debug("VSM parameters \n\t-> phi=%s deg", np.degrees(self.phi))

        # Array which will contain the data cycle (H, ml, mt)
        tab = np.zeros((self.n_H_field * 2 + 1, 3))

        j = 0 # coercitive field index
        k = 0 # remanence index
        #last_Ml = sample.M_f   #last Longitudinal magnetization
        last_ml = self.get_magnetization(sample)[0],        #Actual transverse magnetization
        last_H = self.H_field_max+1                         #Field, higher than the cycle's first point
        H_coer = np.zeros(2, dtype=float)                   #Coercitive fields array initialization
        ml_rem = np.zeros(2, dtype=float)                   #Remanence longitudinal magnetization

        half_cycle = np.linspace(self.H_field_max, -self.H_field_max, self.n_H_field, endpoint=False)
        full_cycle = np.append([half_cycle, half_cycle[::-1]], [self.H_field_max])

        for i, self.H_field in enumerate(full_cycle):
            sample.apply(self)
            sample.relax()

            if export_energy: sample.export_energy(self, i)

            (ml, mt) = self.get_magnetization(sample)
            tab[i] = self.H_field, ml, mt

            # Recording the coercive fields
            if i != 0:                      #not the first point
                if last_ml * ml < 0:
                    H_coer[j] = last_H - last_ml * (self.H_field - last_H) / (ml - last_ml)
                    j += 1

                elif ml == 0:
                    H_coer[j] = self.H_field
                    j +=1
            if j > 2:  #Test in case...
                self.logger.error("there is more than two coercive field values")

            # Recording the remanent moment
            if self.H_field * last_H < 0:
                ml_rem[k] = last_ml + (ml - last_ml) / (self.H_field - last_H) * (0 - last_ml)
                k += 1
            elif self.H_field == 0:
                ml_rem[k] = ml
                k += 1

            #Saving the current values
            last_ml = ml
            last_H = self.H_field

        # Searching the extremum of Mt
        min_mt = np.amin(tab[:, 2])
        max_mt = np.amax(tab[:, 2])
        mt_xtrm = np.array([min_mt, max_mt])

        # Final commands
        # Draw cycle if required
        if display: self.display_cycle(tab, sample)
        # Export the cycle if required
        if export: self.export_cycle(tab, sample, idx)

        return tab, H_coer, mt_xtrm, ml_rem

    def do_rotation(self, sample, phi_step = 5, phi_start = 0, phi_stop = 360, export=True, display=True, export_cycle=True, display_cycle=False):
        """
            Do an azimuthal measurement.
            Export an array with phi, Hc, He
        """
        self.logger.info("Beginning rotation")

        nb_values = int((phi_stop - phi_start)/phi_step) + 1
        #print nb_values

        tab = np.zeros((nb_values, 7))
        tab[:, 0] = np.arange(nb_values)*phi_step

        # Defining id lenght for the cycle files
        lenght = len(str(len(tab)))

        for idx, phi in enumerate(tab[:, 0]):
            self.logger.debug("Set Phi = %s", phi)
            self.set_angle(phi)

            # File Id
            fileId = str(idx).zfill(lenght)
            # Launching cycle
            cycle, H_coer, mt_xtrm, ml_rem = self.do_cycle(sample, idx=fileId, display=display_cycle, export=export_cycle, verbose=False)

            # Recording the coercive fields
            tab[idx, 1:3] = H_coer
            tab[idx, 3:5] = mt_xtrm
            tab[idx, 5:7] = ml_rem

            if export:
                self.export_rotation(tab, sample)

        if display: self.display_rotation(tab, sample)

        return tab

    #Cycle Graph
    def draw_cycle(self, cycle, sample):
        fig = pl.figure()
        ax = fig.add_subplot(111)
        ax.grid(True)
        fact_conv = 1e3 * mu_0      #convert to mT, hence mu_0 * H * 1000
        ax.plot(fact_conv * cycle[:, 0], cycle[:, 1], 'ro-', label='Ml')
        ax.plot(fact_conv * cycle[:, 0], cycle[:, 2], 'bo-', label='Mt')
        ax.plot(fact_conv * cycle[:, 0], np.sqrt(cycle[:, 2]**2 + cycle[:, 1]**2), 'go-', label='Total')
        ax.set_xlabel('Induction B (mT)')
        ax.set_ylabel('Mag. moment (A m**2)')

        # second axis
        #x1, x2 = ax.get_xlim()
        #ax2 = ax.twiny()
        #ax2.set_xlim(convert_field(x1, 'cgs')/fact_conv, convert_field(x2, 'cgs')/fact_conv)
        #ax2.set_xlabel('Field (Oe)')

        #new x axis
        #ax3 = ax.twinx()
        #y1, y2 = ax.get_ylim()
        #ax3.set_ylim(y1 * 1e3 * 1e9, y2 * 1e3 * 1e9)
        #ax3.set_ylabel('Mag. Moment (micro emu)')

        ax.legend()

    def display_cycle(self, cycle, sample):
        # Plot the graph
        self.draw_cycle(cycle, sample)
        # Display it
        pl.show()

    def export_cycle(self, cycle, sample, idx, unit='si'):
        """
            Export the plot in pdf format and the cycle data.
            By default, the data are in SI units.
            arguments:
                cycle : data to export and plot
                sample : the sample
                idx : id number for the pdf
                unit : 'si' by default. Can be changed to 'cgs'
        """
        # Plot the graph
        self.draw_cycle(cycle, sample)

        # Export it
        pl.savefig(self.dos_cycles + '/' + "cycle_n{0}_T{1}_phi{2}.pdf".format(idx, sample.T, round(np.degrees(self.phi)), 2))

        fileName = self.dos_cycles + '/' + "cycle_n{0}_T{1}_phi{2}.dat".format(idx, sample.T, round(np.degrees(self.phi), 2))
        self.logger.info("Exporting cycle graph: {0}".format(fileName))

        if unit == 'cgs':
            cycle[:, 0] = convert_field(cycle[:, 0], 'cgs')
            cycle[:, 1:3] = convert_moment(cycle[:, 1:3], 'cgs')
            columns_unit = 'B field(T)     mu_l(A m**2)       mu_t(A m**2)'
        elif unit == 'si':
            cycle[:, 0] *= mu_0
            columns_unit = 'B field(T)     mu_l(A m**2)       mu_t(A m**2)'
        else:
            self.logger.error('Unit not recognized. Use SI units.')
            columns_unit = 'Error : undefined units'

        header = """Cycle data
Magnetization : Mf = {0} A/m
Ferromagnetic volume : Vf = {1} m**3

{2}""".format(sample.M_f, sample.V_f, columns_unit)

        np.savetxt(fileName, cycle, header=header)
        pl.close()

    #Rotation graph
    def draw_rotation(self, rotation, sample):
        fig = pl.figure(figsize=(13, 13))
        coer = fig.add_subplot(221, polar=True)
        coer.grid(True)
        coer.plot(np.radians(rotation[:, 0]), np.abs((rotation[:, 1]-rotation[:,2])/2) * 1e3 * mu_0, 'ro-', label='Hc (mT)')
        coer.legend()

        ex = fig.add_subplot(222, polar=True)
        ex.grid(True)
        ex.plot(np.radians(rotation[:, 0]), np.abs((rotation[:, 1]+rotation[:,2])/2) * 1e3 * mu_0, 'bo-', label='He (mT)')
        ex.legend()

        trans = fig.add_subplot(223, polar = True)
        trans.grid(True)
        trans.plot(np.radians(rotation[:, 0]), np.abs(rotation[:, 3]), 'go-', label='min(Mt) (A m**2)')
        trans.plot(np.radians(rotation[:, 0]), np.abs(rotation[:, 4]), 'yo-', label='max(Mt) (A m**2)')
        trans.legend()

        rem = fig.add_subplot(224, polar = True)
        rem.grid(True)
        rem.plot(np.radians(rotation[:,0]), np.abs(rotation[:, 5]), 'mo-', label='Mr_1 (A m**2)')
        rem.plot(np.radians(rotation[:,0]), np.abs(rotation[:, 6]), 'co-', label='Mr_2 (A m**2)')
        rem.legend()


    def display_rotation(self, rotation, sample):
        # Plot the graph
        self.draw_rotation(rotation, sample)
        # Show it
        pl.show()

    def export_rotation(self, rotation, sample):
        # Plot the graph
        self.draw_rotation(rotation, sample)

        # Export it
        # Graph and data
        fileNameGraph = self.dos_pdf + '/' + "rotation_T{0}.pdf".format(sample.T)
        fileNameData = self.dos_data + '/' + "rotation_T{0}.dat".format(sample.T)

        pl.savefig(fileNameGraph)
        # And closing
        pl.close()
        # A little chat
        self.logger.info("Exporting rotation \n\t-> graph: {0} \n\t-> data: {1}".format(fileNameGraph, fileNameData))

        #Exporting data
        np.savetxt(fileNameData, rotation, header='theta (deg) \t Hc_l(Oe) \t Hc_r(Oe) \t Mt_min(microemu) \t Mt_max \t Ml_rem1 \t Ml_rem2')
