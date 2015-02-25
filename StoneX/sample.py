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
from StoneX.models import *

##############################
# Defining classes
##############################

def load_class(subSample, model):
    """
        Load the model as well as define the sample.
        Usage : my_sample = load_class(Bilayer, Garcia_otero)
    """
    # We need to make the sample class global, to access it everywhere
    global Sample

    class Sample(subSample, model):
        """
            Main Sample class
        """

        def display_energy(self, nb = 1024):
            # vector theta
            vec_theta = np.linspace(0, 2*np.pi, nb)
            #figure creation
            fig = pl.figure()
            ax = fig.add_subplot(211)
            ax.grid(True)
            # graphs
            ax.set_title("Energy map")
            ax.plot(np.degrees(vec_theta), self.energy(vec_theta), 'r-', label='E')
            ax.plot(np.degrees(vec_theta), self.gradEnergy(vec_theta), 'b-', label='dE')
            ax.plot(np.degrees(vec_theta), self.ggradEnergy(vec_theta), 'g-', label='ddE')
            ax.legend()
            #ax.plot(np.degrees(self.theta_eq), self.energy(self.theta_eq), 'mo')
            #ax.set_xlabel(r'$\theta (degrees)$')
            #ax.set_ylabel(r'Energy (J)')

            #draw the temperature level
            T_sunami = k_B * self.T * np.log(f0 * tau_mes)
            #T_level = np.zeros(vec_theta.size) + self.energy(self.theta_eq) + T_sunami
            #ax.plot(np.degrees(vec_theta), T_level, 'm-')

            axHist = fig.add_subplot(212)
            axHist.grid(True)
            axHist.set_title("Probability")
            axHist.bar(np.degrees(self.theta), self.probM, width=360/(self.theta_n - 1))
            pl.show()

        def export_energy(self, vsm, idx, nb=1024):
            fig = pl.figure(101)
            axis = fig.add_subplot(211)
            vec_theta = np.linspace(0, 2*np.pi, nb)
            axis.plot(np.degrees(vec_theta), self.energy(vec_theta), 'r-', label='E')
            axis.plot(np.degrees(vec_theta), self.gradEnergy(vec_theta), 'b-', label='dE')
            axis.plot(np.degrees(vec_theta), self.ggradEnergy(vec_theta), 'g-', label='ddE')

            ## Plotting temperature level from the most probable M position
            # Thermal energy
            T_sunami = k_B * self.T * np.log(f0 * tau_mes)
            # Most probable position
            theta_eq = self.theta[np.argmax(self.probM)]
            # Relative thermal energy
            T_level = np.zeros(self.theta.size) + self.energy(theta_eq) + T_sunami
            axis.plot(np.degrees(self.theta), T_level, 'm-')


            axis.legend(bbox_to_anchor=(0.9, 0.95), loc=2)
            #axis.plot(np.degrees(self.theta_eq), self.energy(self.theta_eq), 'mo')
            axis.grid(True)
            axis.set_title("Energy")
            #axis.set_xlabel(r'$\theta (degrees)$')
            #axis.set_ylabel(r'Energy (J)')

            axHist = fig.add_subplot(212)
            axHist.grid(True)
            axHist.set_title("M Probability")
            #axHist.bar(np.degrees(self.theta), self.probM, width=360/(self.theta_n - 1))
            try:
                axHist.set_yscale('log')
                axHist.plot(np.degrees(self.theta), self.probM, '-ro')

            except ValueError:
                axHist.set_yscale('linear')
                axHist.plot(np.degrees(self.theta), self.probM, '-ro')


            pl.savefig(vsm.dos_pdf + '/' + "energy_n{0}_T{1}_H{2}_phi{3}.pdf".format(idx, self.T, round(convert_field(vsm.H_field, 'cgs'), 2), round(np.degrees(vsm.phi), 2)))

            tab = np.column_stack((vec_theta, self.energy(vec_theta), self.gradEnergy(vec_theta), self.ggradEnergy(vec_theta)))
            #export the data
            fileName = vsm.dos_energy + '/' + "energy_n{0}_T{1}_H{2}_phi{3}.dat".format(idx, self.T, round(convert_field(vsm.H_field, 'cgs'), 2), round(np.degrees(vsm.phi), 2) )
            print("Export energy : ", fileName)
            np.savetxt(fileName, tab, header="#theta \t E \t dE \t ddE")
            pl.close(101)

    return Sample()




class Domain(object):
    """
        Sub-sample class. Define the energy and paramaters of a F/AF domain.
        The class is define as a box, with two domain : F and AF.
        More informations later.
    """

    def __init__(self):
        self.logger = init_log(__name__, console_level='debug', file_level='info')
        self.logger.info("Sample created.")
        self.model()

        self.logger.warning("Don't forget to update its energy function (sample.apply(vsm))")
        # Definition of all the parameters
        #self.position = (0, 0)

        # Global constants
        self.T = 300                  #temperature
        self.V = (500 * 1e-9)**3    #box volume
        self.S = np.power(self.V, 2/3)               #F/AF surface area m**3
        self.fraction = self.S * 10e-9 / self.V      #Volume fraction of the ferromagnetic part


        #### Ferromagnetic part
        #Geometry, intensive constant
        #self.t_f = 10 * 1e-9                   #thickness
        self.V_f = self.V * self.fraction                # volume m**3, percentage of the total volume
        self.M_f = 400 * 1e-6 * 1e-3 / (1e-4 * 10 * 1e-9)      #magnetization of Py in A/m, 400 μemu / 1cm2 / 10nm
        #self.theta_eq = 0.                      # magnetization orientation
        #Probability of magnetization
        self.theta_n = 720
        self.theta_step = 2 * np.pi / self.theta_n
        self.theta = np.arange(0, 2*np.pi, self.theta_step)
        self.probM = np.zeros(self.theta_n)    #probability array of magnetization direction
        #initialization of the array
        self.probM[0] = 1

        # Magnetization (no need for that anymore, the VSM will measure it)
        #self.M_l = self.M_f
        #self.M_t = 0

        # Anisotropy
        #self.K_f = 5e3 * 1e-7 * 1e6             # anisotropy (from tony)
        self.K_f = convert_field(2, 'si') * mu_0 * self.M_f / 2   #2 Oe uniaxial anisotropy
        self.gamma_f = 0.            #uniaxial anisotropy angle

        self.K_bq = 0.           #biquadratic anisotropy
        self.gamma_bq = 0.       #bq ani. angle

        #### Antiferromagnetic layer
        #Geometry
        #self.t_af = 100e-9
        self.V_af = self.V * (1 - self.fraction)

        #Anisotropy
        self.alpha = 0.              #AF angle moment orientation
        self.K_af = 0.
        self.gamma_af = 0.           #AF anisotropy angle

        ####Interface exchange anisotropy
        self.J_ex = 11e-3 * 1e-7 * 1e4             #J/m**2
        #self.J_ex = mu_0 * self.M_f * self.V_f * convert_field(0.75, 'si') / self.S


    def __str__(self):
        """
            Print the status of the sample
        """
        txt = "\n############\nSample status\n############\n"
        txt += "Intrinsic properties:\n"
        #txt += "  Magnetization angle : {0} deg\n".format(np.degrees(self.theta_eq))
        txt += "  Temperature : {0} K\n".format(self.T)
        txt += "  Magnetization : {0} A/m\n".format(self.M_f)
        txt += "  Mag. Moment : {0} emu\n".format(self.M_f * self.V_f * 1e3)
        #txt += "  Theta Eq : {0} deg\n".format(np.degrees(self.theta_eq))
        txt += "\nGeometry: \n"
        #txt += "  t_f = {0} nm\n".format(1e9 * self.t_f)
        txt += "  V_f = {0} m**3\n".format(self.V_f)
        txt += "  gamma_f = {0} degrees\n".format(np.degrees(self.gamma_f))
        txt += "  V_af = {0} m**3\n".format(self.V_af)
        txt += "  S = {0} m**2\n".format(self.S)
        txt += "  f = {0} \n".format(self.fraction)
        txt += "\nEnergy: \n"
        txt += " J_ex : {0} J/m**2\n".format(self.J_ex)
        txt += " K_f : {0} J/m**3\n".format(self.K_f)
        return txt

    def print_energy(self, vsm):
        txt = "Energy comparison\n"
        txt += " Zeeman {0} J \n".format(mu_0 * vsm.H_field * self.V_f * self.M_f)
        txt += " Anisotropy {0} J \n".format(self.K_f * self.V_f)
        txt += " Exchange {0} J \n".format(self.J_ex * self.S)
        txt += " Temperature {0} J \n".format(k_B * self.T)
        txt += " Energy barrier {0} J\n".format(k_B * self.T * np.log(f0 * tau_mes))
        print(txt)

    """
        def set_theta_eq(self, theta):
        print ("New equilibrium angle : ", theta)
        self.theta_eq = np.radians(theta)
    """
    """
        def get_theta_eq(self):
        print ("Equilibrium angle : ", np.degrees(theta), " deg")
    """
    def set_F(self, frac=None, M_f=None, K_f=None , gamma_f=None):
        """
            Function to modify the ferromagnetic layer parameters.
            t_f in nm
            V_f in m**3
            M_f in SI
            K_f in SI
            gamma_f in degrees
        """


        if frac is not None:
            self.V_f = self.V * frac
            self.V_af = self.V * (1 - frac)

        if M_f is not None:
            self.M_f = M_f

        if K_f is not None:
            self.K_f = K_f

        if gamma_f is not None:
            self.gamma_f = np.radians(gamma_f)



    def set_T(self, T):
        if T < 0:
            print ("No way!!! Are you mad?\n Try again...")
        else:
            print ("New sample temperature : {0} K".format(T))
            self.T = T


    def apply(self, vsm):
        """
            Replace self.update_energy

            Update the energy functions using the parameters of the VSM and the sample.
            Sympy symbolic functions are started with sym_ :
                — sym_zeeman, sym_exchange, ...
            Numpy numeric function :
                — zeeman, exchange,

            The x variable correspond to theta.
        """#


        ### Using analytical functions (more quicker)
        self.define_energy_functions(vsm)
        self.define_gradEnergy_functions(vsm)
        self.define_ggradEnergy_functions(vsm)

        ### Using symbolic function
        #self.define_sym_energy_functions(vsm)
        #self.derivate_sym_energy_functions(vsm)
        #self.convert_sym_energy_functions(vsm)


    def define_energy_functions(self, vsm):
        self.zeeman = lambda x: - mu_0 * vsm.H_field * self.V_f * self.M_f * np.cos(x + vsm.phi)
        self.exchange = lambda x: - self.J_ex * self.S * np.cos(x - self.alpha)
        self.uniaxial = lambda x: self.K_f * self.V_f * np.sin(x - self.gamma_f)**2
        #self.biquadratic = lambda x: self.K_bq * self.V_f * np.sin(x - self.gamma_bq)**2 * np.cos(x - self.gamma_bq)**2
        #self.uniaxialAF = self.K_af * self.V_af * sp.sin(y - self.gamma_af)**2         2 variable function
        self.energy = lambda x: self.zeeman(x) + self.uniaxial(x) + self.exchange(x) #+ self.sym_biquadratic

    def define_gradEnergy_functions(self, vsm):
        self.gradZeeman = lambda x: mu_0 * vsm.H_field * self.V_f * self.M_f * np.sin(x + vsm.phi)
        self.gradExchange = lambda x: self.J_ex * self.S * np.sin(x - self.alpha)
        self.gradUniaxial = lambda x: self.K_f * self.V_f * 2 * np.sin(x - self.gamma_f) * np.cos(x - self.gamma_f)
        self.gradEnergy = lambda x: self.gradZeeman(x) + self.gradExchange(x) + self.gradUniaxial(x)

    def define_ggradEnergy_functions(self, vsm):
        self.ggradZeeman = lambda x: mu_0 * vsm.H_field * self.V_f * self.M_f * np.cos(x + vsm.phi)
        self.ggradExchange = lambda x: self.J_ex * self.S * np.cos(x - self.alpha)
        self.ggradUniaxial = lambda x: self.K_f * self.V_f * 2 * ( -np.sin(x - self.gamma_f)**2 + np.cos(x - self.gamma_f)**2 )
        self.ggradEnergy = lambda x: self.ggradZeeman(x) + self.ggradExchange(x) + self.ggradUniaxial(x)


    def define_sym_energy_functions(self, vsm):
        #print "Updating sample's energy..."
        # Sympy symbolic functions.
        self.sym_zeeman = - mu_0 * vsm.H_field * self.V_f * self.M_f * sp.cos(x + vsm.phi)
        self.sym_exchange = - self.J_ex * self.S * sp.cos(x - self.alpha)
        self.sym_uniaxial = self.K_f * self.V_f * sp.sin(x - self.gamma_f)**2
        self.sym_biquadratic = self.K_bq * self.V_f * sp.sin(x - self.gamma_bq)**2 * sp.cos(x - self.gamma_bq)**2
        #self.sym_uniaxialAF = self.K_af * self.V_af * sp.sin(y - self.gamma_af)**2         2 variable function
        self.sym_energy = self.sym_zeeman + self.sym_uniaxial #+ self.sym_exchange + self.sym_biquadratic

    def derivate_sym_energy_functions(self, vsm):
        # Sympy derivative function
        # First derivative
        self.sym_gradEnergy = self.sym_energy.diff(x)
        # Second derivative
        self.sym_ggradEnergy = self.sym_gradEnergy.diff(x)


    def convert_sym_energy_functions(self, vsm):
        # Conversion of sympy expressions to numpy lambda functions.
        self.energy = sp.lambdify(x, self.sym_energy, 'numpy')
        self.gradEnergy = sp.lambdify(x, self.sym_gradEnergy, 'numpy')
        self.ggradEnergy = sp.lambdify(x, self.sym_ggradEnergy, 'numpy')
