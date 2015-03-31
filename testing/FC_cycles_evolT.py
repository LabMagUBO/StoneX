#!/opt/local/bin/ipython
# -*- coding: utf-8 -*-

"""
    Main program.
    StoneX V0.4 compatible
"""
# Path to the module
module_path = '/Users/zorg/These/Documents/Programmes/StoneX_Project/'

# Define module path
import sys
sys.path.append(module_path)

# Importing StoneX module
from StoneX import *

# Activate logging (cleaning the log on the first logger)
logger = init_log(__name__, console_level='debug', file_level='info', mode='w')
logger.info("Version du programme {}".format(StoneX.__version__))


##########################
#########  VSM  ##########
##########################
# Creating the VSM
vsm = VSM()

# Settings
#vsm.set_output('test')
vsm.set_angle(0)
vsm.H_field = convert_field(20, 'si')
vsm.H_field_max = convert_field(20, 'si')
vsm.n_H_field = 2**7

# Logginq it
logger.info("Etat du VSM : %s", vsm)


##########################
########  SAMPLE  ########
##########################
# Available models : Stoner_2states, Garcia_Otero, Franco_Conde
# Creating the sample
mySample = load_class(Domain, Franco_Conde)

# Settings
mySample.set_discrete(0.5)
mySample.T = 300
mySample.K_f = 200 # Coercivité de ~10 Oe = 1mT à 0K
#mySample.J_ex = sample.J_ex / 10      #réduction de l'échange
mySample.J_ex = 2e-6

logger.info("Echantillon : {0}".format(mySample))


##########################
#########  RUN  ##########
##########################
# Methods
def main(sample, vsm):
    # Create the logger
    logger = init_log(__name__)

    # Set output folder
    vsm.set_output('FC_cycles_evolT')

    #Temperatures
    Ts = np.arange(0, 1500, 200)
    Ts[0] = 1

    # Action
    Hc_evol = swap_T(sample, vsm, Ts)

    # enregistrement des données
    fileName = 'FC_cycles_evolT/Hc_evolT.dat'
    header = 'Temp (K) \t Hc1 (T) \t Hc2 (T)'
    np.savetxt(fileName, Hc_evol, header=header)



def swap_T(sample, vsm, T_tab):
    tab_out = np.zeros((T_tab.size, 3))

    for i, T in enumerate(T_tab):
        #Un peu de blabla
        logger.info('Changement temperature T={0}K'.format(T))

        #Enregistrement température
        sample.T = T
        sample.J_ex = 2e-6 * (600 - T) / 600

        #Cycle à température, n°i
        cycle, H_coer, mt_xtrm, ml_rem = vsm.do_cycle(sample, idx=i, display=False, export=True, export_energy=False, display_energy=False, verbose=True)

        #On enregistre
        tab_out[i, 0] = T
        tab_out[i, 1:3] = H_coer * mu_0


    # et on expulse tout ça après la boucle
    return tab_out



# Do something a least
if __name__ == '__main__':
    main(mySample, vsm)


## Stop the counter and log it
end_time(start_time)
