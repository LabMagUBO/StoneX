# -*- coding: utf-8 -*-
"""
    Fichier de script
"""
from StoneX.sample import *
from StoneX.functions import *
import logging

# Defining the logger
logger = logging.getLogger(__name__)


###########
#Stoner Wohlfart
###########

def stoner_test(sample, vsm):
    """
        Test et debug du modèle de stoner.
    """

    #### VSM ####
    vsm.set_angle(5)
    vsm.H_field = convert_field(5, 'si')
    vsm.H_field_max = convert_field(4, 'si')
    vsm.n_H_field = 2**8

    #### SAMPLE ####
    #sample.J_ex = 0
    sample.T = 0

    sample.apply(vsm)
    sample.relax()
    #sample.export_energy(vsm, 0)

    vsm.H_field = convert_field(-5, 'si')
    sample.apply(vsm)
    sample.relax()
    #sample.export_energy(vsm, 1)

    vsm.H_field = convert_field(0, 'si')
    sample.apply(vsm)
    print(sample.sym_energy)
    sample.relax()
    sample.export_energy(vsm, 2)

    vsm.H_field = convert_field(0.01, 'si')
    sample.apply(vsm)
    print(sample.sym_energy)
    sample.relax()
    sample.export_energy(vsm, 2)

    #vsm.do_cycle(sample, display=False, export=True, export_energy=False, display_energy=False, verbose=True)
    vsm.do_rotation(sample, phi_step = 5, phi_start = 0, phi_stop = 360, export=True, display=False, export_cycle=True, display_cycle=False)

def stoner0K_cycle(sample, vsm):

    #### VSM ####
    vsm.set_angle(0)
    vsm.H_field = convert_field(5, 'si')
    vsm.H_field_max = convert_field(40, 'si')
    vsm.n_H_field = 2**8

    #### SAMPLE ####
    #sample.J_ex = 0
    sample.T = 300

    ### Cycle ###
    vsm.do_cycle(sample, display=False, export=True, export_energy=False, display_energy=False, verbose=True)

    vsm.set_angle(30)
    vsm.do_cycle(sample, display=False, export=True, export_energy=False, display_energy=False, verbose=True)
    vsm.set_angle(60)
    vsm.do_cycle(sample, display=False, export=True, export_energy=False, display_energy=False, verbose=True)
    vsm.set_angle(90)
    vsm.do_cycle(sample, display=False, export=True, export_energy=False, display_energy=False, verbose=True)
    """
    vsm.set_angle(30)
    vsm.set_field(5)
    sample.T = 300
    sample.apply(vsm)

    sample.draw_energy()
    sample.relax()
    sample.draw_energy()
    sample.print_energy(vsm)
    sample.export_energy(vsm,0)

    """
    #sample.T = 4
    #vsm.do_cycle(sample, display=False, export=True, export_energy=True, display_energy=False, verbose=True)
    """
    for t in np.arange(0,10, 4):
        sample.T = t
        sample.apply(vsm)
        vsm.do_cycle(sample, display=False, export=True, export_energy=False, display_energy=False, verbose=True)
    #vsm.do_rotation(sample, phi_stop=10, export=True, display=True, export_cycle=True, display_cycle=False)
    """
    # Make a cycle
    #cycle, coerFields = vsm.do_cycle(sample, export_energy=True)
    #vsm.draw_cycle(cycle, sample)
    #print( convert_field(coerFields, 'cgs') )

def stoner0K_rotation(sample, vsm):

    #### VSM ####
    vsm.set_angle(0)
    vsm.H_field = convert_field(5, 'si')
    vsm.H_field_max = convert_field(14, 'si')
    vsm.n_H_field = 2**8

    #### SAMPLE ####
    sample.J_ex = 1e-6
    sample.T = 300

    ### Rotation ###
    vsm.do_rotation(sample, phi_stop=360, export=True, display=False, export_cycle=True, display_cycle=False)

###########
#Garcia-Otero
###########

def Garcia_Otero_test(sample, vsm):
    #### VSM ####
    vsm.set_angle(5)
    vsm.H_field = convert_field(20, 'si')
    vsm.H_field_max = convert_field(20, 'si')
    vsm.n_H_field = 2**7

    #### SAMPLE ####
    sample.J_ex = 0
    sample.T = 0
    sample.K_f = sample.K_f * 5 # Coercivité de 10 Oe à 0K

    ## Temperature slope
    for t in np.arange(0,310, 100):
        sample.T = t
        sample.apply(vsm)
        print(t)
        sample.print_energy(vsm)
        vsm.do_cycle(sample, display=False, export=True, export_energy=False, display_energy=False, verbose=True)

def Garcia_Otero_cycle(sample, vsm):
    #### VSM ####
    vsm.set_angle(0)
    vsm.H_field = convert_field(20, 'si')
    vsm.H_field_max = convert_field(20, 'si')
    vsm.n_H_field = 2**7

    #### SAMPLE ####
    #sample.J_ex = 0
    sample.T = 300
    sample.K_f = sample.K_f * 5 # Coercivité de 10 Oe à 0K
    sample.J_ex = sample.J_ex / 10      #réduction de l'échange

    ## Cycle
    vsm.do_cycle(sample, display=False, export=True, export_energy=False, display_energy=False, verbose=True)

def Garcia_Otero_rotation(sample, vsm):
    #### VSM ####
    vsm.set_angle(0)
    vsm.H_field = convert_field(20, 'si')
    vsm.H_field_max = convert_field(20, 'si')
    vsm.n_H_field = 2**7

    #### SAMPLE ####
    #sample.J_ex = 0
    sample.T = 300
    sample.K_f = sample.K_f * 5 # Coercivité de 10 Oe à 0K
    sample.J_ex = sample.J_ex / 10      #réduction de l'échange

    ### Rotation ###
    vsm.do_rotation(sample, phi_stop=360, export=True, display=False, export_cycle=True, display_cycle=False)


def Franco_Conde_test(sample, vsm):

    #### VSM ####
    vsm.set_angle(0)
    vsm.H_field = convert_field(20, 'si')
    vsm.H_field_max = convert_field(20, 'si')
    vsm.n_H_field = 2**7

    #### SAMPLE ####
    #sample.J_ex = 0
    sample.T = 300
    sample.K_f = sample.K_f * 5 # Coercivité de 10 Oe à 0K
    sample.J_ex = sample.J_ex / 10      #réduction de l'échange

    logger.info("Sample status : %s", sample)

    #sample.relax()
    #sample.display_energy(0)
    #sample.export_energy(vsm, 0)

    #cycle, H_coer = vsm.do_cycle(sample, display=False, export=True, export_energy=False, display_energy=False, verbose=True)

    #vsm.do_rotation(sample, phi_step = 5, phi_start = 0, phi_stop = 360, export=True, display=False, export_cycle=True, display_cycle=False)

def Franco_Conde_cycle(sample, vsm):
    #### VSM ####
    vsm.set_angle(0)
    vsm.H_field = convert_field(20, 'si')
    vsm.H_field_max = convert_field(20, 'si')
    vsm.n_H_field = 2**7

    logger.info("VSM status\n%s", vsm)

    #### SAMPLE ####
    #sample.J_ex = 0
    sample.T = 300
    sample.K_f = sample.K_f * 5 # Coercivité de 10 Oe à 0K
    sample.J_ex = sample.J_ex / 10      #réduction de l'échange

    logger.info("Sample status : %s", sample)

    vsm.do_cycle(sample, display=False, export=True, export_energy=False, display_energy=False, verbose=True)

def Franco_Conde_rotation(sample, vsm):
    #### VSM ####
    vsm.set_angle(0)
    vsm.H_field = convert_field(20, 'si')
    vsm.H_field_max = convert_field(20, 'si')
    vsm.n_H_field = 2**8

    logger.info("VSM status\n%s", vsm)

    #### SAMPLE ####
    #sample.J_ex = 0
    sample.T = 10
    sample.K_f = sample.K_f * 5 # Coercivité de 10 Oe à 0K
    #sample.J_ex = sample.J_ex / 10      #réduction de l'échange
    sample.J_ex = 0

    logger.info("Sample status : %s", sample)

    vsm.do_rotation(sample, phi_step = 5, phi_start = 0, phi_stop = 360, export=True, display=False, export_cycle=True, display_cycle=False)

def Franco_Conde_evolHc(sample, vsm):
    #### VSM ####
    vsm.set_angle(0)
    vsm.H_field = convert_field(15, 'si')
    vsm.H_field_max = convert_field(10, 'si')
    vsm.n_H_field = 2**8

    logger.info("VSM status\n%s", vsm)

    #### SAMPLE ####
    sample.T = 300
    sample.K_f = sample.K_f * 2  # Coercivité de 10 Oe à 0K
    #sample.J_ex = sample.J_ex / 10      #réduction de l'échange
    sample.J_ex = 0

    logger.info("Sample status: \n{0}".format(sample))

    ## Temperature slope
    T_slope = np.linspace(1, 600, 25)
    coer_fields = np.zeros((T_slope.size, 3))

    # Defining id lenght for the cycle files
    lenght = len(str(len(T_slope)))

    for i, t in enumerate(T_slope):
        print("Temperature = {0} K".format(t))
        sample.T = t
        sample.apply(vsm)

        # File id
        idx = str(i).zfill(lenght)

        cycle, H_coer = vsm.do_cycle(sample, display=False, export=True, export_energy=False, display_energy=False, verbose=True, idx=idx)
        coer_fields[i, 0], coer_fields[i, 1:3] = t, convert_field(H_coer, 'cgs')

        # Export the data
        np.savetxt("output/coer_fields.dat", coer_fields, header='T (K) \t Hc1(Oe) \t Hc2 (O)')

        #plotting
        pl.plot(coer_fields[:, 0], (coer_fields[:, 2] - coer_fields[:, 1])/2, 'ro', label='Hc')
        pl.plot(coer_fields[:, 0], (coer_fields[:, 2] + coer_fields[:, 1])/2, 'bo', label='He')
        pl.grid()
        pl.legend()
        pl.savefig("output/coer_fields.pdf")
        pl.close()
