#!/opt/local/bin/ipython
# -*- coding: utf-8 -*-

"""
    Programme principal
"""


# Import module
from StoneX import *

# Activate logging
logger = logging.getLogger(__name__)

# Define output folder
define_output('output')

# Creating the VSM
vsm = VSM()
logger.info("Etat du VSM : %s", vsm)


# Setting the VSM parameters
#vsm.set_angle(0)
#vsm.H_field = convert_field(+5, 'si')
#vsm.H_field_max = convert_field(4, 'si')
#vsm.n_H_field = 2**7

#####
# Stoner Wohlfart
if False:
    sample = load_class(Domain, Stoner_2states)

    scripts.stoner_test(sample, vsm)
    scripts.stoner0K_cycle(sample, vsm)
    scripts.stoner0K_rotation(sample, vsm)


#####
# Garcia Otero
if False:
    sample = load_class(Domain, Garcia_Otero)
    #scripts.Garcia_Otero_test(sample, vsm)
    #scripts.Garcia_Otero_cycle(sample, vsm)
    scripts.Garcia_Otero_rotation(sample, vsm)

#### Franco-Conde
if True:
    fcsample = load_class(Domain, Franco_Conde)

    # Test Franco-Conde
    #scripts.francoTest(fcsample, vsm)

    # Franco-Conde Cycle
    scripts.Franco_Conde_cycle(fcsample, vsm)

    # Franco-Conde rotation
    #scripts.Franco_Conde_rotation(fcsample, vsm)

    #scripts.Franco_Conde_evolHc(fcsample, vsm)


## Stop the counter and log it
end_time(start_time)
