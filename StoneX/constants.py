#!/opt/local/bin/ipython-2.7
# -*- coding: utf-8 -*-

import numpy as np
#from functions import *

################################################
# Constants and variables definitions
################################################
# Physical constants
mu_0 = 4*np.pi*1e-7 #H/m
k_B = 1.3806488 * 1e-23 #m2 kg s-2 K-1

## Néel relaxation parameters
f0 = 10**9          #attempt frequency, in Hz
tau_mes = 100       #measurement mean time, in s

## Extrinsic variables
H_field = 1e3/4/np.pi * 100                 # Magnetic field (A/m)
phi = np.radians(70)       # Angle (radians), between the field and the main sample's anisotropy
theta_eq = 0                  # Angle (radians), Equilibrium direction of the magnetization


## Intrinsic variables
# Ferromagnetic layer
V_f = 1e-12                     # volume de la couche F (m**3)
M_f = 400e-6/V_f * 1e3          # aimantation de F (A/m)
t_f = 10*1e-9                   # épaiseur de la couche ferro (m)
K_f = 5*1e3 * 1e-7              # anisotropie (J/m**3)

# Antiferromagnetic layer
V_af = 1e-4 * 100e-9        # volume de la couche AF
J_ex = 11*1e-3 * 1e-7 * 1e4        #Énergie surfacique de couplage d'échange

## Modelisation parameters
n_theta = 2**6                                               # nombre de point pour theta
vec_theta = np.linspace(0,2*np.pi,n_theta)

n_H_field = 2**9          # number of step in a H-cycle (make it 2**n)
H_field_max = 1e4/4/np.pi * 200        # maximal field value during a cycle


## Folders
#dos_cycles = "output/cycles"
#dos_data = "output/dat"
#dos_pdf = "output/pdf"
#dos_energy = "output/energy"
