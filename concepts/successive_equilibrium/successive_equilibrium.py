#!/opt/local/bin/ipython-3.4
# -*- coding: utf-8 -*-

import sys
import numpy as np
import numpy.ma as ma
import scipy.ndimage as nd
from scipy import optimize
from matplotlib import pyplot as pl

# Constantes
gamma_f = 0.
gamma_af = 0.
H = -1.5
J = 1
Kaf = 2
T_bas = 0.4
T = 0.4

n = 100
# Discretisation de (alpha, theta)
nb = np.array([n, n])

# Création des variables
theta = np.linspace(0, 360, nb[1], endpoint=False, dtype='float64')
alpha = np.linspace(0, 360, nb[0], endpoint=False, dtype='float64')


# Tableau
Alph, Th = np.meshgrid(alpha, theta, sparse=False, indexing='ij')

# Fonction d'énergie
energy = lambda H, alph, th: - H * np.cos(2*np.pi / 360 * (th - gamma_f)) +\
        np.sin(2*np.pi / 360 * (th - gamma_f))**2 +\
        Kaf * np.sin(2*np.pi/360 * (alph - gamma_af))**2 - J * np.cos(2*np.pi/360 * (alph - th))

# Paysage énergétique
E = np.ma.masked_array(energy(H, Alph, Th))

# Fonction de tracé
def plot_landscape(E, *argv, name='default', contour=True):
    # Nom du fichier
    file = "{}.pdf".format(name)

    # Création de la figure
    fig = pl.figure()
    ax = fig.add_subplot(111, aspect='equal')
    ax.set_xlabel('theta')
    ax.set_ylabel('alpha')

    #ax.contourf(Theta, Alpha, E, cmap=pl.cm.gray, resample=False)
    cax = ax.imshow(E,
        interpolation = 'nearest',
        #origin='upper',
        extent=(0,360,360,0)
    )
    cbar = fig.colorbar(cax)

    if contour:
        C = ax.contour(E, 10,
            colors='black',
            linewidth=.5,
            extent=(0,360,0,360)
        )
        ax.clabel(C, inline=1, fontsize=10)

    pl.grid()

    if argv is not None:
        for k, pos in enumerate(argv):
            if pos.ndim == 1:
                ax.plot(pos[1], pos[0], 'ro')
            else:
                ax.plot(pos[:, 1] % 360, pos[:, 0] % 360, '-go')

    pl.savefig(file)
    pl.close()


# Redefining from search_eq in Stoner_Wohlfarth
def search_eq(E, eq):
    """
        Search the equilibrium state «eq» depending on the initial state «ini» and the energy landscape «E».
        eq is a two-value array [alpha_index, theta_index]
        Return the index of equilibrium
    """
    # stopping criteria
    found = False

    # Overflow security
    k = 0

    # E dimension
    nb = np.shape(E)

    # Path of equilibrium
    N_eq = np.zeros((1, 2))
    N_eq[0] = eq

    while not found:
        k += 1

        #Masking the array
        E = np.ma.array(E)
        E.mask = ma.masked

        # Unmasking all the equilibrium neighbours states
        #   1 1 1
        #   1 0 1
        #   1 1 1
        for i in np.arange(-1, 2) + N_eq[k-1][0]:
            for j in np.arange(-1, 2) + N_eq[k-1][1]:
                if i*j != 1:    #Not the center
                    E.mask[i % nb[0], j % nb[1]] = False

        # Append the minimal state to N_eq
        minIdx = E.argmin()
        N_eq = np.append(N_eq, [[np.floor(minIdx / nb[1]), minIdx % nb[1]]], axis=0)

        # Comparing energy
        barrier = E[N_eq[k][0], N_eq[k][1]] - E[N_eq[k-1][0], N_eq[k-1][1]]

        if barrier < 0:
            # Better state
            # Nothing to do, continue the loop
            continue
        elif barrier == 0:
            # Same energy state. Checking if we have already been there
            if k < 2:
                continue
            elif (N_eq[k] == N_eq[k-2]).all():
                # already been there
                # Stopping the loop
                break
            else:
                continue
        else:
            # no better states
            break

        # If too much loop
        if k == 10000:
            print("search_eq : Not converging.")
            break

    # No forgetting to unmask the array
    E.mask = ma.nomask
    return N_eq


def label_energy(E, eqIdx):
    # No forgetting to unmask the array
    E.mask = False

    # Masking all the values above energy state equilibrium
    mask = E > E[eqIdx[0], eqIdx[1]] + T

    # Labeling the mask (0-> unreachable states ; 1,2,3, i..->ith zone)
    mask_lab, mask_num = nd.measurements.label(np.logical_not(mask))

    #plot_landscape(mask_lab, name='H{}_mask_lab_k{}'.format(H, k), contour=False)

    ## Searching equivalent zones (the array is periodic)
    # First scanning alpha boundaries
    zones_equiv = set()     #set() to contain «unique» equivalent zone label;
    for i in np.arange(nb[0]):
        left = mask_lab[i, 0]
        right = mask_lab[i, -1]
        if left and right and (left != right): #left and right not zero, and different labels
            zones_equiv.add(tuple((min(left, right), max(left, right))))

    # Converting the set to array
    zones_equiv = np.array(list(zones_equiv))

    # Relabeling the upper label to the equivalent lower label
    for i, val in enumerate(zones_equiv):
        # Creating a mask with only the upper label
        mask_rm = mask_lab == val[1]
        #to replace all the zone by le lower label
        mask_lab[mask_rm] = val[0]

    # Then, scanning theta boundaries
    zones_equiv = set()     #set() to contain «unique» equivalent zone label;
    for i in np.arange(nb[1]):
        top = mask_lab[0, i]
        bottom = mask_lab[-1, i]
        if top and bottom and (top != bottom): #left and right not zero, and different labels
            zones_equiv.add( tuple(    (min(top, bottom), max(top, bottom))       )    )

    # Converting the set to array
    zones_equiv = np.array(list(zones_equiv))

    # Relabeling the upper label to the equivalent lower label
    for i, val in enumerate(zones_equiv):
        # Creating a mask with only the upper label
        mask_rm = mask_lab == val[1]
        #to replace all the zone by le lower label
        mask_lab[mask_rm] = val[0]


    # Renaming the last labels (not useful)
    #labels = np.unique(z_lab)
    #z_lab = np.searchsorted(labels, z_lab)

    return mask_lab


def masking_energy(E, eqIdx):
        """
            Method for masking all non-reachable states above the equilibrium.
            After searching all the available states, return the new equilibrium position.
            The non-reachable states are masqued in self.E
        """
        ok = False      # Test result

        # Shape of E[phiIdx, HIdx]
        #nb = np.array([self.alpha.size, self.theta.size])

        # Index of while loop
        k = 0
        # Loop
        while not ok:
            k += 1
            print("loop", k)

            mask_lab = label_energy(E, eqIdx)
            plot_landscape(mask_lab, name='H{}_mask_lab_corr_k{}'.format(H, k), contour=False)

            #Zone where the equilibrium is
            zone_eq = mask_lab[eqIdx[0], eqIdx[1]]
            # Creating a new mask
            new_mask = mask_lab != zone_eq

            #Applying the mask
            E[new_mask] = np.ma.masked

            plot_landscape(E, name='H{}_Landscape_accessible_k{}'.format(H, k))


            # Searching the minimum
            mini = np.ma.min(E)
            # If the minimum over the zone is the same that the eq energy
            if mini == E[eqIdx[0], eqIdx[1]]:
                # finish the loop
                ok = True
            else:
                # Trying to jump over the saddle point
                eqIdx = search_mountainPass(E, E_old, eqIdx)
                # Search stable equilibrium
                eq_trans = search_eq(E, eqIdx) / nb * 360


                plot_landscape(E, eq_trans[0], eq_trans[1:], eq_trans[-1], name='H{}_after_jump'.format(H))

                eqIdx = eq_trans[-1] / 360 * nb

        #return eqIdx
        return eqIdx


def search_mountainPass(E, E_previous, eqIdx):
    """
    Arguments :
        E : accessible energy landscape
        E_old : old accessible energy landscape, for comparison

    Calculated :
        r : searching area radius
    """

    ## Searching the mean radius of the previous reachable domain
    # Copy of the old mask to work with
    E_old = np.ma.copy(E_previous)

    mask_lab = label_energy(E_old, eqIdx)

    # Index where the actual equilibrium is
    domIdx_eq = mask_lab[eqIdx[0], eqIdx[1]]

    # The old mask
    mask_old = mask_lab == domIdx_eq

    #minimum Average radius (case if the domain have a longer axis)
    r2 = np.round(np.sqrt(mask_old.sum()/np.pi), 0)
    r0 = np.max(mask_old.sum(axis=0))/2
    r1 = np.max(mask_old.sum(axis=1))/2
    r = np.min([r0, r1, r2])
    print("rs", r0, r1, r2)
    print("r", r)

    # Successive try positions
    pos = np.array([eqIdx])

    # Trying loop (infinite, with a stop at 10000 tries)
    k = 0   # Loop number
    i = 0   # position number
    while True:
        candidate = pos[i, :] + np.random.randint(-r, r, size=2)
        candidate = candidate % nb

        if E[candidate[0], candidate[1]] is not ma.masked:
            i +=1   # incrementing
            pos = np.append(pos, [candidate], axis=0)   # adding the new position

            # Testing if the try have changed zone
            try_label = mask_lab[candidate[0], candidate[1]]

            if try_label != domIdx_eq and try_label != 0:
                print("find new zone", i)
                plot_landscape(E, pos/n*360, name="Tries")
                break

        k += 1
        if k == 10000:
            print("Not converging")
            plot_landscape(E, pos/n*360, name="Tries_notconv")

            break
    return pos[-1]


## Principal
# Position initiale
eq_ini = np.array([40, 40])

# Indices initiaux
N_eq_ini = eq_ini / 360 * nb

# Recherche du point d'équilibre proche
eq_trans = search_eq(E, N_eq_ini) / nb * 360
eq_fin = eq_trans[-1]

# Tracé
plot_landscape(E, eq_ini, eq_trans[1:], eq_trans[-1], name='H{}_Landscape_ini'.format(H))


# États accessibles selon la température et l'état d'équilibre
Eq_fin = masking_energy(E, eq_fin / 360 * nb)
plot_landscape(E, Eq_fin, name='H{}_Lanscape_final'.format(H))

E_old = E

# Changement du champ
H -= 0.1
# Paysage énergétique
E = np.ma.masked_array(energy(H, Alph, Th))

# États accessibles selon la température et l'état d'équilibre
Eq_fin = masking_energy(E, eq_fin / 360 * nb) *360 / nb
plot_landscape(E, Eq_fin, name='H{}_Lanscape_final'.format(H))
