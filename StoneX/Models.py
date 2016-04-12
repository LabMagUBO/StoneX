#!/opt/local/bin/env python
# -*- coding: utf-8 -*-
"""
    Models Class.

    Attributes :


    Properties :

    Methods :

    Copyright (C) 2016  Jérôme Richy
"""

# General modules
import sys
import numpy as np
# import sympy as sp
import scipy.ndimage as nd
# from sympy.abc import x
from matplotlib import pylab as pl

# My own modules
# from constants import *
from StoneX.Physics import *
from StoneX.Logging import *
from StoneX.Sample import *
# from StoneX.models import *
# from StoneX.constants import *


class Stoner_Wohlfarth(Ferro):
    """
        Stoner_Wohlfarth class.
        Depends on Ferro.
    """
    def __init__(self):
        super().__init__()

        self.logger.debug("Init. Model Class Stoner_Wohlfarth")
        # Naming the model
        self.model = "Stoner_Wohlfarth"

    def calculate_energy(self, vsm):
        """
            Calculate the energy for each parameter's values given by the vsm.
            The parameters are, in order : phi, H, theta
        """

        # sparse=True for saving memory
        # indexing : place the element in the index order
        Hf, Th = np.meshgrid(
            vsm.H, self.theta,
            sparse=True, indexing='ij'
        )

        # Creating a masked array. All values unmasked.
        self.E = np.ma.array(self.energy()(vsm.angle[1], Hf, Th))

    def search_eq(self, E, eq):
        """
            Search the equilibrium state «eq» depending on the initial state
             «ini» and the energy landscape «E».
            Return the index of equilibrium
        """
        # stopping criteria
        found = False

        # Overflow security
        k = 0

        while not found:
            k += 1

            # print("calc barrier")
            barrier = np.array([
                E[(eq - 1) % E.size] - E[eq],
                E[(eq + 1) % E.size] - E[eq]
            ])

            if (barrier >= 0).all():
                # local minimum
                found = True
            elif barrier[0] == barrier[1]:
                self.logger.warn("Symmetrical unstable equilibrium.\
                                Choosing next largeur index (theta[idx+1]).")
                eq = (eq + 1) % E.size
            else:
                idx = np.argmin(barrier)
                eq = (eq + idx * 2 - 1) % E.size

            if k == 10000:
                self.logger.error("search_eq : Not converging.")
                pl.plot(E, '-ro')
                pl.show()
                break

        return eq

    def calculate_magnetization(self, phi, HIdx, eqIdx):
        """
            Calculate the magnetization, only using the local minimum state.
            Will be redefined in Franco_Conde to integrate all the available
            states.
            Arguments : self, phi, HIdx, eqIdx. (to be compatible with Franco-
            Conde)
            Return Mt, Ml.
        """
        return self.Ms * self.V_f * np.sin(self.theta[eqIdx] + phi), \
            self.Ms * self.V_f * np.cos(self.theta[eqIdx] + phi)

    def set_memory(self, vsm):
        """
            Defines the tables for storing calculated data.
        """
        # Creating an array for temperature change.
        # Only one possible value for SW model: T=0K
        self.rotations = np.zeros(1, dtype='object')

        # Data stored in Rotation object.
        self.rotations[0] = Rotation(vsm.phi)
        self.rotations[0].info(self)

        # Loop over phi
        for i, phi in enumerate(vsm.phi):
            # Creating a cycle
            cycle = Cycle(vsm.H, 4)
            cycle.info(self, phi)
            self.rotations[0].cycles[i] = cycle

    def analyse_energy(self, vsm):
        """
            After calculating the energy depending on the parameters, search in
            which state is the magnetization.
            Returns the magnetization path : tab[phi, [H, Mt, Ml, theta_eq]]
        """
        # Direction of the vsm
        phi = vsm.angle[1]
        angle_index = vsm.angle[0]

        # Index of the magnetization equilibrium, to start.
        # Correspond to the global energy minimum
        eq = np.argmin(self.E[0, :])

        # Creating a cycle (H, Mt, Ml, theta_eq)
        cycle = Cycle(vsm.H, 4)
        cycle.info(self, phi)
        # Loop over H (max -> min -> max)
        for j, H in enumerate(vsm.H):
            # self.logger.debug("j = {}, H = {}Oe".format(
            #    j, convert_field(H, 'cgs'))
            # )

            # Finding equilibrium
            eq = self.search_eq(self.E[j, :], eq)

            # Calculating the magnetization
            Mt, Ml = self.calculate_magnetization(phi, j, eq)

            # Storing the results
            cycle.data[j] = np.array([
                H, Mt, Ml, self.theta[eq]
            ])

            # Save energy data, to be plotted later
            cycle.energy[j] = np.ma.copy(self.E[j, :])

        # Saving the cycle
        self.rotations[0].cycles[angle_index] = cycle


# Initial : Stoner_Wohlfarth, AntiFerro, Ferro
class Meiklejohn_Bean(AntiFerro, Stoner_Wohlfarth):
    def __init__(self):
        super().__init__()

        self.logger.debug("Init. Model Class Meiklejohn_Bean")
        self.model = "Meiklejohn_Bean"

        self.S = (200e-9)**2
        self.J_ex = 11e-3 * 1e-7 * 1e4       # J/m**2
        self.J_bq = 0

    # Redefining only the energy
    def energy(self):
        exchange = lambda th: - self.J_ex * self.S * np.cos(th - self.gamma_af)
        biquadratic = lambda th: \
            - self.J_bq * self.S * np.cos(th - self.gamma_af)**2
        return lambda phi, H, th: \
            Ferro.energy(self)(phi, H, th) + exchange(th) + biquadratic(th)


class Garcia_Otero(Meiklejohn_Bean):
    """
        Garcia Otero model.
    """
    def __init__(self):
        # Initiating subclasses
        super().__init__()

        self.logger.debug("Init. Model Class Garcia_Otero")
        # Naming the model
        self.model = "Garcia_Otero"

        # Additionnal parameters
        # Temperature
        self.T = 300            # K

        # Néel relaxation parameters
        self.f0 = 1e9           # attempt frequency, in Hz
        self.tau_mes = 100      # measurement mean time, in s

    # No need to define energy, same as Meiklejohn_Bean

    def masking_energy(self, HIdx, eqIdx):
        """
            Method for masking all non-reachable states above the equilibrium.
            After searching all the available states, return the new
            equilibrium position.
            The non-reachable states are masqued in self.E
        """
        # Loop over local minimum
        while True:
            # No forgetting to unmask the array
            self.E[HIdx, :].mask = False

            # Masking all the values above energy state equilibrium
            mask = self.E[HIdx, :] > self.E[HIdx, eqIdx] \
                + np.log(self.f0 * self.tau_mes) * k_B * self.T

            mask_lab, mask_num = nd.measurements.label(np.logical_not(mask))

            left = mask_lab[0]
            right = mask_lab[-1]
            # if the two extrem are reachable
            if left and right and (left != right):
                # replacement mask
                mask_replace = mask_lab == left
                mask_lab[mask_replace] = right

                # relabelling the mask (not useful, time consuming)
                # labels = np.unique(mask_lab)
                # mask_lab = np.searchsorted(labels, mask_lab)

            # Zone where the equilibrium is
            zone_eq = mask_lab[eqIdx]
            # Creating a new mask
            new_mask = mask_lab != zone_eq

            # Applying the mask
            self.E[HIdx, new_mask] = np.ma.masked

            # Searching the minimum
            mini = np.ma.min(self.E[HIdx, :])

            # If the minimum over the zone is the same that the eq energy
            if mini == self.E[HIdx, eqIdx]:
                # finish the loop
                break
            else:
                # new equilibrium
                self.logger.warn("old eq {}".format(eqIdx))
                eqIdx = self.E[HIdx, :].argmin()
                self.logger.warn("new eq {}".format(eqIdx))

        return eqIdx

    # No need to redefine calculate_magnetization

    # Refining set_memory to store data over T
    def set_memory(self, vsm):
        """
            Defines the tables for storing calculated data.
        """
        # Creating an array for temperature change.
        self.rotations = np.zeros(vsm.T.size, dtype='object')

        # Loop over vsm.T
        for k, T in enumerate(vsm.T):
            # Changing the sample's temperature
            self.T = T

            # Data stored in Rotation object.
            # Cycle containing all the data. Indexes : ([H, Mt, Ml, theta_eq])
            self.rotations[k] = Rotation(vsm.phi)
            self.rotations[k].info(self)

            # Loop over phi
            for i, phi in enumerate(vsm.phi):
                # Creating a cycle
                cycle = Cycle(vsm.H, 4)
                cycle.info(self, phi)
                self.rotations[k].cycles[i] = cycle

    # Redefining analyse_energy from Stoner_Wohlfarth
    def analyse_energy(self, vsm):
        """
            After calculating the energy depending on the parameters, search in
            which state is the magnetization.
            Returns the magnetization path : tab[phi, [H, theta, Mt, Ml]]
        """
        # Direction of the vsm
        phi = vsm.angle[1]
        angle_index = vsm.angle[0]

        # Loop over vsm.T
        for k, T in enumerate(vsm.T):
            # Changing the sample's temperature
            self.T = T
            self.logger.info("Changing temperature : T = {} K".format(T))

            # Index of the magnetization equilibrium, to start.
            # Correspond to the global energy minimum
            eq = np.argmin(self.E[0, :])

            # Creating a cycle
            cycle = Cycle(vsm.H, 4)
            cycle.info(self, phi)

            self.logger.info("Loop over H.")
            # Loop over H (max -> min -> max)
            for j, H in enumerate(vsm.H):
                # Unmasking the energy
                self.E[j, :].mask = False

                # Adjusting equilibrium
                eq = self.search_eq(self.E[j, :], eq)

                # Defining all the accessible states, depending on the
                # temperature (changing equilibrium if necessary)
                eq = self.masking_energy(j, eq)

                # Calculating the magnetization
                Mt, Ml = self.calculate_magnetization(phi, j, eq)

                # Storing the results (H, Mt, Ml, theta)
                cycle.data[j] = np.array([H, Mt, Ml, self.theta[eq]])

                # Storing the energy landscape
                if vsm.plot_energyLandscape or vsm.plot_energyPath:
                    cycle.energy[j] = np.ma.copy(self.E[j, :])

            # Saving cycle
            self.rotations[k].cycles[angle_index] = cycle


class Franco_Conde(Garcia_Otero):
    """
        Franco-Conde model, based on Garcia-Otero.
        Only the calcul of the magnetization differs.
    """
    def __init__(self):
        # Initiating subclasses
        super().__init__()

        self.logger.debug("Init. Model Class Franco_Conde")
        # Naming the model
        self.model = "Franco_Conde"

    # Same energy function as Garcia_Otero

    # Redifining only calculate_magnetization
    def calculate_magnetization(self, phi, HIdx, eqIdx):
        """
            Calculate the magnetization, by integrating all the available
            states.
            Arguments : self, phi, HIdx, eqIdx.
            Return Mt, Ml.
        """
        # View of the energy array
        E = self.E[HIdx, :]

        # Probability array
        P = np.exp(- E / k_B / self.T)

        # Partition function
        Z = np.ma.sum(P)

        # Testing if the sum is inf OR Z==0
        if np.isinf(Z) or Z * self.Ms * self.V_f == 0:
            # Only one state is accessible
            P[:] = 0
            P[eqIdx] = 1
            Z = 1

        # Calculating magnetization
        Ml = self.Ms * self.V_f * np.ma.sum(P * np.cos(self.theta + phi)) \
            / np.ma.sum(P)
        Mt = self.Ms * self.V_f * np.ma.sum(P * np.sin(self.theta + phi)) \
            / np.ma.sum(P)

        return Mt, Ml


class Rotatable_AF(Franco_Conde, AntiFerro_Rotatable):
    def __init__(self):
        # Initiating subclasses
        super().__init__()

        self.logger.debug("Init. Model Class Rotatable_AF")
        # Naming the model
        self.model = "Rotatable_AF"

    # Redefining the energy function
    def energy(self):
        exchange = lambda alph, th: - self.J_ex * self.S * np.cos(th - alph)
        biquadratic = lambda alph, th: \
            - self.J_bq * self.S * np.cos(th - alph)**2

        return lambda phi, H, alph, th: Ferro.energy(self)(phi, H, th) \
            + AntiFerro_Rotatable.energy(self)(alph) \
            + exchange(alph, th) \
            + biquadratic(alph, th)

    # Redefining calculate_energy from Stoner_Wolhfarth, adding the alpha
    # degree of freedom.
    def calculate_energy(self, vsm):
        """
            Calculate the energy for each parameters' values given by the vsm.
            The parameters are, in order : phi, H, alpha, theta
        """

        # sparse=True for saving memory
        # indexing : place the element in the index order
        Hf, Al, Th = np.meshgrid(
            vsm.H, self.alpha, self.theta,
            sparse=True, indexing='ij'
        )

        # Creating a masked array. All values unmasked.
        self.E = np.ma.array(self.energy()(vsm.angle[1], Hf, Al, Th))

    # Redefining from search_eq in Stoner_Wohlfarth
    def search_eq(self, E, eq):
        """
            Search the equilibrium state «eq» depending on the initial state «ini» and the energy landscape «E».
            eq is a two-value array [alpha_index, theta_index]
            Return the index of equilibrium
        """
        # Stopping criteria
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
            E.mask = True

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
                self.logger.error("search_eq : Not converging.")
                break

        # No forgetting to unmask the array
        E.mask = False
        return N_eq[-1], N_eq

    # Debug plot
    def debugplot_landscape(self, E, *argv, name='default', contour=True):
        """
            For debug : plot the energy landscape with successive position
        """
        # Nom du fichier
        file = "{}.pdf".format(name)

        # Création de la figure
        fig = pl.figure()
        ax = fig.add_subplot(111, aspect='equal')
        ax.set_xlabel('theta')
        ax.set_ylabel('alpha')

        cax = ax.imshow(E,
            interpolation = 'nearest',
            origin='upper',
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
                    ax.plot(pos[1]/E.shape[1]*360, pos[0]/E.shape[0]*360, 'ro')
                else:
                    ax.plot(pos[:, 1]/E.shape[1]*360, pos[:, 0]/E.shape[0]*360, '-go')

        pl.savefig(file)
        pl.close()

    def labelling_mask(self, E, eqIdx):
        """
            Method for labelling the reachable zones, taking into account the periodicity of the landscape.
            Returns a labelled mask.
        """
        # Shape of E[phiIdx, HIdx]
        nb = np.array([self.alpha.size, self.theta.size])

        # No forgetting to unmask the array
        E.mask = False

        # Masking all the values above energy state equilibrium
        mask = E > E[eqIdx[0], eqIdx[1]] + np.log(self.f0 * self.tau_mes) * k_B * self.T

        # Labeling the mask (0-> unreachable states ; 1,2,3, i..->ith zone)
        mask_lab, mask_num = nd.measurements.label(np.logical_not(mask))

        if False:
            # Creating figure, with title
            fig = pl.figure()
            fig.set_size_inches(18.5,10.5)
            #fig.suptitle("Model : {}, T = {}K, phi= {}deg, H = {} Oe".format(self.model, self.T, np.degrees(self.phi), np.round(convert_field(H, 'cgs'), 2)))

            # Axis
            ax = fig.add_subplot(111, aspect='equal')

            cax1 = ax.imshow(mask_lab, label="Energy landscape", interpolation = 'nearest', origin='upper', alpha=0.6, extent=(0, 360, 360, 0))
            cbar1 = fig.colorbar(cax1)

            ax.plot(np.degrees(self.theta[eqIdx[1]]), np.degrees(self.alpha[eqIdx[0]]), 'bo', label='old Eq.')

            #ax.set_title("id:{}".format(HIdx))
            ax.set_xlabel("Theta M_f(deg)")
            ax.set_ylabel("Alpha M_af(deg)")

            # Saving the figure
            pl.savefig("debugPlot/test_maskLab_T{}_{}_{}.pdf".format(self.T, eqIdx[0], eqIdx[1] ))
            pl.close()


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

    def search_mountainCol(self, E, mask_lab, eqIdx):
        """
        Arguments :
            E : accessible energy landscape
            E_old : old accessible energy landscape, for comparison

        Calculated :
            r : searching area radius

            Returns all the try indexes
        """
        # Shape of E[phiIdx, HIdx]
        nb = np.array([self.alpha.size, self.theta.size])

        # Initial energy level
        E_ini = E[eqIdx[0], eqIdx[1]]

        if False:
            dos = 'debugPlot/'
            self.debugplot_landscape(mask_lab, name='{}oldmask_eqIdx{}_{}'.format(dos, eqIdx[0], eqIdx[1]))

        # Index where the actual equilibrium is
        zonEq_label = mask_lab[eqIdx[0], eqIdx[1]]

        # The old mask, only containing the reachable domain
        mask_old = mask_lab == zonEq_label

        # Calculate the ration between theta and alpha, if theta and alpha array does not have the same size
        ratio = nb[1] / nb[0]

        # Minimum Average radius (case if the domain have a longer axis) (for theta lenght : self.alpha can be smaller)
        r2 = np.round(np.sqrt(mask_old.sum()/np.pi * ratio), 0)     # circular domain
        r0 = np.max(mask_old.sum(axis=0))/2 * ratio                 # alpha stretch domain
        r1 = np.max(mask_old.sum(axis=1))/2                         # theta stretch domain
        r = np.ceil(np.min([r0, r1, r2]))
        #self.logger.warn("{} {} {} {}".format(r1, r2, r0, r))

        # Successive try index
        tryIdx = np.array([eqIdx])

        # Trying loop (infinite, with a stop at 10000 tries)
        k = 0   # Loop number
        i = 0   # position number
        while True:
            # r and rAlpha are the maximum step for theta and alpha respectively
            rAlpha = np.ceil(r / ratio)
            # Random step for theta and alpha
            thetaStep = np.random.randint(-r, r+1)
            alphaStep = np.random.randint(-rAlpha, rAlpha+1)
            #print("Step", r, rAlpha, thetaStep, alphaStep)

            # New candidate
            candidate = tryIdx[i, :] + np.array([alphaStep, thetaStep])
            candidate = candidate % nb
            #print(candidate, nb)

            if E[candidate[0], candidate[1]] is not np.ma.masked:
                i +=1   # incrementing candidate index
                tryIdx = np.append(tryIdx, [candidate], axis=0)   # adding the new position

                # Testing if the try have changed zone
                try_label = mask_lab[candidate[0], candidate[1]]

                if try_label != zonEq_label and try_label != 0:
                    # Reached a new zone
                    break
                elif E[candidate[0], candidate[1]] < E_ini:
                    # Find a lower energy state
                    break

            k += 1     # incrementing loop index
            if k % 10000 == 0:
                self.logger.warn( "Nb loop : {}".format(k))

                if False:
                    print("plotting")
                    dos = 'debugPlot/'
                    self.debugplot_landscape(mask_old, tryIdx[-100:], name='{}oldmask_eqIdx{}_{}_{}'.format(dos, eqIdx[0], eqIdx[1], k))

                    #sys.exit(0)
                break

            if k == 300000:
                self.logger.warning("Random search not converging. Activate and check the debug map.")

                break

        return tryIdx[-1], tryIdx

    # Redefining from Garcia_Otero
    def masking_energy(self, HIdx, eqIdx):
        """
            Method for masking all non-reachable states above the equilibrium.
            After searching all the available states, return the new
            equilibrium position.
            The non-reachable states are masqued in self.E
        """

        # Shape of E[phiIdx, HIdx]
        nb = np.array([self.alpha.size, self.theta.size])

        # Index of while loop
        k = 0

        # Loop
        while True:
            # Incrementing the counter
            k += 1

            # Labelling the mask from the energy landscape
            mask_lab = self.labelling_mask(self.E[HIdx, :, :], eqIdx)

            # Zone where the equilibrium is
            zone_eq = mask_lab[eqIdx[0], eqIdx[1]]

            # Creating a new mask, corresponding to the reachable zone
            mask_eq = mask_lab == zone_eq

            # Applying the mask, masking all the unreachable states
            self.E[HIdx, ~mask_eq] = np.ma.masked

            # Searching the minimum over the zone
            mini = np.ma.min(self.E[HIdx, :, :])

            # If the minimum over the zone is the same that the eq energy
            if mini == self.E[HIdx, eqIdx[0], eqIdx[1]]:
                # Finish the loop
                break

            else:
                # Getting the previous labelled mask
                oldMask_lab = self.labelling_mask(
                    np.ma.copy(self.E[HIdx - 1, :, :]), eqIdx
                )

                # Calculating the global minimun index
                minIdx = np.ma.argmin(self.E[HIdx, :, :])
                possibleIdx = np.array([int(minIdx / nb[1]), minIdx % nb[1]])

                # If only one zone including the two minima, ie old eq and
                # new eq are in the same old zone
                # Landscape can have infinite zones, but only
                # two are relevant
                if oldMask_lab[possibleIdx[0], possibleIdx[1]] == \
                oldMask_lab[eqIdx[0], eqIdx[1]]:
                    # Calculate the global minimum
                    eqIdx = possibleIdx
                    break

                else:
                    # Going through the landscape col

                    # Debug plot
                    debugPlot = False
                    if debugPlot:
                        dos = 'debugPlot/'
                        self.debugplot_landscape(
                            oldMask_lab,
                            eqIdx,
                            name='{}T{}_H{}_k{}_1oldMasklab'.format(
                                dos, self.T, HIdx, k
                            )
                        )

                        self.debugplot_landscape(
                            mask_lab,
                            name='{}T{}_H{}_k{}_2Masklab'.format(
                                dos, self.T, HIdx, k
                            )
                        )

                        self.debugplot_landscape(
                            self.E[HIdx, :, :],
                            eqIdx,
                            name='{}T{}_H{}_k{}_3Elandscape'.format(
                                dos, self.T, HIdx, k
                            )
                        )

                    # Trying to jump over the saddle point, giving a new
                    # index (unstable)
                    eqIdx, tries = self.search_mountainCol(
                        self.E[HIdx, :, :], oldMask_lab, eqIdx
                    )

                    # Debug plot
                    if debugPlot:
                        self.debugplot_landscape(
                            self.E[HIdx, :, :], tries, tries[-1],
                            name='{}T{}_H{}_k{}_4search_mountainPass'.format(
                                dos, self.T, HIdx, k
                            )
                        )

                    # Search stable equilibrium
                    eqIdx, pathIdx = self.search_eq(self.E[HIdx, :, :], eqIdx)

                    # Debug plot
                    if debugPlot:
                        self.debugplot_landscape(
                            self.E[HIdx, :, :], pathIdx, eqIdx, name='{}T{}_H{}_k{}_5search_eq'.format(
                                dos, self.T, HIdx, k
                            )
                        )

            # If too much loop
            if k > 1000:
                self.logger.error("Too much loop in masking energy. \
                    Unable to find the local minimum.")
                break

        return eqIdx

    # Redefining calculate_magnetization from Franco_Conde
    def calculate_magnetization(self, phi, HIdx, eqIdx):
        """
            Calculate the magnetization, by integrating all the available states.
            Arguments : self, phi, phiIdx, HIdx, eqIdx.
            Return Mt, Ml.
        """
        # Energy array
        E = self.E[HIdx, :, :]

        # Probability array
        P = np.exp(- E / k_B / self.T)

        # Partition function
        Z = np.ma.sum(P)

        # Testing if the sum is inf OR Z * Ms V_f ==0 (ie < 1e-308, depends
        # on the cpu)
        if np.isinf(Z) or Z * self.Ms * self.V_f == 0:
            # Only one state is accessible
            P[:] = 0
            P[eqIdx[0], eqIdx[1]] = 1
            Z = 1

        # Calculating magnetization
        Ml = self.Ms * self.V_f * np.ma.sum(P * np.cos(self.theta + phi)) \
            / np.ma.sum(P)
        Mt = self.Ms * self.V_f * np.ma.sum(P * np.sin(self.theta + phi)) \
            / np.ma.sum(P)

        return Mt, Ml

    # No need to redefine set_memory

    # Redefining analyse_energy from Garcia_Otero
    def analyse_energy(self, vsm):
        """
            After calculating the energy depending on the parameters, search in
            which state is the magnetization.
            Returns the magnetization path :
                tab[phi, [H, Mt, Ml, theta, alpha]]
        """
        # Direction of the vsm
        phi = vsm.angle[1]
        angle_index = vsm.angle[0]

        # Loop over vsm.T
        for k, T in enumerate(vsm.T):
            # Changing the sample's temperature
            self.T = T
            self.logger.info("Changing temperature : T = {} K".format(T))

            # Index of the magnetization equilibrium, to start. Correspond to
            # the global energy minimum
            eq = np.argmin(self.E[0, :, :])
            thIdx = eq % self.theta.size
            alphIdx = np.floor(eq / self.theta.size)
            eq = np.array([alphIdx, thIdx])

            # Creating a cycle
            cycle = Cycle(vsm.H, 5)
            cycle.info(self, phi)

            self.logger.info("Loop over H.")
            # Loop over H (max -> min -> max)
            for j, H in enumerate(vsm.H):
                # Unmasking the energy
                self.E[j, :, :].mask = False

                # Adjusting equilibrium
                eq, pathIdx = self.search_eq(self.E[j, :, :], eq)

                # Defining all the accessible states, depending on the
                # temperature (changing equilibrium if necessary)
                eq = self.masking_energy(j, eq)

                # Calculating the magnetization
                Mt, Ml = self.calculate_magnetization(phi, j, eq)

                # Storing the results (H, Mt, Ml, theta, alpha)
                cycle.data[j] = np.array(
                    [H, Mt, Ml, self.theta[eq[1]], self.alpha[eq[0]]]
                )

                # Storing the energy landscape
                if vsm.plot_energyLandscape or vsm.plot_energyPath:
                    cycle.energy[j] = np.ma.copy(self.E[j, :, :])

            # Saving the cycle
            self.rotations[k].cycles[angle_index] = cycle


class Double_MacroSpin(AntiFerro_Spin, Rotatable_AF):
    def __init__(self):
        # Initiating subclasses
        super().__init__()

        self.logger.debug("Init. Model Class Double_MacroSpin")
        # Naming the model
        self.model = "Double_MacroSpin"

    # Redefining the energy
    def energy(self):
        exchange = lambda alph, th: - self.J_ex * self.S * np.cos(th - alph)
        biquadratic = lambda alph, th: \
            - self.J_bq * self.S * np.cos(th - alph)**2
        return lambda phi, H, alph, th: \
            Ferro.energy(self)(phi, H, th) \
            + AntiFerro_Spin.energy(self)(phi, H, alph) \
            + exchange(alph, th) \
            + biquadratic(alph, th)

    # and the magnetization
    def calculate_magnetization(self, phi, phiIdx, HIdx, eqIdx):
        """
            Calculate the magnetization, by integrating all the available states.
            Arguments : self, phi, phiIdx, HIdx, eqIdx.
            Return Mt, Ml.
        """
        # Column view of alpha
        col_alpha = self.alpha.reshape(self.alpha.size, 1)
        # Energy array
        E = self.E[phiIdx, HIdx]

        # Probability array
        P = np.exp(- E / k_B / self.T)

        # Partition function
        Z = np.ma.sum(P)

        # Testing if the sum is inf OR Z * Ms V_f ==0 (ie < 1e-308, depends on the cpu)
        if np.isinf(Z) or Z * self.Ms * self.V_f == 0:
            # Only one state is accessible
            P[:] = 0
            P[eqIdx[0], eqIdx[1]] = 1
            Z = 1

        # Calculating magnetization
        Ml = ( self.Ms * self.V_f * np.ma.sum(P * np.cos(self.theta - phi))  + self.M_af * self.V_af * np.ma.sum(P * np.cos(col_alpha - phi)) ) / Z
        Mt = ( self.Ms * self.V_f * np.ma.sum(P * np.sin(self.theta - phi))  + self.M_af * self.V_af * np.ma.sum(P * np.sin(col_alpha - phi)) ) / Z

        return Mt, Ml
