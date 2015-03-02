#!/opt/local/bin/ipython-2.7
# -*- coding: utf-8 -*-

#import os
import math
import sys
import numpy as np
#import pylab as pl
#from scipy import optimize, integrate

from StoneX.constants import *
#from classes import *
from StoneX.functions import *
from scipy import optimize, integrate

import logging

#####################
# Model classes
#####################

## General functions
class Model(object):
    """
        Principal class Model.
        Methods:
            — search_eqStates
            — find_0KIndex
            — find_barrier

        Inherited classes:
            — Stoner_conv
            — Stoner_2states
            — Franco_Conde
            — Carey_Respaud

        General method for inherited classes:
            — relax (previous: calculate_magnetization)
    """
    def search_eqStates(self):
        """ Method searching all the equilibrium states at 0K, stable, metastable and instable.
            It returns a (n,3) shape array, n being the index of equilibrium, with :
                eqStates[n, 0] => theta position (in degrees)
                eqStates[n, 1] => energy value
                eqStates[n, 2] => True if stable, False if not

            The search function is using optimize.brentq to find the exact root of the energy derivative (not needed). The nearest discrete theta position is then recorded.

            Before searching the equilibrium states, the energy MUST be updated.

            Parameters :
                — nb[512]: number of points to evaluate the energy function
                — extend[0.1 deg]
        """

        # Parameters

        #nb = 512    # number of points to evaluate the function
        #extend = np.radians(0.1)
        #th = np.radians(np.linspace(0, 360, nb))    # theta vector

        zeros = np.array([], dtype=float)       #Init the array which will contain the zeros of the gradEnergy

        # Check if there is an equilibrium between the edges
        if (self.gradEnergy(self.theta[0]) * self.gradEnergy(self.theta[-1] - 2 * np.pi)) <= 0:
            eq = optimize.brentq(self.gradEnergy, self.theta[-1] - 2 * np.pi, self.theta[0])
            zeros = np.append(zeros, eq)


        # Iterate on theta[i]
        for i, val in enumerate(self.theta[:-1]):
            if ( self.gradEnergy(val) * self.gradEnergy(self.theta[i+1]) <= 0 ):
                eq = optimize.brentq(self.gradEnergy, self.theta[i], self.theta[i+1])

                # If grad(theta[i+1]) = 0, eq will be registered twice
                # Test if the array is not empty first, then if the zero is not the same than the last one
                if  len(zeros) == 0 or zeros[-1] != eq:
                    zeros = np.append(zeros, eq)

        # Placing all the zeros between 0 and 2π, and ordering
        zeros = np.sort(zeros % (2 * np.pi))


        # Creation of the table (theta_eq, Energy, isStable)
        eqStates = np.zeros((0,3))
        mem_val = None
        for i, val in enumerate(zeros):
            if (val != mem_val):
                new_line = np.array([[val, self.energy(val), self.ggradEnergy(val) > 0]])
                eqStates = np.concatenate((eqStates, new_line), axis=0)
                mem_val = val

        # Test the number of equilibrium (in case)
        if ((len(eqStates) % 2) != 0 ) or (len(eqStates) == 0): print("ERREUR search_eqStates : nb d'état impair ou nul.")

        if np.min(eqStates[:,0]) < 0:
            print("minimum", eqStates[:,0], np.min(eqStates))
        return eqStates

    def search_eqStates_manual(self):
        """ Method searching all the equilibrium states at 0K, stable, metastable and instable.
            It returns a (n,3) shape array, n being the index of equilibrium, with :
                eqStates[n, 0] => theta position (in degrees)
                eqStates[n, 1] => energy value
                eqStates[n, 2] => True if stable, False if not

            The search function is using optimize.brentq to find the exact root of the energy derivative (not needed). The nearest discrete theta position is then recorded.

            Before searching the equilibrium states, the energy MUST be updated.

            Parameters :
                — nb[512]: number of points to evaluate the energy function
                — extend[0.1 deg]
        """

        #Init. the array which will contain the zeros of the gradEnergy
        zeros = np.array([], dtype=float)

        # Check if there is an equilibrium between the edges
        gradEnergy_theta = self.gradEnergy(self.theta)
        #if (self.gradEnergy(self.theta[0]) * self.gradEnergy(self.theta[-1] - 2 * np.pi)) <= 0:
        if (gradEnergy_theta[0] * gradEnergy_theta[-1]) <= 0:
            #eq = optimize.brentq(self.gradEnergy, self.theta[-1] - 2 * np.pi, self.theta[0])
            eq = 0
            zeros = np.append(zeros, eq)


        # Iterate on theta[i]
        for i, val in enumerate(self.theta[:-1]):
            if ( (gradEnergy_theta[i] * gradEnergy_theta[i+1]) <= 0 ):
                #eq = optimize.brentq(self.gradEnergy, self.theta[i], self.theta[i+1])
                eq = val

                # If grad(theta[i+1]) = 0, eq will be registered twice
                # Test if the array is not empty first, then if the zero is not the same than the last one
                if  len(zeros) == 0 or zeros[-1] != eq:
                    zeros = np.append(zeros, eq)

        # Placing all the zeros between 0 and 2π, and ordering
        zeros = np.sort(zeros % (2 * np.pi))


        # Creation of the table (theta_eq, Energy, isStable)
        eqStates = np.zeros((0,3))
        mem_val = None
        for i, val in enumerate(zeros):
            if (val != mem_val):
                new_line = np.array([[val, self.energy(val), self.ggradEnergy(val) > 0]])
                eqStates = np.concatenate((eqStates, new_line), axis=0)
                mem_val = val

        # Test the number of equilibrium (in case)
        if ((len(eqStates) % 2) != 0 ) or (len(eqStates) == 0): print("ERREUR search_eqStates : nb d'état impair ou nul.")

        if np.min(eqStates[:,0]) < 0:
            print("minimum", eqStates[:,0], np.min(eqStates))
        return eqStates

    def find_0KIndex(self, eqStates):
        """
            Using the eqStates, calculate in which wheel the moment is, starting with the more probable position as equilibrium.
            It returns the new equilibrium index (0K model).
            Very similar to find_barrier().
        """
        #Finding the magnetization position : maximum of the probability fonction
        theta_eq = self.theta[np.argmax(self.probM)]

        for i, val in enumerate(eqStates):
            # we find between which boundaries the magnetization is
            # knowing the equilibrium are altmost alternatively stable/unstable
            # (still need to take into account this possibility)
            if (theta_eq <= val[0]):
                # which state is stable?
                if val[2]:
                    idx = i
                else:
                    idx = i-1
                break
            elif (i == len(eqStates)-1):  # equilibrium position around 360deg
                if val[2]:
                    idx = i
                else:
                    idx = 0
        # Return the index
        #return eqStates[idx, 0]
        return idx

    def find_barrier(self, eqStates, idx):
        """
            Find the smallest barrier the magnization had to overcome, from the list of equilibrium position given by self.search_eqStates(), and the starting position.
            Return:     array(barrierValue, direction)
                — barrierValue: energy of the barrier
                — direction : 1 or -1, right or left barrier
        """

        # Starting from idx, we search the minimal barrier around it
        tab_eqStates = np.append(eqStates, eqStates)
        leftBarrier = eqStates[(idx-1) % len(eqStates), 1] - eqStates[(idx) % len(eqStates), 1]
        rightBarrier = eqStates[(idx+1) % len(eqStates), 1] - eqStates[(idx) % len(eqStates), 1]
        if leftBarrier <= rightBarrier:
            # M may go to the left
            #print "May go to left"
            return leftBarrier, -1
        else:
            # M may go to the right
            #print "May go to right"
            return rightBarrier, 1

## Stoner-Wolffart
class Stoner_conv(Model):
    """
        TO BE UPDATED
        Define the Stoner-Wohlfarth model (0K model), using the optimize.fmin_bfgs module.
    """
    def relax(self):
        if self.T != 0:
            print ("Warning : the Stoner model is a 0K model. The temperature will not be taken into account.")
        self.theta_eq =  self.search_eq()
        self.M_l = self.M_f * np.cos(self.theta_eq)
        self.M_t = self.M_f * np.sin(self.theta_eq)


    def search_eq(self, theta_error=1e-4):
        """ Function searching the 0K equilibrium from the current position of the magnetization.
            Usage : sample.search_eq()
            Output : the scalar equilibrium angle of the magnetization, in radians
        """
        #print np.degrees(self.theta_eq)
        # First optizimation (it should be enough)
        # Optimization of the volumic energy, the total energy having a too small value
        result_optimize = optimize.fmin_bfgs(lambda x: self.energy(x) / (self.V_f + self.V_af) , self.theta_eq+theta_error, disp=False, full_output=True)
        # In case of error (instable equilibrium)
        while result_optimize[-1] != 0:
            print ("Warning : nouvelle optimisation de l'équilibre")
            result_optimize = optimize.fmin_bfgs(lambda x: self.energy(x) / (self.V_f + self.V_af), result_optimize[0], disp=False, full_output=True)
            #print np.degrees(result_optimize[0])
        eq = result_optimize[0]
        #print np.degrees(eq)
        eq = eq % (2*np.pi)
        #print np.degrees(eq)

        return eq[0]

class Stoner_2states(Model):
    """
        Stoner-Wohlfarth model with manual calculation of the energy barrier (2 states).
    """
    def model(self):
        self.model = 'Stoner-Wohlfarth'
        self.logger.warn("Stoner-Wohlfarth  model loaded.")

    def relax(self):
        """
            Find the new equilibrium position, depending on the new energy function.
            New function, replacing self.calculate_magnetization()
        """

        # First, we search all the equilibrium position
        eqStates = self.search_eqStates_manual()

        # Reading the new equilibrium well
        theta_eq = eqStates[ self.find_0KIndex(eqStates), 0]

        #Clean the probability position array
        self.probM = np.zeros(self.theta_n)
        #The update it
        self.probM[theta_eq / 2 / np.pi * self.theta_n] = 1


## Garcia-Otero
class Garcia_Otero(Model):
    """
        Garcia-Otero model
        À compléter
    """
    def model(self):
        self.model = 'Garcia-Otero'
        self.logger.warn("Garcia-Otero model loaded.")

    ### Main function ####
    def relax(self):
        """ Find the equilibrium position, depending on the temperature."""

        # Search all the equilibrium positions
        eqStates = self.search_eqStates()

        # Search the actual position inside eqStates
        actual_idx = self.find_0KIndex(eqStates)

        # Find the smallest next barrier, and the position where M could go, starting from the actual position (actual_idx)
        barrier, go_to = self.find_barrier(eqStates, actual_idx)

        # If the thermal energy overcome the barrier
        while k_B * self.T * np.log(f0 * tau_mes) >= barrier :
            # We check if the next equilibrium is energetically better
            next_idx = (actual_idx + go_to * 2) % len(eqStates)
            # If the energy state behind the barrier is lower
            if eqStates[next_idx, 1] < eqStates[actual_idx, 1]:
                # the magnetization moves to the next equilibrium
                actual_idx = next_idx
                # searching the smallest barrier, using the new index
                barrier, go_to = self.find_barrier(eqStates, actual_idx)
            else:
                # we stay and update the position
                break

        # Equilibrium angle depending on the index
        theta_eq = eqStates[actual_idx, 0]

        #Clean the probability position array
        self.probM = np.zeros(self.theta_n)
        #Then update it
        self.probM[theta_eq / 2 / np.pi * self.theta_n] = 1


## Franco-Conde
class Franco_Conde(Model):
    """ Franco & Conde model"""
    def model(self):
        self.model = 'Franco-Conde'
        self.logger.warn("Franco-Conde model loaded.")

    def check_newWheel(self, eqStates, pre_index):
        """
            Return the new wheel index, using the different energy wheels and the temperature.
        """
        global k_B, tau_mes, f0


        # Thermal energy :
        T_sunami = k_B * self.T * np.log(f0 * tau_mes)

        # High temperature limit
        if T_sunami >= (np.amax(eqStates[:,1]) - np.amin(eqStates[:, 1])):
            # Return the index of the minimum energy equilibrium
            print("High temp")
            return np.argmin(eqStates[:,1])

        ## Intermediate or low temperature limit
        # Relative index of the right and left accessible barrier
        #sens = 1
        right_relInd_barrier = 1
        left_relInd_barrier = -1

        # searching loop
        while True:
            left_barrier = eqStates[(pre_index + left_relInd_barrier) % len(eqStates), 1] - eqStates[pre_index, 1]
            right_barrier = eqStates[(pre_index + right_relInd_barrier) % len(eqStates), 1] - eqStates[pre_index, 1]

            # if both barrier are too high
            if T_sunami < min(left_barrier, right_barrier):
                #nothing to do
                new_index = pre_index
                return new_index

            else:
                # which direction for the weakest barrier
                if left_barrier < right_barrier:
                    sens = 1
                    dTo_eq = left_relInd_barrier      # distance to the next eq in the left direction
                else:
                    sens = -1
                    dTo_eq = right_relInd_barrier    # distance to the next eq in the right direction

                # If the next equilibrium is lower
                if (eqStates[ (pre_index + sens + dTo_eq) % len(eqStates), 1] - eqStates[pre_index, 1]) < 0:
                    # change of equilibrium
                    pre_index += sens + dTo_eq
                    # reset the energy barrier relative index
                    right_relInd_barrier = 1
                    left_relInd_barrier = -1

                    # redo the loop
                else:
                    # Update the energy barrier
                    if sens == 1:
                        right_relInd_barrier = (right_relInd_barrier + 2)
                    else:
                        left_relInd_barrier = (left_relInd_barrier - 2)


            if (right_relInd_barrier % len(eqStates)) == (left_relInd_barrier % len(eqStates)):
                return np.argmin(eqStates[:,1])

    def find_allowedStates(self, eqStates, new_index):
        """
            Return a list of available states domains.
        """
        # Thermal energy
        T_sunami = k_B * self.T * np.log(f0 * tau_mes)

        # Determine if the high temperature limit is reached
        X = np.linspace(0, 2*np.pi, 256)
        max_energy_theta, = optimize.fmin(lambda x: -self.energy(x), np.argmax(self.energy(X)), disp=False)

        if self.energy(max_energy_theta) <= eqStates[new_index, 1] + T_sunami:
            # high temperature limit
            return np.array([[0, 2*np.pi]])

        # adding the 2π boundary to the eqStates
        eqStates = np.append(eqStates, [[np.pi * 2, self.energy(np.pi * 2), False]], axis=0)
        #print np.degrees(eqStates[:, 0])
        #domains = np.zeros(len(eqStates) + 2, dtype=float)
        #domains = np.array([[0., 0.]], dtype=float)
        domains = np.array([0.], dtype=float)
        # Energy maximum
        E_max = eqStates[new_index, 1] + T_sunami

        # Temperature relative energy function
        f = lambda x: self.energy(x) - E_max

        #print 'eq ', np.degrees(eqStates[:,0])

        #Indexes
        j = 0
        i = 0
        while i < len(eqStates):
            #test if the thermal energy is in between the domain range
            # if higher
            if ( (f(domains[j]) > 0) and (f(eqStates[i,0]) > 0) ) :
                #print "higher", i
                # the domain is not reachable
                # if stable
                if eqStates[i, 2]:
                    # set the boundary to the next unstable position, which can only be higher
                    if (i + 1) < len(eqStates):
                        domains[j] = eqStates[i+1, 0]
                    else:
                        domains[j] = 2*np.pi
                    i += 1
                else:
                    # go to the next domain
                    domains[j] = eqStates[i, 0]
                #print np.degrees(domains)
            #if lower
            elif ( (f(domains[j]) < 0) and (f(eqStates[i,0]) < 0) ):
                # do nothing, domains[j] can still be reached
                # except if we are at the end
                if i == len(eqStates) - 1:
                    domains = np.append(domains, [2*np.pi])
                #print "lower", i
                #print np.degrees(domains)

            else:
                # not lower, so in between
                #searching the zero
                #print "in between", i
                zero = optimize.brentq(f, domains[j], eqStates[i, 0])

                #if stable
                if eqStates[i, 2]:
                    #the zero correspond to the left boundary
                    domains[j] = zero
                    domains = np.append(domains, [eqStates[i, 0], eqStates[i, 0]])

                elif (eqStates[i,0] == 2*np.pi) and (f(eqStates[i, 0]) < 0):
                    #case if the right boundary is 360 degrees and reachable
                    # the zero became the left boundary domain, and 360 the right one
                    domains[j] = zero
                    domains = np.append(domains, [eqStates[i, 0]])

                elif (eqStates[i,0] == 2*np.pi) and (f(eqStates[i, 0]) > 0):
                    #case if the right boundary is 360 degrees and non reachable
                    # the zero became the right boundary domain
                    domains = np.append(domains, [ zero])
                else:
                    # closing the right boundary of the current domain and adding the next one
                    domains = np.append(domains, [ zero, eqStates[i, 0]])
                #print np.degrees(domains)


            j = len(domains) - 1
            i += 1
        #print np.degrees(domains)

        # Joining the domains
        i = 0
        while i < len(domains)-1:
            if domains[i] == domains[i+1]:
                domains = np.delete(domains, [i, i+1])
                i = 0

            i += 1



        # Taking out the extra bound, if odd number and 2*pi boundary
        if (len(domains) % 2 == 1) and domains[-1] == 2*np.pi:
            domains = np.delete(domains, len(domains)-1)
        elif len(domains) % 2 == 1:
            print ("Erreur programmation, nombre de domaine impair. Un cas n'est pas pris en compte.")

        #Reshaping
        domains = domains.reshape(len(domains)/2, 2)
        #print (np.degrees(domains))

        # Determine in which domain the magnetization is
        mask = np.zeros(len(domains), dtype=bool)

        for i in np.arange(len(domains)):
            # if inside the boundary, make it true, else false
            if (eqStates[new_index, 0] <= domains[i,1]) and (eqStates[new_index, 0] >= domains[i,0]):
                mask[i] = True
            else:
                mask[i] = False
        #print mask, (mask[0] + mask[-1]), domains[-1,1] % (2*np.pi)
        #if the left and right boundary are joined ([0,x] and [y, 2π]), and one is reachable, make the other reachable
        if (domains[-1,1] % (2*np.pi) == domains[0,0]) and ((mask[0] + mask[-1]) == True):
            mask[0] = True
            mask[-1] = True


        return domains[mask]

    def is_inside(self, angle, domains):
        for boundaries in domains:
            if (angle <= boundaries[1]) and (angle >= boundaries[0]):
                return True
        return False


    # Main function
    def relax(self):
        # Warning if T=0
        if self.T == 0:
            self.logger.error('T=0 incompatible with Franco-Conde model.')
            sys.exit("End of Program.")
        # First, we need to know the energy extrema
        eqStates = self.search_eqStates()
        #print(np.degrees(eqStates[:,0]), "\n", eqStates[:,1:3])
        # Find the starting range
        pre_index = self.find_0KIndex(eqStates)
        #print(pre_index)
        pre_theta = eqStates[pre_index, 0]
        #print(np.degrees(pre_theta))

        # Determine if the position changes with temperature
        new_index = self.check_newWheel(eqStates, pre_index)
        #print(new_index)

        # updating the theta
        #self.theta_eq = eqStates[new_index, 0]

        # Prévoir température nulle

        # Search the states which can be reached
        domains = self.find_allowedStates(eqStates, new_index)


        #print ("domains", np.degrees(domains))

        ## Now we redefine the magnetization probability array «probM».
        ## The probability is null outside domains, and depends on Boltzman inside
        # Maxwell-Boltzmann function probability (unnormed and without exp to be able to calculate exp(x>700))
        P = lambda x: - self.energy(x) / k_B / self.T
        # Loop on probM elements
        for i, angle in enumerate(self.theta):
            if self.is_inside(angle, domains):
                self.probM[i] = P(angle)
            else:
                self.probM[i] = -float("inf")

        # Applying the exponential to the shifted energy ratio (E/kT)
        self.probM = np.exp(self.probM - np.max(self.probM))

        # Norming probM
        self.probM = self.probM / np.sum(self.probM)


## Franco-Conde
class Franco_Conde_old(Model):
    """ Franco & Conde model"""

    def check_newWheel(self, eqStates, pre_index):
        """
            Return the new wheel index, using the different energy wheels and the temperature.
        """
        global k_B, tau_mes, f0

        # Thermal energy :
        T_sunami = k_B * self.T * np.log(f0 * tau_mes)

        # High temperature limit
        if T_sunami >= (np.amax(eqStates[:,1]) - np.amin(eqStates[:, 1])):
            # Return the index of the minimum energy equilibrium
            print("High temp")
            return np.argmin(eqStates[:,1])

        ## Intermediate or low temperature limit
        # Relative index of the right and left accessible barrier
        #sens = 1
        right_relInd_barrier = 1
        left_relInd_barrier = -1

        # searching loop
        while True:
            left_barrier = eqStates[(pre_index + left_relInd_barrier) % len(eqStates), 1] - eqStates[pre_index, 1]
            right_barrier = eqStates[(pre_index + right_relInd_barrier) % len(eqStates), 1] - eqStates[pre_index, 1]

            # if both barrier are too high
            if T_sunami < min(left_barrier, right_barrier):
                #nothing to do
                new_index = pre_index
                return new_index

            else:
                # which direction for the weakest barrier
                if left_barrier < right_barrier:
                    sens = 1
                    dTo_eq = left_relInd_barrier      # distance to the next eq in the left direction
                else:
                    sens = -1
                    dTo_eq = right_relInd_barrier    # distance to the next eq in the right direction

                # If the next equilibrium is lower
                if (eqStates[ (pre_index + sens + dTo_eq) % len(eqStates), 1] - eqStates[pre_index, 1]) < 0:
                    # change of equilibrium
                    pre_index += sens + dTo_eq
                    # reset the energy barrier relative index
                    right_relInd_barrier = 1
                    left_relInd_barrier = -1

                    # redo the loop
                else:
                    # Update the energy barrier
                    if sens == 1:
                        right_relInd_barrier = (right_relInd_barrier + 2)
                    else:
                        left_relInd_barrier = (left_relInd_barrier - 2)


            if (right_relInd_barrier % len(eqStates)) == (left_relInd_barrier % len(eqStates)):
                return np.argmin(eqStates[:,1])


    def find_allowedStates(self, eqStates, new_index):
        """
            Return a list of available states domains.
        """
        # Thermal energy
        T_sunami = k_B * self.T * np.log(f0 * tau_mes)

        # Determine if the high temperature limit is reached
        X = np.linspace(0, 2*np.pi, 256)
        max_energy_theta, = optimize.fmin(lambda x: -self.energy(x), np.argmax(self.energy(X)), disp=False)

        if self.energy(max_energy_theta) <= eqStates[new_index, 1] + T_sunami:
            # high temperature limit
            return np.array([[0, 2*np.pi]])

        # adding the 2π boundary to the eqStates
        eqStates = np.append(eqStates, [[np.pi * 2, self.energy(np.pi * 2), False]], axis=0)
        #print np.degrees(eqStates[:, 0])
        #domains = np.zeros(len(eqStates) + 2, dtype=float)
        #domains = np.array([[0., 0.]], dtype=float)
        domains = np.array([0.], dtype=float)
        # Energy maximum
        E_max = eqStates[new_index, 1] + T_sunami

        # Temperature relative energy function
        f = lambda x: self.energy(x) - E_max

        #print 'eq ', np.degrees(eqStates[:,0])

        #Indexes
        j = 0
        i = 0
        while i < len(eqStates):
            #test if the thermal energy is in between the domain range
            # if higher
            if ( (f(domains[j]) > 0) and (f(eqStates[i,0]) > 0) ) :
                #print "higher", i
                # the domain is not reachable
                # if stable
                if eqStates[i, 2]:
                    # set the boundary to the next unstable position, which can only be higher
                    if (i + 1) < len(eqStates):
                        domains[j] = eqStates[i+1, 0]
                    else:
                        domains[j] = 2*np.pi
                    i += 1
                else:
                    # go to the next domain
                    domains[j] = eqStates[i, 0]
                #print np.degrees(domains)
            #if lower
            elif ( (f(domains[j]) < 0) and (f(eqStates[i,0]) < 0) ):
                # do nothing, domains[j] can still be reached
                # except if we are at the end
                if i == len(eqStates) - 1:
                    domains = np.append(domains, [2*np.pi])
                #print "lower", i
                #print np.degrees(domains)

            else:
                # not lower, so in between
                #searching the zero
                #print "in between", i
                zero = optimize.brentq(f, domains[j], eqStates[i, 0])

                #if stable
                if eqStates[i, 2]:
                    #the zero correspond to the left boundary
                    domains[j] = zero
                    domains = np.append(domains, [eqStates[i, 0], eqStates[i, 0]])

                elif (eqStates[i,0] == 2*np.pi) and (f(eqStates[i, 0]) < 0):
                    #case if the right boundary is 360 degrees and reachable
                    # the zero became the left boundary domain, and 360 the right one
                    domains[j] = zero
                    domains = np.append(domains, [eqStates[i, 0]])

                elif (eqStates[i,0] == 2*np.pi) and (f(eqStates[i, 0]) > 0):
                    #case if the right boundary is 360 degrees and non reachable
                    # the zero became the right boundary domain
                    domains = np.append(domains, [ zero])
                else:
                    # closing the right boundary of the current domain and adding the next one
                    domains = np.append(domains, [ zero, eqStates[i, 0]])
                #print np.degrees(domains)


            j = len(domains) - 1
            i += 1
        #print np.degrees(domains)

        # Joining the domains
        i = 0
        while i < len(domains)-1:
            if domains[i] == domains[i+1]:
                domains = np.delete(domains, [i, i+1])
                i = 0

            i += 1

        #print np.degrees(domains)

        # Taking out the extra bound, if odd number and 2*pi boundary
        if (len(domains) % 2 == 1) and domains[-1] == 2*np.pi:
            domains = np.delete(domains, len(domains)-1)
        elif len(domains) % 2 == 1:
            print ("Erreur programmation, nombre de domaine impair. Un cas n'est pas pris en compte.")

        #Reshaping
        domains = domains.reshape(len(domains)/2, 2)
        #print (np.degrees(domains))

        # Determine in which domain the magnetization is
        mask = np.zeros(len(domains), dtype=bool)

        for i in np.arange(len(domains)):
            # if inside the boundary, make it true, else false
            if (eqStates[new_index, 0] <= domains[i,1]) and (eqStates[new_index, 0] >= domains[i,0]):
                mask[i] = True
            else:
                mask[i] = False
        #print mask, (mask[0] + mask[-1]), domains[-1,1] % (2*np.pi)
        #if the left and right boundary are joined ([0,x] and [y, 2π]), and one is reachable, make the other reachable
        if (domains[-1,1] % (2*np.pi) == domains[0,0]) and ((mask[0] + mask[-1]) == True):
            mask[0] = True
            mask[-1] = True

        return domains[mask]


    # Main function
    def relax(self):
        # First, we need to know the energy extrema
        eqStates = self.search_eqStates()
        #print(np.degrees(eqStates[:,0]), "\n", eqStates[:,1:3])
        # Find the starting range
        pre_index = self.find_0KIndex(eqStates)
        #print(pre_index)
        pre_theta = eqStates[pre_index, 0]
        #print(np.degrees(pre_theta))

        # Determine if the position changes with temperature
        new_index = self.check_newWheel(eqStates, pre_index)
        #print(new_index)

        # updating the theta
        #self.theta_eq = eqStates[new_index, 0]

        # Search the states which can be reached
        domains = self.find_allowedStates(eqStates, new_index)


        print ("domains", np.degrees(domains))

        # Integrate
        Z_part, int_long, int_trans = 0., 0., 0.
        P = lambda x: np.exp(- self.energy(x) / k_B / self.T)

        for i in np.arange(len(domains)):

            Z, Z_err = integrate.quad(P, domains[i, 0], domains[i,1])
            #print ("Z = ", Z, 'Z_err = ', Z_err )
            Z_part += Z

        # Check if the partition function is not zero
        if Z_part == 0:
            print ("Low temperature limit : 0K model")

            # Uses the theta_eq to calculate the magnetization
            self.M_l = self.M_f * np.cos(self.theta_eq)
            self.M_t = self.M_f * np.sin(self.theta_eq)

        else:
            # Plot the probability function
            """
            fig = pl.figure()
            ax = fig.add_subplot(111)
            ax.grid(True)
            ax.set_title("Probability")
            ax.set_xlim([0,360])
           """

            #print("In-range or high-temperature. Integrates domains...")

            # Calculate the probability, and each magnetization value
            for i in np.arange(len(domains)):
                sampling = np.linspace(domains[i, 0], domains[i,1], 1000)

                """
                ax.plot(np.degrees(sampling), P(sampling)/Z_part, 'r-')
                pl.draw()
                """
                # Magnetization mean in the domain
                L, L_err = integrate.quad(lambda x: P(x) * np.cos(x), domains[i, 0], domains[i,1])
                Tr, Tr_err = integrate.quad(lambda x: P(x) * np.sin(x), domains[i, 0], domains[i,1])

                int_long += L
                int_trans += Tr
                #print(int_long, int_trans, Z_part)

            # Create/Update the magnetization
            self.M_l = self.M_f * int_long / Z_part
            self.M_t = self.M_f * int_trans / Z_part
            #print (int_long / Z_part)


class Carey_Respaud(Model):
    pass
