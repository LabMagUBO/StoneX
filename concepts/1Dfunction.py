#!/usr/bin/ipython-2.7 --pylab
# -*- coding: utf-8 -*-

import numpy as np
import pylab as pl
from scipy import optimize
from time import sleep

#function to plot
def plot_state():
    pl.close()
    pl.plot(x, f(x), 'b-', linewidth=2)
    pl.grid()
    pl.plot(r_pos, f(r_pos), 'ro')
    pl.show()

def print_state():
    x = r_pos[r_pos.size - 1]
    print("Current position : ( %f , %f )" % (x , f(x)) )

def update_plot():
    x = new_pos[0]
    print "updating plot", x, f(x)
    line1.set_ydata(x, f(x))
    fig.canvas.draw()


#Définition de la fonction
def f(x):
    return x**2 + 10*np.sin(x)

# initialisation des variables
x = np.linspace(-5,5, 64)
last_pos = 2
r_pos = np.array([last_pos])

# position initiale
pl.ion()
fig = pl.figure()
ax1 = pl.subplot(2, 1, 1)
ax2 = pl.subplot(2,1,2)
ax1.grid()
ax2.grid()
ax1.set_xlabel('x')
ax1.set_ylabel('E(x)')
ax2.set_xlabel('Iterations')
ax2.set_ylabel('x(i)')
ax1.plot(x, f(x), 'b-', linewidth=3)
ax1.plot(last_pos, f(last_pos), 'ro')
ax2.plot(0, last_pos, 'ro')

pl.draw()

print_state()
sleep(1)


i = 0
while True:
    i = i+1
    new_pos = optimize.fmin_bfgs(f, last_pos, maxiter=1, disp=False, epsilon=0.01)
    r_pos = np.append(r_pos, new_pos)
    ax1.plot(new_pos, f(new_pos), 'ro')
    ax1.plot(r_pos, f(r_pos), 'r-')
    ax2.plot(i, new_pos, 'ro')
    ax2.plot(r_pos, 'r-')
    pl.draw()
    print_state()
    sleep(0.5)


    if new_pos == last_pos:
        break
    else:
        last_pos = new_pos
