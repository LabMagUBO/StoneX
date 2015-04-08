#!/usr/bin/ipython-3.4
# -*- coding: utf-8 -*-

# Modules to import
import numpy as np
import pylab as pl
from scipy import optimize
from time import sleep

import sys

pl.ion()
fig = pl.figure(dpi=150)
ax = pl.subplot(1,1,1)

# Functions
def f(x,y):
    return x*x + 1 +y**2

def F(X):
    x = X[0]
    y = X[1]
    return f(x,y)

n = 1024
x = np.linspace(-10, 10, n)
y = np.linspace(-10, 10, n)
x_grid,y_grid = np.meshgrid(x, y)

#pl.axes([0.025, 0.025, 0.95, 0.95])

ax.contourf(x_grid, y_grid, f(x_grid, y_grid), 10, alpha=.75, cmap=pl.cm.hot)
C = ax.contour(x_grid, y_grid, f(x_grid, y_grid), 10, colors='black', linewidth=.5)
ax.clabel(C, inline=1, fontsize=10)
pl.grid()
#pl.xticks(())
#pl.yticks(())
pl.draw()


# Points de départ
points_depart = np.array([[0,0.2], [9,9], [-7, 9]])
for X_start in points_depart:
    positions = np.array([ X_start, ])

    ax.plot(X_start[0], X_start[1], 'ro')
    pl.draw()

    #Recherche du minimum et convergence
    while True:
        X_eq = optimize.fmin_bfgs(F, X_start, maxiter=1, disp=False)
        print( X_eq, F(X_eq))
        positions = np.append(positions, X_eq)
        ax.plot(X_eq[0], X_eq[1], 'go')
        ax.plot([X_start[0], X_eq[0]], [X_start[1], X_eq[1]], 'g-')
        pl.draw()
        if (X_eq == X_start).all():
            print ("Équilibre atteint")
            sleep(0.1)
            ax.plot(X_start[0], X_start[1], 'bo')
            break
        X_start = X_eq
        sleep(0.1)
