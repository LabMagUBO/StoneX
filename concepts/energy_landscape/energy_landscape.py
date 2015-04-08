#!/opt/local/bin/ipython-3.4
# -*- coding: utf-8 -*-

import sys
import numpy as np
import numpy.ma as ma
import scipy.ndimage as nd
from matplotlib import pyplot as pl

# Fonction
def new_plot():
    fig = pl.figure()
    return fig, fig.add_subplot(111, aspect='equal')
def draw(Z, axe):
    ax.imshow(E, interpolation = 'nearest')
    C = ax.contour(E, 10, colors='black', linewidth=.5)
    ax.clabel(C, inline=1, fontsize=10)

def export(name,fig):
    fig.savefig(name)
    pl.close()



# Fonction de theta, alpha
x_shift = 0.6
y_shift = 0
fx = 2
fy = 1
energy = lambda x, y: (np.cos(2 * np.pi * fx * (x - x_shift)) +  np.cos(2 * np.pi * fy * (y - y_shift)))/2

# Discretization
nb_theta = 2**3
nb_alpha = 2**3
theta = np.linspace(0, 1, nb_theta)
alpha = np.linspace(0, 1, nb_alpha)

#Theta, Alpha = np.meshgrid(theta, alpha)
Alpha, Theta = np.meshgrid(alpha, theta)

E = energy(Theta, Alpha)

#pl.rcParams['image.interpolation'] = 'none'
#pl.rcParams['image.resample'] = False
fig = pl.figure()
ax = fig.add_subplot(111, aspect='equal')
ax.set_xlabel('alpha')
ax.set_ylabel('theta')

#ax.contourf(Theta, Alpha, E, cmap=pl.cm.gray, resample=False)
cax = ax.imshow(E, interpolation = 'nearest',origin='upper')
cbar = fig.colorbar(cax)
C = ax.contour(E, 10, colors='black', linewidth=.5)
ax.clabel(C, inline=1, fontsize=10)

export('landscape.pdf', fig)



## On ajout un niveau
level = 0
mask = E <= level

E_ma = np.ma.array(E, mask=np.logical_not(mask))

fig = pl.figure()
ax = fig.add_subplot(111, aspect='equal')

cax = ax.imshow(E_ma, interpolation = 'nearest', origin='upper')
cbar = fig.colorbar(cax)
#cbar.ax.set_yticklabels(['< -1', '0', '> 1'])# vertically oriented colorbar


#ax.imshow(mask, cmap=pl.cm.gray, alpha=0.7, interpolation = 'nearest')
C = ax.contour(E_ma, 10, colors='black', linewidth=.5)
ax.clabel(C, inline=1, fontsize=10)

#ax.imshow(E_ma, interpolation = 'nearest')
#C = ax.contour(E_ma, 10, colors='black', linewidth=.5)
#ax.clabel(C, inline=1, fontsize=10)

pl.savefig('landscape_flooded.pdf')
pl.close()

### On numérote les zones

z_lab, z_num = nd.measurements.label(mask)
print("Nombre de zones : {}".format(z_num))

fig = pl.figure()
ax = fig.add_subplot(111, aspect='equal')
cax = ax.imshow(z_lab, interpolation = 'nearest', origin='upper')
cbar = fig.colorbar(cax)

pl.savefig('zones.pdf')
pl.close()


#### On recherche les zones communes
#nb est la longueur des tableau
# On scanne le bord theta du table
zones_equiv = set()
for i in np.arange(nb_theta):
    # On regarde si on est sur une zone «à risque» de chaque coté
    left = z_lab[i, 0]
    right = z_lab[i, -1]
    if left and right and (left != right):
        zones_equiv.add(tuple((min(left, right), max(left, right))))

# On regroupe les zones, en les rangeant par ordre croissant :
zones_equiv = np.sort(np.array(list(zones_equiv)), axis=0)

for i, val in enumerate(zones_equiv):
    # On retire l'élément de droite
    z_rm = z_lab == val[1]
    #pour le remplacer par l'indice de gauche
    z_lab[z_rm] = val[0]


# On scanne le bord alpha du table
zones_equiv = set()
for i in np.arange(nb_alpha):
    # On regarde si on est sur une zone «à risque» de chaque coté
    left = z_lab[0, i]
    right = z_lab[-1, i]
    if left and right and (left != right):
        zones_equiv.add(tuple((min(left, right), max(left, right))))

# On regroupe les zones, en les rangeant par ordre croissant :
zones_equiv = np.sort(np.array(list(zones_equiv)), axis=0)


for i, val in enumerate(zones_equiv):
    # On retire l'élément de droite
    z_rm = z_lab == val[1]
    #pour le remplacer par l'indice de gauche
    z_lab[z_rm] = val[0]


labels = np.unique(z_lab)
z_lab = np.searchsorted(labels, z_lab)


fig = pl.figure()
ax = fig.add_subplot(111, aspect='equal')
cax = ax.imshow(z_lab, interpolation = 'nearest', origin='upper')
cbar = fig.colorbar(cax)

pl.savefig('zones_period.pdf')
pl.close()


## selon le départ, on crée un tableau masqué
position = np.array([2/3, 1/2])
pos = position * np.array([ nb_theta, nb_alpha])
print(pos)
lab = z_lab[pos[0], pos[1]]
E_ma = ma.array(E, mask=z_lab != lab)

M = ma.sum((np.cos(theta)*E_ma.T).T) / ma.mean(E_ma)
print(M)
