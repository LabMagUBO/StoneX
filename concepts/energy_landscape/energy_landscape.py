#!/opt/local/bin/ipython-3.4
# -*- coding: utf-8 -*-

import sys
import numpy as np
import numpy.ma as ma
import scipy.ndimage as nd
from scipy import optimize
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

# Expression de la fonction
# X est un vecteur : X[alpha, theta]
# f comprix entre -1 et 1
# X aussi
X_shift = np.array([0, 0])
#fX = np.array([1, 2])

H = -2
J = 1
Kaf = 2
T = 0.15

energy = lambda X: -H * np.cos(2*np.pi * (X[1] - X_shift[1])) +\
    np.sin(2*np.pi * (X[1] - X_shift[1]))**2 +\
    Kaf * np.sin(2*np.pi * (X[0] - X_shift[0]))**2 - J * np.cos(2*np.pi * (X[1] - X[0]))
#energy = lambda X: X[1]
#energy = lambda X: np.cos(2 * np.pi * fX[0] * (X[0] - X_shift[0])) + np.cos(2 * np.pi * (X[0] - X_shift[0]))  np.cos(2 * np.pi * fX[1] * (X[1] - X_shift[1]))

#energy = lambda X: X[1]#np.sin(2 * np.pi * X[0])

# Discretization [theta, alpha]
nb = np.array([2**8, 2**8])#2**6
theta = np.linspace(0, 1, nb[1])
alpha = np.linspace(0, 1, nb[0])

# Départ [theta, alpha]:
X_ini = np.array([0, 0])
#X_ini = np.array([0.5, 0])

# Finding minimum
N_ini = np.floor(X_ini * nb) % nb
#X_eq = optimize.fmin_bfgs(energy, X_ini, disp=True)

print("Position i, j initiale : {}".format(N_ini))



#Theta, Alpha = np.meshgrid(theta, alpha)
Theta, Alpha = np.meshgrid(theta, alpha)

E = energy([Alpha, Theta])

# On recherche l'équilibre
def search_eq(N_ini, E):
    found = False
    k = 0
    nb = np.shape(E)
    print("Shape", nb)
    print("Boucle while")
    N_eq = np.zeros((1, 2))
    N_eq[0] = N_ini
    while not found:
        k = k+1
        print("Boucle n", k)

        E_ma = np.ma.array(E)
        E_ma.mask = True

        for i in np.arange(-1, 2) + N_eq[k-1][0]:
            for j in np.arange(-1, 2) + N_eq[k-1][1]:
                E_ma.mask[i % nb[0], j % nb[1]] = False

        E_ma[N_ini[0], N_ini[1]] = ma.masked

        print("indice du minimum", E_ma.argmin())

        N_eq = np.append(N_eq, [[np.floor(E_ma.argmin() / nb[1]), E_ma.argmin() % nb[1]]], axis=0)
        #print("N_eq :", N_eq)

        print("mouvement ({0},{1})".format(int(N_eq[k][0]-N_eq[k-1][0]), int(N_eq[k][1]-N_eq[k-1][1])))

        #print(E_ma)
        #print(E_ma.mask)
        print("Ini", N_eq[k-1], "eq", N_eq[k])

        E_ini = E[N_eq[k-1][0], N_eq[k-1][1]]
        E_eq = E[N_eq[k][0], N_eq[k][1]]

        print("Énergie", E_ini, E_eq)

        if E_eq < E_ini:
            print("better")

        elif (E_eq == E_ini):
            print("egualité")
            # On regarde si on est pas déjà passé par là
            if k < 2:
                continue
            elif (N_eq[k] == N_eq[k-2]).all():
                print("already here")
                print(E_ma)
                return N_eq
        else:
            print("I stay")
            print(E_eq - E_ini)
            print(E_ma)
            return N_eq


        if k == 1000:
            print("pas de convergence ", k)
            break



    return N_eq

N_eq = search_eq(N_ini, E)
#print("N_eq", N_eq)
X_eq = N_eq / nb
#print(X_eq)

#print("énergie")
#for i, val in enumerate(X_eq):
#    print(energy(val))

#pl.rcParams['image.interpolation'] = 'none'
#pl.rcParams['image.resample'] = False
fig = pl.figure()
ax = fig.add_subplot(111, aspect='equal')
ax.set_xlabel('X[1], theta')
ax.set_ylabel('X[0], alpha')

#ax.contourf(Theta, Alpha, E, cmap=pl.cm.gray, resample=False)
cax = ax.imshow(E,
    interpolation = 'nearest',
    #origin='upper',
    extent=(0,1,1,0)
)
cbar = fig.colorbar(cax)
C = ax.contour(E, 30,
    colors='black',
    linewidth=.5,
    extent=(0,1,0,1)
)
ax.clabel(C, inline=1, fontsize=10)
pl.grid()

ax.plot(X_ini[1] % 1, X_ini[0] % 1, 'ro', label="Ini.")
ax.plot(X_eq[1:-1, 1] % 1, X_eq[1:-1, 0] % 1, 'go', label="Recherche")
ax.plot(X_eq[-1, 1] % 1, X_eq[-1, 0] % 1, 'mo', label="Eq.")

ax.legend()
export('landscape.pdf', fig)


# On oublie les états intermédiaires
X_eq = X_eq[-1]

print("équilibre final", X_eq)
## On ajout un niveau
level = T + energy(X_eq)
mask = E <= level

E_ma = np.ma.array(E, mask=np.logical_not(mask))

fig = pl.figure()
ax = fig.add_subplot(111, aspect='equal')

cax = ax.imshow(E_ma, interpolation = 'nearest', origin='upper', extent=(0,1,1,0))
cbar = fig.colorbar(cax)
#cbar.ax.set_yticklabels(['< -1', '0', '> 1'])# vertically oriented colorbar


#ax.imshow(mask, cmap=pl.cm.gray, alpha=0.7, interpolation = 'nearest')
C = ax.contour(E_ma, 10, colors='black', linewidth=.5, extent=(0,1,0,1))
ax.clabel(C, inline=1, fontsize=10)

#ax.imshow(E_ma, interpolation = 'nearest')
#C = ax.contour(E_ma, 10, colors='black', linewidth=.5)
#ax.clabel(C, inline=1, fontsize=10)

ax.plot(X_eq[1] % 1, X_eq[0] % 1, 'mo', label="Eq.")

pl.savefig('landscape_flooded.pdf')
pl.close()

### On numérote les zones

z_lab, z_num = nd.measurements.label(mask)
print("Nombre de zones : {}".format(z_num))
print(z_lab)

fig = pl.figure()
ax = fig.add_subplot(111, aspect='equal')
cax = ax.imshow(z_lab, interpolation = 'nearest', origin='upper', extent=(0,1,1,0))
cbar = fig.colorbar(cax)

pl.savefig('zones.pdf')
pl.close()


#### On recherche les zones communes
#nb est la longueur des tableau
# On scanne le bord theta du table
zones_equiv = set()
for i in np.arange(nb[0]):
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
for i in np.arange(nb[1]):
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
cax = ax.imshow(z_lab, interpolation = 'nearest', origin='upper', extent=(0,1,1,0))
cbar = fig.colorbar(cax)

pl.savefig('zones_period.pdf')
pl.close()


## selon le départ, on crée un tableau masqué
N_eq = X_eq * nb

lab = z_lab[N_eq[0], N_eq[1]]
E_ma = ma.array(E, mask=z_lab != lab)

fig = pl.figure()
ax = fig.add_subplot(111, aspect='equal')
cax = ax.imshow(E_ma, interpolation = 'nearest', origin='upper', extent=(0,1,1,0))
cbar = fig.colorbar(cax)

pl.savefig('zone_accessible.pdf')
pl.close()

#M = ma.sum((np.cos(theta)*E_ma.T).T) / ma.mean(E_ma)
#print(M)
