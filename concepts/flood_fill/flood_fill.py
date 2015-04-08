import sys
import numpy as np
import scipy as sp
from matplotlib import pyplot as pl

import scipy.ndimage as nd

f = lambda x, y: x**2 + y**2 + 1.1*np.sin(4*y) + 1.1*np.cos(2*(x+y))**2

pos = np.array([0, 1])

n = 50
x = np.linspace(-2., 2., n)
y = np.linspace(-2., 2., n)
X, Y = np.meshgrid(x, y)

Z = f(X, Y)


fig = pl.figure()
ax = fig.add_subplot(111, aspect='equal')
#ax2 = fig.add_subplot(214, aspect='equal')
#ax.pcolormesh(X, Y, Z, cmap = pl.cm.hot)
ax.contourf(X, Y, Z, 15, alpha=0.9, cmap=pl.cm.hot)
C = ax.contour(X, Y, Z, 15, colors='black', linewidth=.5)
ax.clabel(C, inline=1, fontsize=10)


seuil = 0.3

zones_risky = Z < seuil

zones_labeled, zones_num = nd.measurements.label(zones_risky)

ax.contourf(X, Y, zones_risky, 10, alpha=0.5, cmap=pl.cm.gray)

print("Nombre de zones : ",zones_num)

ax.contourf(X, Y, zones_labeled, 15, alpha=0.4, cmap=pl.cm.hot)
pl.savefig('test.pdf')
#actual_label = zones_labeled(pos[0], pos[1])

#zones_actual =
#ax2.contourf(X, Y, zones_labe, 10, alpha=0.5, cmap=pl.cm.gray))
