#!/opt/local/bin/ipython-3.4
# -*- coding: utf-8 -*-

import sys
import numpy as np
import numpy.ma as ma


import numpy as np

np.random.seed(42) # Just to get always the same output

n = 10

a = np.random.randint(-10, 10, (n, n))
x = np.arange(n)
y = np.arange(n)

a_masked = ma.masked_array(a)
a_masked[a >= 0] = ma.masked

if True:
    E_ma = (x * a_masked.T).T + y * a_masked
    m = E_ma.mean()

if False: 
    E = 0

if False:
    print ('a:', a)
    print ('a.argmin():', a.argmin())
    print ('a.min():',a.min())
    print ("")
    print ('a_masked:', a_masked)
    print ('a_masked.argmin():', a_masked.argmin())
    print ('a_masked.min():',a_masked.min())
