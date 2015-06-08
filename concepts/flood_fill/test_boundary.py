import numpy as np
import scipy.ndimage.measurements as img

if False:
    ones = np.ones((10, 10))
    x = np.arange(10)
    y = np.arange(10)

    X, Y = np.meshgrid(x, y)
    Z = (X-3)**2 + (Y-3)**2

    mask = Z < 12

    rolled_mask = np.roll(np.roll(mask, -2, axis=0), -2, axis=1)

    print(rolled_mask)

    labels, nb = img.label(rolled_mask)

    print(labels)

    slices = img.find_objects(labels)
    print(slices)

    print(rolled_mask[slices[3]])


    print(labels[slices[0]])

if True:
    a = np.array([[0,1], [4,3], [2, 2], [0, 1], [1, 0]])

    b = np.array(list(set(tuple(p) for p in a)))
