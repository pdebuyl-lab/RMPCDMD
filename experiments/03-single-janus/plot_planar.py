#!/usr/bin/env python
"""
Display the planar concentration and velocity fields of a RMPCDMD simulation.
"""
import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('file', help="H5MD file")
parser.add_argument('--species', type=int, default=0)
args = parser.parse_args()

import h5py
import matplotlib.pyplot as plt
import numpy as np

with h5py.File(args.file, 'r') as f:
    c = f['fields/planar_concentration']
    v = f['fields/planar_velocity']
    x_min = c.attrs['x_min'][()]
    dx = c.attrs['dx'][()]
    y_min = c.attrs['y_min'][()]
    dy = c.attrs['dy'][()]

    c = c[:]
    v = v[:]

    N_x, N_y = c.shape[:2]

    # x and y must overshoot c.shape by one for pcolormesh
    x = x_min + np.arange(N_x+1)*dx
    y = y_min + np.arange(N_y+1)*dy

    c /= dx*dy

plt.subplot(121, aspect=1)
plt.pcolormesh(x, y, c[:,:,args.species].T, cmap=plt.cm.viridis)
plt.colorbar()

plt.subplot(122, aspect=1)

x, y = np.meshgrid(x[:-1], y[:-1])
plt.quiver(x, y, v[:,:,args.species,0].T, v[:,:,args.species,1].T)

plt.show()
