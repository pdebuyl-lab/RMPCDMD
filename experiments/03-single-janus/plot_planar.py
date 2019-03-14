#!/usr/bin/env python
"""
Display the planar concentration and velocity fields of a RMPCDMD simulation.
"""
import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('file', help="H5MD file", nargs='+')
parser.add_argument('--species', type=int, default=0)
args = parser.parse_args()

import h5py
import matplotlib.pyplot as plt
import numpy as np

c_all = []
v_all = []
v_av_all = []

for filename in args.file:
    with h5py.File(filename, 'r') as f:
        c = f['fields/planar_concentration']
        v = f['fields/planar_velocity']
        x_min = c.attrs['x_min'][()]
        dx = c.attrs['dx'][()]
        y_min = c.attrs['y_min'][()]
        dy = c.attrs['dy'][()]
        thickness = c.attrs['thickness'][()]

        c = c[:]
        v = v[:]

        v_average = c[:,:,0,None]*v[:,:,0,:] + c[:,:,1,None]*v[:,:,1,:]
        mask = c[:,:,0] + c[:,:,1] > 0
        v_average[mask,:] /= c[mask,0,None] + c[mask,1,None]

        N_x, N_y = c.shape[:2]

        # x and y must overshoot c.shape by one for pcolormesh
        x = x_min + np.arange(N_x+1)*dx
        y = y_min + np.arange(N_y+1)*dy

        c /= dx*dy*thickness

        c_all.append(c)
        v_all.append(v)
        v_av_all.append(v_average)

plt.figure()
plt.subplot(111, aspect=1)
c = np.mean(c_all, axis=0)

c_mean = np.mean(c, axis=-1)

c_copy = c[:,:,args.species]
c_copy[c_mean < 0.1] = np.nan

plt.pcolormesh(x, y, c_copy.T, cmap=plt.cm.viridis, rasterized=True)
plt.title(r'$c_{}(x,y)$'.format('ABCD'[args.species]))
plt.colorbar()

plt.figure()
plt.subplot(111, aspect=1)

x, y = np.meshgrid(x[:-1], y[:-1])

v = np.mean(v_all, axis=0)
v = v[:,:,args.species,:]
mask = c[:,:,args.species] < 1
v[mask,:] = 0

plt.quiver(x, y, v[:,:,0].T, v[:,:,1].T, linewidths=0.1)

plt.figure()
plt.subplot(111, aspect=1)

v_average = np.mean(v_av_all, axis=0)

mask = c[:,:,0]+c[:,:,1] < 1
v_average[mask,:] = 0


plt.quiver(x[::2,::2], y[::2,::2], v_average[::2,::2,0].T, v_average[::2,::2,1].T,
           linewidths=0.05, rasterized=True, headwidth=0.1, headlength=0.1, pivot='mid')

plt.show()
