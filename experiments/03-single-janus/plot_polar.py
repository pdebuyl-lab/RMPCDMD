#!/usr/bin/env python
"""
Display the polar concentration and velocity fields of a RMPCDMD simulation.
"""
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('file', help="H5MD file")
parser.add_argument('--species', type=int, default=0)
parser.add_argument('--radius', type=float)
args = parser.parse_args()

import h5py
import matplotlib.pyplot as plt
import numpy as np

with h5py.File(args.file, 'r') as f:
    c = f['fields/polar_concentration']
    v = f['fields/polar_velocity']
    r_min = c.attrs['r_min'][()]
    r_max = c.attrs['r_max'][()]
    dr = c.attrs['dr'][()]
    dtheta = c.attrs['dtheta'][()]

    c = c[:]
    v = v[:]

    N_r, N_theta = c.shape[:2]

    # r and th must overshoot c.shape by one for pcolormesh
    r = r_min + np.arange(N_r+1)*dr
    r = r.reshape((-1,1))
    th = (np.arange(N_theta+1)+0.5)*dtheta
    th = th.reshape((1,-1))

    c[:,:,0] /= (2*np.pi*r**2*np.sin(th)*dr*dtheta)[:-1,:-1]
    c[:,:,1] /= (2*np.pi*r**2*np.sin(th)*dr*dtheta)[:-1,:-1]

    X = r*np.sin(th)
    Z = r*np.cos(th)

plt.subplot(121, aspect=1)
plt.pcolormesh(X, Z, c[:,:,args.species], cmap=plt.cm.viridis)
plt.colorbar()
plt.axhline(0)
if args.radius is not None:
    plt.plot(args.radius*np.sin(th.flatten()),
             args.radius*np.cos(th.flatten()), color='w')

plt.subplot(122, aspect=1)

r = r_min + np.arange(N_r)*dr
th = (np.arange(N_theta)+0.5)*dtheta

r, th = np.meshgrid(r, th)

x = r*np.sin(th)
y = r*np.cos(th)

vx = v[:,:,args.species,0].T*np.sin(th) + v[:,:,args.species,1].T*np.cos(th)
vy = v[:,:,args.species,0].T*np.cos(th) - v[:,:,args.species,1].T*np.sin(th)

plt.quiver(x, y, vx, vy)

plt.show()
