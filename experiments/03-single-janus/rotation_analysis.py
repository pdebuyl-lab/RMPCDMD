#!/usr/bin/env python
"""
Analyze the rotational motion of a L-shaped colloid.
"""
from __future__ import print_function, division

import numpy as np
import h5py
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('file')
parser.add_argument('--plot', action='store_true')
parser.add_argument('--nx', type=int, default=5,
                    help='number of beads along the bottom of the L')
parser.add_argument('--ny', type=int, default=8,
                    help='number of beads along the long arm of the L')
parser.add_argument('--n-tiles', type=int, default=2,
                    help='number of tiles for the thickness of the bottom arm')
parser.add_argument('--arm-width', type=int, default=2,
                    help='width of long arm')
args = parser.parse_args()

with h5py.File(args.file, 'r') as a:
    edges = a['particles/janus/box/edges'][:]
    pos = a['particles/janus/position/value'][:]
    im = a['particles/janus/image/value'][:]
    vel = a['particles/janus/velocity/value'][:]
    pos_com = a['observables/janus_pos/value'][:]

vel_com = vel.mean(axis=1)

pos = pos + im*edges.reshape((1,1,-1))

n_planar = args.nx * args.n_tiles + args.arm_width * (args.ny - args.n_tiles)

print(n_planar)

one_z = pos[:,n_planar,:] - pos[:,0,:]
one_z = one_z / np.sqrt(np.sum(one_z**2, axis=1)).reshape((-1,1))


v12 = vel[:,n_planar-args.arm_width,:] - vel[:,0,:]

v12_inplane = v12 - np.sum(v12*one_z, axis=1).reshape((-1,1))*one_z

off_in = np.sum(v12_inplane*one_z, axis=1)

r12 = pos[:,n_planar-args.arm_width,:] - pos[:,0,:]
dist12 = np.sqrt(np.sum((pos[0,n_planar-args.arm_width,:]-pos[0,0,:])**2))
r12 /= dist12

one_y = np.cross(one_z, r12)

omega_z = np.sum(v12_inplane*one_y, axis=1)/dist12

dir_v = np.sum(r12*vel_com, axis=1)

oz_mean = np.mean(omega_z)

print('mean omega_z', oz_mean)
print('mean directed velocity', dir_v.mean())
print('rotation radius', dir_v.mean()/oz_mean)

if args.plot:
    import matplotlib.pyplot as plt

    plt.figure()

    plt.subplot(221)
    plt.ylabel('directed velocity')
    plt.plot(dir_v)
    plt.subplot(222)
    plt.hist(dir_v, bins=32, normed=True)
    plt.axvline(dir_v.mean(), color='red')

    plt.subplot(223)
    plt.plot(omega_z)
    plt.ylabel('rotational velocity')
    plt.subplot(224)
    plt.hist(omega_z, bins=32, normed=True)
    plt.axvline(oz_mean, color='red')

    plt.figure()

    x, y, z = pos_com.T

    plt.plot(x, y)
    plt.xlabel('x')
    plt.ylabel('y')

    plt.show()
