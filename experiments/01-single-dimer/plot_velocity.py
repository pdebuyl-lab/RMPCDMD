#!/usr/bin/env python
from __future__ import print_function, division

import argparse

description = "Plot the velocity, in x,y,z cartesian coordinates or along the dimer axis."
parser = argparse.ArgumentParser(description=description)
parser.add_argument('file', type=str, help='H5MD datafile')
parser.add_argument('--directed', action='store_true', help="project the velocity along the dimer's axis")
parser.add_argument('--histogram', action='store_true', help="display an histogram of the data")
args = parser.parse_args()

import numpy as np
import h5py
import matplotlib.pyplot as plt

with h5py.File(args.file, 'r') as f:
    r = f['particles/dimer/position/value'][...]
    r_dt = f['particles/dimer/position/time'][()]
    im = f['particles/dimer/image/value'][...]
    v = f['particles/dimer/velocity/value'][...]
    v_dt = f['particles/dimer/velocity/time'][()]
    edges = f['particles/dimer/box/edges'][:].reshape((1,-1))

r += edges*im

assert abs(r_dt-v_dt) < 1e-12
assert r.shape[1]==2
assert r.shape[2]==3
assert v.shape[1]==2
assert v.shape[2]==3

time = np.arange(r.shape[0])*r_dt

v_com = v.mean(axis=1)

if args.directed:
    unit_z = r[:,1,:]-r[:,0,:]
    unit_z /= np.sqrt(np.sum(unit_z**2, axis=1)).reshape((-1,1))
    vz = np.sum(v_com*unit_z, axis=1)
    if args.histogram:
        plt.hist(vz, bins=20)
    else:
        plt.plot(time, vz)
else:
    for i in range(3):
        plt.subplot(3,1,i+1)
        if args.histogram:
            plt.hist(v_com[:,i])
            plt.ylabel(r'$P(v_'+'xyz'[i]+')$')
        else:
            plt.plot(time, v_com[:,i])
            plt.ylabel('xyz'[i])

plt.show()
