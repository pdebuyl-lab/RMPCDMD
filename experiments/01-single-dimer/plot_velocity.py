#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('file', type=str, help='H5MD datafile')
parser.add_argument('--directed', action='store_true')
args = parser.parse_args()

import numpy as np
import h5py
import matplotlib.pyplot as plt

with h5py.File(args.file, 'r') as f:
    r = f['particles/dimer/position/value'][...]
    r_dt = f['particles/dimer/position/time'][()]
    v = f['particles/dimer/velocity/value'][...]
    v_dt = f['particles/dimer/velocity/time'][()]

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
    plt.plot(time, vz)
else:
    plt.plot(time, v_com)

plt.show()
