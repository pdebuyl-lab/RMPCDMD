#!/usr/bin/env python
from __future__ import print_function, division

import argparse

description = "Plot the mean square displacement of the dimer's center of mass."
parser = argparse.ArgumentParser(description=description)
parser.add_argument('file', type=str, help='H5MD datafile', nargs='+')
args = parser.parse_args()

import numpy as np
import h5py
import matplotlib.pyplot as plt

msd_data = []
vz_data = []
for filename in args.file:
    with h5py.File(filename, 'r') as f:
        r = f['particles/dimer/position/value'][:,:,:]
        r0 = f['particles/dimer/position/value'][0,:,:]
        r_dt = f['particles/dimer/position/time'][()]
        im = f['particles/dimer/image/value'][:,:,:]
        v = f['particles/dimer/velocity/value'][...]
        v_dt = f['particles/dimer/velocity/time'][()]
        edges = f['particles/dimer/box/edges'][:].reshape((1,-1))
        assert abs(r_dt-v_dt) < 1e-12
        assert r.shape[1]==2
        assert r.shape[2]==3
        r += im*edges
        v_com = v.mean(axis=1)
        unit_z = r[:,1,:]-r[:,0,:]
        unit_z /= np.sqrt(np.sum(unit_z**2, axis=1)).reshape((-1,1))
        vz = np.sum(v_com*unit_z, axis=1)
        vz_data.append(vz.mean())
        r -= r0.reshape((1,2,3))
        msd_data.append(np.sum(np.mean(r,axis=1)**2, axis=1))

msd_data = np.array(msd_data)
m = msd_data.mean(axis=0)
s = msd_data.std(axis=0)
vz = np.mean(vz_data, axis=0)
time = np.arange(r.shape[0])*r_dt

plt.plot(time, m, 'k-')
plt.plot(time, m+s, 'k--')
plt.plot(time, m-s, 'k--')

D = np.mean(m[1:]/time[1:])/6
plt.plot(time, 6*D*time, 'k:')

print("Estimated D_eff", D)

plt.show()
