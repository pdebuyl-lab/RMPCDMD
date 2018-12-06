#!/usr/bin/env python
"Plot the mean square displacement of the dimer's center of mass."
from __future__ import print_function, division

import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--index', type=int, help='index of colloid', default=0)
parser.add_argument('file', type=str, help='H5MD datafile', nargs='+')
parser.add_argument('--dimer', action='store_true')
parser.add_argument('--slicer', type=int, nargs=2)
parser.add_argument('--group')
args = parser.parse_args()
print(args)

import numpy as np
import matplotlib.pyplot as plt
import h5py
import tidynamics

msd_data = []
vz_data = []

if args.slicer:
    slicer = slice(args.slicer[0], None, args.slicer[1])
else:
    if args.dimer:
        slicer = slice(2)
    else:
        slicer = args.index

for filename in args.file:
    with h5py.File(filename, 'r') as f:
        particles = f['particles']
        if args.group:
            group = particles[args.group]
        else:
            group = particles[list(particles.keys())[0]]
        r = group['position/value'][:,slicer,:]
        r_dt = group['position/time'][()]
        im = group['image/value'][:,slicer,:]
        v = group['velocity/value'][:,slicer,:]
        v_dt = group['velocity/time'][()]
        edges = group['box/edges'][:].reshape((1,-1))
        r += im*edges
        assert abs(r_dt-v_dt) < 1e-12
        if args.dimer:
            assert r.shape[1]==2
            assert r.shape[2]==3
            r = r.mean(axis=1)
            v_com = v.mean(axis=1)

        if r.ndim == 3:
            for i in range(r.shape[1]):
                msd_data.append(tidynamics.msd(r[:,i,:]))
        else:
            msd_data.append(tidynamics.msd(r))

msd_data = np.array(msd_data)
m = msd_data.mean(axis=0)
s = msd_data.std(axis=0)
time = np.arange(r.shape[0])*r_dt

plt.plot(time, m, 'k-')
plt.plot(time, m+s, 'k--')
plt.plot(time, m-s, 'k--')

fit = np.polyfit(time[:r.shape[0]//4], m[:r.shape[0]//4], 1)
D = fit[0] / 6
plt.plot(time, 6*D*time, 'k:')

print("Estimated D_eff", D)

plt.show()
