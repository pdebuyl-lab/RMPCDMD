#!/usr/bin/env python3
from __future__ import print_function, division

import argparse

"""
Display the number of particles of each available species as a function of time.
"""

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('file', type=str, help='H5MD datafile')
parser.add_argument('--tau', type=float, required=True,
    help='MPCD timestep (used to scale time when not present in H5MD file).')
parser.add_argument('--species', type=int, default=1,
                    help='species to fit (starts at zero)')
parser.add_argument('--rate', type=float,
                    help='Rate of the exponential')
args = parser.parse_args()

import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.optimize import leastsq

with h5py.File(args.file, 'r') as f:
    n_solvent = f['observables/n_solvent/value'][:]
    n_solvent_dt = f['observables/n_solvent/step'][()]*args.tau

time = np.arange(n_solvent.shape[0])*n_solvent_dt

def fitfunc(p, t):
    return p[0]*np.exp(-p[1]*t) + p[2]

def errfunc(p, t, y):
    return fitfunc(p, t) - y

p, success = leastsq(errfunc, [n_solvent[0,args.species], 1e-3, n_solvent[:,args.species].mean()], args=(time, n_solvent[:,args.species]))

assert success in (1, 2, 3, 4)

n_plots = n_solvent.shape[1]
for i in range(n_plots):
    plt.subplot(n_plots, 1, i+1)
    plt.plot(time, n_solvent[:,i])
    if i==args.species:
        plt.plot(time, fitfunc(p, time))
        if args.rate:
            plt.plot(time, n_solvent[0,i]*np.exp(-args.rate*time))

print("rate:", p[1])

plt.show()
