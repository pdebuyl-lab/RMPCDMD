#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('file', type=str, help='H5MD datafile')
parser.add_argument('--obs', type=str, help='Observables to plot, e.g. \'temperature\'', nargs='+')
parser.add_argument('--traj', type=str)
args = parser.parse_args()

import numpy as np
import h5py
import matplotlib.pyplot as plt

with h5py.File(args.file, 'r') as f:
    if args.obs:
        for obs in args.obs:
            g = f['observables'][obs]
            if 'time' in g:
                time = g['time']
            else:
                time = g['step']
            if len(time.shape)==1:
                plt.plot(time, g['value'])
            else:
                plt.plot(np.arange(g['value'].shape[0])*time, g['value'])

    if args.traj:
        g = f['particles'][args.traj]
        plt.plot(g['value'][:,0,:])

plt.show()

