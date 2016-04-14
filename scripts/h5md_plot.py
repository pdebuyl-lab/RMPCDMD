#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('file', type=str, help='H5MD datafile')
parser.add_argument('--obs', type=str, help='Observables to plot, e.g. \'temperature\'', nargs='+')
parser.add_argument('--traj', type=str)
parser.add_argument('--com', action='store_true')
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
        group, traj = args.traj.split('/')
        edges = f['particles'][group]['box/edges'][:]
        if args.com:
            if traj=='position':
                data = f['particles'][group][traj]['value'][:,:,:] + f['particles'][group]['image']['value'][:,:,:]*edges.reshape((1,1,3))
            else:
                data = f['particles'][group][traj]['value'][:,:,:]
            data = data.mean(axis=1)
        else:
            data = f['particles'][group][traj]['value'][:,0,:]
        for i in range(3):
            plt.subplot(3,1,i+1)
            plt.plot(data[:,i])
            plt.ylabel('xyz'[i])

plt.show()

