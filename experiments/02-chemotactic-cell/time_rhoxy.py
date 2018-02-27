#!/usr/bin/env python
"Plot the concentration field rho(x,y) in the course of time"
from __future__ import print_function, division

import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('file', type=str, help='H5MD datafile')
parser.add_argument('--species', type=int, default=1)
parser.add_argument('--interval', type=int, default=20)
args = parser.parse_args()

import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.animation as animation

f = plt.figure()
c = None

with h5py.File(args.file, 'r') as a: 

    rho_xy = a['fields/rho_xy/value']
    traj = a['particles/dimer/position/value'][:]
    edges = a['particles/dimer/box/edges'][:]
    traj += a['particles/dimer/image/value'][:]*edges.reshape((1,1,-1))

    x = np.arange(rho_xy.shape[1]+1)
    y = np.arange(rho_xy.shape[2]+1)
    X, Y = np.meshgrid(x, y)

    def update_plot(i):
        global c
        plt.title(str(i))
        l = plt.pcolormesh(X, Y, rho_xy[i*args.interval,:,:,args.species].T)
        if c is None:
            c = plt.colorbar()
        plt.plot(traj[:i*args.interval,0,0],traj[:i*args.interval,0,1], c='r', lw=4)
        plt.plot(traj[:i*args.interval,1,0],traj[:i*args.interval,1,1], c='k', lw=4)
        return l

    ani = animation.FuncAnimation(f, update_plot, rho_xy.shape[0]//args.interval, interval=100, repeat=False)

    plt.show()
