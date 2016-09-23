#!/usr/bin/env python
from __future__ import print_function, division

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('file', type=str, help='H5MD datafile')
parser.add_argument('--obs', type=str,
                    help='Observables to plot, e.g. \'temperature\'',
                    nargs='+')
parser.add_argument('--traj', type=str)
parser.add_argument('--field', type=str)
parser.add_argument('--com', action='store_true')
parser.add_argument('--hist', action='store_true')
parser.add_argument('--mean', action='store_true')
parser.add_argument('--index', default=0, type=int)
args = parser.parse_args()

import numpy as np
import h5py
import matplotlib.pyplot as plt

def get_edges(particles_group):
    box_group = particles_group['box']
    if isinstance(box_group['edges'],h5py.Group):
        edges = box_group['edges/value'][:]
        assert np.allclose(edges - edges[0].reshape((1,-1)), np.zeros_like(edges))
        return edges[0]
    else:
        return box_group['edges'][:]

def get_time(element):
    if 'time' in element:
        time = element['time']
    else:
        time = element['step']
    if 'offset' in time.attrs:
        offset = time.attrs['offset'][()]
    else:
        offset = 0
    if time.shape==():
        return offset + np.arange(element['value'].shape[0])*time[()]
    else:
        return time[:]

with h5py.File(args.file, 'r') as f:
    if args.obs:
        for obs in args.obs:
            g = f['observables'][obs]
            plt.plot(get_time(g), g['value'])

    elif args.traj:
        group, traj = args.traj.split('/')
        edges = get_edges(f['particles'][group])
        if args.com:
            if traj == 'position':
                pos = f['particles'][group][traj]['value']
                im = f['particles'][group]['image']['value']
                data = pos[:,:,:] + im[:,:,:]*edges.reshape((1, 1, 3))
            else:
                data = f['particles'][group][traj]['value'][:,:,:]
            data = data.mean(axis=1)
        else:
            data = f['particles'][group][traj]['value'][:,args.index,:]
        for i in range(3):
            plt.subplot(3, 1, i+1)
            if args.hist:
                plt.hist(data[:,i], bins=32)
                print('xyz'[i], 'mean', data[:,i].mean(),
                      'std', data[:,i].std())
            else:
                plt.plot(get_time(f['particles'][group][traj]), data[:,i])
                plt.ylabel('xyz'[i])

    elif args.field:
        g = f['fields'][args.field]
        value = g['value']
        if len(value.shape)==3:
            value = value[:,:,args.index]
        else:
            value = value[:]
        if args.mean:
            plt.plot(value.mean(axis=0))
        else:
            N = value.shape[0]
            plt.plot(value[::N//4].T)
        plt.title(args.field)

plt.show()
