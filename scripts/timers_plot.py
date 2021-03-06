#!/usr/bin/env python
from __future__ import print_function

import argparse

parser = argparse.ArgumentParser('')
parser.add_argument('file', type=str, help='H5MD datafile')
parser.add_argument('--plot', action='store_true', help='display the timers as a bar graph')
args = parser.parse_args()

import numpy as np
import h5py

timers_data = None
with h5py.File(args.file, 'r') as f:
    if 'timers' not in f:
        raise Exception('No timers group found')
    timers_group = f['timers']
    timers_names = timers_group.keys()
    timers = [(name, timers_group[name][()]) for name in timers_names]

timers.sort(key=lambda x: x[1])
if timers:
    if args.plot:
        import matplotlib.pyplot as plt
        y_pos = np.arange(len(timers))
        plt.barh(y_pos, [t[1] for t in timers], align='center')
        plt.yticks(y_pos, [t[0] for t in timers])
        plt.show()
    else:
        for name, value in timers:
            print(name, value)
