#!/usr/bin/env python3
from __future__ import print_function, division
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('file', type=str, help='H5MD datafile')
parser.add_argument('--plot', action='store_true')
args = parser.parse_args()

import numpy as np
import h5py
from scipy.optimize import leastsq

with h5py.File(args.file, 'r') as a:
    tz = a['/observables/tz/value'][:]
    tz_count = a['/observables/tz_count/value'][:]
    vx = a['/observables/vx/value'][:]
    vx_count = a['/observables/vx_count/value'][:]

vx = vx[-1]
Lz = len(vx)
z = np.arange(Lz)+0.5

# Fit the parabolic profile
def fitfunc(p, t):
    return p[0] * t*(Lz-t) / (Lz/2)**2

def errfunc(p, t, y):
    return fitfunc(p, t) - y

p, success = leastsq(errfunc, [0.1,], args=(z, vx))
assert success in (1, 2, 3, 4)

print("max_speed = ", p[0])

if args.plot:
    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(z, vx, 'k-', lw=2)
    plt.plot(z, fitfunc(p, z), 'k--')
    plt.xlabel(r'$z$', fontsize=20)
    plt.ylabel(r'$v_x(z)$', fontsize=20)
    plt.title('velocity')
    plt.show()
