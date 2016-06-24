#!/usr/bin/env python
from __future__ import print_function, division

import argparse

description = "Plot the histogram of B particles in a cylindrical shell around the dimer."
parser = argparse.ArgumentParser(description=description)
parser.add_argument('file', type=str, help='H5MD datafile')
args = parser.parse_args()

import numpy as np
import h5py
import matplotlib.pyplot as plt

with h5py.File(args.file) as a:
    r_min = a['fields/cylindrical_shell_histogram'].attrs['r_min'][()]
    r_max = a['fields/cylindrical_shell_histogram'].attrs['r_max'][()]
    xmin = a['fields/cylindrical_shell_histogram'].attrs['xmin'][()]
    dx = a['fields/cylindrical_shell_histogram'].attrs['dx'][()]
    csh = a['fields/cylindrical_shell_histogram/value'][:]
    n_steps = csh.shape[0]
    z = xmin + np.arange(csh.shape[1])*dx

plt.plot(z, csh[n_steps//2:].mean(axis=0) / (np.pi*(r_max**2-r_min**2)*dx) )
plt.xlabel(r'$z$', fontsize=20)
plt.ylabel(r'$c_B(z)$', fontsize=20)
plt.show()
