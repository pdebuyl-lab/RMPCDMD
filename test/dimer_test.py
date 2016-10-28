#!/usr/bin/env python
from __future__ import print_function, division

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('file', type=str, help='H5MD datafile')
args = parser.parse_args()

import numpy as np
import h5py

with h5py.File(args.file, 'r') as f:
    internal_energy = f['observables/internal_energy/value'][:]
    internal_energy -= internal_energy.mean()
    print(internal_energy.mean(), internal_energy.std())

    center_of_mass_velocity = f['observables/center_of_mass_velocity/value'][:]
    center_of_mass_velocity -= center_of_mass_velocity.mean(axis=0).reshape((1,3))
    print(np.max(np.abs(center_of_mass_velocity)))
    assert internal_energy.std() < 0.1
    assert np.max(np.abs(center_of_mass_velocity)) < 1e-15
