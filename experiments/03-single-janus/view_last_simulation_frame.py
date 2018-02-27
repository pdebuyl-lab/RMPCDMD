#!/usr/bin/env python
# Author: Pierre de Buyl http://pdebuyl.be/
# License: BSD 3-clause
"""Display the last frame of a dimer simulation, optionally with solvent
species B
"""
from __future__ import print_function, division

import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('file', type=str, help='H5MD datafile')
parser.add_argument('--show-B', action='store_true',
                    help='show solvent species B')
parser.add_argument('--unwrap', action='store_true', help='unwrap coordinates')
parser.add_argument('--skip', type=int, help='skip initial points', default=0)
args = parser.parse_args()

if args.show_B and args.unwrap:
    raise ValueError('Cannot use show-B and unwrap options together')

import numpy as np
import h5py
from mayavi import mlab

with h5py.File(args.file, 'r') as f:
    edges = f['particles/janus/box/edges'][:]
    so_pos = f['particles/solvent/position/value'][-1]
    so_species = f['particles/solvent/species/value'][-1]
    mask = so_species==2
    pos = so_pos[mask]
    mask = so_species==1
    other_pos = so_pos[mask]
    del so_pos
    del so_species
    del mask
    dimer_pos = f['particles/janus/position/value'][-1]
    dimer_all_pos = f['particles/janus/position/value'][args.skip:]
    dimer_all_im = f['particles/janus/image/value'][args.skip:]
    sigma_C = f['parameters/sigma'][()]
    sigma_N = sigma_C
    dimer_species = f['particles/janus/species'][:]

if args.unwrap:
    dimer_all_pos += dimer_all_im * edges.reshape((1,-1))

sigma = np.ones(dimer_species.shape[0], dtype=float)
sigma[dimer_species==1] = sigma_C
sigma[dimer_species==2] = sigma_N

mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(600, 600))

if args.show_B:
    mlab.points3d(pos[:,0], pos[:,1], pos[:,2], scale_mode='none',
                  scale_factor=0.2, color=(0.9, 0.6, 0))

# scaled value of sigma such that the norm of the vector (scaled_sigma,
# scaled_sigma, scaled_sigma) is 2*sigma, that is the diameter we wish for the
# sphere glyph
scaled_sigma = 2*sigma/np.sqrt(3)

# From the mayavi documentation, this is how to size and color spheres
# independently, here size is given by 2*sigma (diameter) and color by
# -dimer_species, so that the catalytic bead (species=1) is red and the
# non-catalytic one is blue (species=2) with the default colormap in mayavi
dimer_plot = mlab.quiver3d(dimer_all_pos[-1,:,0], dimer_all_pos[-1,:,1],
                           dimer_all_pos[-1,:,2], scaled_sigma, scaled_sigma,
                           scaled_sigma, scalars=-dimer_species, mode='sphere',
                           scale_factor=1, scale_mode='vector', resolution=30)
dimer_plot.glyph.color_mode = 'color_by_scalar'
dimer_plot.glyph.glyph_source.glyph_source.center = [0, 0, 0]

mlab.axes(extent=[0, edges[0], 0, edges[1], 0, edges[2]])

if args.unwrap:
    mlab.plot3d(dimer_all_pos[:,0,0], dimer_all_pos[:,0,1],
                dimer_all_pos[:,0,2], line_width=10, color=(0.8, 0.4, 0))

mlab.show()
