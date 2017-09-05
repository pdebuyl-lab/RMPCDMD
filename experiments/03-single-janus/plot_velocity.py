#!/usr/bin/env python
"""
Program to display the directed velocity of a self-propelled rigid-body colloid
along its axis.
"""

import argparse


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('file', type=str, help='H5MD datafile')
parser.add_argument('--plot', action='store_true',
                    help='Display the result graphically')
parser.add_argument('--no-q', action='store_true',
                    help='Do not use quaternion data')
args = parser.parse_args()

import numpy as np
import matplotlib.pyplot as plt
from transforms3d import quaternions
import pyh5md

with pyh5md.File(args.file, 'r') as f:
    obs = f['observables']
    u = pyh5md.element(obs, 'u').value[:]
    if 'q' in obs:
        q = pyh5md.element(obs, 'q').value[:]
    else:
        q = None
    janus_vel = pyh5md.element(obs, 'janus_vel')
    dt = janus_vel.time
    janus_vel = janus_vel.value[:]

one_z = np.array([0., 0., 1.])

# Change of quaternion storage convention
if q is not None and not args.no_q:
    tmp_s = q[:,3].copy()
    tmp_v = q[:,:3].copy()
    q[:,1:4] = tmp_v
    q[:,0] = tmp_s
    del tmp_s, tmp_v

    # Obtain u by rotation of 1_z with quaternion
    u_q = np.array([quaternions.rotate_vector(one_z, q_var) for q_var in q])

    dir_vel = np.sum(janus_vel*u_q, axis=1)
else:
    print('no q')
    dir_vel = np.sum(janus_vel*u, axis=1)

m = dir_vel.mean()

print(args.file, m)
if args.plot:
    plt.subplot(121)
    plt.title('directed velocity')
    plt.plot(np.arange(len(janus_vel))*dt, dir_vel)
    plt.subplot(122)
    plt.title('histogram of directed velocity')
    plt.hist(dir_vel, bins=32)
    plt.axvline(m, c='red')

    plt.show()
