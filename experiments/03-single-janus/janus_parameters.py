#!/usr/bin/env python
from __future__ import print_function, division

import sys
import argparse
import math

"""Write a parameter file for the program single_janus_pbc of RMPCDMD"""

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--out', type=str, help='Name of output file')
parser.add_argument('--N-loop', type=int, help='Number of MPCD loops',
                    default=1024)
parser.add_argument('--N-MD', type=int, help='Number of MD loops', default=100)
parser.add_argument('-L', type=int, nargs='+', help='Length of the box',
                    default=[32])
parser.add_argument('-T', type=float, help='Temperature', default=1)
parser.add_argument('--sigma', type=float,
                    help='LJ sigma for colloid solvent interaction',
                    default=3)
parser.add_argument('--sigma-colloid', type=float,
                    help='LJ sigma for colloid colloid interaction and '
                    'for setting the mass of the colloids',
                    default=2)
parser.add_argument('--alpha', type=float, help='MPCD collision angle',
                    default=math.pi/2)
parser.add_argument('--tau', type=float, help='MPCD collision time', default=1)
parser.add_argument('--datafile', type=str, help='Datafile for janus particle',
                    default='janus_structure.h5')
parser.add_argument('--data-group', type=str,
                    help='Name of particles group in input file',
                    default='janus')
parser.add_argument('--prob', type=float,
                    help='Probability of surface reaction',
                    default=1)
parser.add_argument('--epsilon-C', type=float, help='Interaction for C bead',
                    default=(1, 1), nargs=2)
parser.add_argument('--epsilon-N', type=float, help='Interaction for N bead',
                    default=(1, 1), nargs=2)
parser.add_argument('--bulk-rate', type=float, help='Rate of bulk reaction',
                    default=0.01)
parser.add_argument('--colloid-sampling', type=int,
                    help='Interval for colloid sampling', default=50)
parser.add_argument('--reaction-radius',
                    help='reaction radius around the c.o.m. of the colloid')
parser.add_argument('--ywall', action='store_true',
                    help='enable confinement of the colloid'
                    ' in the y direction')
parser.add_argument('--ywall-bc', help='boundary condition for the fluid wall',
                    choices=['BOUNCE_BACK', 'SPECULAR', 'PERIODIC'],
                    default='PERIODIC')
parser.add_argument('--ywall-shift',
                    help='shift of the wall colloid potential',
                    type=float, default=1)
parser.add_argument('--ywall-epsilon',
                    help='magnitude of the wall colloid potential',
                    type=float, default=1)
parser.add_argument('--polar-r-max', type=float,
                    help='radius of polar histogram', default=10)

args = parser.parse_args()

if args.reaction_radius is not None:
    r_radius = args.reaction_radius
else:
    r_radius = 7.3*args.sigma/3

link_treshold = 2.7*args.sigma/3

ywall_logical = 'T' if args.ywall else 'F'

if len(args.L)==1:
    box_L = '{L} {L} {L}'.format(L=args.L[0])
elif len(args.L)==3:
    box_L = '{L[0]} {L[1]} {L[2]}'.format(L=args.L)
else:
    raise ValueError('L must be 1 or 3 integers')

output = """# physical parameters
T = {T}
L = {box_L}
rho = 10
tau = {tau}
alpha = {alpha}
probability = {prob}

# simulation parameters
N_MD = {N_MD}
N_loop = {N_loop}
colloid_sampling = {colloid_sampling}
do_solvent_io = F
equilibration_loops = 50
data_filename = {datafile}
data_group = {data_group}
reaction_radius = {r_radius}
link_treshold = {link_treshold}
do_read_links = F
polar_r_max = {polar_r_max}
bulk_rate = {bulk_rate}

# wall parameters
do_ywall = {ywall_logical}
wall_sigma = {sigma}
wall_shift = {ywall_shift}
wall_epsilon = {ywall_epsilon}
fluid_wall = {ywall_bc}


# interaction parameters
sigma_colloid = {sigma_colloid}
epsilon_colloid = 2
do_lennard_jones = F
do_elastic = F
do_rattle = F
rattle_pos_tolerance = 1d-8
rattle_vel_tolerance = 1d-8
do_quaternion = T
quaternion_treshold = 1d-13

sigma = {sigma}
epsilon_N = {epsilon_N[0]} {epsilon_N[1]}
epsilon_C = {epsilon_C[0]} {epsilon_C[1]}"""

output = output.format(r_radius=r_radius, link_treshold=link_treshold,
                       box_L=box_L,
                       ywall_logical=ywall_logical, **args.__dict__)


if args.out:
    with open(args.out, 'w') as out_f:
        print(output, file=out_f)
else:
    print(output, file=sys.stdout)
