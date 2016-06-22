#!/usr/bin/env python
from __future__ import print_function, division

import sys
import argparse
import numpy as np

USAGE = """Compute typical properties for a bulk MPCD fluid."""

def v_of_eta(rho, g, Lz, eta):
    return rho*g*Lz**2/(8*eta)

def bulk_command(args):
    parser_bulk = argparse.ArgumentParser(usage=USAGE)
    parser_bulk.add_argument('--tau', type=float, required=True, help='MPCD collision time')
    parser_bulk.add_argument('-T', type=float, required=True, help='Temperature')
    parser_bulk.add_argument('--rho', type=float, required=True, help='Average number density')
    parser_bulk.add_argument('--alpha', type=float, default=np.pi/2, help='Collision angle')
    parser_bulk.add_argument('--length', type=float, help='Typical length-scale for diffusive time')
    parser_bulk.add_argument('--force', type=float, help='Force for a Poiseuille flow')
    parser_bulk.add_argument('--heigth', type=float, help='Height of the channel for a Poiseuille flow')
    args = parser_bulk.parse_args(args)
    tau = args.tau
    T = args.T
    rho = args.rho
    gamma = rho
    alpha = args.alpha
    m = 1
    a = 1
    # Kapral review Eq. 55
    eta_kin = T * tau * rho / (2*m) * \
        (5*gamma-(gamma - 1 + np.exp(-gamma))*(2 - np.cos(alpha)-np.cos(2*alpha)))/ \
        ((gamma - 1 + np.exp(-gamma))*(2 - np.cos(alpha)-np.cos(2*alpha)))
    # Kapral review Eq. 56
    eta_coll = m / (18 * a * tau) * (gamma - 1 + np.exp(-gamma))*(1-np.cos(alpha))
    eta = eta_kin + eta_coll
    print("Viscosity", eta)
    D_fluid = T*tau/(2*m) * (3*gamma/((gamma - 1 + np.exp(-gamma))*(1-np.cos(alpha))) - 1)
    print("Self-diffusion D", D_fluid)
    if args.length:
        print("Diffusion time", args.length**2/D_fluid)
    if args.force and args.heigth:
        v_max = v_of_eta(rho, args.force, args.heigth, eta)
        v_av = 2/3*v_max
        print("Flow maximum ", v_max)
        print("Flow average ", v_av)
        print("Poiseuille flow Peclet number", v_av*args.heigth/D_fluid)

if __name__=='__main__':
    bulk_command(sys.argv[1:])
