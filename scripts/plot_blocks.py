"""
Plot the logarithmic block-averaged data stored in the H5MD simulation file.
"""
from __future__ import print_function, division
import h5py
import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('file', type=str, help='H5MD datafile')
args = parser.parse_args()


def get_block_data(block, dt):
    t_data = []
    data = []
    for i in range(block.shape[0]):
        t = np.arange(block.shape[1])*block.shape[1]**i
        t_data.append(t[1:])
        data.append(block[i,1:,:,:].reshape((-1,3)))

    return np.concatenate(t_data), np.concatenate(data)


def read_data(group):
    value = group['value'][:]
    count = group['count'][:]
    value /= count.reshape((-1, 1, 1, 1))
    tau = group['time'][()]
    return value, count, tau


a = h5py.File(args.file, 'r')

do_msd = 'mean_square_displacement' in a['block_correlators']
if do_msd:
    g = a['block_correlators/mean_square_displacement']
    msd, msd_count, msd_tau = read_data(g)

do_oacf = 'orientation_autocorrelation' in a['block_correlators']
if do_oacf:
    g = a['block_correlators/orientation_autocorrelation']
    oacf, oacf_count, oacf_tau = read_data(g)

do_omega = 'omega_body_autocorrelation' in a['block_correlators']
if do_omega:
    g = a['block_correlators/omega_body_autocorrelation']
    omega, omega_count, omega_tau = read_data(g)

do_vacf = 'velocity_autocorrelation' in a['block_correlators']
if do_vacf:
    g = a['block_correlators/velocity_autocorrelation']
    vacf, vacf_count, vacf_tau = read_data(g)

do_pvacf = 'parallel_velocity_autocorrelation' in a['block_correlators']
if do_pvacf:
    g = a['block_correlators/parallel_velocity_autocorrelation']
    pvacf, pvacf_count, pvacf_tau = read_data(g)

do_tvacf = 'transverse_velocity_autocorrelation' in a['block_correlators']
if do_tvacf:
    g = a['block_correlators/transverse_velocity_autocorrelation']
    tvacf, tvacf_count, tvacf_tau = read_data(g)

if do_msd:
    plt.figure()
    plt.title('msd')

    t_data, msd_data = get_block_data(msd, msd_tau)
    plt.plot(t_data, msd_data, marker='o')
    plt.loglog()

if do_oacf:
    plt.figure()
    plt.title('oacf')

    t_data, oacf_data = get_block_data(oacf, oacf_tau)
    oacf_data = np.sum(oacf_data, axis=-1)
    plt.plot(t_data, oacf_data, marker='o')
    plt.xscale('log')

if do_omega:
    plt.figure()
    plt.title('omega body')

    t_data, omega_data = get_block_data(omega, omega_tau)
    plt.plot(t_data, omega_data, marker='o')
    plt.xscale('log')

if do_vacf or do_pvacf or do_tvacf:
    plt.figure()

if do_vacf:
    plt.subplot(311)
    plt.title('vacf')

    t_data, vacf_data = get_block_data(vacf, vacf_tau)
    plt.plot(t_data, vacf_data, marker='o')
    plt.xscale('log')

if do_pvacf:
    plt.subplot(312)
    plt.title('parallel vacf')

    t_data, pvacf_data = get_block_data(pvacf, pvacf_tau)
    plt.plot(t_data, pvacf_data, marker='o')
    plt.xscale('log')

if do_tvacf:
    plt.subplot(313)
    plt.title('transverse vacf')

    t_data, tvacf_data = get_block_data(tvacf, tvacf_tau)
    plt.plot(t_data, tvacf_data, marker='o')
    plt.xscale('log')

plt.show()
