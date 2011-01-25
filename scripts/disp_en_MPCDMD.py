#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

a = np.loadtxt('energy')

plt.subplot(311)
plt.plot(a[:,6]-a[0,6])

plt.subplot(312)
for i in range(7):
    plt.plot(a[:,i]-a[0,i])

plt.legend(['at_sol_en', 'at_at_en', 'sol_kin', 'at_kin', 'excess','total (no ex)', 'total'])

plt.subplot(313)
for i in range(2):
    plt.plot(a[:,i]-a[0,i])

plt.show()
