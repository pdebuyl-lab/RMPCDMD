#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
fig = plt.figure()
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)

a = np.loadtxt('at_x')
for i in range(a.shape[1]/3):
    ax1.plot(a[:,3*i])
    ax2.plot(a[:,1+3*i])
    ax3.plot(a[:,2+3*i])

plt.show()
