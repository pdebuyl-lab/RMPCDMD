#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

a = np.loadtxt('at_x')
for i in range(a.shape[1]/3):
    ax.plot(a[:,3*i],a[:,1+3*i],a[:,2+3*i])
    ax.plot(a[-2:,3*i],a[-2:,1+3*i],a[-2:,2+3*i],ls='',marker='o')

ax.set_xlim3d(0,24)
ax.set_ylim3d(0,24)
ax.set_zlim3d(0,24)

plt.show()
