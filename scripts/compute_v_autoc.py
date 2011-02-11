import numpy as np
import matplotlib.pyplot as plt

v = np.loadtxt('at_v')
v = v[:,:]

data = []
for dt in range(1,500):
    vf = (v[500:,:]*v[500-dt:-dt,:]).sum()/3.
    data.append(vf)

data = np.array(data)

plt.plot(data)

plt.show()
