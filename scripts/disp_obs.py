import numpy as np
import matplotlib.pyplot as plt
import h5py

a = h5py.File('data.h5', 'r')

all_obs = a['observables'].keys()

for o in all_obs:
    try:
        d_set = a['observables'][o]['samples']
        t_set = a['observables'][o]['time']

        if (len(d_set.shape)==1):
            plt.plot(t_set, d_set[:] - d_set[0], label=o)
    except:
        print o, 'not a good observable'

a.close()

plt.legend()
plt.show()
