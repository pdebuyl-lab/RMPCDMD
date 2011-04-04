import numpy as np
import matplotlib.pyplot as plt
import h5py
from sys import argv

a = h5py.File('data.h5', 'r')

all_obs = a['observables'].keys()

if (len(argv) > 1):
    o = argv[1]
    d_set = a['observables'][o]['samples']
    t_set = a['observables'][o]['time']
    plt.plot(t_set, d_set, label=o)

else:

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
