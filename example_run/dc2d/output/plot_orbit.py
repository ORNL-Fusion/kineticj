#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#import sys
#sys.path.append('~/code/netcdf4-pre/')
#from netCDF4 import Dataset


orbit = np.loadtxt("orbit000.txt",skiprows=2)
x=orbit[:,1]
y=orbit[:,2]
z=orbit[:,3]

fig=plt.figure()
ax = fig.gca(projection='3d')
ax.plot(x,y,z)
#plt.show()


plt.figure()
plt.plot(orbit[:,0], orbit[:,7])
plt.show()

#ar2 = '../data/ar2_kj.nc'
#fh = Dataset(ar2, mode='r')

#r = fh.variables['r'][:]
#er = fh.variables['er'][:]

#plt.figure()
#plt.plot(r,er)
#plt.show()
