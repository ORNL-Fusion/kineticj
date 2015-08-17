#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pylab import *

ic = np.loadtxt("particles.txt")

num_orbits_full = 20
num_orbits_sample = 5
nPerRF = 1
iVar = 1

fig=plt.figure()
for iOrb in range(0,num_orbits_sample):
    iOrbstr = "{0:03d}".format((num_orbits_full/num_orbits_sample)*iOrb)
    orbit = np.loadtxt("unperturbed/orbit_"  + str(iVar) + "_" + iOrbstr + ".txt",skiprows=2)
    x=orbit[:,1]
    y=orbit[:,2]
    z=orbit[:,3]

    orbit_pert = np.loadtxt("perturbed/orbit_"  + str(iVar) + "_" + iOrbstr + ".txt",skiprows=2)
    xp=orbit_pert[:,1]
    yp=orbit_pert[:,2]
    zp=orbit_pert[:,3]

    ax = fig.gca(projection='3d')
    ax.plot(x,y,z, 'k')
    ax.plot([x[0],x[0]],[y[0],y[0]], [z[0],z[0]], 'ok')
    ax.plot(xp,yp,zp)

ax.plot(ic[:,0], ic[:,1], ic[:,2], '.r')

savefig("orbits_on_larmor_circle.png")
plt.show()
