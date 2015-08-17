#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pylab import *

from readGEQDSK3 import eqdskClass

eqdsk = eqdskClass()
filename = '/Users/dvk/scratch/aorsa2d/nstx_March30_2015/nstx_dfactor0_25/input/NSTX_130608A03_trxpl_0.410_plasma_state.geq'
eqdsk.readGEQDSK3(filename)

plt.contourf(eqdsk.r, eqdsk.z, eqdsk.psizr, 20)
plt.plot(eqdsk.rbbbs, eqdsk.zbbbs, 'k', linewidth = 2)
plt.plot(eqdsk.rlim, eqdsk.zlim, 'k', linewidth = 1)

nPerRF = 20
nVar = 1
nOrb = 64

var_num = 0
orb_num = 0

#for var_num in range (0,nVar):
#    for orb_num  in range(0,nOrb):
iVarstr = "{0:03d}".format(var_num)
iOrbstr = "{0:03d}".format(orb_num)

orbit = np.loadtxt("unperturbed/orbit_" + iVarstr + "_" + iOrbstr + ".txt",skiprows=2, usecols=(0,1,2,3,4,5,6))
x=orbit[:,1]
y=orbit[:,2]
z=orbit[:,3]

plt.plot(sqrt(x*x + y*y),z, 'r')
#plt.plot(x,z, 'r')

"""
orbit_pert = np.loadtxt("perturbed/orbit_" + iVarstr + "_" + iOrbstr + ".txt",skiprows=2)
xp=orbit_pert[:,1]
yp=orbit_pert[:,2]
zp=orbit_pert[:,3]
plt.plot(sqrt(xp*xp + yp*yp),zp, 'k')
"""
plt.show()
