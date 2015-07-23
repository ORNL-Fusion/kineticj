#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pylab import *

#ic = np.loadtxt("particles.txt")
nPerRF = 20
nVar = 1
nOrb = 10

for var_num in range (0,nVar):
    #for orb_num in [0]:
    for orb_num  in range(0,nOrb):
        #var_num = 1
        #orb_num = 7
        iVarstr = "{0:03d}".format(var_num)
        iOrbstr = "{0:03d}".format(orb_num)

        orbit = np.loadtxt("unperturbed/orbit_" + iVarstr + "_" + iOrbstr + ".txt",skiprows=2, usecols=(0,1,2,3,4,5,6))
        x=orbit[:,1]
        y=orbit[:,2]
        z=orbit[:,3]

        tfull = orbit[:,0]
        t = tfull[0::nPerRF]

        orbit_pert = np.loadtxt("perturbed/orbit_" + iVarstr + "_" + iOrbstr + ".txt",skiprows=2)
        xp=orbit_pert[:,1]
        yp=orbit_pert[:,2]
        zp=orbit_pert[:,3]

        fig=plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot(x,y,z)
        ax.plot(xp,yp,zp)
        ax.set_zlim3d(-0.0008, 0.0008)
        savefig("pert_vs_unpert_orbit" + iVarstr  + "_" + iOrbstr + ".png")
        plt.close()
        
        run_int = np.zeros(len(tfull))
        dt = diff(tfull)[0]
        run_int[0] = 0.0*orbit[0,5]
        for i in range(1, len(tfull)):
            run_int[i] = run_int[i-1] + 1.0*(orbit_pert[i,5] - orbit_pert[i-1,5])

        vperp_diff = orbit_pert[:,5]# - orbit[:,5]
        plt.figure()
        #plt.plot(-tfull, run_int, '--o')
        #plt.plot(-t, np.mean(vperp_diff.reshape(-1, nPerRF), axis=1), '--o')
        plt.plot(-tfull, vperp_diff)
        plt.ylim( [10000, 400000] )
        savefig("vell_diff" + iVarstr + "_" + iOrbstr + ".png")
        plt.close()
plt.show()

#plt.figure()
#plt.plot(-orbit_pert[0:19,0], orbit[0:19,9])
#plt.show()

#plt.figure()
#plt.plot(orbit[:,0], orbit[:,7])
#plt.show()

#ar2 = '../data/ar2_kj.nc'
#fh = Dataset(ar2, mode='r')

#r = fh.variables['r'][:]
#er = fh.variables['er'][:]

#plt.figure()
#plt.plot(r,er)
#plt.show()
