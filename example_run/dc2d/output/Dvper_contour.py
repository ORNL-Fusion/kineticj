#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from scipy import stats
from mpl_toolkits.mplot3d import Axes3D

unpert_orbit_dir = "unperturbed/"
pert_orbit_dir = "perturbed/"

num_orbits = 10
numVar = 5
nPerRF = 20

var = np.linspace(1.2,1.4,numVar)
iVarstr = "{0:03d}".format(0)
tfull = np.loadtxt(unpert_orbit_dir + "orbit_" + iVarstr + "_000.txt", skiprows = 2, usecols=(0,))
nsteps = tfull.shape[0]
t = tfull[0::nPerRF]

Vper_unpert_RFavg = np.zeros((nsteps/nPerRF,num_orbits))
Vper_pert_RFavg = np.zeros((nsteps/nPerRF,num_orbits))

Dvper = np.zeros( (nsteps/nPerRF -1,numVar))

plt.figure()
for iVar in range(0,numVar):
    for iOrb in range(0,num_orbits):
        iOrbstr = "{0:03d}".format(iOrb)
        iVarstr = "{0:03d}".format(iVar)

        Vper_unpert_tmp = np.loadtxt(unpert_orbit_dir + "orbit_" + iVarstr + "_" + iOrbstr + ".txt", skiprows = 2, usecols=(5,) )
        Vper_pert_tmp = np.loadtxt(pert_orbit_dir + "orbit_" + iVarstr + "_" + iOrbstr + ".txt", skiprows = 2, usecols=(5,))
        
        Vper_unpert_RFavg[:,iOrb] = np.mean(Vper_unpert_tmp.reshape(-1, nPerRF), axis=1)
        Vper_pert_RFavg[:,iOrb] = np.mean(Vper_pert_tmp.reshape(-1, nPerRF), axis=1)

    nBins = 1001
    xlim = np.max( [ np.abs(np.min(Vper_pert_RFavg - Vper_unpert_RFavg)), np.max(Vper_pert_RFavg - Vper_unpert_RFavg)])
    for it in range(1,nsteps/nPerRF):
        histVar, binVar = np.histogram(Vper_pert_RFavg[it,:] - Vper_unpert_RFavg[it,:], density=True, bins=np.linspace(-xlim,xlim,nBins))
        midPoints = binVar[:-1] + np.diff(binVar)/2
        Dvper[it-1,iVar] = np.dot(diff(binVar), histVar*(midPoints*midPoints))#/(-2*t[it])
#    print np.min(Vper_pert - Vper_unpert)
#    print np.max(Vper_pert - Vper_unpert)

plt.figure()
for i in range(0,numVar):
    plt.plot(-t[1:],Dvper[:,i], '--o')
plt.savefig("Dvper_vs_t_4_traj_1orb.png")
plt.show()

for it in range(1,nsteps/nPerRF):
    plt.plot(var, Dvper[it-1,:], label="RFcycle = " + str(it+1))
plt.legend()
plt.savefig("Dvper_vs_x.png")

fig = plt.figure()
Varmg, tmg = np.meshgrid(var, t[1:])
plt.contourf(Varmg, -tmg, log10(Dvper))
plt.plot(np.zeros( (2,1)) + 1.388, [-t[1], -t[-1]])
plt.colorbar()
plt.savefig("ContourDvper_no_t.png")
plt.show()
