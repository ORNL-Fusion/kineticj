#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from scipy import stats
from mpl_toolkits.mplot3d import Axes3D

unpert_orbit_dir = "unperturbed/"
pert_orbit_dir = "perturbed/"

num_orbits_full = 20
num_orbits_sample = 20
numVar = 4
nPerRF = 20

tfull = np.loadtxt(unpert_orbit_dir + "orbit_" + str(0) + "_000.txt", skiprows = 2, usecols=(0,))
nsteps = tfull.shape[0]
t = tfull[0::nPerRF]

Vper_unpert_RFavg = np.zeros((nsteps/nPerRF,num_orbits_sample))
Vper_pert_RFavg = np.zeros((nsteps/nPerRF,num_orbits_sample))

Dvper = np.zeros( (nsteps/nPerRF -1,numVar))

plt.figure()

#for iVar in range(0,numVar):
iVar = 2
for iOrb in range(0,num_orbits_sample):
    iOrbstr = "{0:03d}".format((num_orbits_full/num_orbits_sample)*iOrb)
    print iOrbstr
    Vper_unpert_tmp = np.loadtxt(unpert_orbit_dir + "orbit_" + str(iVar) + "_" + iOrbstr + ".txt", skiprows = 2, usecols=(5,) )
    Vper_pert_tmp = np.loadtxt(pert_orbit_dir + "orbit_" + str(iVar) + "_" + iOrbstr + ".txt", skiprows = 2, usecols=(5,))
    
    Vper_unpert_RFavg[:,iOrb] = np.mean(Vper_unpert_tmp.reshape(-1, nPerRF), axis=1)
    Vper_pert_RFavg[:,iOrb] = np.mean(Vper_pert_tmp.reshape(-1, nPerRF), axis=1)

    nBins = 300
    xlim = np.max( [ np.abs(np.min(Vper_pert_RFavg - Vper_unpert_RFavg)), np.max(Vper_pert_RFavg - Vper_unpert_RFavg)])
    for it in range(1,nsteps/nPerRF):
        histVar, binVar = np.histogram(Vper_pert_RFavg[it,:] - Vper_unpert_RFavg[it,:], density=True, bins=np.linspace(-xlim,xlim,nBins))
        midPoints = binVar[:-1] + np.diff(binVar)/2
        ############################
        Dvper[it-1,iVar] = np.dot(diff(binVar), histVar*(midPoints))
        ##########################
        #Dvper[it-1,iVar] = np.dot(diff(binVar), histVar*(midPoints*midPoints))

    for i in range(0,numVar):
        plt.plot(-t[1:],Dvper[:,i], '-o')

plt.savefig("Dvper_vs_t_single_orb_comp.png")
plt.show()
