#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from scipy import stats

unpert_orbit_dir = "unperturbed/"
pert_orbit_dir = "perturbed/"

num_orbits = 200
iVar =0

nPerRF = 1

tfull = np.loadtxt(unpert_orbit_dir + "orbit_" + str(0) + "_000.txt", skiprows = 2, usecols=(0,))
nsteps = tfull.shape[0]
t = tfull[0::nPerRF]

Vper_unpert = np.zeros((nsteps/nPerRF,num_orbits))
Vper_pert = np.zeros((nsteps/nPerRF,num_orbits))

hist_fn = np.zeros( (num_orbits,) )
hist_fn_old = np.zeros( (num_orbits,) )

for iOrb in range(0,num_orbits):
    iOrbstr = "{0:03d}".format(iOrb)
    Vper_unpert_tmp = np.loadtxt(unpert_orbit_dir + "orbit_" + str(iVar) + "_" + iOrbstr + ".txt", skiprows = 2,usecols=(5,))
    Vper_pert_tmp = np.loadtxt(pert_orbit_dir + "orbit_" + str(iVar) + "_" + iOrbstr + ".txt", skiprows = 2,usecols=(5,))

    Vper_unpert[:,iOrb] = np.mean(Vper_unpert_tmp.reshape(-1, nPerRF), axis=1)
    Vper_pert[:,iOrb] = np.mean(Vper_pert_tmp.reshape(-1, nPerRF), axis=1)

#MinMax = [np.min(Vper_pert-Vper_unpert)]#, np.max(Vper_pert-Vper_unpert)]
xlim = np.max(Vper_pert)

#xlim = np.max(np.abs(MinMax))

nBins = 31
for i in range(0,nsteps/nPerRF-1):
 #   a = ((i -1.0)/i)*hist_fn_old
 #   b = (1.0/i)*(Vper_pert[i,:] - Vper_unpert[i,:])
    hist_fn = Vper_pert[i,:]# - Vper_unpert[i,:]

    plt.figure()
    plt.hist(hist_fn,bins=np.linspace(-xlim,xlim, nBins))
    plt.xlim(-xlim,xlim)
    plt.ylim(0,num_orbits)
    istr = "{0:03d}".format(i)
    savefig("plots_unpert_vs_pert/" + "hist" + istr + ".png")
    plt.close()
