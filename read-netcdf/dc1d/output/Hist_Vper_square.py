#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from scipy import stats

unpert_orbit_dir = "unperturbed/"
pert_orbit_dir = "perturbed/"

num_orbits = 1000
iVar =0

nPerRF = 1

iVarstr = "{0:03d}".format(0)
tfull = np.loadtxt(unpert_orbit_dir + "orbit_" + iVarstr + "_000.txt", skiprows = 2, usecols=(0,))
nsteps = tfull.shape[0]
t = tfull[0::nPerRF]

Vper_pert = np.zeros((nsteps/nPerRF,num_orbits))

hist_fn = np.zeros( (num_orbits,) )
hist_fn_old = np.zeros( (num_orbits,) )

iVar = 0
for iOrb in range(0,num_orbits):
    iOrbstr = "{0:03d}".format(iOrb)
    iVarstr = "{0:03d}".format(iVar)

    Vper_pert_tmp = np.loadtxt(pert_orbit_dir + "orbit_" + iVarstr + "_" + iOrbstr + ".txt", skiprows = 2, usecols=(5,))
 
    Vper_pert[:,iOrb] = np.mean(Vper_pert_tmp.reshape(-1, nPerRF), axis=1)

xlim = 2000
nBins = 50
for i in range(0,nsteps/nPerRF-1):
    hist_fn = Vper_pert[i,:] - Vper_pert[0,:]

    plt.figure()
    #plt.hist(hist_fn, nBins)
    plt.hist(hist_fn,bins=np.linspace(-xlim,xlim, nBins))
    plt.xlim(-xlim,xlim)
    plt.ylim(0,num_orbits)
    istr = "{0:03d}".format(i)
    savefig("plots_unpert_vs_pert/" + "hist" + istr + ".png")
    plt.close()
