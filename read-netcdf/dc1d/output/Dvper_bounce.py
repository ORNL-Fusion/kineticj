#!/usr/bin/env python

### compute change in DV after 1-D particle in periodic domain
### returns where it started

import numpy as np
import matplotlib.pyplot as plt

nR = 25
nPhase = 10

rminD = 0.92
rmaxD = 1.72
dr = rmaxD - rminD

rGridMin = 0.925
rGridMax = 1.7
r = np.linspace(rGridMin, rGridMax, nR)

Dvperp = np.zeros( (nR,nPhase))

orb_num = 0
for iR in range(0, nR):
    for iPhase in range(0, nPhase):
    
        var_num = iPhase*nR + iR
        
        iVarstr = "{0:03d}".format(var_num)
        iOrbstr = "{0:03d}".format(orb_num)

        orbit = np.loadtxt("perturbed/orbit_" + iVarstr + "_" + iOrbstr + ".txt",skiprows=2, usecols=(0,1,2,3,4,5,6))

        ### detect when particle returns in r
        iCross = np.argmin( abs(orbit[:,1] - orbit[0,1] + dr))
        print orbit[iCross,1] - orbit[0,1]
        Dvperp[iR, iPhase] = ( (orbit[iCross,5] - orbit[0,5])**2)/(-2.0*(orbit[iCross,0] - orbit[0,0]))

plt.figure()
plt.plot(r, np.mean(Dvperp[:,], axis=1), '--o')
plt.show()
