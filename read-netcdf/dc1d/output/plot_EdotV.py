#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pylab import *

var_num = 3
num_points = 19

orbit = np.loadtxt("unperturbed/orbit_" + str(var_num) + "_000.txt",skiprows=2)
x=orbit[:,1]
y=orbit[:,2]
z=orbit[:,3]

vx = orbit[:,6]
vy = orbit[:,7]
vz = orbit[:,8]

E_x_re = orbit[:,9]
E_x_im = orbit[:,10]

E_y_re = orbit[:,11]
E_y_im = orbit[:,12]

E_z_re = orbit[:,13]
E_z_im = orbit[:,14]

Ex_a = np.sqrt(E_x_re*E_x_re + E_x_im*E_x_im)
Ey_a = np.sqrt(E_y_re*E_y_re + E_y_im*E_y_im)
Ez_a = np.sqrt(E_z_re*E_z_re + E_z_im*E_z_im)

Ex_ph = np.arctan2(E_x_im,E_x_re)
Ey_ph = np.arctan2(E_y_im,E_y_re)
Ez_ph = np.arctan2(E_z_im,E_z_re)

EdotV = Ex_a*Ex_ph*vx +  Ey_a*Ey_ph*vy +  Ez_a*Ez_ph*vz

plt.figure()
plt.plot(-orbit[0:num_points,0], EdotV[0:num_points])

orbit_pert = np.loadtxt("perturbed/orbit_" + str(var_num) + "_000.txt",skiprows=2)
xp=orbit_pert[:,1]
yp=orbit_pert[:,2]
zp=orbit_pert[:,3]

vx_pert = orbit_pert[:,6]
vy_pert = orbit_pert[:,7]
vz_pert = orbit_pert[:,8]

E_x_re_pert = orbit_pert[:,9]
E_x_im_pert = orbit_pert[:,10]

E_y_re_pert = orbit_pert[:,11]
E_y_im_pert = orbit_pert[:,12]

E_z_re_pert = orbit_pert[:,13]
E_z_im_pert = orbit_pert[:,14]

Ex_a_pert = np.sqrt(E_x_re_pert*E_x_re_pert + E_x_im_pert*E_x_im_pert)
Ey_a_pert = np.sqrt(E_y_re_pert*E_y_re_pert + E_y_im_pert*E_y_im_pert)
Ez_a_pert = np.sqrt(E_z_re_pert*E_z_re_pert + E_z_im_pert*E_z_im_pert)

#Ex_ph_pert = np.arctan2(E_x_im_pert,E_x_re_pert)
#Ey_ph_pert = np.arctan2(E_y_im_pert,E_y_re_pert)
#Ez_ph_pert = np.arctan2(E_z_im_pert,E_z_re_pert)

Ex_ph_pert = np.arctan2(E_x_re_pert,E_x_im_pert)
Ey_ph_pert = np.arctan2(E_y_re_pert,E_y_im_pert)
Ez_ph_pert = np.arctan2(E_z_re_pert,E_z_im_pert)

EdotV_pert = Ex_a_pert*np.cos(Ex_ph_pert)*vx_pert +  Ey_a_pert*np.cos(Ey_ph_pert)*vy_pert +  Ez_a_pert*np.cos(Ez_ph_pert)*vz_pert

plt.figure()
plt.plot(-orbit[0:num_points,0], EdotV_pert[0:num_points])


fig=plt.figure()
ax = fig.gca(projection='3d')
ax.plot(x,y,z)
ax.plot(xp,yp,zp)
savefig("pert_vs_unpert_orbit.png")

plt.show()
