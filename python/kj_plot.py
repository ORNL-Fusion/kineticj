from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

ncIn = Dataset('input/input-data.nc')

f = 13e6
e = 1.60217662e-19
q = 1.0*e 
m = 1.6726219e-27

xIn = ncIn.variables['r'][:]
brIn = ncIn.variables['B0_r'][:]
btIn = ncIn.variables['B0_p'][:]
bzIn = ncIn.variables['B0_z'][:]
erReIn = ncIn.variables['e_r_re'][:]
erImIn = ncIn.variables['e_r_im'][:]

B = np.sqrt(np.power(brIn,2)+np.power(btIn,2)+np.power(bzIn,2))
wc = q*B/m
w = 2*np.pi*f
w_wc = w/wc
mod = np.mod(w_wc,1.0)
kPar = 0
T_eV = 5e3
vth = np.sqrt(2*T_eV*e/m)


nc = Dataset('output/jP2.nc')

x = nc.variables['x'][:]

jx_re = nc.variables['j1xc_re'][:]
jx_im = nc.variables['j1xc_im'][:]
jy_re = nc.variables['j1yc_re'][:]
jy_im = nc.variables['j1yc_im'][:]
jz_re = nc.variables['j1zc_re'][:]
jz_im = nc.variables['j1zc_im'][:]

#print jx_re

xRng = [1,5]

fig, (ax0,ax1,ax2,ax3,ax4) = plt.subplots(nrows=5,figsize=(6, 10))

ax0.plot(x,jx_re)
ax0.plot(x,jx_im)
ax0.set_title('Jx')
ax0.set_xlim(xRng)

ax1.plot(x,jy_re)
ax1.plot(x,jy_im)
ax1.set_title('Jy')
ax1.set_xlim(xRng)

ax2.plot(x,jz_re)
ax2.plot(x,jz_im)
ax2.set_title('Jz')
ax2.set_xlim(xRng)

res = w-kPar*vth-1*wc
ax3.semilogy(xIn,1/np.abs(res))
res = w-kPar*vth-2*wc
ax3.semilogy(xIn,1/np.abs(res))
ax3.set_xlim(xRng)
ax3.set_title('1/(|w-vth*kpar-n*wc|) for n=1,2')

ax4.plot(xIn,erReIn)
ax4.plot(xIn,erImIn)
ax4.set_xlim(xRng)
ax4.set_title('Ex(Re,Im) with kx = pi/(2*rhoL) for rhoL@2nd harmonic')

fig.subplots_adjust(hspace=0.8)
plt.show()

