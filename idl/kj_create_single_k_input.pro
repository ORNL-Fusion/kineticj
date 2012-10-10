pro kj_create_single_k_input

@constants

f_Hz = 1.15d8
E_eV = 0.5d3
n_e = 1.1d14

w = 2*!pi*f_Hz
vTh = sqrt(2d0/3d0*E_eV*_e/me)
wpe = sqrt(n_e*_e^2/(me*e0))
kPar=sqrt((w^2-wpe^2)/(2*vTh^2))
vPhs = w/kPar
lambdaPar = 2*!Pi/kPar

print, 'Parallel k: ', kPar
print, 'Parallel Wavelength: ', lambdaPar
print, 'Parallel Phase Velocity: ', vPhs

kPer = 0
lambda=0
I0 = 1
zeta0=w/(kPar*vTh)
print, 'Zeta_0: ', zeta0

Zp = complex(0.240566,-0.0201775) ; From mathematica Zp worksheet function

K3 = 1d0 - wpe^2 * exp(-lambda) / (w*kPar*vTh) * I0 * zeta0 * Zp
II = complex(0,1)
sig33 = -(K3 - 1d0)*II*w*e0

print, 'K3: ', K3

print, 'Sigma(3,3): ', sig33

xOffset = 100
nPts = 1001
nCycles = 29 
xRange = lambdaPar*nCycles
dx = xRange / (nPts-1)
x = fIndGen(nPts)*dx+xOffSet
E0 = 5000 
E = E0*complex(cos(kPar*x),sin(kPar*x))
Jp = sig33*E

p = plot(x,E,layout=[1,2,1])
!null = plot(x,imaginary(E),/over,color='b')
p = plot(x,Jp,layout=[1,2,2],/current)
!null = plot(x,imaginary(Jp),/over,color='b')

br = fltArr(nPts)+0.01 
bt = br*0
bz = br*0

Er = E
Et = E*0
Ez = E*0

Jpr = Jp
Jpt = Jp*0
Jpz = Jp*0


; Write netCDF file

nc_id = nCdf_create ('kj_single_1d.nc', /clobber )

	nCdf_control, nc_id, /fill
	
	nr_id = nCdf_dimDef ( nc_id, 'nR', n_elements(x) )
	scalar_id = nCdf_dimDef ( nc_id, 'scalar', 1 )

	freq_id = nCdf_varDef ( nc_id, 'freq', scalar_id, /float )
	r_id = nCdf_varDef ( nc_id, 'r', nr_id, /float )

	B0_r_id = nCdf_varDef ( nc_id, 'B0_r', nr_id, /float )
	B0_p_id = nCdf_varDef ( nc_id, 'B0_p', nr_id, /float )
	B0_z_id = nCdf_varDef ( nc_id, 'B0_z', nr_id, /float )

	e_r_re_id = nCdf_varDef ( nc_id, 'e_r_re', nr_id, /float )
	e_r_im_id = nCdf_varDef ( nc_id, 'e_r_im', nr_id, /float )
	e_p_re_id = nCdf_varDef ( nc_id, 'e_p_re', nr_id, /float )
	e_p_im_id = nCdf_varDef ( nc_id, 'e_p_im', nr_id, /float )
	e_z_re_id = nCdf_varDef ( nc_id, 'e_z_re', nr_id, /float )
	e_z_im_id = nCdf_varDef ( nc_id, 'e_z_im', nr_id, /float )

	jP_r_re_id = nCdf_varDef ( nc_id, 'jP_r_re', nr_id, /float )
	jP_r_im_id = nCdf_varDef ( nc_id, 'jP_r_im', nr_id, /float )
	jP_p_re_id = nCdf_varDef ( nc_id, 'jP_p_re', nr_id, /float )
	jP_p_im_id = nCdf_varDef ( nc_id, 'jP_p_im', nr_id, /float )
	jP_z_re_id = nCdf_varDef ( nc_id, 'jP_z_re', nr_id, /float )
	jP_z_im_id = nCdf_varDef ( nc_id, 'jP_z_im', nr_id, /float )

	nCdf_control, nc_id, /enDef

	nCdf_varPut, nc_id, freq_id, f_Hz

	nCdf_varPut, nc_id, r_id, x

	nCdf_varPut, nc_id, B0_r_id, br
	nCdf_varPut, nc_id, B0_p_id, bt 
	nCdf_varPut, nc_id, B0_z_id, bz

	nCdf_varPut, nc_id, e_r_re_id,real_part(Er) 
	nCdf_varPut, nc_id, e_r_im_id,imaginary(Er) 
	nCdf_varPut, nc_id, e_p_re_id,real_part(Et)
	nCdf_varPut, nc_id, e_p_im_id,imaginary(Et)
	nCdf_varPut, nc_id, e_z_re_id,real_part(Ez)
	nCdf_varPut, nc_id, e_z_im_id,imaginary(Ez)

	nCdf_varPut, nc_id, jP_r_re_id,real_part(Jpr) 
	nCdf_varPut, nc_id, jP_r_im_id,imaginary(Jpr) 
	nCdf_varPut, nc_id, jP_p_re_id,real_part(Jpt) 
	nCdf_varPut, nc_id, jP_p_im_id,imaginary(Jpt) 
	nCdf_varPut, nc_id, jP_z_re_id,real_part(Jpz) 
	nCdf_varPut, nc_id, jP_z_im_id,imaginary(Jpz) 

nCdf_close, nc_id


stop

end
