pro kj_create_single_k_input, $
        b0=_b0, kPar=_kPar, kPer=_kPer, f_Hz=_f_Hz, n_m3=_n_m3, $
        Er=Er, Et=Et, Ez=Ez, x=x, writeOutput=writeOutput, $
        E1Multiplier=_E1Multiplier, $
        E2Multiplier=_E2Multiplier, $
        E3Multiplier=_E3Multiplier

@constants

if keyword_set(_b0) then b0 = _b0 else b0 = 1
if keyword_set(_kPar) then kPar = _kPar else kPar = 1
if keyword_set(_kPer) then kPer = _kPer else kPer = 1
if keyword_set(_f_Hz) then f_Hz = _f_Hz else f_Hz = 13.56e6 
if keyword_set(_n_m3) then n_m3 = _n_m3 else n_m3 = 1e19 
if keyword_set(_E1Multiplier) then E1Multiplier = _E1Multiplier else E1Multiplier = 0 
if keyword_set(_E2Multiplier) then E2Multiplier = _E2Multiplier else E2Multiplier = 0 
if keyword_set(_E3Multiplier) then E3Multiplier = _E3Multiplier else E3Multiplier = 0 

lambdaPar = 2*!Pi/kPar
lambdaPer = 2*!Pi/kPer

xOffset = 1e5
nPts = 301L
nCycles = 5 
xRange = lambdaPer*nCycles
dx = xRange / (nPts-1)
x = fIndGen(nPts)*dx+xOffSet

print, 'Parallel k: ', kPar
print, 'Perp k: ', kPer
print, 'Parallel Wavelength: ', lambdaPar
print, 'nPhi to use: ', 2*!pi*xOffset/lambdaPar
print, 'xGridMin: ', x[0]
print, 'xGridMax: ', x[-1]

EmagR = 1
EmagT = 0.1 
EmagZ = 2

Er = EmagR*exp(-_II*kPer*x)
Et = EmagT*exp(-_II*kPer*x)
Ez = EmagZ*exp(-_II*kPer*x)

Er = Er * E1Multiplier
Et = Et * E2Multiplier
Ez = Ez * E3Multiplier

br = fltArr(nPts)
bt = fltArr(nPts)+b0 
bz = fltArr(nPts)

Jpr = fltArr(nPts)
Jpt = fltArr(nPts)
Jpz = fltArr(nPts)

density = fltArr(nPts,1)+n_m3

;p = plot(x-xOffSet,Er)
;!null = plot(x-xOffSet,imaginary(Er),/over,color='r')
;p = plot([x,x+(x[-1]-x[0])]-xOffSet,[Er,Er])
;!null = plot([x,x+(x[-1]-x[0])]-xOffSet,imaginary([Er,Er]),/over,color='r')

; Write netCDF file

if keyword_set(writeOutput) then begin

    nc_id = nCdf_create ('data/kj_single_1d.nc', /clobber )
    
    	nCdf_control, nc_id, /fill
    	
    	nr_id = nCdf_dimDef ( nc_id, 'nR', n_elements(x) )
    	nS_id = nCdf_dimDef ( nc_id, 'nSpec', 1 )
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
    
    	b_r_rb_id = nCdf_varDef ( nc_id, 'b_r_re', nr_id, /float )
    	b_r_im_id = nCdf_varDef ( nc_id, 'b_r_im', nr_id, /float )
    	b_p_rb_id = nCdf_varDef ( nc_id, 'b_p_re', nr_id, /float )
    	b_p_im_id = nCdf_varDef ( nc_id, 'b_p_im', nr_id, /float )
    	b_z_rb_id = nCdf_varDef ( nc_id, 'b_z_re', nr_id, /float )
    	b_z_im_id = nCdf_varDef ( nc_id, 'b_z_im', nr_id, /float )
    
    	jP_r_re_id = nCdf_varDef ( nc_id, 'jP_r_re', nr_id, /float )
    	jP_r_im_id = nCdf_varDef ( nc_id, 'jP_r_im', nr_id, /float )
    	jP_p_re_id = nCdf_varDef ( nc_id, 'jP_p_re', nr_id, /float )
    	jP_p_im_id = nCdf_varDef ( nc_id, 'jP_p_im', nr_id, /float )
    	jP_z_re_id = nCdf_varDef ( nc_id, 'jP_z_re', nr_id, /float )
    	jP_z_im_id = nCdf_varDef ( nc_id, 'jP_z_im', nr_id, /float )
    
    	density_id = nCdf_varDef ( nc_id, 'density_m3', [nr_id,nS_id], /float )
    
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
    
    	nCdf_varPut, nc_id, b_r_rb_id,real_part(Er)*0 
    	nCdf_varPut, nc_id, b_r_im_id,imaginary(Er)*0 
    	nCdf_varPut, nc_id, b_p_rb_id,real_part(Et)*0
    	nCdf_varPut, nc_id, b_p_im_id,imaginary(Et)*0
    	nCdf_varPut, nc_id, b_z_rb_id,real_part(Ez)*0
    	nCdf_varPut, nc_id, b_z_im_id,imaginary(Ez)*0
    
    	nCdf_varPut, nc_id, jP_r_re_id,real_part(Jpr) 
    	nCdf_varPut, nc_id, jP_r_im_id,imaginary(Jpr) 
    	nCdf_varPut, nc_id, jP_p_re_id,real_part(Jpt) 
    	nCdf_varPut, nc_id, jP_p_im_id,imaginary(Jpt) 
    	nCdf_varPut, nc_id, jP_z_re_id,real_part(Jpz) 
    	nCdf_varPut, nc_id, jP_z_im_id,imaginary(Jpz) 
    
    	nCdf_varPut, nc_id, density_id, density 
    
    nCdf_close, nc_id

endif

end
