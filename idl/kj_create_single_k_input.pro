pro kj_create_single_k_input, $
        b0=_b0, bUnit=_bUnit, kx=_kx, f_Hz=_f_Hz, n_m3=_n_m3, $
        Er=Er, Et=Et, Ez=Ez, x=x, writeOutput=writeOutput, $
        E1Multiplier=_E1Multiplier, $
        E2Multiplier=_E2Multiplier, $
        E3Multiplier=_E3Multiplier, $
        fileName=_fileName, $
        nPts = _nPts

@constants

if keyword_set(_nPts) then nPts = _nPts else nPts = 301
if keyword_set(_b0) then b0 = _b0 else b0 = 1
if keyword_set(_bUnit) then bUnit = _bUnit else bUnit = [0,0,1]
if keyword_set(_kx) then kx = _kx else kx = 1
if keyword_set(_f_Hz) then f_Hz = _f_Hz else f_Hz = 13.56e6 
if keyword_set(_n_m3) then n_m3 = _n_m3 else n_m3 = 1e19 
if keyword_set(_E1Multiplier) then E1Multiplier = _E1Multiplier else E1Multiplier = 0 
if keyword_set(_E2Multiplier) then E2Multiplier = _E2Multiplier else E2Multiplier = 0 
if keyword_set(_E3Multiplier) then E3Multiplier = _E3Multiplier else E3Multiplier = 0 
if keyword_set(_fileName) then fileName = _fileName else fileName = 'input/input-data.nc' 

lambda = 2*!Pi/kx

xOffset = 0 
nCycles = 5 
xRange = lambda*nCycles
dx = xRange / (nPts-1)
x = fIndGen(nPts)*dx+xOffSet

print, 'Parallel k: ', kx
print, 'Parallel Wavelength: ', lambda
print, 'nPhi to use: ', 2*!pi*xOffset/lambda
print, 'xGridMin: ', x[0]
print, 'xGridMax: ', x[-1]

EmagR = 1.0
EmagT = 1.0 
EmagZ = 1.0 

Er = EmagR*exp(_II*kx*x)
Et = EmagT*exp(_II*kx*x)
Ez = EmagZ*exp(_II*kx*x)

Er = Er * E1Multiplier
Et = Et * E2Multiplier
Ez = Ez * E3Multiplier

if n_elements(b0) gt 1 then begin

    br = b0*bUnit[0]
    bt = b0*bUnit[1]
    bz = b0*bUnit[2]

endif else begin

    br = fltArr(nPts)+b0*bUnit[0]
    bt = fltArr(nPts)+b0*bUnit[1]
    bz = fltArr(nPts)+b0*bUnit[2]

endelse

Jpr = fltArr(nPts)
Jpt = fltArr(nPts)
Jpz = fltArr(nPts)

density = fltArr(nPts,1)+n_m3

; Write netCDF file

if keyword_set(writeOutput) then begin

    nc_id = nCdf_create (fileName, /clobber )
    
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
