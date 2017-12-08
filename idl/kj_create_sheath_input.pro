pro kj_create_sheath_input

	@dlg_constants

	freq = 1.15d8
	e_density_m3 = 1.1d14
	T_eV = 1d0
	T_keV = T_eV * 1d-3
	n20 = e_density_m3 / 10d20
	bMag = 0.001
	w_ce = e * bMag / me

	lambda_D = 2.35d-5 * sqrt(T_keV/n20)

	nR = 1000
	rMin = 0
	rMax = 1000*lambda_D

	print, 'rMin: ', rMin
	print, 'rMax: ', rMax
	print, 'w_ce: ', w_ce
	print, 'w_rf: ', freq * 2 * !pi

	r = fIndGen(nR)*(rMax-rMin)/(nR-1)+rMin

	br = fltArr(nR) + bMag
	bt = fltArr(nR) + 0d0
	bz = fltArr(nR) + 0d0

	er_re = fltArr(nR) + 10.1
	er_im = fltArr(nR) + 0.0

	; Write netCDF file

	nc_id = nCdf_create ('kj_sheath_1d.nc', /clobber )

	nCdf_control, nc_id, /fill
	
	nr_id = nCdf_dimDef ( nc_id, 'nR', n_elements(r) )
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

	jA_r_re_id = nCdf_varDef ( nc_id, 'jA_r_re', nr_id, /float )
	jA_r_im_id = nCdf_varDef ( nc_id, 'jA_r_im', nr_id, /float )
	jA_p_re_id = nCdf_varDef ( nc_id, 'jA_p_re', nr_id, /float )
	jA_p_im_id = nCdf_varDef ( nc_id, 'jA_p_im', nr_id, /float )
	jA_z_re_id = nCdf_varDef ( nc_id, 'jA_z_re', nr_id, /float )
	jA_z_im_id = nCdf_varDef ( nc_id, 'jA_z_im', nr_id, /float )

	nCdf_control, nc_id, /enDef

	_empty = fltArr(nR)

	nCdf_varPut, nc_id, freq_id, freq

	nCdf_varPut, nc_id, r_id, r 

	nCdf_varPut, nc_id, B0_r_id, br
	nCdf_varPut, nc_id, B0_p_id, bt 
	nCdf_varPut, nc_id, B0_z_id, bz

	nCdf_varPut, nc_id, e_r_re_id, er_re
	nCdf_varPut, nc_id, e_r_im_id, er_im
	nCdf_varPut, nc_id, e_p_re_id, _empty
	nCdf_varPut, nc_id, e_p_im_id, _empty 
	nCdf_varPut, nc_id, e_z_re_id, _empty  
	nCdf_varPut, nc_id, e_z_im_id, _empty 

	nCdf_varPut, nc_id, jP_r_re_id, _empty 
	nCdf_varPut, nc_id, jP_r_im_id, _empty 
	nCdf_varPut, nc_id, jP_p_re_id, _empty 
	nCdf_varPut, nc_id, jP_p_im_id, _empty 
	nCdf_varPut, nc_id, jP_z_re_id, _empty 
	nCdf_varPut, nc_id, jP_z_im_id, _empty 

	nCdf_varPut, nc_id, jA_r_re_id, _empty 
	nCdf_varPut, nc_id, jA_r_im_id, _empty 
	nCdf_varPut, nc_id, jA_p_re_id, _empty 
	nCdf_varPut, nc_id, jA_p_im_id, _empty 
	nCdf_varPut, nc_id, jA_z_re_id, _empty 
	nCdf_varPut, nc_id, jA_z_im_id, _empty 

	nCdf_close, nc_id
end
