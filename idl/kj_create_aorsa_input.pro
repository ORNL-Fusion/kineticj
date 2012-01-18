pro kj_create_aorsa_input

	solutionFileList = file_search ( 'solution*.nc' )
	dataFileList = file_search ( 'runData*.nc' )

	listNo = 0

	cdfId = ncdf_open ( dataFileList[listNo], /noWrite ) 
		nCdf_varGet, cdfId, 'freq', freq 
		nCdf_varGet, cdfId, 'capR', x 
		nCdf_varGet, cdfId, 'bmod', bmod
		nCdf_varGet, cdfId, 'brU', brU
		nCdf_varGet, cdfId, 'btU', btU
		nCdf_varGet, cdfId, 'bzU', bzU
	ncdf_close, cdfId

	cdfId = ncdf_open ( solutionFileList[listNo], /noWrite ) 

		nCdf_varGet, cdfId, 'er_re', er_re 
		nCdf_varGet, cdfId, 'et_re', et_re 
		nCdf_varGet, cdfId, 'ez_re', ez_re 
		nCdf_varGet, cdfId, 'er_im', er_im 
		nCdf_varGet, cdfId, 'et_im', et_im 
		nCdf_varGet, cdfId, 'ez_im', ez_im 

		nCdf_varGet, cdfId, 'jP_r_re', jP_r_re 
		nCdf_varGet, cdfId, 'jP_t_re', jP_t_re 
		nCdf_varGet, cdfId, 'jP_z_re', jP_z_re 
		nCdf_varGet, cdfId, 'jP_r_im', jP_r_im 
		nCdf_varGet, cdfId, 'jP_t_im', jP_t_im 
		nCdf_varGet, cdfId, 'jP_z_im', jP_z_im 

	ncdf_close, cdfId


	; Write netCDF file

	nc_id = nCdf_create ('kj_aorsa_1d.nc', /clobber )

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

	j_r_re_id = nCdf_varDef ( nc_id, 'j_r_re', nr_id, /float )
	j_r_im_id = nCdf_varDef ( nc_id, 'j_r_im', nr_id, /float )
	j_p_re_id = nCdf_varDef ( nc_id, 'j_p_re', nr_id, /float )
	j_p_im_id = nCdf_varDef ( nc_id, 'j_p_im', nr_id, /float )
	j_z_re_id = nCdf_varDef ( nc_id, 'j_z_re', nr_id, /float )
	j_z_im_id = nCdf_varDef ( nc_id, 'j_z_im', nr_id, /float )

	nCdf_control, nc_id, /enDef

	nCdf_varPut, nc_id, freq_id, freq

	nCdf_varPut, nc_id, r_id, x 

	nCdf_varPut, nc_id, B0_r_id, brU*bMod
	nCdf_varPut, nc_id, B0_p_id, btU*bMod 
	nCdf_varPut, nc_id, B0_z_id, bzU*bMod

	nCdf_varPut, nc_id, e_r_re_id, er_re
	nCdf_varPut, nc_id, e_r_im_id, er_im
	nCdf_varPut, nc_id, e_p_re_id, et_re
	nCdf_varPut, nc_id, e_p_im_id, et_im
	nCdf_varPut, nc_id, e_z_re_id, ez_re 
	nCdf_varPut, nc_id, e_z_im_id, ez_im

	nCdf_varPut, nc_id, j_r_re_id, (total(jP_r_re,3))[*,0] 
	nCdf_varPut, nc_id, j_r_im_id, (total(jP_r_im,3))[*,0] 
	nCdf_varPut, nc_id, j_p_re_id, (total(jP_t_re,3))[*,0] 
	nCdf_varPut, nc_id, j_p_im_id, (total(jP_t_im,3))[*,0] 
	nCdf_varPut, nc_id, j_z_re_id, (total(jP_z_re,3))[*,0] 
	nCdf_varPut, nc_id, j_z_im_id, (total(jP_z_im,3))[*,0] 

	nCdf_close, nc_id

end
