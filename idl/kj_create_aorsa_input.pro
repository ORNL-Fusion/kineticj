pro kj_create_aorsa_input;, rsfwcPath=rsfwcPath

;; rsfwcPath is the path to the rsfwc output file to which we will match
;; the grid to, so the aorsa jP can be used as input to rsfwc.
;
;	cdfId = ncdf_open(rsfwcPath)
;
;		ncdf_varget, cdfId, 'freq', freq 
;		ncdf_varget, cdfId, 'r', r 
;		ncdf_varget, cdfId, 'r_', r_
;
;	ncdf_close, cdfId
;
	solutionFileList = file_search ( 'solution*.nc' )
	dataFileList = file_search ( 'runData*.nc' )

	listNo = 0

	cdfId = ncdf_open ( dataFileList[listNo], /noWrite ) 

		nCdf_varGet, cdfId, 'freq', freq 
		nCdf_varGet, cdfId, 'capR', a_r 
		nCdf_varGet, cdfId, 'bmod', bmod
		nCdf_varGet, cdfId, 'brU', brU
		nCdf_varGet, cdfId, 'btU', btU
		nCdf_varGet, cdfId, 'bzU', bzU

		nCdf_varGet, cdfId, 'jr_re', jA_r_re 
		nCdf_varGet, cdfId, 'jt_re', jA_t_re 
		nCdf_varGet, cdfId, 'jz_re', jA_z_re 
		nCdf_varGet, cdfId, 'jr_im', jA_r_im 
		nCdf_varGet, cdfId, 'jt_im', jA_t_im 
		nCdf_varGet, cdfId, 'jz_im', jA_z_im 

	ncdf_close, cdfId

	cdfId = ncdf_open ( solutionFileList[listNo], /noWrite ) 

		nCdf_varGet, cdfId, 'er_re', er_re 
		nCdf_varGet, cdfId, 'et_re', et_re 
		nCdf_varGet, cdfId, 'ez_re', ez_re 
		nCdf_varGet, cdfId, 'er_im', er_im 
		nCdf_varGet, cdfId, 'et_im', et_im 
		nCdf_varGet, cdfId, 'ez_im', ez_im 

		nCdf_varGet, cdfId, 'jP_r_re', a_jP_r_re 
		nCdf_varGet, cdfId, 'jP_t_re', a_jP_t_re 
		nCdf_varGet, cdfId, 'jP_z_re', a_jP_z_re 
		nCdf_varGet, cdfId, 'jP_r_im', a_jP_r_im 
		nCdf_varGet, cdfId, 'jP_t_im', a_jP_t_im 
		nCdf_varGet, cdfId, 'jP_z_im', a_jP_z_im 

	ncdf_close, cdfId

	;; Optionally remove the cold antenna piece of the e field.
	;remove_cold = 1
	;	if remove_cold then begin
	;	print, 'REMOVING COLD PIECE'
	;	cdfId = ncdf_open ( '../langmuir_cold/solution001.nc', /noWrite ) 
	;		nCdf_varGet, cdfId, 'er_re', er_re_c 
	;		nCdf_varGet, cdfId, 'et_re', et_re_c 
	;		nCdf_varGet, cdfId, 'ez_re', ez_re_c 
	;		nCdf_varGet, cdfId, 'er_im', er_im_c 
	;		nCdf_varGet, cdfId, 'et_im', et_im_c 
	;		nCdf_varGet, cdfId, 'ez_im', ez_im_c 
	;	ncdf_close, cdfId
	;	er_re =  er_re - er_re_c 
	;	et_re =  et_re - et_re_c 
	;	ez_re =  ez_re - ez_re_c 
	;	er_im =  er_im - er_im_c 
	;	et_im =  et_im - et_im_c 
	;	ez_im =  ez_im - ez_im_c 
	;	stop
	;endif

	a_jP_r_re = (total(a_jP_r_re,3))[*,0]
	a_jP_r_im = (total(a_jP_r_im,3))[*,0]
	a_jP_t_re = (total(a_jP_t_re,3))[*,0]
	a_jP_t_im = (total(a_jP_t_im,3))[*,0]
	a_jP_z_re = (total(a_jP_z_re,3))[*,0]
	a_jP_z_im = (total(a_jP_z_im,3))[*,0]

	;jP_r_re = interpol(a_jP_r_re,a_r,r,/spline)
	;jP_r_im = interpol(a_jP_r_im,a_r,r,/spline)
	;jP_t_re = interpol(a_jP_t_re,a_r,r,/spline)
	;jP_t_im = interpol(a_jP_t_im,a_r,r,/spline)
	;jP_z_re = interpol(a_jP_z_re,a_r,r,/spline)
	;jP_z_im = interpol(a_jP_z_im,a_r,r,/spline)

	;jP_r_re_ = interpol(a_jP_r_re,a_r,r_,/spline)
	;jP_r_im_ = interpol(a_jP_r_im,a_r,r_,/spline)
	;jP_t_re_ = interpol(a_jP_t_re,a_r,r_,/spline)
	;jP_t_im_ = interpol(a_jP_t_im,a_r,r_,/spline)
	;jP_z_re_ = interpol(a_jP_z_re,a_r,r_,/spline)
	;jP_z_im_ = interpol(a_jP_z_im,a_r,r_,/spline)


	; Write netCDF file

	nc_id = nCdf_create ('kj_aorsa_1d.nc', /clobber )

	nCdf_control, nc_id, /fill
	
	nr_id = nCdf_dimDef ( nc_id, 'nR', n_elements(r) )
	;nr_id_ = nCdf_dimDef ( nc_id, 'nR_', n_elements(r_) )
	scalar_id = nCdf_dimDef ( nc_id, 'scalar', 1 )

	freq_id = nCdf_varDef ( nc_id, 'freq', scalar_id, /float )
	r_id = nCdf_varDef ( nc_id, 'r', nr_id, /float )
	;r_id_ = nCdf_varDef ( nc_id, 'r_', nr_id_, /float )

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

	;jP_r_re_id_ = nCdf_varDef ( nc_id, 'jP_r_re_', nr_id_, /float )
	;jP_r_im_id_ = nCdf_varDef ( nc_id, 'jP_r_im_', nr_id_, /float )
	;jP_p_re_id_ = nCdf_varDef ( nc_id, 'jP_p_re_', nr_id_, /float )
	;jP_p_im_id_ = nCdf_varDef ( nc_id, 'jP_p_im_', nr_id_, /float )
	;jP_z_re_id_ = nCdf_varDef ( nc_id, 'jP_z_re_', nr_id_, /float )
	;jP_z_im_id_ = nCdf_varDef ( nc_id, 'jP_z_im_', nr_id_, /float )

	jA_r_re_id = nCdf_varDef ( nc_id, 'jA_r_re', nr_id, /float )
	jA_r_im_id = nCdf_varDef ( nc_id, 'jA_r_im', nr_id, /float )
	jA_p_re_id = nCdf_varDef ( nc_id, 'jA_p_re', nr_id, /float )
	jA_p_im_id = nCdf_varDef ( nc_id, 'jA_p_im', nr_id, /float )
	jA_z_re_id = nCdf_varDef ( nc_id, 'jA_z_re', nr_id, /float )
	jA_z_im_id = nCdf_varDef ( nc_id, 'jA_z_im', nr_id, /float )

	nCdf_control, nc_id, /enDef

	nCdf_varPut, nc_id, freq_id, freq

	nCdf_varPut, nc_id, r_id, a_r 
	;nCdf_varPut, nc_id, r_id_, r_ 

	nCdf_varPut, nc_id, B0_r_id, brU*bMod
	nCdf_varPut, nc_id, B0_p_id, btU*bMod 
	nCdf_varPut, nc_id, B0_z_id, bzU*bMod

	nCdf_varPut, nc_id, e_r_re_id, er_re
	nCdf_varPut, nc_id, e_r_im_id, er_im
	nCdf_varPut, nc_id, e_p_re_id, et_re
	nCdf_varPut, nc_id, e_p_im_id, et_im
	nCdf_varPut, nc_id, e_z_re_id, ez_re 
	nCdf_varPut, nc_id, e_z_im_id, ez_im

	nCdf_varPut, nc_id, jP_r_re_id, a_jP_r_re 
	nCdf_varPut, nc_id, jP_r_im_id, a_jP_r_im 
	nCdf_varPut, nc_id, jP_p_re_id, a_jP_t_re 
	nCdf_varPut, nc_id, jP_p_im_id, a_jP_t_im 
	nCdf_varPut, nc_id, jP_z_re_id, a_jP_z_re 
	nCdf_varPut, nc_id, jP_z_im_id, a_jP_z_im 

	;nCdf_varPut, nc_id, jP_r_re_id_, jP_r_re_ 
	;nCdf_varPut, nc_id, jP_r_im_id_, jP_r_im_ 
	;nCdf_varPut, nc_id, jP_p_re_id_, jP_t_re_ 
	;nCdf_varPut, nc_id, jP_p_im_id_, jP_t_im_ 
	;nCdf_varPut, nc_id, jP_z_re_id_, jP_z_re_ 
	;nCdf_varPut, nc_id, jP_z_im_id_, jP_z_im_ 

	nCdf_varPut, nc_id, jA_r_re_id, jA_r_re[*,0] 
	nCdf_varPut, nc_id, jA_r_im_id, jA_r_im[*,0] 
	nCdf_varPut, nc_id, jA_p_re_id, jA_t_re[*,0] 
	nCdf_varPut, nc_id, jA_p_im_id, jA_t_im[*,0] 
	nCdf_varPut, nc_id, jA_z_re_id, jA_z_re[*,0] 
	nCdf_varPut, nc_id, jA_z_im_id, jA_z_im[*,0] 

	nCdf_close, nc_id

end
