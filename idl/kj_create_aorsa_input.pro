pro kj_create_aorsa_input, rsfwcPath=rsfwcPath

	match_rsfwc = 0
	if(keyword_set(rsfwcPath))then begin
		; rsfwcPath is the path to the rsfwc output file to which we will match
		; the grid to, so the aorsa jP can be used as input to rsfwc.

		cdfId = ncdf_open(rsfwcPath)

			ncdf_varget, cdfId, 'freq', freq 
			ncdf_varget, cdfId, 'r', r 
			ncdf_varget, cdfId, 'r_', r_

		ncdf_close, cdfId

		match_rsfwc = 1
	endif

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

	a_br = brU*bMod
	a_bt = brU*bMod
	a_bz = brU*bMod

	a_er_re = er_re
	a_er_im = er_im
	a_et_re = et_re
	a_et_im = et_im
	a_ez_re = ez_re
	a_ez_im = ez_im

	a_jA_r_re = jA_r_re[*,0]
	a_jA_r_im = jA_r_im[*,0]
	a_jA_t_re = jA_t_re[*,0]
	a_jA_t_im = jA_t_im[*,0]
	a_jA_z_re = jA_z_re[*,0]
	a_jA_z_im = jA_z_im[*,0]

	if match_rsfwc eq 1 then begin

		_r = r
		_r_ = r_

		_br = interpol(a_br,a_r,r,/spline)
		_bt = interpol(a_bt,a_r,r,/spline)
		_bz = interpol(a_bz,a_r,r,/spline)

		_br_ = interpol(a_br,a_r,r_,/spline)
		_bt_ = interpol(a_bt,a_r,r_,/spline)
		_bz_ = interpol(a_bz,a_r,r_,/spline)

		_er_re = interpol(a_er_re,a_r,r,/spline)
		_er_im = interpol(a_er_im,a_r,r,/spline)
		_et_re = interpol(a_et_re,a_r,r,/spline)
		_et_im = interpol(a_et_im,a_r,r,/spline)
		_ez_re = interpol(a_ez_re,a_r,r,/spline)
		_ez_im = interpol(a_ez_im,a_r,r,/spline)

		_er_re_ = interpol(a_er_re,a_r,r_,/spline)
		_er_im_ = interpol(a_er_im,a_r,r_,/spline)
		_et_re_ = interpol(a_et_re,a_r,r_,/spline)
		_et_im_ = interpol(a_et_im,a_r,r_,/spline)
		_ez_re_ = interpol(a_ez_re,a_r,r_,/spline)
		_ez_im_ = interpol(a_ez_im,a_r,r_,/spline)

		_jP_r_re = interpol(a_jP_r_re,a_r,r,/spline)
		_jP_r_im = interpol(a_jP_r_im,a_r,r,/spline)
		_jP_t_re = interpol(a_jP_t_re,a_r,r,/spline)
		_jP_t_im = interpol(a_jP_t_im,a_r,r,/spline)
		_jP_z_re = interpol(a_jP_z_re,a_r,r,/spline)
		_jP_z_im = interpol(a_jP_z_im,a_r,r,/spline)

		_jP_r_re_ = interpol(a_jP_r_re,a_r,r_,/spline)
		_jP_r_im_ = interpol(a_jP_r_im,a_r,r_,/spline)
		_jP_t_re_ = interpol(a_jP_t_re,a_r,r_,/spline)
		_jP_t_im_ = interpol(a_jP_t_im,a_r,r_,/spline)
		_jP_z_re_ = interpol(a_jP_z_re,a_r,r_,/spline)
		_jP_z_im_ = interpol(a_jP_z_im,a_r,r_,/spline)

		_jA_r_re = interpol(a_jA_r_re,a_r,r,/spline)
		_jA_r_im = interpol(a_jA_r_im,a_r,r,/spline)
		_jA_t_re = interpol(a_jA_t_re,a_r,r,/spline)
		_jA_t_im = interpol(a_jA_t_im,a_r,r,/spline)
		_jA_z_re = interpol(a_jA_z_re,a_r,r,/spline)
		_jA_z_im = interpol(a_jA_z_im,a_r,r,/spline)

		_jA_r_re_ = interpol(a_jA_r_re,a_r,r_,/spline)
		_jA_r_im_ = interpol(a_jA_r_im,a_r,r_,/spline)
		_jA_t_re_ = interpol(a_jA_t_re,a_r,r_,/spline)
		_jA_t_im_ = interpol(a_jA_t_im,a_r,r_,/spline)
		_jA_z_re_ = interpol(a_jA_z_re,a_r,r_,/spline)
		_jA_z_im_ = interpol(a_jA_z_im,a_r,r_,/spline)

	endif else begin

		r = a_r
		_r = r

		_br = interpol(a_br,a_r,r,/spline)
		_bt = interpol(a_bt,a_r,r,/spline)
		_bz = interpol(a_bz,a_r,r,/spline)

		_er_re = interpol(a_er_re,a_r,r,/spline)
		_er_im = interpol(a_er_im,a_r,r,/spline)
		_et_re = interpol(a_et_re,a_r,r,/spline)
		_et_im = interpol(a_et_im,a_r,r,/spline)
		_ez_re = interpol(a_ez_re,a_r,r,/spline)
		_ez_im = interpol(a_ez_im,a_r,r,/spline)

		_jP_r_re = interpol(a_jP_r_re,a_r,r,/spline)
		_jP_r_im = interpol(a_jP_r_im,a_r,r,/spline)
		_jP_t_re = interpol(a_jP_t_re,a_r,r,/spline)
		_jP_t_im = interpol(a_jP_t_im,a_r,r,/spline)
		_jP_z_re = interpol(a_jP_z_re,a_r,r,/spline)
		_jP_z_im = interpol(a_jP_z_im,a_r,r,/spline)

		_jA_r_re = interpol(a_jA_r_re,a_r,r,/spline)
		_jA_r_im = interpol(a_jA_r_im,a_r,r,/spline)
		_jA_t_re = interpol(a_jA_t_re,a_r,r,/spline)
		_jA_t_im = interpol(a_jA_t_im,a_r,r,/spline)
		_jA_z_re = interpol(a_jA_z_re,a_r,r,/spline)
		_jA_z_im = interpol(a_jA_z_im,a_r,r,/spline)


	endelse


	; Write netCDF file

	nc_id = nCdf_create ('kj_aorsa_1d.nc', /clobber )

	nCdf_control, nc_id, /fill
	nCdf_control, nc_id, /verbose
	
	nr_id = nCdf_dimDef ( nc_id, 'nR', n_elements(_r) )
	scalar_id = nCdf_dimDef ( nc_id, 'scalar', 1 )
	if match_rsfwc eq 1 then nr_id_ = nCdf_dimDef ( nc_id, 'nR_', n_elements(_r_) )

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


	if match_rsfwc eq 1 then begin

		r_id_ = nCdf_varDef ( nc_id, 'r_', nr_id_, /float )

		B0_r_id_ = nCdf_varDef ( nc_id, 'B0_r_', nr_id_, /float )
		B0_p_id_ = nCdf_varDef ( nc_id, 'B0_p_', nr_id_, /float )
		B0_z_id_ = nCdf_varDef ( nc_id, 'B0_z_', nr_id_, /float )

		e_r_re_id_ = nCdf_varDef ( nc_id, 'e_r_re_', nr_id_, /float )
		e_r_im_id_ = nCdf_varDef ( nc_id, 'e_r_im_', nr_id_, /float )
		e_p_re_id_ = nCdf_varDef ( nc_id, 'e_p_re_', nr_id_, /float )
		e_p_im_id_ = nCdf_varDef ( nc_id, 'e_p_im_', nr_id_, /float )
		e_z_re_id_ = nCdf_varDef ( nc_id, 'e_z_re_', nr_id_, /float )
		e_z_im_id_ = nCdf_varDef ( nc_id, 'e_z_im_', nr_id_, /float )

		jP_r_re_id_ = nCdf_varDef ( nc_id, 'jP_r_re_', nr_id_, /float )
		jP_r_im_id_ = nCdf_varDef ( nc_id, 'jP_r_im_', nr_id_, /float )
		jP_p_re_id_ = nCdf_varDef ( nc_id, 'jP_p_re_', nr_id_, /float )
		jP_p_im_id_ = nCdf_varDef ( nc_id, 'jP_p_im_', nr_id_, /float )
		jP_z_re_id_ = nCdf_varDef ( nc_id, 'jP_z_re_', nr_id_, /float )
		jP_z_im_id_ = nCdf_varDef ( nc_id, 'jP_z_im_', nr_id_, /float )

		jA_r_re_id_ = nCdf_varDef ( nc_id, 'jA_r_re_', nr_id_, /float )
		jA_r_im_id_ = nCdf_varDef ( nc_id, 'jA_r_im_', nr_id_, /float )
		jA_p_re_id_ = nCdf_varDef ( nc_id, 'jA_p_re_', nr_id_, /float )
		jA_p_im_id_ = nCdf_varDef ( nc_id, 'jA_p_im_', nr_id_, /float )
		jA_z_re_id_ = nCdf_varDef ( nc_id, 'jA_z_re_', nr_id_, /float )
		jA_z_im_id_ = nCdf_varDef ( nc_id, 'jA_z_im_', nr_id_, /float )

	endif

	nCdf_control, nc_id, /enDef

	nCdf_varPut, nc_id, freq_id, freq

	nCdf_varPut, nc_id, r_id, _r 

	nCdf_varPut, nc_id, B0_r_id, _br 
	nCdf_varPut, nc_id, B0_p_id, _bt 
	nCdf_varPut, nc_id, B0_z_id, _bz

	nCdf_varPut, nc_id, e_r_re_id, _er_re
	nCdf_varPut, nc_id, e_r_im_id, _er_im
	nCdf_varPut, nc_id, e_p_re_id, _et_re
	nCdf_varPut, nc_id, e_p_im_id, _et_im
	nCdf_varPut, nc_id, e_z_re_id, _ez_re 
	nCdf_varPut, nc_id, e_z_im_id, _ez_im

	nCdf_varPut, nc_id, jP_r_re_id, _jP_r_re 
	nCdf_varPut, nc_id, jP_r_im_id, _jP_r_im 
	nCdf_varPut, nc_id, jP_p_re_id, _jP_t_re 
	nCdf_varPut, nc_id, jP_p_im_id, _jP_t_im 
	nCdf_varPut, nc_id, jP_z_re_id, _jP_z_re 
	nCdf_varPut, nc_id, jP_z_im_id, _jP_z_im 

	nCdf_varPut, nc_id, jA_r_re_id, _jA_r_re 
	nCdf_varPut, nc_id, jA_r_im_id, _jA_r_im 
	nCdf_varPut, nc_id, jA_p_re_id, _jA_t_re 
	nCdf_varPut, nc_id, jA_p_im_id, _jA_t_im 
	nCdf_varPut, nc_id, jA_z_re_id, _jA_z_re 
	nCdf_varPut, nc_id, jA_z_im_id, _jA_z_im 

	if match_rsfwc then begin

		nCdf_varPut, nc_id, r_id_, _r_

		nCdf_varPut, nc_id, B0_r_id_, _br_ 
		nCdf_varPut, nc_id, B0_p_id_, _bt_ 
		nCdf_varPut, nc_id, B0_z_id_, _bz_

		nCdf_varPut, nc_id, e_r_re_id_, _er_re_
		nCdf_varPut, nc_id, e_r_im_id_, _er_im_
		nCdf_varPut, nc_id, e_p_re_id_, _et_re_
		nCdf_varPut, nc_id, e_p_im_id_, _et_im_
		nCdf_varPut, nc_id, e_z_re_id_, _ez_re_ 
		nCdf_varPut, nc_id, e_z_im_id_, _ez_im_

		nCdf_varPut, nc_id, jP_r_re_id_, _jP_r_re_ 
		nCdf_varPut, nc_id, jP_r_im_id_, _jP_r_im_ 
		nCdf_varPut, nc_id, jP_p_re_id_, _jP_t_re_ 
		nCdf_varPut, nc_id, jP_p_im_id_, _jP_t_im_ 
		nCdf_varPut, nc_id, jP_z_re_id_, _jP_z_re_ 
		nCdf_varPut, nc_id, jP_z_im_id_, _jP_z_im_ 

		nCdf_varPut, nc_id, jA_r_re_id_, _jA_r_re_ 
		nCdf_varPut, nc_id, jA_r_im_id_, _jA_r_im_ 
		nCdf_varPut, nc_id, jA_p_re_id_, _jA_t_re_ 
		nCdf_varPut, nc_id, jA_p_im_id_, _jA_t_im_ 
		nCdf_varPut, nc_id, jA_z_re_id_, _jA_z_re_ 
		nCdf_varPut, nc_id, jA_z_im_id_, _jA_z_im_ 

	endif

	nCdf_close, nc_id
stop
end
