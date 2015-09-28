pro kj_combine_species

    baseName = 'colestock-kashuba-'
    dirList = [ baseName+'spec0', baseName+'spec1', baseName+'spec2' ]

    nS = n_elements(dirList)
    j = kj_read_jp_nc( runDirName = dirList[0] )

    nR = n_elements(j.r)
    nZ = n_elements(j.z)

    jP_r_s = complexArr(nR,nZ,nS)
    jP_t_s = complexArr(nR,nZ,nS)
    jP_z_s = complexArr(nR,nZ,nS)

    for s=0,nS-1 do begin

        j = kj_read_jp_nc( runDirName = dirList[s] )

        jP_r_s[*,*,s] = j.jP_r
        jP_t_s[*,*,s] = j.jP_t
        jP_z_s[*,*,s] = j.jP_z

    endfor

	; Write kj_jP in file

	nc_id = nCdf_create ('kj_jP_combined.nc', /clobber )

	nCdf_control, nc_id, /fill
	
	nr_id = nCdf_dimDef ( nc_id, 'nR', n_elements(j.r) )
	nx_id = nCdf_dimDef ( nc_id, 'nX', n_elements(j.r) )
	ny_id = nCdf_dimDef ( nc_id, 'nY', n_elements(j.z) )
	ns_id = nCdf_dimDef ( nc_id, 'nSpec', nS )
	scalar_id = nCdf_dimDef ( nc_id, 'scalar', 1 )

	freq_id = nCdf_varDef ( nc_id, 'freq', scalar_id, /float )
	r_id = nCdf_varDef ( nc_id, 'r', nr_id, /float )
	z_id = nCdf_varDef ( nc_id, 'z', ny_id, /float )

	jP_r_re_id = nCdf_varDef ( nc_id, 'jP_r_re', [nx_id,ny_id,ns_id], /float )
	jP_r_im_id = nCdf_varDef ( nc_id, 'jP_r_im', [nx_id,ny_id,ns_id], /float )
	jP_p_re_id = nCdf_varDef ( nc_id, 'jP_t_re', [nx_id,ny_id,ns_id], /float )
	jP_p_im_id = nCdf_varDef ( nc_id, 'jP_t_im', [nx_id,ny_id,ns_id], /float )
	jP_z_re_id = nCdf_varDef ( nc_id, 'jP_z_re', [nx_id,ny_id,ns_id], /float )
	jP_z_im_id = nCdf_varDef ( nc_id, 'jP_z_im', [nx_id,ny_id,ns_id], /float )

	nCdf_control, nc_id, /enDef

	nCdf_varPut, nc_id, freq_id, j.freq

	nCdf_varPut, nc_id, r_id, j.r
	nCdf_varPut, nc_id, z_id, j.z

	nCdf_varPut, nc_id, jP_r_re_id, real_part(jP_r_s)
	nCdf_varPut, nc_id, jP_r_im_id, imaginary(jP_r_s) 
	nCdf_varPut, nc_id, jP_p_re_id, real_part(jP_t_s)
	nCdf_varPut, nc_id, jP_p_im_id, imaginary(jP_t_s) 
	nCdf_varPut, nc_id, jP_z_re_id, real_part(jP_z_s)
	nCdf_varPut, nc_id, jP_z_im_id, imaginary(jP_z_s) 

	nCdf_close, nc_id


end
