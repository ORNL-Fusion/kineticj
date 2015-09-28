function kj_read_jp_nc, runDirName = runDirName

    if keyword_set(runDirName) then fPath = runDirName else fPath = './'
    fName = fPath+'/output/kj_jP.nc'

	cdfId = ncdf_open(fName)

		ncdf_varget, cdfId, 'r', r
		ncdf_varget, cdfId, 'z', z
		ncdf_varget, cdfId, 'freq', freq

		ncdf_varget, cdfId, 'jP_r_re', jP_r_re 
		ncdf_varget, cdfId, 'jP_r_im', jP_r_im

	    ncdf_varget, cdfId, 'jP_t_re', jP_t_re 
		ncdf_varget, cdfId, 'jP_t_im', jP_t_im

	    ncdf_varget, cdfId, 'jP_z_re', jP_z_re 
		ncdf_varget, cdfId, 'jP_z_im', jP_z_im

	nCdf_close,	cdfId 

    nR = n_elements(r)
    nZ = n_elements(z)
    nSpec = n_elements(jP_r_re[0,0,*]) 

    data = { r:r, $
            z:z, $
            freq:freq, $
            jP_r : complex(jP_r_re,jP_r_im), $
            jP_t : complex(jP_t_re,jP_t_im), $
            jP_z : complex(jP_z_re,jP_z_im) }
    
    return, data 

end
