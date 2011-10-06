pro kj_plot_current

	fileList = file_search ( 'output/jP*' )

	cdfId = ncdf_open(fileList[0])
	ncdf_varget, cdfId, 't', t
	nCdf_close, cdfId
	
	nT = n_elements(t)
	nF = n_elements(fileList)

	j1x = fltArr ( nT, nF )
	xF = fltArr ( nF )

	for f=0,nF-1 do begin

		cdfId = ncdf_open(fileList[f])

			ncdf_varget, cdfId, 'x', x
			ncdf_varget, cdfId, 't', t

			ncdf_varget, cdfId, 'j1x', j1x_0 
			ncdf_varget, cdfId, 'j1y', j1y_0 
			ncdf_varget, cdfId, 'j1z', j1z_0 

		nCdf_close,	cdfId 

		j1x[*,f] = j1x_0
		xF[f] = x

	endfor

	cdfId = ncdf_open('data/rsfwc_1d.nc')

		ncdf_varget, cdfId, 'wrf', wrf
		ncdf_varget, cdfId, 'r', r 

		ncdf_varget, cdfId, 'e_r_re', er_re
		ncdf_varget, cdfId, 'e_r_im', er_im
		ncdf_varget, cdfId, 'e_p_re', ep_re
		ncdf_varget, cdfId, 'e_p_im', ep_im
		ncdf_varget, cdfId, 'e_z_re', ez_re
		ncdf_varget, cdfId, 'e_z_im', ez_im

		ncdf_varget, cdfId, 'j_r_re', jr_re
		ncdf_varget, cdfId, 'j_r_im', jr_im
		ncdf_varget, cdfId, 'j_p_re', jp_re
		ncdf_varget, cdfId, 'j_p_im', jp_im
		ncdf_varget, cdfId, 'j_z_re', jz_re
		ncdf_varget, cdfId, 'j_z_im', jz_im

		ncdf_varget, cdfId, 'jA_r_re', jAr_re
		ncdf_varget, cdfId, 'jA_r_im', jAr_im
		ncdf_varget, cdfId, 'jA_p_re', jAp_re
		ncdf_varget, cdfId, 'jA_p_im', jAp_im
		ncdf_varget, cdfId, 'jA_z_re', jAz_re
		ncdf_varget, cdfId, 'jA_z_im', jAz_im

	ncdf_close, cdfId


	for i=0,n_elements(t)-1 do begin	
		plot, r,jr_re*cos(wrf*t[i])+jr_im*sin(wrf*t[i]),yRange=[-100,100], xRange=[9.5,10.5]
		;plots, r,jAr_re*cos(wrf*t[i])+jAr_im*sin(wrf*t[i])
		for f=0,nF-1 do begin
			plots, xF[f], j1x[i,f], psym=4
		endfor
		wait, 0.2
	endfor


stop
end
