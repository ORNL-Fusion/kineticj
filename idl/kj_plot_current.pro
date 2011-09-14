pro kj_plot_current

	fileList = file_search ( 'output/orbits*' )

	cdfId = ncdf_open(fileList[0])
	ncdf_varget, cdfId, 't', t
	nCdf_close, cdfId
	
	nT = n_elements(t)
	nF = n_elements(fileList)

	j1x = fltArr ( nT, nF )

	for f=0,nF-1 do begin

		cdfId = ncdf_open(fileList[f])

			ncdf_varget, cdfId, 'x', x
			ncdf_varget, cdfId, 'y', y
			ncdf_varget, cdfId, 'z', z
			ncdf_varget, cdfId, 't', t
			ncdf_varget, cdfId, 'e1_x', e1_x 
			ncdf_varget, cdfId, 'e1_y', e1_y 
			ncdf_varget, cdfId, 'e1_z', e1_z 

			ncdf_varget, cdfId, 'j1x', j1x_0 
			ncdf_varget, cdfId, 'j1y', j1y_0 
			ncdf_varget, cdfId, 'j1z', j1z_0 

			ncdf_varget, cdfId, 'v1x', v1x 
			ncdf_varget, cdfId, 'v1y', v1y 
			ncdf_varget, cdfId, 'v1z', v1z 

		nCdf_close,	cdfId 

		j1x[*,f] = j1x_0

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

	ncdf_close, cdfId


	for i=0,n_elements(t)-1 do begin	
		plot, r,jr_re*cos(wrf*t[i])+jr_im*sin(wrf*t[i]),yRange=[-500,500], xRange=[9.5,10.5]
		for f=0,nF-1 do begin
			plots, 9.95+f*0.1/10, j1x[i,f], psym=4
		endfor
		wait, 0.05
	endfor


stop
end
