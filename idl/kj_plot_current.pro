pro kj_plot_current

	cdfId = ncdf_open('output/orbits.nc')

		ncdf_varget, cdfId, 'x', x
		ncdf_varget, cdfId, 'y', y
		ncdf_varget, cdfId, 'z', z
		ncdf_varget, cdfId, 't', t
		ncdf_varget, cdfId, 'e1_x', e1_x 
		ncdf_varget, cdfId, 'e1_y', e1_y 
		ncdf_varget, cdfId, 'e1_z', e1_z 

		ncdf_varget, cdfId, 'j1x', e1_x 
		ncdf_varget, cdfId, 'j1y', e1_y 
		ncdf_varget, cdfId, 'j1z', e1_z 

		ncdf_varget, cdfId, 'v1x', e1_x 
		ncdf_varget, cdfId, 'v1y', e1_y 
		ncdf_varget, cdfId, 'v1z', e1_z 

	nCdf_close,	cdfId 


	cdfId = ncdf_open('data/rsfwc_1d.nc')

		ncdf_varget, cdfId, 'wrf', wrf
		ncdf_varget, cdfId, 'r', r 

		ncdf_varget, cdfId, 'e_r_re', er_re
		ncdf_varget, cdfId, 'e_r_im', er_im
		ncdf_varget, cdfId, 'e_p_re', ep_re
		ncdf_varget, cdfId, 'e_p_im', ep_im
		ncdf_varget, cdfId, 'e_z_re', ez_re
		ncdf_varget, cdfId, 'e_z_im', ez_im


	ncdf_close, cdfId

end
