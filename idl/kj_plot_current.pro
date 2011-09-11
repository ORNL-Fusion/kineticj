pro kj_plot_current

	cdfId = ncdf_open('output/orbits_1.0.nc')

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

	cdfId = ncdf_open('output/orbits_1.0A.nc')

		ncdf_varget, cdfId, 'x', x
		ncdf_varget, cdfId, 'y', y
		ncdf_varget, cdfId, 'z', z
		ncdf_varget, cdfId, 't', t
		ncdf_varget, cdfId, 'e1_x', e1_x 
		ncdf_varget, cdfId, 'e1_y', e1_y 
		ncdf_varget, cdfId, 'e1_z', e1_z 

		ncdf_varget, cdfId, 'j1x', j1x_0A 
		ncdf_varget, cdfId, 'j1y', j1y_0A 
		ncdf_varget, cdfId, 'j1z', j1z_0A 

		ncdf_varget, cdfId, 'v1x', v1x 
		ncdf_varget, cdfId, 'v1y', v1y 
		ncdf_varget, cdfId, 'v1z', v1z 

	nCdf_close,	cdfId 


	cdfId = ncdf_open('output/orbits_1.01.nc')

		ncdf_varget, cdfId, 'x', x
		ncdf_varget, cdfId, 'y', y
		ncdf_varget, cdfId, 'z', z
		ncdf_varget, cdfId, 't', t
		ncdf_varget, cdfId, 'e1_x', e1_x 
		ncdf_varget, cdfId, 'e1_y', e1_y 
		ncdf_varget, cdfId, 'e1_z', e1_z 

		ncdf_varget, cdfId, 'j1x', j1x_1
		ncdf_varget, cdfId, 'j1y', j1y_1
		ncdf_varget, cdfId, 'j1z', j1z_1

		ncdf_varget, cdfId, 'v1x', v1x 
		ncdf_varget, cdfId, 'v1y', v1y 
		ncdf_varget, cdfId, 'v1z', v1z 

	nCdf_close,	cdfId 

	cdfId = ncdf_open('output/orbits_1.02.nc')

		ncdf_varget, cdfId, 'x', x
		ncdf_varget, cdfId, 'y', y
		ncdf_varget, cdfId, 'z', z
		ncdf_varget, cdfId, 't', t
		ncdf_varget, cdfId, 'e1_x', e1_x 
		ncdf_varget, cdfId, 'e1_y', e1_y 
		ncdf_varget, cdfId, 'e1_z', e1_z 

		ncdf_varget, cdfId, 'j1x', j1x_2
		ncdf_varget, cdfId, 'j1y', j1y_2
		ncdf_varget, cdfId, 'j1z', j1z_2

		ncdf_varget, cdfId, 'v1x', v1x 
		ncdf_varget, cdfId, 'v1y', v1y 
		ncdf_varget, cdfId, 'v1z', v1z 

	nCdf_close,	cdfId 

	cdfId = ncdf_open('output/orbits_1.03.nc')

		ncdf_varget, cdfId, 'x', x
		ncdf_varget, cdfId, 'y', y
		ncdf_varget, cdfId, 'z', z
		ncdf_varget, cdfId, 't', t
		ncdf_varget, cdfId, 'e1_x', e1_x 
		ncdf_varget, cdfId, 'e1_y', e1_y 
		ncdf_varget, cdfId, 'e1_z', e1_z 

		ncdf_varget, cdfId, 'j1x', j1x_3
		ncdf_varget, cdfId, 'j1y', j1y_3
		ncdf_varget, cdfId, 'j1z', j1z_3

		ncdf_varget, cdfId, 'v1x', v1x 
		ncdf_varget, cdfId, 'v1y', v1y 
		ncdf_varget, cdfId, 'v1z', v1z 

	nCdf_close,	cdfId 

	cdfId = ncdf_open('output/orbits_1.04.nc')

		ncdf_varget, cdfId, 'x', x
		ncdf_varget, cdfId, 'y', y
		ncdf_varget, cdfId, 'z', z
		ncdf_varget, cdfId, 't', t
		ncdf_varget, cdfId, 'e1_x', e1_x 
		ncdf_varget, cdfId, 'e1_y', e1_y 
		ncdf_varget, cdfId, 'e1_z', e1_z 

		ncdf_varget, cdfId, 'j1x', j1x_4
		ncdf_varget, cdfId, 'j1y', j1y_4
		ncdf_varget, cdfId, 'j1z', j1z_4

		ncdf_varget, cdfId, 'v1x', v1x 
		ncdf_varget, cdfId, 'v1y', v1y 
		ncdf_varget, cdfId, 'v1z', v1z 

	nCdf_close,	cdfId 

	cdfId = ncdf_open('output/orbits_1.05.nc')

		ncdf_varget, cdfId, 'x', x
		ncdf_varget, cdfId, 'y', y
		ncdf_varget, cdfId, 'z', z
		ncdf_varget, cdfId, 't', t
		ncdf_varget, cdfId, 'e1_x', e1_x 
		ncdf_varget, cdfId, 'e1_y', e1_y 
		ncdf_varget, cdfId, 'e1_z', e1_z 

		ncdf_varget, cdfId, 'j1x', j1x_5
		ncdf_varget, cdfId, 'j1y', j1y_5
		ncdf_varget, cdfId, 'j1z', j1z_5

		ncdf_varget, cdfId, 'v1x', v1x 
		ncdf_varget, cdfId, 'v1y', v1y 
		ncdf_varget, cdfId, 'v1z', v1z 

	nCdf_close,	cdfId 

	cdfId = ncdf_open('output/orbits_1.05A.nc')

		ncdf_varget, cdfId, 'x', x
		ncdf_varget, cdfId, 'y', y
		ncdf_varget, cdfId, 'z', z
		ncdf_varget, cdfId, 't', t
		ncdf_varget, cdfId, 'e1_x', e1_x 
		ncdf_varget, cdfId, 'e1_y', e1_y 
		ncdf_varget, cdfId, 'e1_z', e1_z 

		ncdf_varget, cdfId, 'j1x', j1x_5A
		ncdf_varget, cdfId, 'j1y', j1y_5A
		ncdf_varget, cdfId, 'j1z', j1z_5A

		ncdf_varget, cdfId, 'v1x', v1x 
		ncdf_varget, cdfId, 'v1y', v1y 
		ncdf_varget, cdfId, 'v1z', v1z 

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

		ncdf_varget, cdfId, 'j_r_re', jr_re
		ncdf_varget, cdfId, 'j_r_im', jr_im
		ncdf_varget, cdfId, 'j_p_re', jp_re
		ncdf_varget, cdfId, 'j_p_im', jp_im
		ncdf_varget, cdfId, 'j_z_re', jz_re
		ncdf_varget, cdfId, 'j_z_im', jz_im

	ncdf_close, cdfId

	for i=0,n_elements(t)-1 do begin	
		plot, r,jr_re*cos(wrf*t[i])+jr_im*sin(wrf*t[i]),yRange=[-500,500]
		plots, 1.00, j1x_0[i], psym=4
		plots, 1.00, j1x_0A[i], psym=4
		plots, 1.01, j1x_1[i], psym=4
		plots, 1.02, j1x_2[i], psym=4
		plots, 1.03, j1x_3[i], psym=4
		plots, 1.04, j1x_4[i], psym=4
		plots, 1.05, j1x_5[i], psym=4
		plots, 1.05, j1x_5A[i], psym=4
		wait, 0.05
	endfor


stop
end
