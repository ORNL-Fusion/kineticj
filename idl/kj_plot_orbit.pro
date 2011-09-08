pro kj_plot_orbit

	cdfId = ncdf_open('output/orbits.nc')

		ncdf_varget, cdfId, 'x', x
		ncdf_varget, cdfId, 'y', y
		ncdf_varget, cdfId, 'z', z
		ncdf_varget, cdfId, 't', t
		ncdf_varget, cdfId, 'e1_x', e1_x 
		ncdf_varget, cdfId, 'e1_y', e1_y 
		ncdf_varget, cdfId, 'e1_z', e1_z 

	nCdf_close,	cdfId 

	p1 = plot3d ( x, y, z )
	p2 = plot3d ( x, y, z, aspect_ratio = 1.0, aspect_z = 1.0 )

stop
end
