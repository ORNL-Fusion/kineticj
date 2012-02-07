pro kj_plot_orbit

	fileList = file_search ( 'output/jP*' )

	cdfId = ncdf_open(fileList[0])
		ncdf_varget, cdfId, 'freq', freq 
		ncdf_varget, cdfId, 't', tJ 
	nCdf_close, cdfId
	
	fileList = file_search ( 'output/orbits*' )

	cdfId = ncdf_open(fileList[0])
		ncdf_varget, cdfId, 't', t
	nCdf_close, cdfId

	nT = n_elements(t)
	nF = n_elements(fileList)

	maxE = 0
	for f=0,nF-1 do begin

		cdfId = ncdf_open(fileList[f])

			ncdf_varget, cdfId, 'x', x_0
			ncdf_varget, cdfId, 'y', y_0
			ncdf_varget, cdfId, 'z', z_0

			ncdf_varget, cdfId, 't', t_0

			ncdf_varget, cdfId, 'v1x', v1x_0 
			ncdf_varget, cdfId, 'v1y', v1y_0 
			ncdf_varget, cdfId, 'v1z', v1z_0 

			ncdf_varget, cdfId, 'e1_x', e1x_0 
			ncdf_varget, cdfId, 'e1_y', e1y_0 
			ncdf_varget, cdfId, 'e1_z', e1z_0 

		nCdf_close,	cdfId 

		if(max(e1x_0) gt maxE) then maxE = max(e1x_0)

	endfor

	pNum = 20
	;p=plot(t_0*freq,v1x_0[*,pNum])
	;p=plot(t_0*freq,e1x_0[*,pNum]

stop
end
