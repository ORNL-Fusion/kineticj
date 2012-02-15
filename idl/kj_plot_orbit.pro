pro kj_plot_orbit

	runFile = 'kj.cfg'
	cfg = kj_read_cfg (runFile)	

	; Read particle list

	particleList = cfg.particleList_fName
	cdfId = ncdf_open(particleList)
		ncdf_varget, cdfId, 'vx', p_vx 
		ncdf_varget, cdfId, 'weight', p_weight
	nCdf_close, cdfId

	binSize = 0.3e6
	vMax = 0.5e8
	vMin = -vMax 
	nBins = (vMax-vMin)/binSize
	histX = fIndGen(nBins)*binSize+vMin+binSize/2.0
	histY = fltArr(nBins)
	for n=0,n_elements(p_vx)-1 do begin
		ii = (p_vx[n]-vMin)/(vMax-vMin)*nBins
		histY[ii] += p_weight[n]
	endfor

	p=plot(histX,histY)

	fileList = file_search ( 'output/'+cfg.runIdent+'/jP*' )

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


	binSize = 0.1
	vMax = 910
	vMin = 890 
	nBins = (vMax-vMin)/binSize
	histXa = fIndGen(nBins)*binSize+vMin+binSize/2.0
	histYa = fltArr(nBins)
	for n=0,n_elements(p_vx)-1 do begin
		ii = (v1x_0[0,5,n]-vMin)/(vMax-vMin)*nBins
		histYa[ii] += p_weight[n]
	endfor

	p=plot(histXa,histYa)


stop
end
