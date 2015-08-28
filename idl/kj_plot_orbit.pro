pro kj_plot_orbit

	runFile = './'
	cfg = kj_read_cfg (runFile)	

	; Read particle list

	particleList = cfg.particleList_fName
	cdfId = ncdf_open(particleList)
		ncdf_varget, cdfId, 'vx', p_vx 
		ncdf_varget, cdfId, 'weight', p_weight
		ncdf_varget, cdfId, 'vTh', p_vTh
	nCdf_close, cdfId
	df0_dv = deriv(p_vx,p_weight)

	;vMax = 5e7
	;vMin = -vMax
	;nBins = 1000000L 
	;binSize = (vMax-vMin)/nBins
	;nBins = (vMax-vMin)/binSize
	;histX = fIndGen(nBins)*binSize+vMin+binSize/2.0
	;histY = fltArr(nBins)
	;for n=0,n_elements(p_vx)-1 do begin
	;	ii = (p_vx[n]-vMin)/(vMax-vMin)*nBins
	;	histY[ii] += p_weight[n]
	;endfor

	;iiPlot = where(histY gt 0)

	;histX = histX[iiPlot]
	;histY = histY[iiPlot]
	;plotHist=plot(histX[iiPlot],histY[iiPlot],thick=3,transparency=50)

	fileList = file_search ( 'output/'+cfg.runIdent+'/jP*' )

	spatialPoint = 1 

	cdfId = ncdf_open(fileList[spatialPoint])
		ncdf_varget, cdfId, 'freq', freq 
		ncdf_varget, cdfId, 't', tJ 
		ncdf_varget, cdfId, 'x', x
		ncdf_varget, cdfId, 'j1x', j1x_0 
	nCdf_close, cdfId

	wrf = freq * 2 * !pi

	fileList = file_search ( 'output/orbits*' )

	cdfId = ncdf_open(fileList[0])
		ncdf_varget, cdfId, 't', t
	nCdf_close, cdfId

	nT = n_elements(t)
	nF = n_elements(fileList)

	maxE = 0

	for f=spatialPoint,spatialPoint do begin

		cdfId = ncdf_open(fileList[f])
			print, fileList[f]

			ncdf_varget, cdfId, 'x', x_0
			ncdf_varget, cdfId, 'y', y_0
			ncdf_varget, cdfId, 'z', z_0

			ncdf_varget, cdfId, 'vx', vx_0
			ncdf_varget, cdfId, 'vy', vy_0
			ncdf_varget, cdfId, 'vz', vz_0

			ncdf_varget, cdfId, 't', t_0

			ncdf_varget, cdfId, 'v1x', v1x_0 
			ncdf_varget, cdfId, 'v1y', v1y_0 
			ncdf_varget, cdfId, 'v1z', v1z_0 
stop
			ncdf_varget, cdfId, 'v1x_re', v1x_0_re 
			ncdf_varget, cdfId, 'v1y_re', v1y_0_re 
			ncdf_varget, cdfId, 'v1z_re', v1z_0_re 

			ncdf_varget, cdfId, 'v1x_im', v1x_0_im 
			ncdf_varget, cdfId, 'v1y_im', v1y_0_im 
			ncdf_varget, cdfId, 'v1z_im', v1z_0_im 

			ncdf_varget, cdfId, 'e1_x', e1x_0 
			ncdf_varget, cdfId, 'e1_y', e1y_0 
			ncdf_varget, cdfId, 'e1_z', e1z_0 

			ncdf_varget, cdfId, 'e1_x_re', e1x_0_re 
			ncdf_varget, cdfId, 'e1_y_re', e1y_0_re 
			ncdf_varget, cdfId, 'e1_z_re', e1z_0_re 

			ncdf_varget, cdfId, 'e1_x_im', e1x_0_im 
			ncdf_varget, cdfId, 'e1_y_im', e1y_0_im 
			ncdf_varget, cdfId, 'e1_z_im', e1z_0_im 

			e1x_0 = complex(e1x_0_re,e1x_0_im)
			e1y_0 = complex(e1y_0_re,e1y_0_im)
			e1z_0 = complex(e1z_0_re,e1z_0_im)

			v1x_0 = complex(v1x_0_re,v1x_0_im)
			v1y_0 = complex(v1y_0_re,v1y_0_im)
			v1z_0 = complex(v1z_0_re,v1z_0_im)


		nCdf_close,	cdfId 

		if(max(e1x_0) gt maxE) then maxE = max(e1x_0)

	endfor

	nSteps = n_elements(v1x_0[*,0,0])
	nJp = n_elements(v1x_0[0,*,0])
	nP = n_elements(v1x_0[0,0,*])

	pNum = 100 

	phaseNum = 0
	p=plot(t_0*freq,v1x_0[*,phaseNum,pNum],title='v1x_0')
	p=plot(t_0*freq,e1x_0[*,pNum],title='e1x_0')
	kPar= 12;22.11
	vPhs = wrf/kPar

    p = plot3d(x_0[*,pNum],y_0[*,pNum],z_0[*,pNum])

    SubSampleFac = 5
    nT = n_elements(t_0)/SubSampleFac
    nV = n_elements(vx_0[0,*])/SubSampleFac
    ThisX = congrid(t_0*freq,nT)
    ThisY = congrid(reform(vx_0[0,*]/vPhs),nV)

    ThisSurf = congrid(e1x_0,nT,nV)
	s=surface(ThisSurf, ThisX, ThisY,depth_cue = [0,2], /perspective)

    ThisSurf = congrid(reform(v1x_0[*,0,*],nSteps,nP),nT,nV)
	s=surface(ThisSurf, ThisX, ThisY,depth_cue = [0,2], /perspective)

    ThisSurf = congrid(-reform(v1x_0[*,0,*]*transpose(rebin(df0_dv,nP,nSteps)),nSteps,nP),nT,nV)
	s=surface(ThisSurf, ThisX, ThisY,depth_cue = [0,2], /perspective, title='Re[f1]')
	
    ThisSurf = imaginary(congrid(-reform(v1x_0[*,0,*]*transpose(rebin(df0_dv,nP,nSteps)),nSteps,nP),nT,nV))
	s=surface(ThisSurf, ThisX, ThisY,depth_cue = [0,2], /perspective, title='Im[f1]')


	p=plot( vx_0[0,*]/vPhs, v1x_0[0,0,*],symbol="s")

	print, 'This particle v/vPhs: ', vx_0[0,pNum]/vPhs

	;p=plot(p_vx,p_weight,thick=3,transparency=50,layout=[9,2,0+1],axis_style=0)
	;p=plot(p_vx+v1x_0[nSteps-2,0*nJp/9,*],p_weight,/over,color='blue',thick=2,transparency=30)
	;for pp=1,8 do begin
	;	p=plot(p_vx,p_weight,thick=3,transparency=50,layout=[9,2,pp+1],/current,axis_style=0)
	;	p=plot(p_vx+v1x_0[nSteps-2,pp*nJp/9,*],p_weight,/over,color='blue',thick=2,transparency=30)
	;endfor

	p=plot(tJ,j1x_0,thick=3,color='red',transparency=20,$
			position=[0.02,0.05,0.98,0.45],font_size=8)
stop
end
