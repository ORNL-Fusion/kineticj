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

	phaseOffSet = !pi/2

	yr = max(j1x)*1.5
	p=plot(t,j1x[*,0],/noData,yRange=[-yr,yr])
	for f=0,n_elements(xF)-1 do begin
		p=plot(t,j1x[*,f],/over)
		iiHere = where(abs(r-xF[f]) eq min(abs(r-xF[f])))
		tmp = (jr_re[iiHere])[0]*cos(wrf*t+phaseOffSet)+(jr_im[iiHere])[0]*sin(wrf*t+phaseOffSet)
		p=plot(t,tmp, /over,color='blue')
	endfor

	p=plot(t,j1x[*,0],/noData,yRange=[0,3])
	for f=0,n_elements(xF)-1 do begin
		iiHere = where(abs(r-xF[f]) eq min(abs(r-xF[f])))
		tmp = (jr_re[iiHere])[0]*cos(wrf*t+phaseOffSet)+(jr_im[iiHere])[0]*sin(wrf*t+phaseOffSet)
		p=plot(t,tmp/j1x[*,f],/over,lineStyle='none',symbol='Plus')
	endfor

	for i=0,n_elements(t)-1 do begin	
		plot, r,(jr_re*cos(wrf*t[i]+phaseOffSet)+jr_im*sin(wrf*t[i]+phaseOffSet)),$
				yRange=[-100,100], xRange=[9.5,10.5]
		;plots, r,(jAr_re*cos(wrf*t[i])+jAr_im*sin(wrf*t[i]))/10
		;for f=0,nF-1 do begin
			plots, xF, j1x[i,*]/(2*!pi),psym=-4
		;endfor
		wait, 0.2
	endfor


stop
end
