pro kj_plot_jp_iteration_range, noOverlay=noOverlay, noIterations=noIterations

	itNo = 5
	iterationDir = 'output/mpe_it_'+string(itNo,format='(i3.3)')

	overlayFile = !null
	overlayColors = !null

	overlayFile = [iterationDir+'/mpe_extrapolated_jP.nc']
	overlayColors = [overlayColors,'blue']

	overlayFile = [overlayFile,'data/kj_aorsa_1d.nc']
	overlayColors = [overlayColors,'red']

	;overlayFile = [overlayFile,'data/rsfwc_1d_0.5keV_this_001.nc']
	;overlayColors = [overlayColors,'black']


	fileList = file_search ( iterationDir+'/kj_jP*')
	nFiles = n_elements(fileList)

	cdfId = ncdf_open(fileList[0])
	ncdf_varget, cdfId, 'r', r 
	ncdf_close, cdfId
	nX = n_elements(r)

	jP_r_all = complexarr(nX,nFiles)

	for f=0,nFiles-1 do begin

		cdfId = ncdf_open(fileList[f])
		ncdf_varget, cdfId, 'r', r 
		ncdf_varget, cdfId, 'r_', r_ 
		ncdf_varget, cdfId, 'jP_r_re', jP_r_re
		ncdf_varget, cdfId, 'jP_r_im', jP_r_im
		ncdf_varget, cdfId, 'jP_p_re', jP_p_re
		ncdf_varget, cdfId, 'jP_p_im', jP_p_im
		ncdf_varget, cdfId, 'jP_z_re', jP_z_re
		ncdf_varget, cdfId, 'jP_z_im', jP_z_im
		ncdf_close, cdfId

		jP_r_all[*,f] = complex(jP_r_re,jP_r_im)

	endfor

	xMin = 97.0
	xMax = 103.0
	pr=plot(r,real_part(jP_r_all[*,0]),$
			color='black',thick=1,buffer=1, dim=[1200,400],transparency=50,$
			xTitle='x [m]', yTitle=' jP_1 [A/m^-1]',xRange=[xMin,xMax],$
			title='KineticJ jP_1 Iterations',noData=noIterations)
	for f=1,nFiles-1,1 do !null=plot(r,real_part(jP_r_all[*,f]),/over,color='black',thick=1,transparency=50,nodata=noIterations)

	if not keyword_set(noOverlay) then begin
	nOverlays = n_elements(overlayFile)
	for n=0,nOverlays-1 do begin
		cdfId = ncdf_open(overlayFile[n])
		ncdf_varget, cdfId, 'r', r 
		ncdf_varget, cdfId, 'jP_r_re', jP_r_re
		ncdf_varget, cdfId, 'jP_r_im', jP_r_im
		ncdf_varget, cdfId, 'jP_p_re', jP_p_re
		ncdf_varget, cdfId, 'jP_p_im', jP_p_im
		ncdf_varget, cdfId, 'jP_z_re', jP_z_re
		ncdf_varget, cdfId, 'jP_z_im', jP_z_im
		ncdf_close, cdfId

		jP_r_overlay = complex(jP_r_re,jP_r_im)

		!null=plot(r,real_part(jP_r_overlay),/over,$
				transparency=50,color=overlayColors[n],thick=3)
	endfor
	endif
	pr.save, 'jPr_mpe.eps'
	stop
end
