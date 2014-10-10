; Iterate kj with rsfwc with file based communiction

pro kj_iterate, jPFile=jPFile, itStartNo=itStartNo, nIterations=nIterations

	if keyword_set(itStartNo) then itStart=itStartNo else itStart=0
	if keyword_set(nIterations) then nIt=nIterations else nIt=2

    KJ_BINRARY = '~/code/kineticj/bin/kineticj'
    KJ_BINARY_GC = '~/code/kineticj/bin/kineticj-gc'
 
	cd, current=RootRunDir
    RootRunDir = RootRunDir+'/'
    RunDir0 = 'mpe_0/k_0/'
    RsDir = 'rsfwc/'
    RsRunDir0 = RootRunDir+RunDir0+RsDir
	RsCfg = kj_read_rsfwc_cfg(RsRunDir0)

    kj_jP_FileName = 'kj_jp.nc'
    rs_FileName = 'rsfwc_1d.nc'

    kjSpecies = ['spec_e','spec_D','spec_H']
    ElectronSpecStr = kjSpecies[0]

    kjConfigs = []
    nS = n_elements(kjSpecies)
    for s=0,nS-1 do begin
	    kjCfg = kj_read_cfg(RunDir0+kjSpecies[s]+'/')
        kjConfigs = [kjConfigs,kjCfg]
    endfor

	jAmpMax = 1.0
	jAmpStep = 1.0 

	nk = 10 
	jGuessFileList = strArr(nk)

    if nk gt 10 then stop ; Filenames cannot handle this
    if nIt+itStart gt 10 then stop ; Filenames cannot handle this either

    ThisRunFolder = ''

	for it=itStart,itStart+nIt-1 do begin

		for k=0,nk-1 do begin

			print, 'MPE Iteration: ', string(it,format='(i3.3)'),' of ', $
					string(itStart+nIt-1,format='(i3.3)'), ' and Picard iteration: ', $
					string(k,format='(i3.3)'), ' of ', string(nk-1,format='(i3.3)')

            LastRunFoler = ThisRunFolder
            ThisRunFolder = RootRunDir+'/mpe_'+string(it,format='(i1.1)')+'/k_'+string(k,format='(i1.1)')+'/'
            ThisRsRunFolder = ThisRunFolder+RsDir

            if k gt 0 then begin
                stop 
                file_delete, ThisRunFolder, /recursive, /allow_nonexistent
                file_copy, LastRunFolder, ThisRunFolder, /recursive 
            endif

			RsCfg['jAmp'] = ((k+1)*jAmpStep)<jAmpMax

			if(k eq 0 and not keyword_set(jPFile) ) then begin
				RsCfg['kjInput']=0 
				RsCfg['kj_jP_fileName'] = ''
			endif else if(k eq 0 and keyword_set(jPFile) ) then begin
				print, 'Continuing withh file ... ', jPFile
				RsCfg['kjInput']=1 
				RsCfg['kj_jP_fileName'] = jPFile
			endif else begin
				RsCfg['kjInput']=1
				RsCfg['kj_jP_fileName'] = 'kj_jP.nc'
			endelse

			kj_write_rsfwc_cfg, RsCfg, ThisRsRunFolder

			cd, ThisRsRunFolder
			spawn, 'idl -quiet run_rsfwc'
			cd, RootRunDir

            for s=0,nS-1 do begin

                kjCfg = kjConfigs[s]        
                ThisKJRunFolder = ThisRunFolder+kjSpecies[s]+'/'

			    kjCfg['eField_fName'] = 'data/'+rs_FileName

			    kj_write_kj_cfg, kjCfg, ThisKJRunFolder 

                file_copy, ThisRsRunFolder+rs_FileName, ThisKJRunFolder+'data/'

                cd, ThisKJRunFolder
                if StrCmp(kjSpecies[s],ElectronSpecStr) then begin
			        spawn, KJ_BINARY_GC
                endif else begin
			        spawn, KJ_BINARY
               endelse
			    spawn, 'idl -quiet run_kj_plot_current'
            endfor

            ; Sum over the kJ species and then create the list below

			;jGuessFileList[k] = ThisKJRunFolder+'output/kj_jP_'+kjCfg.runIdent+'.nc'
stop
		endfor

		; Read the previous n guesses and apply vector extrapolation

		mpe_it_dir = 'output/mpe_it_'+string(it,format='(i3.3)')
		spawn, 'mkdir  ' + mpe_it_dir

		jGuess = !null

		for k=0,nk-1 do begin

			cdfId = ncdf_open(jGuessFileList[k])

				print, jGuessFileList[k]
				spawn, 'cp '+jGuessFileList[k]+' '+mpe_it_dir+'/'

				ncdf_varget, cdfId, 'r', r 
				ncdf_varget, cdfId, 'r_', r_ 

				ncdf_varget, cdfId, 'jP_r_re', jP_r_re
				ncdf_varget, cdfId, 'jP_r_im', jP_r_im
				ncdf_varget, cdfId, 'jP_p_re', jP_p_re
				ncdf_varget, cdfId, 'jP_p_im', jP_p_im
				ncdf_varget, cdfId, 'jP_z_re', jP_z_re
				ncdf_varget, cdfId, 'jP_z_im', jP_z_im

			ncdf_close, cdfId

			nX = n_elements(r)

			jGuess = [[jGuess],[complex(jP_r_re,jP_r_im)]]

		endfor

		x = jGuess
		_k = n_elements(x[0,*])

		;s_re = kj_mpe(real_part(x))
		;s_im = kj_mpe(imaginary(x))

		;s = complex(s_re,s_im)
	
		s = kj_mpe(x)
		s_re = real_part(s)
		s_im = imaginary(s)

		spline_sigma = 0.01
		s_ = complex(spline(r,s_re,r_,spline_sigma),spline(r,s_im,r_,spline_sigma))
		s_re_ = real_part(s_)
		s_im_ = imaginary(s_)


		print, 'Writing vector extrapolated jP to file ... ', jGuessFileList[0]
		cdfId = ncdf_open(jGuessFileList[0],/write)

			jP_r_re_id = nCdf_varid(cdfId, 'jP_r_re')
			jP_r_im_id = nCdf_varid(cdfId, 'jP_r_im')
			jP_r_re_id_ = nCdf_varid(cdfId, 'jP_r_re_')
			jP_r_im_id_ = nCdf_varid(cdfId, 'jP_r_im_')
	
			nCdf_varPut, cdfId, jP_r_re_id, s_re
			nCdf_varPut, cdfId, jP_r_im_id, s_im
			nCdf_varPut, cdfId, jP_r_re_id_, s_re_
			nCdf_varPut, cdfId, jP_r_im_id_, s_im_
	
		nCdf_close, cdfId

		spawn, 'cp '+jGuessFileList[0]+' '+mpe_it_dir+'/mpe_extrapolated_jP.nc'

		pr=plot(s_re,color='b',thick=6,buffer=1, dim=[1200,400],transparency=50)
		for k=0,nk-1 do !null=plot(real_part(jGuess[*,k]),/over,transparency=50)

		pi=plot(s_im,color='b',thick=6,buffer=1, dim=[1200,400],transparency=50)
		for k=0,nk-1 do !null=plot(imaginary(jGuess[*,k]),/over,transparency=50)

		pr.save, 'jPr.png'
		pi.save, 'jPi.png'

		pr.save, 'jPr.eps'
		pi.save, 'jPi.eps'

		jPFile = file_baseName(jGuessFileList[0])

	endfor

	stop

end


