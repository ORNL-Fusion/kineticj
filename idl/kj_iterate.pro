; Iterate kj with rsfwc with file based communiction

pro kj_iterate, jPFile=jPFile, itStartNo=itStartNo, nIterations=nIterations, useAR2=useAR2

    if keyword_set(useAR2) then useAORSA = 1 else useAORSA = 0
	if keyword_set(itStartNo) then itStart=itStartNo else itStart=0
	if keyword_set(nIterations) then nIt=nIterations else nIt=2

    cartesian_offset = 0.0

    KJ_BINARY = '~/code/kineticj/bin/kineticj'
    KJ_BINARY_GC = '~/code/kineticj/bin/kineticj'
 
	cd, current=RootRunDir
    RootRunDir = RootRunDir+'/'
    RunDir0 = 'mpe_0/k_0/'

    if useAORSA then begin
        ar2Dir = 'ar2/'
        ar2RunDir0 = RootRunDir+RunDir0+ar2Dir
        ar2nml = ar2_read_namelist(ar2RunDir0)
        ar2_read_ar2input, ar2RunDir0+'input/ar2Input.nc', ar2=ar2input
        ar2_FileName = 'ar2_kj.nc'
    endif else begin
        RsDir = 'rsfwc/'
        RsRunDir0 = RootRunDir+RunDir0+RsDir
	    RsCfg = kj_read_rsfwc_cfg(RsRunDir0)
        rs_FileName = 'rsfwc_1d.nc'
    endelse

    kj_jP_FileName = 'kj-jp-on-rs-grid.nc'

    kjSpecies = ['spec_D','spec_H','spec_e'] ; The order here MUST match the RSFWC spec order (i.e., electrons last)
    ElectronSpecStr = kjSpecies[2]

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

            LastRunFolder = ThisRunFolder
            ThisMPEDir = 'mpe_'+string(it,format='(i1.1)')+'/'
            ThisPicardDir = 'k_'+string(k,format='(i1.1)')+'/'
            ThisRunFolder = RootRunDir+ThisMPEDir+ThisPicardDir
            ThisRsRunFolder = ThisRunFolder+RsDir

            if k gt 0 then begin
                stop 
                file_delete, ThisRunFolder, /recursive, /allow_nonexistent
                file_copy, LastRunFolder, ThisRunFolder, /recursive 
                file_copy, jGuessFileList[k-1], ThisRsRunFolder
            endif

            ;if not useAORSA then RsCfg['jAmp'] = ((k+1)*jAmpStep)<jAmpMax

			if(k eq 0 and not keyword_set(jPFile) ) then begin
                if useAORSA then begin
                endif else begin
				    RsCfg['kjInput']=0 
				    RsCfg['kj_jP_fileName'] = ''
                endelse
			endif else if(k eq 0 and keyword_set(jPFile) ) then begin
				print, 'Continuing with file ... ', jPFile
                if useAORSA then begin
                endif else begin
				    RsCfg['kjInput']=1 
				    RsCfg['kj_jP_fileName'] = jPFile
                endelse
			endif else begin
                if useAORSA then begin
                endif else begin
				    RsCfg['kjInput']=1
				    RsCfg['kj_jP_fileName'] = kj_jP_fileName
                endelse
			endelse

            if useAORSA then begin
            endif else begin
			    kj_write_rsfwc_cfg, RsCfg, ThisRsRunFolder
            endelse

			cd, ThisRsRunFolder
			spawn, 'idl -quiet run_rsfwc'
			cd, RootRunDir

            for s=0,nS-1 do begin

                kjCfg = kjConfigs[s]        
                ThisKJRunFolder = ThisRunFolder+kjSpecies[s]+'/'

			    kjCfg['input_fName'] = 'data/'+rs_FileName

			    kj_write_kj_cfg, kjCfg, ThisKJRunFolder 

                file_copy, ThisRsRunFolder+rs_FileName, ThisKJRunFolder+'data/', /overwrite

                cd, ThisKJRunFolder
                if StrCmp(kjSpecies[s],ElectronSpecStr) then begin
			        spawn, KJ_BINARY_GC
                endif else begin
			        spawn, KJ_BINARY
                endelse
			    spawn, 'idl -quiet run_kj_plot_current'
            endfor

            cd, RootRunDir+ThisMPEDir+ThisPicardDir

            ; Sum over the kJ species and then create the list below

            print, kjSpecies+'/output/'+kj_jP_FileName
            kj_combine_spec_jp, kjSpecies+'/output/'+kj_jP_FileName, SumFileName = kj_jP_FileName, cartesian_offset = cartesian_offset

			jGuessFileList[k] = RootRunDir+ThisMPEDir+ThisPicardDir+kj_jP_FileName 
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


