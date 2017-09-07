; Iterate kj with rsfwc with file based communiction

pro kj_iterate, $
        kjDeltaFile=kjDeltaFile, $
        continueFromLastMPE=_continueFromLastMPE, $
        nIterations=nIterations, $
        useAR2=useAR2, $
        useKJStix = _useKJStix

    kjIterationStateFileName = 'kj-iteration-state.sav'

    if keyword_set(useAR2) then useAORSA = 1 else useAORSA = 0
	if keyword_set(itStartNo) then itStart=itStartNo else itStart=0
	if keyword_set(nIterations) then nIt=nIterations else nIt=2
    if keyword_set(_useKJStix) then useKJStix = _useKJStix else useKJStix = 1
    if keyword_set(_continueFromLastMPE) then begin
        restore, kjIterationStateFileName
        itStart = it+1 
    endif else begin
        itStart = 0
    endelse

    useKJFull = 0
    if not useKJStix then useKJFull = 1

    cartesian_offset = 0.0

    KJ_BINARY = '~/code/kineticj/bin/kineticj'
    KJ_BINARY_GC = '~/code/kineticj/bin/kineticj'
 
	cd, current=RootPath
    RootPath = RootPath+'/'
    RunDir0 = 'mpe_0/k_0/'

    if useAORSA then begin
        ar2Dir = 'ar2/'
        ar2RunDir0 = RootPath+RunDir0+ar2Dir
        ar2nml = ar2_read_namelist(ar2RunDir0)
        ar2_read_ar2input, ar2RunDir0+'input/ar2Input.nc', ar2=ar2input
        ar2_FileName = 'ar2_kj.nc'
    endif else begin
        RsDir = 'rs/'
        RsRunDir0 = RootPath+RunDir0+RsDir
	    RsCfg = kj_read_rsfwc_cfg(RsRunDir0)
        rs_FileName = 'rs_solution.nc'
    endelse

    kjDeltaFileName = 'kj-delta.nc'

    if useKJFull then begin

        kjSpecies = ['spec_D','spec_H','spec_e'] 
        ; The order here MUST match the RSFWC spec order (i.e., electrons last)

        ElectronSpecStr = kjSpecies[2]

        kjConfigs = []
        nS = n_elements(kjSpecies)
        for s=0,nS-1 do begin
	        kjCfg = kj_read_cfg(RunDir0+kjSpecies[s]+'/')
            kjConfigs = [kjConfigs,kjCfg]
        endfor

    endif

	jAmpMax = 1.0
	jAmpStep = 1.0 

	nk = 4 
	picardDeltaFileList = strArr(nk)

    if nk gt 10 then stop ; Filenames cannot handle this
    if nIt+itStart gt 10 then stop ; Filenames cannot handle this either

    ThisMPEPath = ''
    ThisPicardPath = ''

	for it=itStart,itStart+nIt-1 do begin

        LastMPEPath = ThisMPEPath 
        ThisMPEDir  = 'mpe_'+string(it,format='(i1.1)')+'/'
        ThisMPEPath = RootPath+ThisMPEDir

        if it gt 0 then begin
            file_copy, LastMPEPath, ThisMPEPath, /recursive 
        endif

		for k=0,nk-1 do begin

			print, 'MPE Iteration: ', string(it,format='(i3.3)'),' of ', $
					string(itStart+nIt-1,format='(i3.3)'), ' and Picard iteration: ', $
					string(k,format='(i3.3)'), ' of ', string(nk-1,format='(i3.3)')

            LastPicardPath = ThisPicardPath
            ThisPicardDir = 'k_'+string(k,format='(i1.1)')+'/'
            ThisPicardPath = ThisMPEPath+ThisPicardDir
            ThisRSPath = ThisPicardPath+RsDir

            if k gt 0 then begin
                file_delete, ThisPicardPath, /recursive, /allow_nonexistent
                file_copy, LastPicardPath, ThisPicardPath, /recursive 
                file_copy, picardDeltaFileList[k-1], ThisRSPath, /overWrite
            endif

			if(k eq 0 and not keyword_set(kjDeltaFile) ) then begin
                if useAORSA then begin
                endif else begin
				    RsCfg['kjInput']=0 
				    RsCfg['kjDeltaFileName'] = ''
                endelse
			endif else if(k eq 0 and keyword_set(kjDeltaFile) ) then begin
				print, 'Continuing with file ... ', kjDeltaFile
                if useAORSA then begin
                endif else begin
				    RsCfg['kjInput']=1 
				    RsCfg['kjDeltaFileName'] = kjDeltaFile
                endelse
			endif else begin
                if useAORSA then begin
                endif else begin
				    RsCfg['kjInput']=1
				    RsCfg['kjDeltaFileName'] = kjDeltaFileName
                endelse
			endelse

            if useAORSA then begin
            endif else begin
			    kj_write_rsfwc_cfg, RsCfg, ThisRSPath
            endelse

			cd, ThisRSPath
			spawn, 'idl -quiet run_rs'
			cd, RootPath


            if UseKJStix then begin

                ; Evaluate kinetic-j using the stix method

                cd, ThisRSPath

                kj_stix_current, overPlotSolution = 1, useRS=1, hot=1, $
                        jr = this_jr, jt = this_jt, jz = this_jz, $
                        kjDeltaFileName = kjDeltaFileName, $
                        referenceSolutionDir = rsRunDir0
                
                file_copy, kjDeltaFileName, ThisPicardPath, /overWrite

                cd, RootPath

            endif else begin
            
                ; Evaluate kinetic-j using the numeric integrals

                for s=0,nS-1 do begin

                    kjCfg = kjConfigs[s]        
                    ThisKJRunFolder = ThisPicardPath+kjSpecies[s]+'/'

			        kjCfg['input_fName'] = 'data/'+rs_FileName

			        kj_write_kj_cfg, kjCfg, ThisKJRunFolder 

                    file_copy, ThisRSPath+rs_FileName, ThisKJRunFolder+'data/', /overwrite

                    cd, ThisKJRunFolder
                    if StrCmp(kjSpecies[s],ElectronSpecStr) then begin
			            spawn, KJ_BINARY_GC
                    endif else begin
			            spawn, KJ_BINARY
                    endelse
			        spawn, 'idl -quiet run_kj_plot_current'
                endfor

                cd, RootPath+ThisMPEDir+ThisPicardDir

                ; Sum over the kJ species and then create the list below

                print, kjSpecies+'/output/'+kjDeltaFileName
                kj_combine_spec_jp, kjSpecies+'/output/'+kjDeltaFileName, SumFileName = kjDeltaFileName, cartesian_offset = cartesian_offset

            endelse

			picardDeltaFileList[k] = RootPath+ThisMPEDir+ThisPicardDir+kjDeltaFileName 

		endfor

		; Read the previous n guesses and apply vector extrapolation

		jrDelta_picard = !null
        jtDelta_picard = !null
        jzDelta_picard = !null

		for k=0,nk-1 do begin

			cdfId = ncdf_open(picardDeltaFileList[k])

				print, picardDeltaFileList[k]

				ncdf_varget, cdfId, 'r', r 

				ncdf_varget, cdfId, 'jP_r_re', jP_r_re
				ncdf_varget, cdfId, 'jP_r_im', jP_r_im
				ncdf_varget, cdfId, 'jP_t_re', jP_t_re
				ncdf_varget, cdfId, 'jP_t_im', jP_t_im
				ncdf_varget, cdfId, 'jP_z_re', jP_z_re
				ncdf_varget, cdfId, 'jP_z_im', jP_z_im

			ncdf_close, cdfId

			nX = n_elements(r)

			jrDelta_picard = [[jrDelta_picard],[complex(jP_r_re,jP_r_im)]]
			jtDelta_picard = [[jtDelta_picard],[complex(jP_t_re,jP_t_im)]]
			jzDelta_picard = [[jzDelta_picard],[complex(jP_z_re,jP_z_im)]]

		endfor

		jrDelta_mpe = kj_mpe(jrDelta_picard)
		jrDelta_mpe_re = real_part(jrDelta_mpe)
		jrDelta_mpe_im = imaginary(jrDelta_mpe)

		jtDelta_mpe = kj_mpe(jtDelta_picard)
		jtDelta_mpe_re = real_part(jtDelta_mpe)
		jtDelta_mpe_im = imaginary(jtDelta_mpe)

		jzDelta_mpe = kj_mpe(jzDelta_picard)
		jzDelta_mpe_re = real_part(jzDelta_mpe)
		jzDelta_mpe_im = imaginary(jzDelta_mpe)

		print, 'Writing vector extrapolated jP to file ... ', picardDeltaFileList[0]

        jP_MPE_FileName = ThisMPEDir + 'kjDelta_mpe.nc'

		cdfId = ncdf_open(picardDeltaFileList[0],/write)

			jP_r_re_id = nCdf_varid(cdfId, 'jP_r_re')
	        jP_r_im_id = nCdf_varid(cdfId, 'jP_r_im')
		    jP_t_re_id = nCdf_varid(cdfId, 'jP_t_re')
			jP_t_im_id = nCdf_varid(cdfId, 'jP_t_im')
			jP_z_re_id = nCdf_varid(cdfId, 'jP_z_re')
			jP_z_im_id = nCdf_varid(cdfId, 'jP_z_im')
	
			nCdf_varPut, cdfId, jP_r_re_id, jrDelta_mpe_re
			nCdf_varPut, cdfId, jP_r_im_id, jrDelta_mpe_im
	    	nCdf_varPut, cdfId, jP_t_re_id, jtDelta_mpe_re
			nCdf_varPut, cdfId, jP_t_im_id, jtDelta_mpe_im
		 	nCdf_varPut, cdfId, jP_z_re_id, jzDelta_mpe_re
			nCdf_varPut, cdfId, jP_z_im_id, jzDelta_mpe_im
		
		nCdf_close, cdfId

        file_copy, picardDeltaFileList[0], jP_MPE_FileName, /overWrite

        range = max(abs([jrDelta_mpe_re]))
		pr=plot(jrDelta_mpe_re,color='b',thick=6,buffer=1, dim=[1200,1200],transparency=50,yRange=[-1,1]*range,layout=[1,3,1])
		for k=0,nk-1 do !null=plot(real_part(jrDelta_Picard[*,k]),/over,transparency=50)
        range = max(abs([jtDelta_mpe_re]))
		pr=plot(jtDelta_mpe_re,color='b',thick=6,transparency=50,yRange=[-1,1]*range,layout=[1,3,2],/current)
		for k=0,nk-1 do !null=plot(real_part(jtDelta_Picard[*,k]),/over,transparency=50)
        range = max(abs([jzDelta_mpe_re]))
		pr=plot(jzDelta_mpe_re,color='b',thick=6,transparency=50,yRange=[-1,1]*range,layout=[1,3,3],/current)
		for k=0,nk-1 do !null=plot(real_part(jzDelta_Picard[*,k]),/over,transparency=50)

        range = max(abs([jrDelta_mpe_im]))
		pi=plot(jrDelta_mpe_im,color='b',thick=6,buffer=1, dim=[1200,1200],transparency=50,yRange=[-1,1]*range,layout=[1,3,1])
		for k=0,nk-1 do !null=plot(real_part(jrDelta_Picard[*,k]),/over,transparency=50)
        range = max(abs([jtDelta_mpe_im]))
		pi=plot(jtDelta_mpe_im,color='b',thick=6,transparency=50,yRange=[-1,1]*range,layout=[1,3,2],/current)
		for k=0,nk-1 do !null=plot(real_part(jtDelta_Picard[*,k]),/over,transparency=50)
        range = max(abs([jzDelta_mpe_im]))
		pi=plot(jzDelta_mpe_im,color='b',thick=6,transparency=50,yRange=[-1,1]*range,layout=[1,3,3],/current)
        for k=0,nk-1 do !null=plot(real_part(jzDelta_Picard[*,k]),/over,transparency=50)

		pr.save, ThisMPEPath+'jPr.png'
		pi.save, ThisMPEPath+'jPi.png'

		pr.save, ThisMPEPath+'jPr.eps'
		pi.save, ThisMPEPath+'jPi.eps'

		kjDeltaFile = file_baseName(picardDeltaFileList[0])

        cd, RootPath

        save, ThisMPEPath, ThisPicardPath, kjDeltaFile, it, fileName = kjIterationStateFileName

	endfor

	stop

end


