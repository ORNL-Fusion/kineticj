function kj_read_rsfwc_cfg, runFile

	;runFile = 'data/rsfwc_input.pro'
	openr, lun, runFile, /get_lun
	runFileArray = ''
	line = ''
	while not eof(lun) do begin
		readf, lun, line
		runFileArray = [runFileArray,line]
	endwhile
	free_lun, lun

	for f=0,n_elements(runFileArray)-1 do begin
		if(strMatch(runFileArray[f],'*kj_jP_fileName*'))then begin
				tmp = strSplit(runFileArray[f],"'",/extract)
				if n_elements(tmp) ge 2 then begin
						kj_jP_fileName=(tmp)[1]
				endif else begin
						kj_jP_fileName=''
				endelse

		endif
		if(strMatch(runFileArray[f],'*runIdent*'))then runIdent=(strSplit(runFileArray[f],"'",/extract))[1]
		if(strMatch(runFileArray[f],'*kjInput*'))then kjInput=fix((strSplit(runFileArray[f],'=',/extract))[1])
		if(strMatch(runFileArray[f],'*jAmp*'))then jAmp=float((strSplit(runFileArray[f],'=',/extract))[1])
	endfor

	cfg = create_struct ( name='rsfwcCfg', $
			'kj_jP_fileName', kj_jP_fileName, $
		   	'runIdent', runIdent, $
			'kjInput', kjInput, $
		    'jAmp', jAmp )

	return, cfg


end



