function kj_read_cfg, runFile

	;runFile = 'kj.cfg'
	openr, lun, runFile, /get_lun
	runFileArray = ''
	line = ''
	while not eof(lun) do begin
		readf, lun, line
		runFileArray = [runFileArray,line]
	endwhile
	free_lun, lun

	for f=0,n_elements(runFileArray)-1 do begin
		if(strMatch(runFileArray[f],'*eField_fName*'))then eField_fName=(strSplit(runFileArray[f],'"',/extract))[1]
		;if(strMatch(runFileArray[f],'*particleList_fName*'))then particleList_fName=(strSplit(runFileArray[f],'"',/extract))[1]
		if(strMatch(runFileArray[f],'*runIdent*'))then runIdent=(strSplit(runFileArray[f],'"',/extract))[1] 
		if(strMatch(runFileArray[f],'*xGridMin*'))then xGridMin=float((strSplit(runFileArray[f],'=,;',/extract))[1])
		if(strMatch(runFileArray[f],'*xGridMax*'))then xGridMax=float((strSplit(runFileArray[f],'=,;',/extract))[1])
		if(strMatch(runFileArray[f],'*nXGrid*'))then nXGrid=fix((strSplit(runFileArray[f],'=,;',/extract))[1])
		if(strMatch(runFileArray[f],'*nRFCycles*'))then nRFCycles=fix((strSplit(runFileArray[f],'=,;',/extract))[1])
		if(strMatch(runFileArray[f],'*nStepsPerCycle*'))then nStepsPerCycle=fix((strSplit(runFileArray[f],'=,;',/extract))[1])
		if(strMatch(runFileArray[f],'*nJpCycles*'))then nJpCycles=fix((strSplit(runFileArray[f],'=,;',/extract))[1])
		if(strMatch(runFileArray[f],'*nJpPerCycle*'))then nJpPerCycle=fix((strSplit(runFileArray[f],'=,;',/extract))[1])
	endfor

	cfg = create_struct ( name='kjCfg', $
			'xGridMin', xGridMin, $
			'xGridMax', xGridMax, $
			'nXGrid', nXGrid, $
			'nRFCycles', nRFCycles, $
			'nStepsPerCycle', nStepsPerCycle, $
			'nJpCycles', nJpCycles, $
			'nJpPerCycle', nJpPerCycle, $
			;'particleList_fName', particleList_fName, $
			'runIdent', runIdent, $
			'eField_fName', eField_fName )

	return, cfg
end 


