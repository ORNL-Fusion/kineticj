; Iterate kj with rsfwc with file based communiction

pro kj_iterate, jPFile=jPFile, itStartNo=itStartNo, nIterations=nIterations

	if keyword_set(itStartNo) then itStart=itStartNo else itStart=0
	if keyword_set(nIterations) then nIt=nIterations else nIt=20

	cd, current=runDir
	runIdent = file_baseName(runDir)
	rsfwcCfg = kj_read_rsfwc_cfg('data/rsfwc_input.pro')
	kjCfg = kj_read_cfg('kj.cfg')

	for it=itStart,itStart+nIt do begin

		thisIdent = runIdent+'_'+string(it+1,format='(i3.3)')
		lastIdent = runIdent+'_'+string(it+1-1,format='(i3.3)')

		rsfwcCfg.runIdent = thisIdent 

		if(it eq 0 and not keyword_set(jPFile) ) then begin
			rsfwcCfg.kjInput=0 
			rsfwcCfg.kj_jP_fileName = ''
		endif else if(it eq 0 and keyword_set(jPFile) ) then begin
			print, 'Continuing with file ... ', jPFile
			rsfwcCfg.kjInput=1 
			rsfwcCfg.kj_jP_fileName = jPFile
		endif else begin
			rsfwcCfg.kjInput=1
			rsfwcCfg.kj_jP_fileName = 'kj_jP_'+lastIdent+'.nc'
		endelse

		kj_write_rsfwc_cfg, rsfwcCfg, it

		cd, 'data'
		spawn, 'idl run_rsfwc'
		cd, runDir

		kjCfg.eField_fName = 'data/rsfwc_1d_'+rsfwcCfg.runIdent+'.nc'
		kjCfg.runIdent = thisIdent 

		kj_write_kj_cfg, kjCfg, it

		spawn, '~/code/kineticj/bin/kineticj'
		spawn, 'idl run_kj_plot_current'
		spawn, 'cp output/kj_jP_'+thisIdent+'.nc data/'

	endfor

	stop

end


