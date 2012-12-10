pro kj_write_kj_cfg, cfg, it

	itStr = string(it+1-1,format='(i+4.3)')
	;if(it gt 0) then $
		spawn, 'mv kj.cfg kj'+itStr+'.cfg'
	openw, lun, 'kj.cfg', /get_lun 	
		printf, lun, 'xGridMin = ', strTrim(string(cfg.xGridMin),2),';'
		printf, lun, 'xGridMax = ', strTrim(string(cfg.xGridMax),2),';'
		printf, lun, 'nXGrid = ', strTrim(string(cfg.nXGrid),2),';'
		printf, lun, 'nRFCycles = ', strTrim(string(cfg.nRFCycles),2),';' 
		printf, lun, 'nStepsPerCycle = ', strTrim(string(cfg.nStepsPerCycle),2),';' 
		printf, lun, 'nJpCycles = ', strTrim(string(cfg.nJpCycles),2),';' 
		printf, lun, 'nJpPerCycle = ', strTrim(string(cfg.nJpPerCycle),2),';' 
		printf, lun, 'particleList_fName = ', '"',cfg.particleList_fName,'"',';'
		printf, lun, 'eField_fName = ', '"',cfg.eField_fName,'"',';'
		printf, lun, 'runIdent = ', '"',cfg.runIdent,'"',';'
	free_lun, lun

end


