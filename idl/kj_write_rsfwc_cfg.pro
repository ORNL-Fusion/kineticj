pro kj_write_rsfwc_cfg, cfg, it

	itStr = string(it-1,format='(i+4.3)')
	;if(it gt 0) then $
		spawn, 'mv data/rsfwc_input.pro data/rsfwc_input'+itStr+'.pro'
	openw, lun, 'data/rsfwc_input.pro', /get_lun 	

		printf, lun, 'kjInput = ', strTrim(string(cfg.kjInput),2)
		printf, lun, 'kj_jP_fileName = ', "'",cfg.kj_jP_fileName,"'"
		printf, lun, 'runIdent = ', "'",cfg.runIdent,"'"

	free_lun, lun

end
	
