pro kj_write_kj_cfg, cfg, RunDir

    fName = 'kj.cfg'
    file_delete, RunDir+fName
	openw, lun, RunDir+fName, /get_lun 	

        keys = cfg.keys()
        for f=0,n_elements(cfg)-1 do begin
                key = keys[f]
                value = cfg[keys[f]]
                if size(value,/type) eq 7 then begin ; just wrap " " around string types
                    printf, lun, key, ' = ', '"',value,'";'
                endif else if size(value,/type) eq 4 then begin ; floats 
                    printf, lun, key, ' = ', strTrim( string(value,format='(f24.12)') ,2),';'
                endif else if size(value,/type) eq 5 then begin ; doubles
                    printf, lun, key, ' = ', strTrim( string(value,format='(f24.12)') ,2),';'
                endif else if size(value,/type) eq 2 then begin ; integers 
                    printf, lun, key, ' = ', strTrim( string(value) ,2),';'
                endif else begin
                    printf, lun, key, ' = ', value,';'
                endelse
        endfor

	free_lun, lun

end


