function kj_read_cfg, RunDir

    cd, RunDir, current = OldDir
    RunFile = 'kj.cfg'

	openr, lun, runFile, /get_lun
	runFileArray = ''
	line = ''
	while not eof(lun) do begin
		readf, lun, line
		runFileArray = [runFileArray,line]
	endwhile
	free_lun, lun
    runFileArray = runFileArray[1:-1]

    h = orderedHash()
    StrData = StrSplit(RunFileArray,'=,;',/extract)
    for f=0,n_elements(StrData)-1 do begin
        ThisStr = StrData[f]
        if not (strMatch(thisStr,'//*'))[0] then begin
            ThisKey = StrTrim(ThisStr[0],2)
            ThisValue = StrTrim(ThisStr[1],2)
            ; Determine type of value
            i = StRegEx(ThisValue,'^[-+]?[0-9]+$',/bool) ; integer
            d = StRegEx(ThisValue,'^[-+]?[0-9]+\.?[0-9]+$',/bool) ; decimal

            if i eq 1 then begin
                    h = h + orderedHash(ThisKey,fix(ThisValue))
            endif else if d eq 1 then begin
                    h = h + orderedHash(ThisKey,float(ThisValue))
            endif else begin
                    len = StrLen(ThisValue)
                    str = StrMid(ThisValue,1,len-2)
                    h = h + orderedHash(ThisKey,str)
            endelse
        endif

    endfor

    cd, OldDir

	return, h

end 


