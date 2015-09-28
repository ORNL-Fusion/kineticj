pro kj_read_f1

    fileName = 'output/f1.txt'
    nLines = file_lines(fileName)
    openr, lun, fileName, /get_lun

    offset = 1

    vx = fltArr(nLines-offset)
    vy = fltArr(nLines-offset)
    vz = fltArr(nLines-offset)

    f1 = complexArr(nLines-offset)

    skip_lun, lun, offset, /lines
    for l=0,nLines-1-offset do begin

        readf, lun, _vx, _vy, _vz, _f1re, _f1im 

        vx[l] = _vx
        vy[l] = _vy
        vz[l] = _vz

        f1[l] = complex(_f1re,_f1im)

    endfor
	free_lun, lun
stop
end
