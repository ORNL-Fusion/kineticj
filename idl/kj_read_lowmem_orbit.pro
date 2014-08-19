pro kj_read_lowmem_orbit

    fileName = 'output/orbit.txt'
    nLines = file_lines(fileName)
    openr, lun, fileName, /get_lun

    offset = 1

    t = fltArr(nLines-offset)
    x = fltArr(nLines-offset)
    y = fltArr(nLines-offset)
    z = fltArr(nLines-offset)

    e1 = complexArr(nLines-offset)
    e2 = complexArr(nLines-offset)
    e3 = complexArr(nLines-offset)

    b1 = complexArr(nLines-offset)
    b2 = complexArr(nLines-offset)
    b3 = complexArr(nLines-offset)

    skip_lun, lun, offset, /lines
    for l=0,nLines-1-offset do begin

        readf, lun, _t, _x, _y, _z, $
                _e1_re, _e1_im, _e2_re, _e2_im, _e3_re, _e3_im, $
                _b1_rb, _b1_im, _b2_rb, _b2_im, _b3_rb, _b3_im

        t[l] = _t
        x[l] = _x
        y[l] = _y
        z[l] = _z

        e1[l] = complex(_e1_re,_e1_im)
        e2[l] = complex(_e2_re,_e2_im)
        e3[l] = complex(_e3_re,_e3_im)

        b1[l] = complex(_b1_rb,_b1_im)
        b2[l] = complex(_b2_rb,_b2_im)
        b3[l] = complex(_b3_rb,_b3_im)


    endfor

    fileName = 'output/orbit_v1.txt'
    nLines = file_lines(fileName)
    openr, lun, fileName, /get_lun

    offset = 1

    v1 = complexArr(nLines-offset)
    v2 = complexArr(nLines-offset)
    v3 = complexArr(nLines-offset)

    skip_lun, lun, offset, /lines
    for l=0,nLines-1-offset do begin

        readf, lun, _t, $
                _v1_re, _v1_im, _v2_re, _v2_im, _v3_re, _v3_im

        v1[l] = complex(_v1_re,_v1_im)
        v2[l] = complex(_v2_re,_v2_im)
        v3[l] = complex(_v3_re,_v3_im)

    endfor

    free_lun, lun
    p=plot3d(x,y,z,aspect_ratio=1.0,aspect_z=1.0)

    p=plot(t,e1,layout=[1,3,1])
    p=plot(t,e2,layout=[1,3,2],/current)
    p=plot(t,e3,layout=[1,3,3],/current)

    p=plot(t,v1,layout=[1,3,1])
    p=plot(t,v2,layout=[1,3,2],/current)
    p=plot(t,v3,layout=[1,3,3],/current)
 
stop
end
