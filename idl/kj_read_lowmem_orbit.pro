pro kj_read_lowmem_orbit

    fileName = 'output/orbit.txt'
    nLines = file_lines(fileName)
    openr, lun, fileName, /get_lun

    offset = 2

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
	free_lun, lun

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

    fileName = 'output/orbit_e1_dot_grad_df0_dv.txt'
    nLines = file_lines(fileName)
    openr, lun, fileName, /get_lun

    offset = 1

    e1_dot_grad_1 = complexArr(nLines-offset)
    e1_dot_grad_2 = complexArr(nLines-offset)
    e1_dot_grad_3 = complexArr(nLines-offset)

    skip_lun, lun, offset, /lines
    for l=0,nLines-1-offset do begin

        readf, lun, _t, $
                _v1_re, _v1_im, _v2_re, _v2_im, _v3_re, _v3_im

        e1_dot_grad_1[l] = complex(_v1_re,_v1_im)
        e1_dot_grad_2[l] = complex(_v2_re,_v2_im)
        e1_dot_grad_3[l] = complex(_v3_re,_v3_im)

    endfor
	free_lun, lun

    fileName = 'output/df0dv.txt'
    nLines = file_lines(fileName)
    openr, lun, fileName, /get_lun

    offset = 1

    vx = fltArr(nLines-offset)
    vy = fltArr(nLines-offset)
    vz = fltArr(nLines-offset)

    vAlp = fltArr(nLines-offset)
    vBet = fltArr(nLines-offset)
    vPar = fltArr(nLines-offset)
    vPer = fltArr(nLines-offset)
    vPhs = fltArr(nLines-offset)

    df0dv_x = fltArr(nLines-offset)
    df0dv_y = fltArr(nLines-offset)
    df0dv_z = fltArr(nLines-offset)

    skip_lun, lun, offset, /lines
    for l=0,nLines-1-offset do begin

        readf, lun, _t, $
                _vx, _vy, _vz, _vAlp, _vBet, _vPar, _vPer, _vPhs, _df0dv_x, _df0dv_y, _df0dv_z

        vx[l] = _vx
        vy[l] = _vy
        vz[l] = _vz

        vAlp[l] = _vAlp
        vBet[l] = _vBet
        vPar[l] = _vPar
        vPer[l] = _vPer
        vPhs[l] = _vPhs

        df0dv_x[l] = _df0dv_x
        df0dv_y[l] = _df0dv_y
        df0dv_z[l] = _df0dv_z

    endfor
	free_lun, lun


    p=plot3d(x,y,z,aspect_ratio=1.0,aspect_z=1.0)

	fs = 12
    p=plot(t,e1,layout=[1,3,1], title="e1(t')", font_size=fs)
    p=plot(t,imaginary(e1),/over,color='r')
    p=plot(t,e2,layout=[1,3,2],/current, font_size=fs)
    p=plot(t,imaginary(e2),/over,color='r')
    p=plot(t,e3,layout=[1,3,3],/current, font_size=fs)
    p=plot(t,imaginary(e3),/over,color='r')

    p=plot(t,e1_dot_grad_1,layout=[1,3,1], title="e1 . gradv f0 (t')", font_size=fs)
    p=plot(t,e1_dot_grad_2,layout=[1,3,2],/current, font_size=fs)
    p=plot(t,e1_dot_grad_3,layout=[1,3,3],/current, font_size=fs)

    p=plot(t,v1,layout=[1,3,1], title='dv1(t)', font_size=fs)
    p=plot(t,v2,layout=[1,3,2],/current, font_size=fs)
    p=plot(t,v3,layout=[1,3,3],/current, font_size=fs)

    p=plot(t,vx,layout=[1,3,1], title='vx,vy,vz(t)', font_size=fs)
    p=plot(t,vy,layout=[1,3,2],/current, font_size=fs)
    p=plot(t,vz,layout=[1,3,3],/current, font_size=fs)

    p=plot(t,vAlp,layout=[1,3,1], title='vAlp,vBet,vPar(t)', font_size=fs)
    p=plot(t,vBet,layout=[1,3,2],/current, font_size=fs)
    p=plot(t,vPar,layout=[1,3,3],/current, font_size=fs)

    p=plot(t,vPar,layout=[1,3,1], title='vPar,vPer,vPhs(t)', font_size=fs)
    p=plot(t,vPer,layout=[1,3,2],/current, font_size=fs)
    p=plot(t,vPhs,layout=[1,3,3],/current, font_size=fs)

    p=plot(t,df0dv_x,layout=[1,3,1], title='df0dv_x,df0dv_y,df0dv_z(t)', font_size=fs)
    p=plot(t,df0dv_y,layout=[1,3,2],/current, font_size=fs)
    p=plot(t,df0dv_z,layout=[1,3,3],/current, font_size=fs)


stop
end
