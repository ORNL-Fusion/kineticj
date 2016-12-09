pro kj_plot_orbit, overPlotLorentz=_overPlotLorentz, runDir=_runDir

    if keyword_set(_overPlotLorentz) then overPlotLorentz=_overPlotLorentz else overPlotLorentz=0
    if keyword_set(_runDir) then runDir=_runDir else runDir='./'

    lorentzDir = '~/scratch/kineticj/colestock-kashuba-lorentz/'

	cfg = kj_read_cfg (runDir)	

    fileName = runDir+'/output/orbit_e1_dot_grad_df0_dv.txt'
    nHeader = 1
    N = file_lines(fileName)-nHeader
    e1_dot_data = replicate({t:0.0,re:0.0,im:0.0,rePer:0.0,imPer:0.0,rePar:0.0,imPar:0.0},N)

    openr, lun, fileName, /get_Lun
        skip_lun, lun, nHeader, /lines
        readf, lun, e1_dot_data
    close, lun

    if overPlotLorentz then begin

        fileName = lorentzDir + 'output/orbit_e1_dot_grad_df0_dv.txt'
        nHeader = 1
        N = file_lines(fileName)-nHeader
        e1_dot_data_lorentz = replicate({t:0.0,re:0.0,im:0.0,rePer:0.0,imPer:0.0,rePar:0.0,imPar:0.0},N)
        openr, lun, fileName, /get_Lun
            skip_lun, lun, nHeader, /lines
            readf, lun, e1_dot_data_lorentz
        close, lun

        t = e1_dot_data_lorentz.t
        wc = 5.19163e11
        f = wc/(2*!pi)
        reA = 8.08e-14
        reP = 31.0*!dtor ; 0 gives int = 0 
        imA = 5.0e-14
        imP = 128.5*!dtor ; 90 gives int = 0 
        re = reA*cos(2*!pi*f*t+reP)
        im = -imA*sin(2*!pi*f*t+imP)
        h = hanning(n_elements(re))
        h[0:N/2]=1
        reAnalytic = re*h
        imAnalytic = im*h

        ;e1_dot_data_lorentz.re = reAnalytic
        ;e1_dot_data_lorentz.im = imAnalytic

        re = e1_dot_data_lorentz.re
        im = e1_dot_data_lorentz.im

        rePer = e1_dot_data_lorentz.rePer
        imPer = e1_dot_data_lorentz.imPer

        rePar = e1_dot_data_lorentz.rePar
        imPar = e1_dot_data_lorentz.imPar

        p=plot(t,re,thick=2)
        p=plot(t,reAnalytic,/over)
        p=plot(t,im,color='r',thick=2,/over)
        p=plot(t,imAnalytic,color='r',/over)
        p=plot(t,rePer,/over,color='b')
        p=plot(t,imPer,color='b',thick=2,/over)
 
    endif

    margin = [0.2,0.2,0.2,0.2]

    sW = 100
    order = 100
    dt = e1_dot_data[1].t-e1_dot_data[0].t
    if overPlotLorentz then begin
        dtL = t[1]-t[0]
    endif

    print, 'Int [re]: ', total(e1_dot_data.re)
    print, 'Int [im]: ', total(e1_dot_data.im)

    p=plot(e1_dot_data.t,e1_dot_data.re,layout=[1,2,1],margin=margin,title='e1 . \/v f0')
    if overPlotLorentz then begin
        print, 'Int Lorentz [re]: ', total(re)
        print, 'Int Lorentz [im]: ', total(im)
        p=plot(t,re,/over)
        p=plot(t,rePer,/over,color='b')
        p=plot(t,rePar,/over,color='g')
    endif
    p=plot(e1_dot_data.t,e1_dot_data.im,/over,color='r')
    if overPlotLorentz then begin
        p=plot(t,im,/over,color='r')
        p=plot(t,imPer,/over,color='b',lineStyle='2')
        p=plot(t,imPar,/over,color='g',lineStyle='2')
    endif

    p=plot(e1_dot_data.t,total((e1_dot_data.re),/cum)*dt,thick=2,layout=[1,2,2],/current,margin=margin,title="int(e1 . \/vf0)")
    if overPlotLorentz then begin
        p=plot(t,total((re),/cum)*dtL,/over)
    endif
    p=plot(e1_dot_data.t,total((e1_dot_data.im),/cum)*dt,/over,thick=2,color='r')
    if overPlotLorentz then begin
        p=plot(t,total((im),/cum)*dtL,/over,color='r')
    endif

    fileName = runDir+'/output/orbit.txt'
    nHeader = 2
    N = file_lines(fileName)-nHeader
    orbit = replicate({t:0.0, x:0.0,y:0.0,z:0.0,$
            exr:0.0,exc:0.0,eyr:0.0,eyc:0.0,ezr:0.0,ezc:0.0,$
            bxr:0.0,bxc:0.0,byr:0.0,byc:0.0,bzr:0.0,bzc:0.0, $
            vCrossBxr:0.0,vCrossBxi:0.0,vCrossByr:0.0,vCrossByi:0.0,vCrossBzr:0.0,vCrossBzi:0.0,$
            status:0},N)
    openr, lun, fileName, /get_Lun
        skip_lun, lun, nHeader, /lines
        readf, lun, orbit 
    close, lun

    if overPlotLorentz then begin
        fileName = lorentzDir+'output/orbit.txt'
        nHeader = 2
        N = file_lines(fileName)-nHeader
        orbit_lorentz = replicate({t:0.0, x:0.0,y:0.0,z:0.0,$
                exr:0.0,exc:0.0,eyr:0.0,eyc:0.0,ezr:0.0,ezc:0.0,$
                bxr:0.0,bxc:0.0,byr:0.0,byc:0.0,bzr:0.0,bzc:0.0, $
                vCrossBxr:0.0,vCrossBxi:0.0,vCrossByr:0.0,vCrossByi:0.0,vCrossBzr:0.0,vCrossBzi:0.0,$
                status:0},N)
        openr, lun, fileName, /get_Lun
            skip_lun, lun, nHeader, /lines
            readf, lun, orbit_lorentz 
        close, lun
    endif

    p=plot(orbit[*].t, orbit[*].x,layout=[1,3,1],margin=margin,title='Position (x,y,z)' )
    if overPlotLorentz then p=plot(orbit_lorentz[*].t, orbit_lorentz[*].x,/over )

    p=plot(orbit[*].t, orbit[*].y,layout=[1,3,2],/current,margin=margin)
    if overPlotLorentz then p=plot(orbit_lorentz[*].t, orbit_lorentz[*].y,/over )

    p=plot(orbit[*].t, orbit[*].z,layout=[1,3,3],/current,margin=margin)
    if overPlotLorentz then p=plot(orbit_lorentz[*].t, orbit_lorentz[*].z,/over )



    p=plot(orbit[*].t, orbit[*].exr,layout=[1,3,1],margin=margin,title='E field' )
    p=plot(orbit[*].t, orbit[*].exc,/over,color='r' )
    if overPlotLorentz then p=plot(orbit_lorentz[*].t, orbit_lorentz[*].x,/over )

    p=plot(orbit[*].t, orbit[*].eyr,layout=[1,3,2],/current,margin=margin)
    p=plot(orbit[*].t, orbit[*].eyc,/over,color='r' )
    if overPlotLorentz then p=plot(orbit_lorentz[*].t, orbit_lorentz[*].y,/over )

    p=plot(orbit[*].t, orbit[*].ezr,layout=[1,3,3],/current,margin=margin)
    p=plot(orbit[*].t, orbit[*].ezc,/over,color='r' )
    if overPlotLorentz then p=plot(orbit_lorentz[*].t, orbit_lorentz[*].z,/over )


   if overPlotLorentz then begin

        meanx = mean(orbit_lorentz.x)
        range = 10e-3
        p3 = plot3d(orbit[*].x-meanx,orbit[*].y,orbit[*].z,xrange=[-1,+1]*range,yrange=[min(orbit.y),max(orbit.y)],zrange=[-1,1]*range,$
            xtitle='X',ytitle='Y',ztitle='Z')
        p3=plot3d(orbit_lorentz[*].x-meanx,orbit_lorentz[*].y,orbit_lorentz[*].z,/over )
    endif


    if overPlotLorentz then begin

        print, 'freq: ', f
        period = 1/f
        print, 'period: ', period

        angle1 = reP 
        angle2 = imP 
        angle1frac = angle1/(2*!pi)
        angle2frac = angle2/(2*!pi)

        nPts1 = (angle1frac*period)/dtL
        nPts2 = (angle2frac*period)/dtL

        print, 'angle1: ', angle1*!radeg
        print, 'angle2: ', angle2*!radeg

        print, 'nPts1: ', nPts1
        print, 'angle1frac: ', angle1frac
        print, 'angle2frac: ', angle2frac
        print, 'nPts2: ', nPts2

        if angle1*!radeg ge 0 then angle1end = 180 
        if angle1*!radeg ge 180 then angle1end = 180+180 

        if angle2*!radeg ge 0 then angle2end = 90
        if angle2*!radeg ge 90 then angle2end = 90+180
        if angle2*!radeg ge 270 then angle2end = 90+180+180

        print, 'Angle1End: ', angle1end
        print, 'Angle2End: ', angle2end

        print, 'Offest[re] analytic: ', ( sin(angle1end*!dtor) - sin(angle1) ) * reA * dtL * 2 * !pi
        print, 'Offest[im] analytic: ', -( cos(angle2) - cos(angle2end*!dtor) ) * imA * dtL * 2 * !pi

        print, 'Int all[re]: ', int_tabulated(t,re,/double)
        print, 'Int all[im]: ', int_tabulated(t,im,/double)

        if nPts1 gt 1 and nPts2 gt 1 then begin
            print, 'Offset[re]: ', int_tabulated(t[0:nPts1],re[0:nPts1])
            print, 'Offset[im]: ', int_tabulated(t[0:nPts2],im[0:nPts2])

            print, 'Int of the rest[re]: ', int_tabulated(t[nPts1+1:-1],re[nPts1+1:-1])
            print, 'Int of the rest[im]: ', int_tabulated(t[nPts2+1:-1],im[nPts2+1:-1])

            p=plot(t[0:nPts1],re[0:nPts1],symbol='D')
            p=plot(t[0:nPts2],im[0:nPts2],/over,color='r',symbol='D')
        endif

    endif 
stop
end
