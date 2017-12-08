function int_tabulated_2d, x, y, f

    nx = n_elements(f[*,0])
    ny = n_elements(f[0,*])

    resA = fltArr(nx)
    for i=0,nx-1 do begin

           resA[i] = int_tabulated (y,f[i,*]) 

    endfor

    return, int_tabulated (x,resA)
end

function int_tabulated_3d, x, y, z, f

    nx = n_elements(f[*,0,0])
    ny = n_elements(f[0,*,0])
    nz = n_elements(f[0,0,*])

    resB = fltArr(nx,ny)

    for i=0,nx-1 do begin
        for j=0,ny-1 do begin
    
           resB[i,j] = int_tabulated (z,f[i,j,*]) 

        endfor
    endfor

    return, int_tabulated_2d (x,y,resB)
end


function trapezoidal, x, y, f

    nx = n_elements(f[*,0])
    ny = n_elements(f[0,*])

    result = 0

    resA = fltArr(nx)
    for i=0,nx-1 do begin

           resA[i] = int_tabulated (y,f[i,*]) 

    endfor

    return, int_tabulated (x,resA)
end

pro kj_read_f1

    @dlg_constants

    fileName = 'output/f1.txt'
    nLines = file_lines(fileName)
    openr, lun, fileName, /get_lun

	cfg = kj_read_cfg('./')

	nx = cfg['nP_Vx'] 
	ny = cfg['nP_Vy']
	nz = cfg['nP_Vz']


    offset = 1

    vx = fltArr(nX)
    vy = fltArr(nY)
    vz = fltArr(nZ)

    vx3 = fltArr(nX,nY,nZ)
    vy3 = fltArr(nY,nY,nZ)
    vz3 = fltArr(nZ,nY,nZ)


    f1 = complexArr(nX,nY,nZ)

    	skip_lun, lun, offset, /lines

    	for i=0,nx-1 do begin
    		for j=0,ny-1 do begin
    			for k=0,nz-1 do begin
        			readf, lun, _vx, _vy, _vz, _f1re, _f1im 

        			vx[i] = _vx
        			vy[j] = _vy
        			vz[k] = _vz

        			vx3[i,j,k] = _vx
        			vy3[i,j,k] = _vy
        			vz3[i,j,k] = _vz

        			f1[i,j,k] = complex(_f1re,_f1im)
			endfor
		endfor
    	endfor

	free_lun, lun
	;dx = (max(vx)-min(vx))/(nx-1)
	;vxGrid = min(vx)+fIndGen(nx)*dx
	;dy = (max(vy)-min(vy))/(ny-1)
	;vyGrid = min(vy)+fIndGen(ny)*dy
	;dz = (max(vz)-min(vz))/(nz-1)
	;vzGrid = min(vz)+fIndGen(nz)*dz

	;rxy = GridData(vx,vy,f1,/grid,xout=vxGrid,yout=vyGrid,title='f1_xy')
	;rxz = GridData(vx,vz,f1,/grid,xout=vxGrid,yout=vzGrid,title='f1_xz')
	;ryz = GridData(vy,vz,f1,/grid,xout=vyGrid,yout=vyGrid,title='f1_yz')

	rxy = reform(real_part(f1[*,*,nz/2]))
	rxz = reform(real_part(f1[*,ny/2,*]))
	ryz = reform(real_part(f1[nx/2,*,*]))

	n_lev = 15L
	maxL = max(abs([max(abs(rxy)),max(abs(rxz)),max(abs(ryz))]))/1.0
	power = 3L
	dLev = maxL/n_lev^power
	levels = findgen(n_lev)^power*dLev
	colors = 255-(bytscl(findgen(n_lev),top=254)+1)

	c=contour(rxy,vx,vy,/fill,$
		layout=[2,2,1],rgb_table=1,c_value=levels,c_color=colors,title='f1_xy')
	c=contour(-rxy,vx,vy,/fill,$
		layout=[2,2,1],rgb_table=7,c_value=levels,c_color=colors,/current)
	
	c=contour(rxz,vx,vz,/fill,$
		layout=[2,2,2],/current,rgb_table=1,c_value=levels,c_color=colors,title='f1_xz')
	c=contour(-rxz,vx,vz,/fill,$
		layout=[2,2,2],/current,rgb_table=7,c_value=levels,c_color=colors)
	
	c=contour(ryz,vy,vz,/fill,$
		layout=[2,2,3],rgb_table=1,c_value=levels,c_color=colors,/current,title='f1_yz')
	c=contour(-ryz,vy,vz,/fill,$
		layout=[2,2,3],rgb_table=7,c_value=levels,c_color=colors,/current)
stop
    zSlice = 2
	rxy = reform(imaginary(f1[*,*,zSlice]))
    ySlice = 2 
	rxz = reform(imaginary(f1[*,ySlice,*]))
    xSlice = 2 
	ryz = reform(imaginary(f1[xSlice,*,*]))

	c=contour(rxy,vx,vy,/fill,$
		layout=[2,2,1],rgb_table=1,c_value=levels,c_color=colors)
	c=contour(-rxy,vx,vy,/fill,$
		layout=[2,2,1],rgb_table=7,c_value=levels,c_color=colors,/current)
	
	c=contour(rxz,vx,vz,/fill,$
		layout=[2,2,2],/current,rgb_table=1,c_value=levels,c_color=colors)
	c=contour(-rxz,vx,vz,/fill,$
		layout=[2,2,2],/current,rgb_table=7,c_value=levels,c_color=colors)
	
	c=contour(rxy,vx,vy,/fill,$
		layout=[2,2,3],rgb_table=1,c_value=levels,c_color=colors,/current)
	c=contour(-rxy,vx,vy,/fill,$
		layout=[2,2,3],rgb_table=7,c_value=levels,c_color=colors,/current)


	nPer = 50
	nPar = 100
	nPhs = 90

	vPerMin = 0
	vPerMax = max(vx)
	dvPer = (vPerMax-vPerMin)/(nPer-1)
	vPer = fIndGen(nPer)*dvPer + vPerMin 

	vParMin = -max(vx) 
	vParMax = max(vx)
	dvPar = (vParMax-vParMin)/(nPar-1)
	vPar = fIndGen(nPar)*dvPar + vParMin 

	vPhsMin = 0
	vPhsMax = 2*!pi
	dvPhs = (vPhsMax-vPhsMin)/(nPhs-1)
	vPhs = fIndGen(nPhs)*dvPhs + vPhsMin 

	f1_CYL = fltArr(nPer,nPar,nPhs)	

	for i=0,nPer-1 do begin
		for j=0,nPar-1 do begin
			for k=0,nPhs-1 do begin
				thisX = vPer[i]*cos(vPhs[k])
				thisY = vPer[i]*sin(vPhs[k])
				thisZ = vPar[j]

				_x = (thisX-min(vx))/(vx[-1]-vx[0])*(nx-1)
				_y = (thisY-min(vy))/(vy[-1]-vy[0])*(ny-1)
				_z = (thisZ-min(vz))/(vz[-1]-vz[0])*(nz-1)

				f1_CYL[i,j,k] = interpolate(f1,_x,_y,_z)
			endfor
		endfor
	endfor

	f1pp = total(f1_CYL,3)/nPhs

	c=contour(transpose(f1pp),vPar,vPer,/fill,$
		rgb_table=1,c_value=levels,c_color=colors,title='f1',$
		xTitle='vPar',ytitle='vPer')
	c=contour(-transpose(f1pp),vPar,vPer,/fill,$
		rgb_table=7,c_value=levels,c_color=colors,/over)

	this = rxz*rebin(vx,nx,nz)

	maxL = max(abs([max(abs(this))]))/1.0
	power = 3L
	dLev = maxL/n_lev^power
	levels = findgen(n_lev)^power*dLev
	colors = 255-(bytscl(findgen(n_lev),top=254)+1)

	c=contour(this,vx,vz,/fill,$
		layout=[2,2,1],rgb_table=1,c_value=levels,c_color=colors,$
                    title='vx * f1_xz', xtitle='vx',ytitle='vz')
	c=contour(-this,vx,vz,/fill,$
		layout=[2,2,1],rgb_table=7,c_value=levels,c_color=colors,/over)

        print, 'vx * f1_xz [total]: ', total(this)*(vx[1]-vx[0])*(vz[1]-vz[0])
        print, 'vx * f1_xz [int tab]: ', int_tabulated_2d(vx,vz,this)  



	this = rxy*transpose(rebin(vy,ny,nx))
	c=contour(this,vx,vy,/fill,$
		layout=[2,2,2],/current,rgb_table=1,c_value=levels,c_color=colors,$
                    title='vy * f1_xy',xtitle='vx',ytitle='vy')
	c=contour(-this,vx,vy,/fill,$
		layout=[2,2,2],rgb_table=7,c_value=levels,c_color=colors,/over)

        print, 'vy * f1_xy [total]: '  ,total(this)*(vx[1]-vx[0])*(vy[1]-vy[0])
        print, 'vy * f1_xy [int tab]: ',int_tabulated_2d(vx,vy,this)  


	this = ryz*rebin(vy,ny,nz)
	c=contour(this,vy,vz,/fill,$
		layout=[2,2,3],/current,rgb_table=1,c_value=levels,c_color=colors,$
                title='vy * f1_yz',xtitle='vy',ytitle='vz')
	c=contour(-this,vy,vz,/fill,$
		layout=[2,2,3],rgb_table=7,c_value=levels,c_color=colors,/over)

        print, 'vy * f1_yz [total]: '  ,total(this)*(vy[1]-vy[0])*(vz[1]-vz[0])
        print, 'vy * f1_yz [int tab]: ',int_tabulated_2d(vy,vz,this)  


    this = ryz*transpose(rebin(vz,nz,ny))
	c=contour(this,vy,vz,/fill,$
		layout=[2,2,4],/current,rgb_table=1,c_value=levels,c_color=colors,$
                title='vz * f1_yz',xtitle='vy',ytitle='vz')
	c=contour(-this,vy,vz,/fill,$
		layout=[2,2,4],rgb_table=7,c_value=levels,c_color=colors,/over)

        print, 'vz * f1_yz [total]: '  , total(this)*(vy[1]-vy[0])*(vz[1]-vz[0])
        print, 'vz * f1_yz [int tab]: ', int_tabulated_2d(vy,vz,this)  



        dV = (vx[1]-vx[0])*(vy[1]-vy[0])*(vz[1]-vz[0])

        print, 'jy [total]: ', total((vx3*f1))*dV*_e
        print, 'jx [int tab]: ', int_tabulated_3d( vx, vy, vz, real_part(vx3*f1))*_e

        print, 'jy [total]: ', total((vy3*f1))*dV*_e
        print, 'jy [int tab]: ', int_tabulated_3d( vx, vy, vz, real_part(vy3*f1))*_e

        print, 'jz [total]: ', total((vz3*f1))*dV*_e
        print, 'jz [int tab]: ', int_tabulated_3d( vx, vy, vz, real_part(vz3*f1))*_e



stop
end
