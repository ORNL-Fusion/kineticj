pro kj_read_f1

    fileName = 'output/f1.txt'
    nLines = file_lines(fileName)
    openr, lun, fileName, /get_lun

	cfg = kj_read_cfg('./')

	nx = cfg['nPx'] 
	ny = cfg['nPy']
	nz = cfg['nPz']


    offset = 1

		

    vx = fltArr(nX)
    vy = fltArr(nY)
    vz = fltArr(nZ)

    f1 = complexArr(nX,nY,nZ)

    	skip_lun, lun, offset, /lines

    	for i=0,nx-1 do begin
    		for j=0,ny-1 do begin
    			for k=0,nz-1 do begin
        			readf, lun, _vx, _vy, _vz, _f1re, _f1im 

        			vx[i] = _vx
        			vy[j] = _vy
        			vz[k] = _vz

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
stop
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

	rxy = reform(imaginary(f1[*,*,nz/2]))
	rxz = reform(imaginary(f1[*,ny/2,*]))
	ryz = reform(imaginary(f1[nx/2,*,*]))

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
		layout=[2,2,1],rgb_table=1,c_value=levels,c_color=colors,title='vx moment')
	c=contour(-this,vx,vz,/fill,$
		layout=[2,2,1],rgb_table=7,c_value=levels,c_color=colors,/current)
	
	;this = rxz*transpose(rebin(vz,nz,nx))
	;c=contour(this,vz,vx,/fill,$
	;	layout=[2,2,2],/current,rgb_table=1,c_value=levels,c_color=colors,title='vz moment',$
    ;    xtitle='vz',ytitle='vx')
	;c=contour(-this,vz,vx,/fill,$
	;	layout=[2,2,2],rgb_table=7,c_value=levels,c_color=colors,/current)
	
stop
end
