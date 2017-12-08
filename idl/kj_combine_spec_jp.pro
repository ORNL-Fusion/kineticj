pro kj_combine_spec_jp, FilesToSum, SumFileName = SumFileName, $
        cartesian_offset = _cartesian_offset, $
        doPlots = _doPlots
  
    if keyword_set(_cartesian_offset) then cartesian_offset = _cartesian_offset else cartesian_offset = 0
    if keyword_set(_doPlots) then doPlots = _doPlots else doPlots = 0

    if keyword_set(SumFileName) then begin
            OutFileName = SumFileName
    endif else begin
            OutFileName = 'kj_jp.nc'
    endelse

    nS = n_elements(FilesToSum)

    kjTmp = kj_read_kj_jp(FilesToSum[0])

    nR = n_elements(kjTmp.jPr)
    nR_ = n_elements(kjTmp.jPr_)

    jPSum_r = complexArr(nR) 
    jPSum_t = complexArr(nR) 
    jPSum_z = complexArr(nR) 

    jPSpec_r = complexArr(nR,1,nS) 
    jPSpec_t = complexArr(nR,1,nS) 
    jPSpec_z = complexArr(nR,1,nS) 

    for s=0,nS-1 do begin

        kj = kj_read_kj_jp(FilesToSum[s])

        jPSum_r = jPSum_r + kj.jPr
        jPSum_t = jPSum_t + kj.jPt
        jPSum_z = jPSum_z + kj.jPz

        jPSpec_r[*,0,s] = kj.jPr
        jPSpec_t[*,0,s] = kj.jPt
        jPSpec_z[*,0,s] = kj.jPz

        print, n_elements(kj.jpr)

        if doPlots then begin

            p=plot(kj.r, kj.jPr, yRange=yRange, layout=[nS+1,3,s+1], /current, title=FilesToSum[s])
            p=plot(kj.r, imaginary(kj.jPr), color='r',/over,yRange=yRange)

            p=plot(kj.r, kj.jPt, yRange=yRange, layout=[nS+1,3,s+nS+2], /current)
            p=plot(kj.r, imaginary(kj.jPt), color='r',/over,yRange=yRange)

            p=plot(kj.r, kj.jPz, yRange=yRange, layout=[nS+1,3,s+2*nS+3], /current)
            p=plot(kj.r, imaginary(kj.jPz), color='r',/over,yRange=yRange )

        endif

    endfor

    if doPlots then begin

        p=plot(kj.r, jPsum_r, yRange=yRange, layout=[nS+1,3,s+1], /current, title='Sum')
        p=plot(kj.r, imaginary(jPsum_r), color='r',/over,yRange=yRange)

        p=plot(kj.r, jPsum_t, yRange=yRange, layout=[nS+1,3,s+nS+2], /current)
        p=plot(kj.r, imaginary(jPsum_t), color='r',/over,yRange=yRange)

        p=plot(kj.r, jPsum_z, yRange=yRange, layout=[nS+1,3,s+2*nS+3], /current)
        p=plot(kj.r, imaginary(jPsum_z), color='r',/over,yRange=yRange )

    endif

    jPSpec_r[0] = 0
    jPSpec_r[-1] = 0
    jPSpec_t[0] = 0
    jPSpec_t[-1] = 0
    jPSpec_z[0] = 0
    jPSpec_z[-1] = 0

    ;; Try smoothing the kj current

    ;sWidth = 60
    ;for s=0,nS-1 do begin
    ;    jPSpec_r[*,0,s] = smooth(jPSpec_r[*,0,s],sWidth,/edge_truncate)   
    ;    jPSpec_t[*,0,s] = smooth(jPSpec_t[*,0,s],sWidth,/edge_truncate)   
    ;    jPSpec_z[*,0,s] = smooth(jPSpec_z[*,0,s],sWidth,/edge_truncate)   
    ;endfor

    ; Try mixing with the real AORSA jP to test sensitivity to the iteration being not quite correct.

    ;ar2 = ar2_read_solution('/Users/dg6/scratch/aorsa2d/colestock-kashuba-reference',1)
    rs = rsfwc_read_solution('rsfwc')

    ; Use RS
    jP_r_prev = rs.jP_r_spec
    jP_t_prev = rs.jP_t_spec
    jP_z_prev = rs.jP_z_spec

    ;; Use AR2
    ;jP_r_prev = ar2.jP_r
    ;jP_t_prev = ar2.jP_t
    ;jP_z_prev = ar2.jP_z

    kjFrac = 1.0
    prevFrac = 1 - kjFrac

    ;; Test if it's the boundary by weighting those preferentially 

    ;h0 = hanning(n_elements(jPSpec_r[*,0,0]))
    ;h = jPSpec_r*0
    ;for s=0,nS-1 do h[*,0,s] = complex(h0,h0)
    ;prevFrac = (1-h*kjFrac)

    jPSpec_r = jPSpec_r*(1-prevFrac) + jP_r_prev*prevFrac
    jPSpec_t = jPSpec_t*(1-prevFrac) + jP_t_prev*prevFrac
    jPSpec_z = jPSpec_z*(1-prevFrac) + jP_z_prev*prevFrac

    ; Try using only the delta on the RHS

    jPSpec_r = jPSpec_r - jP_r_prev 
    jPSpec_t = jPSpec_t - jP_t_prev 
    jPSpec_z = jPSpec_z - jP_z_prev 

    w = kj.freq * 2 * !pi

    @dlg_constants

	; Write kj_jP in file for next iterate

	nc_id = nCdf_create (OutFileName, /clobber )

	nCdf_control, nc_id, /fill
	
	nr_id = nCdf_dimDef ( nc_id, 'nR', n_elements(kj.r) )
	nX_id = nCdf_dimDef ( nc_id, 'nX', n_elements(kj.r) )

	nz_id = nCdf_dimDef ( nc_id, 'nZ', 1 )
	nrH_id = nCdf_dimDef ( nc_id, 'nR_', n_elements(kj.r_) )
	ny_id = nCdf_dimDef ( nc_id, 'nY', 1 )
	scalar_id = nCdf_dimDef ( nc_id, 'scalar', 1 )
    nS_id = nCdf_dimDef( nc_id, 'nSpec', nS)

	freq_id = nCdf_varDef ( nc_id, 'freq', scalar_id, /float )
	r_id = nCdf_varDef ( nc_id, 'r', nX_id, /float )
	rH_id = nCdf_varDef ( nc_id, 'r_', nrH_id, /float )

	z_id = nCdf_varDef ( nc_id, 'z', nY_id, /float )

	jP_r_re_id = nCdf_varDef ( nc_id, 'jP_r_re', [nX_id,nY_id,nS_id], /float )
	jP_r_im_id = nCdf_varDef ( nc_id, 'jP_r_im', [nX_id,nY_id,nS_id], /float )
	jP_p_re_id = nCdf_varDef ( nc_id, 'jP_t_re', [nX_id,nY_id,nS_id], /float )
	jP_p_im_id = nCdf_varDef ( nc_id, 'jP_t_im', [nX_id,nY_id,nS_id], /float )
	jP_z_re_id = nCdf_varDef ( nc_id, 'jP_z_re', [nX_id,nY_id,nS_id], /float )
	jP_z_im_id = nCdf_varDef ( nc_id, 'jP_z_im', [nX_id,nY_id,nS_id], /float )

	nCdf_control, nc_id, /enDef

	nCdf_varPut, nc_id, freq_id, kj.freq

	nCdf_varPut, nc_id, z_id, 0.0

	nCdf_varPut, nc_id, r_id, kj.r + cartesian_offset
	nCdf_varPut, nc_id, rH_id, kj.r_ + cartesian_offset

	nCdf_varPut, nc_id, jP_r_re_id, real_part(jPSpec_r)
	nCdf_varPut, nc_id, jP_r_im_id, imaginary(jPSpec_r) 
	nCdf_varPut, nc_id, jP_p_re_id, real_part(jPSpec_t)
	nCdf_varPut, nc_id, jP_p_im_id, imaginary(jPSpec_t) 
	nCdf_varPut, nc_id, jP_z_re_id, real_part(jPSpec_z)
	nCdf_varPut, nc_id, jP_z_im_id, imaginary(jPSpec_z) 

	nCdf_close, nc_id

    if doPlots then begin

        p=plot(total(jPSpec_r,3), layout=[2,3,1])
        p=plot(total(jP_r_prev,3),/over,color='red')
        p=plot(imaginary(total(jPSpec_r,3)), layout=[2,3,2],/current)
        p=plot(imaginary(total(jP_r_prev,3)),/over,color='red')

        p=plot(total(jPSpec_t,3), layout=[2,3,3],/current)
        p=plot(total(jP_t_prev,3),/over,color='red')
        p=plot(imaginary(total(jPSpec_t,3)), layout=[2,3,4],/current)
        p=plot(imaginary(total(jP_t_prev,3)),/over,color='red')

        p=plot(total(jPSpec_z,3), layout=[2,3,5],/current)
        p=plot(total(jP_z_prev,3),/over,color='red')
        p=plot(imaginary(total(jPSpec_z,3)), layout=[2,3,6],/current)
        p=plot(imaginary(total(jP_z_prev,3)),/over,color='red')

    endif

    print, norm(jPSpec_r[*]), norm(jP_r_prev[*])
    print, norm(jPSpec_t[*]), norm(jP_t_prev[*])
    print, norm(jPSpec_z[*]), norm(jP_z_prev[*])

end
