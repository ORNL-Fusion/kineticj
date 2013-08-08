pro kj_plot_current, noInterp = noInterp, sig33 = sig33, noTimeDep = noTimeDep, noIterate=noIterate

	@constants

	cd, current=currentDir

	runFile = 'kj.cfg'
	cfg = kj_read_cfg (runFile)	

	openr, lun, runFile, /get_lun
	runFileArray = ''
	line = ''
	while not eof(lun) do begin
		readf, lun, line
		runFileArray = [runFileArray,line]
	endwhile
	free_lun, lun
	for f=0,n_elements(runFileArray)-1 do begin
		if(strMatch(runFileArray[f],'*eField_fName*'))then $
				eField_fName=(strSplit(runFileArray[f],'"',/extract))[1]
	endfor

	print, 'eField_fName: ', eField_fName

	cdfId = ncdf_open(eField_fName)

		ncdf_varget, cdfId, 'freq', freq 
		ncdf_varget, cdfId, 'r', r 
		sheath = 0
		if(strMatch(eField_fName,'*sheath*') $
				or strMatch(eField_fName,'*mets*') $
                or strMatch(eField_fName,'*single*') $
                or strMatch(eField_fName,'*upshift*') )then begin
			if(strMatch(eField_fName,'*sheath*'))then sheath = 1
			r_ = r[0:-2]+(r[1]-r[0])/2.0
		endif else begin
			ncdf_varget, cdfId, 'r_', r_
		endelse

		ncdf_varget, cdfId, 'e_r_re', er_re
		ncdf_varget, cdfId, 'e_r_im', er_im
		ncdf_varget, cdfId, 'e_p_re', ep_re
		ncdf_varget, cdfId, 'e_p_im', ep_im
		ncdf_varget, cdfId, 'e_z_re', ez_re
		ncdf_varget, cdfId, 'e_z_im', ez_im

		ncdf_varget, cdfId, 'jP_r_re', jPr_re
		ncdf_varget, cdfId, 'jP_r_im', jPr_im
		ncdf_varget, cdfId, 'jP_p_re', jPt_re
		ncdf_varget, cdfId, 'jP_p_im', jPt_im
		ncdf_varget, cdfId, 'jP_z_re', jPz_re
		ncdf_varget, cdfId, 'jP_z_im', jPz_im

	ncdf_close, cdfId

	wrf = freq * 2 * !pi

	r_prevIterate = r

	jPr_prevIterate = complex ( jPr_re, jPr_im ) 
	jPt_prevIterate = complex ( jPt_re, jPt_im ) 
	jPz_prevIterate = complex ( jPz_re, jPz_im ) 

	spline_sigma = 0.01

	jPr_prevIterate_ = complex(spline(r,jPr_re,r_,spline_sigma),spline(r,jPr_im,r_,spline_sigma))
	jPt_prevIterate_ = complex(spline(r,jPt_re,r_,spline_sigma),spline(r,jPt_im,r_,spline_sigma))
	jPz_prevIterate_ = complex(spline(r,jPz_re,r_,spline_sigma),spline(r,jPz_im,r_,spline_sigma))

	fileList = file_search ( 'output/'+cfg.runIdent+'/jP*' )

	cdfId = ncdf_open(fileList[0])
	ncdf_varget, cdfId, 't', t
	nCdf_close, cdfId
	
	nT = n_elements(t)
	nF = n_elements(fileList)

	j1x = fltArr ( nT, nF )
	j1xc = complexArr ( nF )
	j1 = complexArr ( nF )
	xF = fltArr ( nF )

	;hanWindow = hanning (nT, alpha=0.5 )

	for f=0,nF-1 do begin

		cdfId = ncdf_open(fileList[f])

			ncdf_varget, cdfId, 'x', x
			ncdf_varget, cdfId, 't', t

			if not keyword_set(noTimeDep) then begin
				ncdf_varget, cdfId, 'j1x', j1x_0 
				ncdf_varget, cdfId, 'j1y', j1y_0 
				ncdf_varget, cdfId, 'j1z', j1z_0 
			endif

			ncdf_varget, cdfId, 'j1xc_re', j1xc_re 
			ncdf_varget, cdfId, 'j1xc_im', j1xc_im

		nCdf_close,	cdfId 

		xF[f] = x

		j1xc[f] = complex(j1xc_re,j1xc_im)

		if not keyword_set(noTimeDep) then begin

			j1x[*,f] = j1x_0;-mean(j1x_0) ; not sure if there is a nicer way to do this

			dt = t[1]-t[0]

			;; Test code
			;ampRe = 1000.0
			;ampIm = -300.0
			;amp = sqrt(ampRe^2+ampIm^2)
			;phs = atan(ampIm,ampRe)
			;j1x = (amp * exp (-II*(wrf*t+phs)))

			jrFFT = fft ( j1x[0:-1,f], /center )
			freqAxis = (fIndGen(nT)-(nT/2)) / (nT*dt)
			freqNorm = freqAxis / freq

			;p=plot(freqNorm,abs(jrFFT)^2);,xRange=[0,5])
			;p=plot(freqNorm,imaginary(jrFFT),color='blue');,xRange=[0,5])

			; Positive (right) frequency 
			iiAntFreq = where(abs(freqNorm-1) eq min(abs(freqNorm-1)),iiAntFreqCnt)
			rpL = real_part ( jrFFT[iiAntFreq[0]] )
			ipL = imaginary ( jrFFT[iiAntFreq[0]] )
			; Negative (left) frequency 
			iiAntFreq = where(abs(freqNorm+1) eq min(abs(freqNorm+1)),iiAntFreqCnt)
			rpR = real_part ( jrFFT[iiAntFreq[0]] )
			ipR = -imaginary ( jrFFT[iiAntFreq[0]] )

			print, 'Re: ', rpL, rpR, rpL+rpR
			print, 'Im: ', ipL, ipR, ipL+ipR

			j1[f] = complex ( rpL+rpR, ipL+ipR )

		endif else begin
            j1[f] = j1xc[f]
        endelse
	
	endfor

	; Create Debye length axis

	n = 1.1d14
	n_20 = n/10d0^20
	T_keV = 0.0001
   	lambda_D = 2.35d-5*sqrt(T_keV/n_20)
	print, "Debye Length: ", lambda_D


if not keyword_set(noIterate) then begin

	if not keyword_set(noTimeDep) then begin

		; Fudge Factoring
		j1 = conj(j1) ; Not sure why we neet a conj here. Most likely due to the way I'm doing the fft

		; Create a jP for rsfcw_1d

		jROut  = complex(spline(xf,real_part(j1),r,spline_sigma),spline(xf,imaginary(j1),r,spline_sigma))
		jROut_ = complex(spline(xf,real_part(j1),r_,spline_sigma),spline(xf,imaginary(j1),r_,spline_sigma))

		jTOut = jROut*0
		jTOut_ = jROut_*0

		jZOut = jROut*0
		jZOut_ = jROut_*0

		; Average the iterations to test stability and convergence.

	endif else begin


		jROut  = complex(spline(xf,real_part(j1),r,spline_sigma),spline(xf,imaginary(j1),r,spline_sigma))
		jROut_ = complex(spline(xf,real_part(j1),r_,spline_sigma),spline(xf,imaginary(j1),r_,spline_sigma))

		jTOut = jROut*0
		jTOut_ = jROut_*0

		jZOut = jROut*0
		jZOut_ = jROut_*0

	endelse

	relaxFactor = 0.0

	jROut = jROut*(1-relaxFactor) + jPr_prevIterate*relaxFactor
	jTOut = jTOut*(1-relaxFactor) + jPt_prevIterate*relaxFactor
	jZOut = jZOut*(1-relaxFactor) + jPz_prevIterate*relaxFactor
	
	jROut_ = jROut_*(1-relaxFactor) + jPr_prevIterate_*relaxFactor
	jTOut_ = jTOut_*(1-relaxFactor) + jPt_prevIterate_*relaxFactor
	jZOut_ = jZOut_*(1-relaxFactor) + jPz_prevIterate_*relaxFactor

	; Write the iteration error to a file for analysis

	nc_id = nCdf_create ('output/kj_itErr_'+cfg.runIdent+'.nc', /clobber )

		nCdf_control, nc_id, /verbose

		nr_id = nCdf_dimDef ( nc_id, 'nR', n_elements(r) )
		r_id = nCdf_varDef ( nc_id, 'r', nr_id, /float )

		jP_err_re_id = nCdf_varDef ( nc_id, 'jP_err_re', nr_id, /float )
		jP_err_im_id = nCdf_varDef ( nc_id, 'jP_err_im', nr_id, /float )

		nCdf_control, nc_id, /enDef

		nCdf_varPut, nc_id, r_id, r
		jP_err = jROut - jPr_prevIterate
		nCdf_varPut, nc_id, jP_err_re_id, real_part(jP_err)
		nCdf_varPut, nc_id, jP_err_im_id, imaginary(jP_err) 

	nCdf_close, nc_id


	; Plot comparison with previous iterate

	xrange = [min(r),max(r)]

	if not keyword_set(sig33) then begin
	if(not sheath)then begin
		c_pb_re=plot(r_prevIterate,jPr_prevIterate,$
				thick=3.0,$
				xrange=xRange,$
				name='prevIterate_re',$
				color='b',$
				window_title='kj',$
				buffer=1,$
				dimensions=[1200,400])
		c_pb_im=plot(r_prevIterate,imaginary(jPr_prevIterate),thick=2.0,xrange=xRange,/over,name='prevIterate_im',color='b')
	endif	
	if(sheath)then begin
		pk_re=plot(xF/lambda_D,j1,thick=3.0,name='kj_re',color='black')
		pk_im=plot(xF/lambda_D,imaginary(j1),/over,color='black',thick=2.0,name='kj_im')
	endif else begin
		pk_re=plot(xF,j1,thick=3.0,name='kj_re',color='black',/over)
		pk_im=plot(xF,imaginary(j1),/over,color='black',thick=2.0,name='kj_im')
	endelse

	if(not sheath)then begin
	l=legend(target=[c_pb_re,c_pb_im,pk_re,pk_im],$
			position=[0.98,0.9],/norm,font_size=10,horizontal_alignment='RIGHT')
	endif
	c_pb_re.save, 'kj_jP.png', resolution=72
	endif

	; Write kj_jP in file for next iterate

	nc_id = nCdf_create ('output/kj_jP_'+cfg.runIdent+'.nc', /clobber )

	nCdf_control, nc_id, /fill
	
	nr_id = nCdf_dimDef ( nc_id, 'nR', n_elements(r) )
	nrH_id = nCdf_dimDef ( nc_id, 'nR_', n_elements(r_) )
	scalar_id = nCdf_dimDef ( nc_id, 'scalar', 1 )

	freq_id = nCdf_varDef ( nc_id, 'freq', scalar_id, /float )
	r_id = nCdf_varDef ( nc_id, 'r', nr_id, /float )
	rH_id = nCdf_varDef ( nc_id, 'r_', nrH_id, /float )

	jP_r_re_id = nCdf_varDef ( nc_id, 'jP_r_re', nr_id, /float )
	jP_r_im_id = nCdf_varDef ( nc_id, 'jP_r_im', nr_id, /float )
	jP_p_re_id = nCdf_varDef ( nc_id, 'jP_p_re', nr_id, /float )
	jP_p_im_id = nCdf_varDef ( nc_id, 'jP_p_im', nr_id, /float )
	jP_z_re_id = nCdf_varDef ( nc_id, 'jP_z_re', nr_id, /float )
	jP_z_im_id = nCdf_varDef ( nc_id, 'jP_z_im', nr_id, /float )

	jP_r_re_id_ = nCdf_varDef ( nc_id, 'jP_r_re_', nrH_id, /float )
	jP_r_im_id_ = nCdf_varDef ( nc_id, 'jP_r_im_', nrH_id, /float )
	jP_p_re_id_ = nCdf_varDef ( nc_id, 'jP_p_re_', nrH_id, /float )
	jP_p_im_id_ = nCdf_varDef ( nc_id, 'jP_p_im_', nrH_id, /float )
	jP_z_re_id_ = nCdf_varDef ( nc_id, 'jP_z_re_', nrH_id, /float )
	jP_z_im_id_ = nCdf_varDef ( nc_id, 'jP_z_im_', nrH_id, /float )

	nCdf_control, nc_id, /enDef

	nCdf_varPut, nc_id, freq_id, freq

	nCdf_varPut, nc_id, r_id, r
	nCdf_varPut, nc_id, rH_id, r_

	nCdf_varPut, nc_id, jP_r_re_id, real_part(jROut)
	nCdf_varPut, nc_id, jP_r_im_id, imaginary(jROut) 
	nCdf_varPut, nc_id, jP_p_re_id, real_part(jTOut)
	nCdf_varPut, nc_id, jP_p_im_id, imaginary(jTOut) 
	nCdf_varPut, nc_id, jP_z_re_id, real_part(jZOut)
	nCdf_varPut, nc_id, jP_z_im_id, imaginary(jZOut) 

	nCdf_varPut, nc_id, jP_r_re_id_, real_part(jROut_)
	nCdf_varPut, nc_id, jP_r_im_id_, imaginary(jROut_) 
	nCdf_varPut, nc_id, jP_p_re_id_, real_part(jTOut_)
	nCdf_varPut, nc_id, jP_p_im_id_, imaginary(jTOut_) 
	nCdf_varPut, nc_id, jP_z_re_id_, real_part(jZOut_)
	nCdf_varPut, nc_id, jP_z_im_id_, imaginary(jZOut_) 

	nCdf_close, nc_id

endif

	; Interpolate the E field to the Jp locations to calculated sig33
	E_at_Jp = complex(interpol(er_re,r,xF ,/spline),interpol(er_im,r,xF ,/spline)) 

	sig33 = j1/E_at_Jp

    print, 'Sig33: ', sig33

	restore, 'AnalyticSig33.sav'
	this_margin = 0.2
	p=plot(SPoints, SPoints_sig33,Layout=[1,3,1],title='sig33',$
			xRange=[min(SPoints),max(SPoints)],$
			margin=this_margin, color='grey',thick=2)
	p=plot(SPoints, imaginary(SPoints_sig33),/over,color='r',thick=2)

	SPoints_sig33_FApprox = SPoints_sig33_FApprox/1e6
	p=plot(SPoints, SPoints_sig33_FApprox,/over,thick=2)
	p=plot(SPoints, imaginary(SPoints_sig33_FApprox),/over,color='r',thick=2)

	p=plot(xf, sig33,/over,thick=4,color='dark slate grey')
	p=plot(xf, imaginary(sig33),/over,thick=4,color='orange red')
	p=plot(s_Coord, kb, Layout=[1,3,2], /current, thick=2, color='g',title='local kPrl',$
			xRange=[min(SPoints),max(SPoints)],margin=this_margin)
	p=plot(s_Coord, Eb, Layout=[1,3,3], /current, thick=2, title='ePrl',$
			xRange=[min(SPoints),max(SPoints)],$
			margin = this_margin)
	p=plot(s_Coord, imaginary(Eb), /over, thick=2, color='r',margin=this_margin)

	p.save, 'kj_sig33_analytic_comparison.png', resolution=128
	p.save, 'kj_sig33_analytic_comparison.pdf', /close


	; try to reconstruct some Z-function data

   	wpe  = sqrt(n_e*_e^2/(me*e0))
	k_ = interpol(kb,s_coord,xF,/spline) 
	K3 = -sig33/(II*wrf*e0)+1	
	Zeta_Zp = -(K3-1)/(wpe^2)*wrf*k_*vTh
 	zeta = wrf/(k_*vTh)
	Zp = Zeta_Zp / zeta

	Z_ = ComplexArr(n_elements(zeta))
	Zp_ = ComplexArr(n_elements(zeta))	
	for i=0,n_elements(zeta)-1 do begin
		Z_[i] = kj_zfunction(zeta[i],k_[i]/abs(k_[i]),Zp=tmp)
		Zp_[i] = tmp
	endfor

	xrange=[-5,5]
	p=plot(zeta,zP,LineStyle=6,Symbol="o",xrange=xrange)
	p=plot(zeta,imaginary(zP),LineStyle=6,Symbol="o",sym_color="r",/over)

	iiSorted = sort(zeta)

	p=plot(zeta[iiSorted],Zp_[iiSorted],/over,thick=2,transparency=50)
	p=plot(zeta[iiSorted],imaginary(Zp_[iiSorted]),/over,thick=2,transparency=50,color='r')
stop

end
