pro kj_plot_current

	@constants

	runFile = 'kj.cfg'
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

	cdfId = ncdf_open(eField_fName)

		ncdf_varget, cdfId, 'freq', freq 
		ncdf_varget, cdfId, 'r', r 
		ncdf_varget, cdfId, 'r_', r_

		ncdf_varget, cdfId, 'e_r_re', er_re
		ncdf_varget, cdfId, 'e_r_im', er_im
		ncdf_varget, cdfId, 'e_p_re', ep_re
		ncdf_varget, cdfId, 'e_p_im', ep_im
		ncdf_varget, cdfId, 'e_z_re', ez_re
		ncdf_varget, cdfId, 'e_z_im', ez_im

		ncdf_varget, cdfId, 'j_r_re', jr_re
		ncdf_varget, cdfId, 'j_r_im', jr_im
		ncdf_varget, cdfId, 'j_p_re', jp_re
		ncdf_varget, cdfId, 'j_p_im', jp_im
		ncdf_varget, cdfId, 'j_z_re', jz_re
		ncdf_varget, cdfId, 'j_z_im', jz_im

		ncdf_varget, cdfId, 'jA_r_re', jAr_re
		ncdf_varget, cdfId, 'jA_r_im', jAr_im
		ncdf_varget, cdfId, 'jA_p_re', jAp_re
		ncdf_varget, cdfId, 'jA_p_im', jAp_im
		ncdf_varget, cdfId, 'jA_z_re', jAz_re
		ncdf_varget, cdfId, 'jA_z_im', jAz_im

	ncdf_close, cdfId

	cdfId = ncdf_open('data/kj_aorsa_1d_0.5keV.nc')

		ncdf_varget, cdfId, 'freq', ao_freq 
		ncdf_varget, cdfId, 'r', ao_r 

		ncdf_varget, cdfId, 'e_r_re', ao_er_re
		ncdf_varget, cdfId, 'e_r_im', ao_er_im
		ncdf_varget, cdfId, 'e_p_re', ao_ep_re
		ncdf_varget, cdfId, 'e_p_im', ao_ep_im
		ncdf_varget, cdfId, 'e_z_re', ao_ez_re
		ncdf_varget, cdfId, 'e_z_im', ao_ez_im

		ncdf_varget, cdfId, 'j_r_re', ao_jr_re
		ncdf_varget, cdfId, 'j_r_im', ao_jr_im
		ncdf_varget, cdfId, 'j_p_re', ao_jp_re
		ncdf_varget, cdfId, 'j_p_im', ao_jp_im
		ncdf_varget, cdfId, 'j_z_re', ao_jz_re
		ncdf_varget, cdfId, 'j_z_im', ao_jz_im

		ncdf_varget, cdfId, 'jA_r_re', ao_jAr_re
		ncdf_varget, cdfId, 'jA_r_im', ao_jAr_im
		ncdf_varget, cdfId, 'jA_p_re', ao_jAp_re
		ncdf_varget, cdfId, 'jA_p_im', ao_jAp_im
		ncdf_varget, cdfId, 'jA_z_re', ao_jAz_re
		ncdf_varget, cdfId, 'jA_z_im', ao_jAz_im

	ncdf_close, cdfId


	wrf = freq * 2 * !pi

	r_hot = ao_r
	j1_hot = complex ( ao_jr_re, ao_jr_im ) ;+ complex ( ao_jAr_re, ao_jAr_im ) 

	r_cold = r
	j1_cold = complex ( jr_re, jr_im ) ;+ complex ( jAr_re, jAr_im ) 

	fileList = file_search ( 'output/jP*' )

	cdfId = ncdf_open(fileList[0])
	ncdf_varget, cdfId, 't', t
	nCdf_close, cdfId
	
	nT = n_elements(t)
	nF = n_elements(fileList)

	j1x = fltArr ( nT, nF )
	j1 = complexArr ( nF )
	xF = fltArr ( nF )

	;hanWindow = hanning (nT, alpha=0.5 )

	for f=0,nF-1 do begin

		cdfId = ncdf_open(fileList[f])

			ncdf_varget, cdfId, 'x', x
			ncdf_varget, cdfId, 't', t

			ncdf_varget, cdfId, 'j1x', j1x_0 
			ncdf_varget, cdfId, 'j1y', j1y_0 
			ncdf_varget, cdfId, 'j1z', j1z_0 

		nCdf_close,	cdfId 

		j1x[*,f] = j1x_0-mean(j1x_0) ; not sure if there is a nicer way to do this
		xF[f] = x

		dt = t[1]-t[0]
		;jrFFT = fft ( hanWindow * j1x[*,f] )
		;ampRe = 1000.0
		;ampIm = 0.0
		;if x gt 99.95 then stop
		;j1x_0 = ampRe * cos ( wrf * t ) + ampIm * ii * sin ( wrf * t )
		jrFFT = fft ( j1x[0:-2,f] )
		freqAxis = fIndGen(nT-1) / ((nT-1)*dt)
		freqNorm = freqAxis / freq

		;p=plot(freqNorm,abs(jrFFT)^2,xRange=[0,5])
		;p=plot(freqNorm,imaginary(jrFFT),color='blue',xRange=[0,5])

		iiAntFreq = where(abs(freqNorm-1) eq min(abs(freqNorm-1)),iiAntFreqCnt)
		rp = real_part ( jrFFT[iiAntFreq[0]] )
		ip = imaginary ( jrFFT[iiAntFreq[0]] )

		j1[f] = complex ( rp, ip )

	endfor

	fudgeFac = 1.0*!pi ; Not sure why we need a pi here, most likely IDLs fft.

	; Create a jP for rsfcw_1d

	jAR = complex(jAr_re,jAr_im)
	jAT = complex(jAp_re,jAp_im)
	jAZ = complex(jAz_re,jAz_im)

	jAR_ = interpol ( jAR, r, r_ ) 
	jAT_ = interpol ( jAT, r, r_ ) 
	jAZ_ = interpol ( jAZ, r, r_ ) 

	jROut = interpol ( j1*fudgeFac, xF, r ) ;- jAR
	jROut_ = interpol ( j1*fudgeFac, xF, r_ ) ;- jAR_

	jTOut = jROut*0
	jTOut_ = jROut_*0

	jZOut = jROut*0
	jZOut_ = jROut_*0

	; Plot comparison with previous iterate

	xrange = [min(r),max(r)]

	c_pb_re=plot(r_cold,j1_cold,thick=3.0,xrange=xRange,name='cold_re',transp=50,color='b',$
			window_title='kj')
	c_pb_im=plot(r_cold,imaginary(j1_cold),thick=2.0,xrange=xRange,/over,name='cold_im',transp=50,color='b')

	h_pb_re=plot(r_hot,j1_hot,thick=3.0,name='hot_re',transp=50,color='r',/over)
	h_pb_im=plot(r_hot,imaginary(j1_hot),thick=2.0,/over,name='hot_im',color='r',transp=50)
	
	;pk_re=plot(xF,j1*fudgeFac,/over,thick=3.0,name='kj_re',color='black')
	;pk_im=plot(xF,imaginary(j1*fudgeFac),/over,color='black',thick=2.0,name='kj_im',transp=50)
	pk_re=plot(r,jROut,/over,thick=3.0,name='kj_re',color='black')
	pk_im=plot(r,imaginary(jROut),/over,color='black',thick=2.0,name='kj_im',transp=50)


	l=legend(target=[c_pb_re,c_pb_im,h_pb_re,h_pb_im,pk_re,pk_im],$
			position=[0.98,0.9],/norm,font_size=10,horizontal_align='RIGHT')

	; Write kj_jP in file for next iterate

	nc_id = nCdf_create ('data/kjInput.nc', /clobber )

	nCdf_control, nc_id, /fill
	
	nr_id = nCdf_dimDef ( nc_id, 'nR', n_elements(r) )
	nrH_id = nCdf_dimDef ( nc_id, 'nR_', n_elements(r_) )
	scalar_id = nCdf_dimDef ( nc_id, 'scalar', 1 )

	freq_id = nCdf_varDef ( nc_id, 'freq', scalar_id, /float )
	r_id = nCdf_varDef ( nc_id, 'r', nr_id, /float )
	rH_id = nCdf_varDef ( nc_id, 'r_', nrH_id, /float )

	jP_r_re_id = nCdf_varDef ( nc_id, 'kj_jP_r_re', nr_id, /float )
	jP_r_im_id = nCdf_varDef ( nc_id, 'kj_jP_r_im', nr_id, /float )
	jP_p_re_id = nCdf_varDef ( nc_id, 'kj_jP_p_re', nr_id, /float )
	jP_p_im_id = nCdf_varDef ( nc_id, 'kj_jP_p_im', nr_id, /float )
	jP_z_re_id = nCdf_varDef ( nc_id, 'kj_jP_z_re', nr_id, /float )
	jP_z_im_id = nCdf_varDef ( nc_id, 'kj_jP_z_im', nr_id, /float )

	jP_r_re_id_ = nCdf_varDef ( nc_id, 'kj_jP_r_re_', nrH_id, /float )
	jP_r_im_id_ = nCdf_varDef ( nc_id, 'kj_jP_r_im_', nrH_id, /float )
	jP_p_re_id_ = nCdf_varDef ( nc_id, 'kj_jP_p_re_', nrH_id, /float )
	jP_p_im_id_ = nCdf_varDef ( nc_id, 'kj_jP_p_im_', nrH_id, /float )
	jP_z_re_id_ = nCdf_varDef ( nc_id, 'kj_jP_z_re_', nrH_id, /float )
	jP_z_im_id_ = nCdf_varDef ( nc_id, 'kj_jP_z_im_', nrH_id, /float )

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

stop
end