pro kj_plot_current

	@constants

	;cdfId = ncdf_open('data/rsfwc_1d.nc')
	cdfId = ncdf_open('data/kj_aorsa_1d.nc')

		ncdf_varget, cdfId, 'freq', freq 
		ncdf_varget, cdfId, 'r', r 

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

	ncdf_close, cdfId

	wrf = freq * 2 * !pi

	j1_cold = complex ( jr_re, jr_im ) 

	fileList = file_search ( 'output/jP*' )

	cdfId = ncdf_open(fileList[0])
	ncdf_varget, cdfId, 't', t
	nCdf_close, cdfId
	
	nT = n_elements(t)
	nF = n_elements(fileList)

	j1x = fltArr ( nT, nF )
	j1 = complexArr ( nF )
	xF = fltArr ( nF )

	hanWindow = hanning (nT, alpha=0.5 )

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
		rp = -!pi*real_part ( jrFFT[iiAntFreq[0]] ) ; not sure exactly where this -!pi comes from :s
		ip = !pi*imaginary ( jrFFT[iiAntFreq[0]] )

		j1[f] = complex ( rp, ip )

	endfor

	fudgeFac = 1.0 
	xrange = [min(r),max(r)]
	pb_re=plot(r,j1_cold,thick=3.0,xrange=xRange,name='base_re',transp=50)
	pb_im=plot(r,imaginary(j1_cold),thick=3.0,xrange=xRange,/over,name='base_im',color='r',transp=50)
	pk_re=plot(xF,j1*fudgeFac,/over,thick=1.0,name='kj_re')
	pk_im=plot(xF,imaginary(j1)*fudgeFac,/over,color='r',thick=1.0,name='kj_im')

	l=legend(target=[pb_re,pb_im,pk_re,pk_im],position=[0.98,0.9],/norm,font_size=10,horizontal_align='RIGHT')


stop
end
