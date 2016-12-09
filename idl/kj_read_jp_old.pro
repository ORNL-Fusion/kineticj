pro kj_read_jp_old, x=xF, j1x=j1x, j1y=j1y, j1z=j1z
	
	fileList = file_search ( 'output/jP*' )

	cdfId = ncdf_open(fileList[0])
	ncdf_varget, cdfId, 't', t
	nCdf_close, cdfId
	
	nT = n_elements(t)
	nF = n_elements(fileList)

	j1x = fltArr ( nT, nF )

	j1xc = complexArr ( nF )
	j1yc = complexArr ( nF )
	j1zc = complexArr ( nF )

	j1x = complexArr ( nF )
	j1y = complexArr ( nF )
	j1z = complexArr ( nF )

	xF = fltArr ( nF )

	for f=0,nF-1 do begin

		cdfId = ncdf_open(fileList[f])

			ncdf_varget, cdfId, 'x', x
			ncdf_varget, cdfId, 't', t

			ncdf_varget, cdfId, 'j1xc_re', j1xc_re 
			ncdf_varget, cdfId, 'j1xc_im', j1xc_im

			ncdf_varget, cdfId, 'j1yc_re', j1yc_re 
			ncdf_varget, cdfId, 'j1yc_im', j1yc_im

			ncdf_varget, cdfId, 'j1zc_re', j1zc_re 
			ncdf_varget, cdfId, 'j1zc_im', j1zc_im

		nCdf_close,	cdfId 

		xF[f] = x

		j1xc[f] = complex(j1xc_re,j1xc_im)
		j1yc[f] = complex(j1yc_re,j1yc_im)
		j1zc[f] = complex(j1zc_re,j1zc_im)

        j1x[f] = j1xc[f]
        j1y[f] = j1yc[f]
        j1z[f] = j1zc[f]
	
	endfor

end
