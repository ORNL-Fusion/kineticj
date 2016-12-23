pro kj_read_jp_old, x=xF, j1x=j1x, j1y=j1y, j1z=j1z, oldFormat=_oldFormat, fileName=_fileName

if keyword_set(_oldFormat) then begin

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

endif else begin

    if keyword_set(_fileName) then fileName = _fileName else fileName = 'jP2.nc' 

	cdfId = ncdf_open(fileName)

		ncdf_varget, cdfId, 'x', x2 

		ncdf_varget, cdfId, 'j1xc_re', jPr_re2
		ncdf_varget, cdfId, 'j1xc_im', jPr_im2
		ncdf_varget, cdfId, 'j1yc_re', jPt_re2
		ncdf_varget, cdfId, 'j1yc_im', jPt_im2
		ncdf_varget, cdfId, 'j1zc_re', jPz_re2
		ncdf_varget, cdfId, 'j1zc_im', jPz_im2

	ncdf_close, cdfId

    xF = x2

    j1x=complex(jPr_re2,jPr_im2)
    j1y=complex(jPt_re2,jPt_im2)
    j1z=complex(jPz_re2,jPz_im2)

endelse

end
