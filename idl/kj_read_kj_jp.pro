function kj_read_kj_jp, FileName

; Read in jP on some grid

    print, fileName

    cdfId = ncdf_open(FileName)

    nCdf_varGet, cdfId, 'freq', freq
    nCdf_varGet, cdfId, 'r', kj_r

	nCdf_varGet, cdfId, 'jP_r_re', kj_jPr_re
	nCdf_varGet, cdfId, 'jP_r_im', kj_jPr_im
	nCdf_varGet, cdfId, 'jP_t_re', kj_jPt_re
	nCdf_varGet, cdfId, 'jP_t_im', kj_jPt_im
	nCdf_varGet, cdfId, 'jP_z_re', kj_jPz_re
	nCdf_varGet, cdfId, 'jP_z_im', kj_jPz_im

	ncdf_close, cdfId

    jpR = complex(kj_jPr_re,kj_jPr_im)
    jpT = complex(kj_jPt_re,kj_jPt_im)
    jpZ = complex(kj_jPz_re,kj_jPz_im)

    kj_r_ = kj_r[0:-2] + (kj_r[1]-kj_r[0])/2

    nR = n_elements(kj_r_)
    nS = n_elements(kj_jPr_re[0,0,*])
    nZ = n_elements(kj_jPr_re[0,*,0])

    jPr_ = complexArr(nR,nZ,nS)
    jPt_ = complexArr(nR,nZ,nS)
    jPz_ = complexArr(nR,nZ,nS)

    spline_sigma = 0.01

    ; Interpolate to the half grid

    for s=0,nS-1 do begin

        jPr_[*,0,s] = complex(spline(kj_r,kj_jpR_re[*,0,s],kj_r_,spline_sigma),spline(kj_r,kj_jpR_im[*,0,s],kj_r_,spline_sigma))
        jPt_[*,0,s] = complex(spline(kj_r,kj_jpT_re[*,0,s],kj_r_,spline_sigma),spline(kj_r,kj_jpT_im[*,0,s],kj_r_,spline_sigma))
        jPz_[*,0,s] = complex(spline(kj_r,kj_jpZ_re[*,0,s],kj_r_,spline_sigma),spline(kj_r,kj_jpZ_im[*,0,s],kj_r_,spline_sigma))

    endfor

    kjIn = { $
        freq : freq, $
        r : kj_r, $
        r_ : kj_r_, $
        jPr : jPr, $
        jPt : jPt, $
        jPz : jPz, $
        jPr_ : jPr_, $
        jPt_ : jPt_, $
        jPz_ : jPz_ $
   }

    return, kjIn

end
