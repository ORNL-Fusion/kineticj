function kj_read_kj_jp, FileName, r, rH=r_;, PerSpecies = PerSpecies

; Read in jP on some grid

    cdfId = ncdf_open(FileName)

    nCdf_varGet, cdfId, 'freq', freq
    nCdf_varGet, cdfId, 'r', kj_r

	nCdf_varGet, cdfId, 'jP_r_re', kj_jpR_re
	nCdf_varGet, cdfId, 'jP_r_im', kj_jpR_im
	nCdf_varGet, cdfId, 'jP_p_re', kj_jpT_re
	nCdf_varGet, cdfId, 'jP_p_im', kj_jpT_im
	nCdf_varGet, cdfId, 'jP_z_re', kj_jpZ_re
	nCdf_varGet, cdfId, 'jP_z_im', kj_jpZ_im

	nCdf_varGet, cdfId, 'jP_r_re_spec',  kj_jpR_spec_re
	nCdf_varGet, cdfId, 'jP_r_im_spec',  kj_jpR_spec_im
	nCdf_varGet, cdfId, 'jP_p_re_spec',  kj_jpT_spec_re
	nCdf_varGet, cdfId, 'jP_p_im_spec',  kj_jpT_spec_im
	nCdf_varGet, cdfId, 'jP_z_re_spec',  kj_jpZ_spec_re
	nCdf_varGet, cdfId, 'jP_z_im_spec',  kj_jpZ_spec_im

	ncdf_close, cdfId

    nR = n_elements(r)
    nS = n_elements(kj_jpR_spec_re[0,*])

	kj_jpR = complex(kj_jpR_re,kj_jpR_im)
	kj_jpT = complex(kj_jpT_re,kj_jpT_im)
	kj_jpZ = complex(kj_jpZ_re,kj_jpZ_im)

	kj_jpR_spec =  complex(kj_jpR_spec_re, kj_jpR_spec_im)
	kj_jpT_spec =  complex(kj_jpT_spec_re, kj_jpT_spec_im)
	kj_jpZ_spec =  complex(kj_jpZ_spec_re, kj_jpZ_spec_im)

	kj_jpR_spec =  complex(kj_jpR_spec_re, kj_jpR_spec_im)
	kj_jpT_spec =  complex(kj_jpT_spec_re, kj_jpT_spec_im)
	kj_jpZ_spec =  complex(kj_jpZ_spec_re, kj_jpZ_spec_im)

; Now interpolate to the desired grid(s)

    spline_sigma = 0.01

    ; Interpolate to the full grid

    jpR = complex(spline(kj_r,real_part(kj_jpR),r,spline_sigma),spline(kj_r,imaginary(kj_jpR),r,spline_sigma))
    jpT = complex(spline(kj_r,real_part(kj_jpT),r,spline_sigma),spline(kj_r,imaginary(kj_jpT),r,spline_sigma))
    jpZ = complex(spline(kj_r,real_part(kj_jpZ),r,spline_sigma),spline(kj_r,imaginary(kj_jpZ),r,spline_sigma))

    jpR_spec = complexArr(nR,nS) 
	jpT_spec = complexArr(nR,nS)
	jpZ_spec = complexArr(nR,nS)

    for s=0,nS-1 do begin
        jpR_spec[*,s] = complex(spline(kj_r,real_part(kj_jpR_spec[*,s]),r,spline_sigma),spline(kj_r,imaginary(kj_jpR_spec[*,s]),r,spline_sigma))
        jpT_spec[*,s] = complex(spline(kj_r,real_part(kj_jpT_spec[*,s]),r,spline_sigma),spline(kj_r,imaginary(kj_jpT_spec[*,s]),r,spline_sigma))
        jpZ_spec[*,s] = complex(spline(kj_r,real_part(kj_jpZ_spec[*,s]),r,spline_sigma),spline(kj_r,imaginary(kj_jpZ_spec[*,s]),r,spline_sigma))
    endfor
 
    ; Interpolate to the half grid
    if keyword_set(r_) then begin

		nR_ = n_elements(r_)

        jpR_ = complex(spline(kj_r,real_part(kj_jpR),r_,spline_sigma),spline(kj_r,imaginary(kj_jpR),r_,spline_sigma))
        jpT_ = complex(spline(kj_r,real_part(kj_jpT),r_,spline_sigma),spline(kj_r,imaginary(kj_jpT),r_,spline_sigma))
        jpZ_ = complex(spline(kj_r,real_part(kj_jpZ),r_,spline_sigma),spline(kj_r,imaginary(kj_jpZ),r_,spline_sigma))

        jpR_spec_ = complexArr(nR_,nS) 
	    jpT_spec_ = complexArr(nR_,nS)
	    jpZ_spec_ = complexArr(nR_,nS)

        for s=0,nS-1 do begin
            jpR_spec_[*,s] = complex(spline(kj_r,real_part(kj_jpR_spec[*,s]),r_,spline_sigma),spline(kj_r,imaginary(kj_jpR_spec[*,s]),r_,spline_sigma))
            jpT_spec_[*,s] = complex(spline(kj_r,real_part(kj_jpT_spec[*,s]),r_,spline_sigma),spline(kj_r,imaginary(kj_jpT_spec[*,s]),r_,spline_sigma))
            jpZ_spec_[*,s] = complex(spline(kj_r,real_part(kj_jpZ_spec[*,s]),r_,spline_sigma),spline(kj_r,imaginary(kj_jpZ_spec[*,s]),r_,spline_sigma))
        endfor

	endif


	if not keyword_set(r_)then begin

    	kjIn = { $
    	    freq : freq, $
    	    r : r, $
    	    r_ : r_, $
    	    jpR : jpR, $
    	    jpT : jpT, $
    	    jpZ : jpZ, $
    	    jpR_spec : jpR_spec, $
    	    jpT_spec : jpT_spec, $
    	    jpZ_spec : jpZ_spec $
    	}

	endif else begin

    	kjIn = { $
    	    freq : freq, $
    	    r : r, $
    	    r_ : r_, $
    	    jpR : jpR, $
    	    jpT : jpT, $
    	    jpZ : jpZ, $
    	    jpR_ : jpR_, $
    	    jpT_ : jpT_, $
    	    jpZ_ : jpZ_, $
    	    jpR_spec : jpR_spec, $
    	    jpT_spec : jpT_spec, $
    	    jpZ_spec : jpZ_spec, $
    	    jpR_spec_ : jpR_spec_, $
    	    jpT_spec_ : jpT_spec_, $
    	    jpZ_spec_ : jpZ_spec_ $
    	}

	endelse

    return, kjIn

end
