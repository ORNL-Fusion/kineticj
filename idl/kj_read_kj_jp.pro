function kj_read_kj_jp, FileName;, PerSpecies = PerSpecies

    cdfId = ncdf_open(FileName)

    nCdf_varGet, cdfId, 'freq', freq
    nCdf_varGet, cdfId, 'r', r
    nCdf_varGet, cdfId, 'r_', r_

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

	kj_jpR = complex(kj_jpR_re,kj_jpR_im)
	kj_jpT = complex(kj_jpT_re,kj_jpT_im)
	kj_jpZ = complex(kj_jpZ_re,kj_jpZ_im)

    ; Interpolate to the half grid
    spline_sigma = 0.01
    kj_jpR_ = complex(spline(r,real_part(kj_jpR),r_,spline_sigma),spline(r,imaginary(kj_jpR),r_,spline_sigma))
    kj_jpT_ = complex(spline(r,real_part(kj_jpT),r_,spline_sigma),spline(r,imaginary(kj_jpT),r_,spline_sigma))
    kj_jpZ_ = complex(spline(r,real_part(kj_jpZ),r_,spline_sigma),spline(r,imaginary(kj_jpZ),r_,spline_sigma))


    nS = n_elements(kj_jpR_spec_re[0,*])

	kj_jpR_spec =  complex(kj_jpR_spec_re, kj_jpR_spec_im)
	kj_jpT_spec =  complex(kj_jpT_spec_re, kj_jpT_spec_im)
	kj_jpZ_spec =  complex(kj_jpZ_spec_re, kj_jpZ_spec_im)

	kj_jpR_spec =  complex(kj_jpR_spec_re, kj_jpR_spec_im)
	kj_jpT_spec =  complex(kj_jpT_spec_re, kj_jpT_spec_im)
	kj_jpZ_spec =  complex(kj_jpZ_spec_re, kj_jpZ_spec_im)

    kj_jpR_spec_ = complexArr(nR-1,nS) 
	kj_jpT_spec_ = complexArr(nR-1,nS)
	kj_jpZ_spec_ = complexArr(nR-1,nS)

    for s=0,nS-1 do begin
        kj_jpR_spec_[*,s] = complex(spline(r,real_part(kj_jpR_spec[*,s]),r_,spline_sigma),spline(r,imaginary(kj_jpR_spec[*,s]),r_,spline_sigma))
        kj_jpT_spec_[*,s] = complex(spline(r,real_part(kj_jpT_spec[*,s]),r_,spline_sigma),spline(r,imaginary(kj_jpT_spec[*,s]),r_,spline_sigma))
        kj_jpZ_spec_[*,s] = complex(spline(r,real_part(kj_jpZ_spec[*,s]),r_,spline_sigma),spline(r,imaginary(kj_jpZ_spec[*,s]),r_,spline_sigma))
    endfor

    kjIn = { $
            freq : freq, $
            r : r, $
            r_ : r_, $
        jpR : kj_jpR, $
        jpT : kj_jpT, $
        jpZ : kj_jpZ, $
        jpR_ : kj_jpR_, $
        jpT_ : kj_jpT_, $
        jpZ_ : kj_jpZ_, $
        jpR_spec : kj_jpR_spec, $
        jpT_spec : kj_jpT_spec, $
        jpZ_spec : kj_jpZ_spec, $
        jpR_spec_ : kj_jpR_spec_, $
        jpT_spec_ : kj_jpT_spec_, $
        jpZ_spec_ : kj_jpZ_spec_ $
    }

    return, kjIn

end
