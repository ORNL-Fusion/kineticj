function kj_read_kj_jp, FileName, r, rH = r_

; Read in jP on some grid

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

    nR = n_elements(r)
    if keyword_set(nH) then nR_ = n_elements(r_) else nR_ = !null

	kj_jpR_spec = complex(kj_jpR_re,kj_jpR_im)
	kj_jpT_spec = complex(kj_jpT_re,kj_jpT_im)
	kj_jpZ_spec = complex(kj_jpZ_re,kj_jpZ_im)

    kj_nR = n_elements(kj_jPr_spec[*,0,0])
    kj_nZ = n_elements(kj_jPr_spec[0,*,0])
    kj_nS = n_elements(kj_jPr_spec[0,0,*])

    ;_shift = [0,0,2]
    ;kj_jPr_spec = shift(kj_jPr_spec,_shift)
    ;kj_jPt_spec = shift(kj_jPt_spec,_shift)
    ;kj_jPz_spec = shift(kj_jPz_spec,_shift)

; Now interpolate to the desired grid(s)

    iiBad = where(r lt min(kj_r)*0.9999 or r gt max(kj_r)*1.0001,iiBadCnt)
    if iiBadCnt gt 0 then stop

    spline_sigma = 0.01

    ; Interpolate to the full grid

    jpR_spec = complexArr(nR,kj_nS) 
	jpT_spec = complexArr(nR,kj_nS)
	jpZ_spec = complexArr(nR,kj_nS)

    if kj_nZ eq 1 then begin
    for s=0,kj_nS-1 do begin
        jpR_spec[*,s] = complex(spline(kj_r,real_part(kj_jpR_spec[*,0,s]),r,spline_sigma),spline(kj_r,imaginary(kj_jpR_spec[*,0,s]),r,spline_sigma))
        jpT_spec[*,s] = complex(spline(kj_r,real_part(kj_jpT_spec[*,0,s]),r,spline_sigma),spline(kj_r,imaginary(kj_jpT_spec[*,0,s]),r,spline_sigma))
        jpZ_spec[*,s] = complex(spline(kj_r,real_part(kj_jpZ_spec[*,0,s]),r,spline_sigma),spline(kj_r,imaginary(kj_jpZ_spec[*,0,s]),r,spline_sigma))
    endfor
    endif else begin
            print, 'Have not implemented 2D interp yet'
            stop
    endelse 
    ; Interpolate to the half grid
    if keyword_set(r_) then begin

        iiBad_ = where(r_ lt min(kj_r) or r_ gt max(kj_r),iiBadCnt_)
        if iiBadCnt_ gt 0 then stop

		nR_ = n_elements(r_)

        jpR_spec_ = complexArr(nR_,kj_nS) 
	    jpT_spec_ = complexArr(nR_,kj_nS)
	    jpZ_spec_ = complexArr(nR_,kj_nS)

        if kj_nZ eq 1 then begin
        for s=0,kj_nS-1 do begin
            jpR_spec_[*,s] = complex(spline(kj_r,real_part(kj_jpR_spec[*,0,s]),r_,spline_sigma),spline(kj_r,imaginary(kj_jpR_spec[*,0,s]),r_,spline_sigma))
            jpT_spec_[*,s] = complex(spline(kj_r,real_part(kj_jpT_spec[*,0,s]),r_,spline_sigma),spline(kj_r,imaginary(kj_jpT_spec[*,0,s]),r_,spline_sigma))
            jpZ_spec_[*,s] = complex(spline(kj_r,real_part(kj_jpZ_spec[*,0,s]),r_,spline_sigma),spline(kj_r,imaginary(kj_jpZ_spec[*,0,s]),r_,spline_sigma))
        endfor
        endif else begin
            print, 'Have not implemented 2D interp yet'
            stop
        endelse

	endif


	if not keyword_set(r_)then begin
        r_=!null
        jPr_spec_ = !null
        jPt_spec_ = !null
        jPz_spec_ = !null
     
    endif

    if size(jPr_spec,/n_dim) le 1 then begin
        jPr = jPr_spec
        jPt = jPt_spec
        jPz = jPz_spec
        jPr_ = jPr_spec_
        jPt_ = jPt_spec_
        jPz_ = jPz_spec_
    endif else begin
        sumDim = 2
        jPr = total(jPr_spec,sumDim)
        jPt = total(jPt_spec,sumDim)
        jPz = total(jPz_spec,sumDim)
        jPr_ = total(jPr_spec_,sumDim)
        jPt_ = total(jPt_spec_,sumDim)
        jPz_ = total(jPz_spec_,sumDim)
    endelse

    kjIn = { $
        freq : freq, $
        r : r, $
        r_ : r_, $
        jPr_spec : jPr_spec, $
        jPt_spec : jPt_spec, $
        jPz_spec : jPz_spec, $
        jPr_spec_ : jPr_spec_, $
        jPt_spec_ : jPt_spec_, $
        jPz_spec_ : jPz_spec_, $
        jPr : jPr, $
        jPt : jPt, $
        jPz : jPz, $
        jPr_ : jPr_, $
        jPt_ : jPt_, $
        jPz_ : jPz_ $
   }

    return, kjIn

end
