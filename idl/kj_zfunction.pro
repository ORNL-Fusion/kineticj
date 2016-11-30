function kj_dn, n ; Used in asymptotic series for arg >> 1

    d = dblArr(n+1)
    d[0] = 1

    if n eq 0 then begin
        return, d[0]
    endif else begin
       
        for _n=1,n do begin
            d[_n] = (2d0*_n+1d0)/2d0 * d[_n-1]
        endfor

    endelse
    return, d[n]

end

function kj_zfunction, arg, Zp=Zp

    ; Plasma dispersion function, either table lookup using 
    ; mathematica created table (see zFunction_table_creator.nb in kinetic-j repo)
    ; or using the asymptotic series for large argument where it's
    ; off the end of the table

    @constants

    zfunFileName = expand_path('~/code/kineticj/idl/zFunction.nc')
	cdfId = ncdf_open(zFunFileName)

		ncdf_varget, cdfId, 'arg_re', arg_re
		ncdf_varget, cdfId, 'Z_re', Z_re
		ncdf_varget, cdfId, 'Z_im', Z_im
		ncdf_varget, cdfId, 'Zp_re', Zp_re
		ncdf_varget, cdfId, 'Zp_im', Zp_im

	ncdf_close, cdfId

    ; First interpolate from table

	this_Z_re = interpol(Z_re,arg_re,real_part(arg),/spline)
	this_Z_im = interpol(Z_im,arg_re,real_part(arg),/spline)
	this_Zp_re = interpol(Zp_re,arg_re,real_part(arg),/spline)
	this_Zp_im = interpol(Zp_im,arg_re,real_part(arg),/spline)

	Z = complex(this_Z_re,this_Z_im)
	Zp = complex(this_Zp_re,this_Zp_im)

    ; Then overwrite elements off table with asymptotic series for arg >> 1

    iiLarge = where(abs(arg) gt max(arg_re), iiLargeCnt)

    ; For large argument (off the end of the table)

	if iiLargeCnt gt 0 then begin
    
        nMax = 40

        for i=0,iiLargeCnt-1 do begin

            ; Asymptotic series for x >> 1

            print, arg[iiLarge[i]]
            stop

            sum = 0
            for n=0,nMax do begin
                sum += kj_dn(n) / (arg[iiLarge[i]]^(2*n))
            endfor
            sig = 1
            if imaginary(arg[iiLarge[i]]) gt 0 then sig = 0
            if imaginary(arg[iiLarge[i]]) lt 0 then sig = 2

            Z[iiLarge[i]] = _ii * sig * sqrt(!dpi) * exp(-arg[iiLarge[i]]^2) - sum / arg[iiLarge[i]]
            Zp[iiLarge[i]] = -2d0 * ( 1d0 + arg[iiLarge[i]] * Z[iiLarge[i]] )

        endfor

	endif 

	return, Z

end


