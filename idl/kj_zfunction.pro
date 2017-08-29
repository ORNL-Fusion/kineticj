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

    ; Plasma dispersion function using  
    ; mathematica created table 
    ; (see zFunction_table_creator.nb in kinetic-j repo)
    ; Don't use the asymptoptic expansion anymore, since 
    ; we have the table out to +/- 1e5 in the Re(arg) and
    ; could even go larger with the non-uniform (log10)
    ; spacing of table. 

    @constants

    zfunFileName = expand_path('~/code/kineticj/mathematica/zFunction.nc')
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

	return, Z

end


