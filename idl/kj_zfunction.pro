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

function kj_zfunction, zeta, Zp=Zp, zetaC=_zetaC

    ; zeta  = (w +/- wc) /(|kparalllel|*vThermal)
    ; zetaC = nuCollision/(|kparalllel|*vThermal)

    if keyword_set(_zetaC) then zetaC = _zetaC else zetaC = 0

    common zFunctionTable, Z_re, Z_im, Zp_re, Zp_im, zeta_re

    ; Plasma dispersion function using  
    ; mathematica created table 
    ; (see zFunction_table_creator.nb in kinetic-j repo)
    ; Don't use the asymptoptic expansion anymore, since 
    ; we have the table out to +/- 1e5 in the Re(zeta) and
    ; could even go lzetaer with the non-uniform (log10)
    ; spacing of table. 

    @constants

    ; If first call then read file

    if size(zeta_re,/type) eq 0 then begin

        print, 'OPENING ZFUNCTION FILE'

        zfunFileName = expand_path('~/code/kineticj/mathematica/zFunction.nc')
	    cdfId = ncdf_open(zFunFileName)

	    	ncdf_varget, cdfId, 'arg_re', zeta_re
	    	ncdf_varget, cdfId, 'Z_re', Z_re
	    	ncdf_varget, cdfId, 'Z_im', Z_im
	    	ncdf_varget, cdfId, 'Zp_re', Zp_re
	    	ncdf_varget, cdfId, 'Zp_im', Zp_im

	    ncdf_close, cdfId

    endif

    if zeta lt min(zeta_re) or zeta gt max(zeta_re) then begin
     
        ; Analytic limits off the end of the table

        Z  = complex(-1/zeta,0)
        Zp = complex(1/zeta^2,0)

    endif else begin

        ; First interpolate from table

	    this_Z_re = interpol(Z_re,zeta_re,real_part(zeta),/spline)
	    this_Z_im = interpol(Z_im,zeta_re,real_part(zeta),/spline)
	    this_Zp_re = interpol(Zp_re,zeta_re,real_part(zeta),/spline)
	    this_Zp_im = interpol(Zp_im,zeta_re,real_part(zeta),/spline)

	    Z  = complex(this_Z_re,this_Z_im)
	    Zp = complex(this_Zp_re,this_Zp_im)

    endelse

    ; From David Smithe to add collisional damping.
    ; zetaC = nuCollision/(|kparalllel|*vThermal)

    factor = complex(1.0,0.0)-complex(0.0,zetaC)*Z
    Z  = Z / factor
    Zp = Zp/ (factor*factor) 

	return, Z

end


