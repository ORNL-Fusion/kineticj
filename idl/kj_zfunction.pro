function kj_zfunction, arg, Sgn_kPrl, Zp=Zp

	if abs(imaginary(arg)) gt 0 then stop

	if Sgn_kPrl ge 0 then begin
		zfunFileName = 'zFunction_kPrl_p.nc'
	endif else begin
		zfunFileName = 'zFunction_kPrl_n.nc'
	endelse
	cdfId = ncdf_open(zFunFileName)

		ncdf_varget, cdfId, 'arg_re', arg_re
		ncdf_varget, cdfId, 'Z_re', Z_re
		ncdf_varget, cdfId, 'Z_im', Z_im
		ncdf_varget, cdfId, 'Zp_re', Zp_re
		ncdf_varget, cdfId, 'Zp_im', Zp_im

	ncdf_close, cdfId

	if arg gt max(arg_re) or arg lt min(arg_re) then begin
			print, 'ERROR: Z-function argument off the end of the interpolation table.'
			print, 'ERROR: Please re-run the mathematica for appropriate range.'
			stop
	endif

	this_Z_re = interpol(Z_re,arg_re,real_part(arg),/spline)
	this_Z_im = interpol(Z_im,arg_re,real_part(arg),/spline)
	this_Zp_re = interpol(Zp_re,arg_re,real_part(arg),/spline)
	this_Zp_im = interpol(Zp_im,arg_re,real_part(arg),/spline)

	Z = complex(this_Z_re,this_Z_im)
	Zp = complex(this_Zp_re,this_Zp_im)

	return, Z

end


