function kj_zfunction, arg, Zp=Zp

	if abs(imaginary(arg)) gt 0 then stop

	zfunFileName = 'zFunction.nc'
	cdfId = ncdf_open(zFunFileName)

		ncdf_varget, cdfId, 'arg_re', arg_re
		ncdf_varget, cdfId, 'Z_re', Z_re
		ncdf_varget, cdfId, 'Z_im', Z_im
		ncdf_varget, cdfId, 'Zp_re', Zp_re
		ncdf_varget, cdfId, 'Zp_im', Zp_im

	ncdf_close, cdfId

	this_Z_re = interpol(Z_re,arg_re,real_part(arg),/spline)
	this_Z_im = interpol(Z_im,arg_re,real_part(arg),/spline)
	this_Zp_re = interpol(Zp_re,arg_re,real_part(arg),/spline)
	this_Zp_im = interpol(Zp_im,arg_re,real_part(arg),/spline)

	Z = complex(this_Z_re,this_Z_im)
	Zp = complex(this_Zp_re,this_Zp_im)

	return, Z

end


