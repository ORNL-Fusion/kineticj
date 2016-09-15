pro kj_look_at_input_efield

    kjInput = kj_read_cfg('./')

	cdfId = ncdf_open ( kjInput['eField_fName'], /noWrite ) 

        nCdf_varGet, cdfId, 'r', r

        nCdf_varGet, cdfId, 'er_re', er_re
        nCdf_varGet, cdfId, 'er_im', er_im
        nCdf_varGet, cdfId, 'et_re', et_re
        nCdf_varGet, cdfId, 'et_im', et_im
	    nCdf_varGet, cdfId, 'ez_re', ez_re
        nCdf_varGet, cdfId, 'ez_im', ez_im
	ncdf_close, cdfId

    er = complex(er_re,er_im)
    et = complex(et_re,et_im)
    ez = complex(ez_re,ez_im)

    p=plot(r,er,layout=[1,3,1])
    p=plot(r,imaginary(er),/over,color='r')
    p=plot(r,et,layout=[1,3,2],/current)
    p=plot(r,imaginary(et),/over,color='r')
    p=plot(r,ez,layout=[1,3,3],/current)
    p=plot(r,imaginary(ez),/over,color='r')

stop

end
