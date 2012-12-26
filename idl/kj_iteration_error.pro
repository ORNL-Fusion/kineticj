pro kj_iteration_error

	itErrFiles = file_search ('output/kj_itErr*.nc') 	

	nIt = n_elements(itErrFiles)

	jPerr = !null

	for i=0,nIt-1 do begin

		cdfId = ncdf_open(itErrFiles[i])

			ncdf_varget, cdfId, 'r', r
			ncdf_varget, cdfId, 'jP_err_re', jP_err_re
			ncdf_varget, cdfId, 'jP_err_im', jP_err_im

		ncdf_close, cdfId

		thisjPerr = complex(jP_err_re,jP_err_im)

		jPerr = [[jPerr],[thisjPerr]]

	endfor

	;s1=surface(real_part(jPerr))
	;s2=surface(imaginary(jPerr))

	scale = 10
	nlevels=20
	levels = (fIndGen(nLevels)+1)/(nLevels)*scale
	colors = 255-(bytScl(levels,top=253)+1)
	c1=contour(real_part(jPerr),r,fIndGen(nIt),c_value=levels,rgb_indices=colors,rgb_table=7,/fil,/buffer)
	!null=contour(-real_part(jPerr),r,fIndGen(nIt),c_value=levels,rgb_indices=colors,rgb_table=1,/fil,/over)

	c2=contour(imaginary(jPerr),r,fIndGen(nIt),c_value=levels,rgb_indices=colors,rgb_table=7,/fil,/buffer)
	!null=contour(-imaginary(jPerr),r,fIndGen(nIt),c_value=levels,rgb_indices=colors,rgb_table=1,/fil,/over)

	c1.save, 'itErr_r.png', res='100'
	c2.save, 'itErr_i.png', res='100'

	stop
end
