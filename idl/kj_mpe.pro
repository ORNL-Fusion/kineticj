function kj_mpe, x

	type = size(x,/type)

	N = n_elements(x[*,0])
	_k = n_elements(x[0,*])

	U = x[*,1:_k-2]-x[*,0:_k-3]
	uk = x[*,_k-1]-x[*,_k-2]
	k = n_elements(U[0,*])

	c = make_array(k+1,type=type)

	c[0:k-1] = la_least_squares(-transpose(U),-uk,method=3,status=stat,/double)
	; Check least squares solution
	;p=plot(-uk)
	;!null=plot(transpose(U)##c[0:k-1],/over,color='b')

	if stat ne 0 then stop

	c[k] = -1
	alpha = 0

	for i=0,k do alpha = alpha + c[i]

	_gamma = c / alpha

	s = make_array(N,type=type)

	if alpha eq 0 then stop

	for j=0,k do s = s + _gamma[j] * x[*,j]

	return, s

end
