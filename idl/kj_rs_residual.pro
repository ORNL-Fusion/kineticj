function kj_rs_residual

    rs = rs_read_solution('./') 
    
    r  = rs('r')
    dR = r[1]-r[0]
    r_ = r[1:-1]+dR/2

    nR = n_elements(r)

    LHS = rs_lhs()
    RHS = rs_rhs()

    res = RHS - LHS

    ; Put the residual back on the full grid

	ii_R	= 2 + lIndGen(nR)*3
	ii_T	= ii_R[0:nR-2]+1
	ii_Z	= ii_T+1

    res_r  = res[ii_R]
    res_t_ = [0,res[ii_T],0]
    res_z_ = [0,res[ii_Z],0]

    r__ = [r_[0]-dR,r_,r_[-1]+dR]

    res_t_re = interpol(real_part(res_t_),r__,r,/spline)
    res_z_re = interpol(real_part(res_z_),r__,r,/spline)
    res_t_im = interpol(imaginary(res_t_),r__,r,/spline)
    res_z_im = interpol(imaginary(res_z_),r__,r,/spline)

    res_t = dcomplex(res_t_re,res_t_im)
    res_z = dcomplex(res_z_re,res_z_im)

    res1 = [res_r,res_t,res_z]

    solHash = HASH() 

    solHash = solHash + HASH('res_re',real_part(res1));
    solHash = solHash + HASH('res_im',imaginary(res1));

    solHash = solHash + HASH('res_r_re',real_part(res_r));
    solHash = solHash + HASH('res_r_im',imaginary(res_r));
    solHash = solHash + HASH('res_t_re',real_part(res_t));
    solHash = solHash + HASH('res_t_im',imaginary(res_t));
    solHash = solHash + HASH('res_z_re',real_part(res_z));
    solHash = solHash + HASH('res_z_im',imaginary(res_z));

    NCDF_PUT, 'output/kj-rs-res.nc', /NEW, VARIABLES=solHash

    print, 'norm(residual): ', norm(res1)

    return, res1 

end
