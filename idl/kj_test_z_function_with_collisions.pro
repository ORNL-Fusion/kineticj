pro kj_test_z_function_with_collisions

    @dlg_constants

    f = 42e6
    amu = 1.0;_me_amu 
    atomicZ = +1
    B = 0.6
    density = 1e19
    harmonicNumber = 3
    kPar = 100
    kPer = 10
    T_eV = 1e-3

    N = 35
    nu_omg_min = 0
    nu_omg_max =3 
    nu_omg = findgen(N)/(N-1)*(nu_omg_max-nu_omg_min)+nu_omg_min

    eps_zz_c = complexarr(N) 
    eps_zz_h = complexarr(N) 
    eps_xx_c = complexarr(N) 
    eps_xx_h = complexarr(N) 

    sig_zz_c = complexarr(N) 
    sig_zz_h = complexarr(N) 
    sig_xx_c = complexarr(N) 
    sig_xx_h = complexarr(N) 

    for n=0,N-1 do begin

        eps_c = kj_epsilon_cold( f, amu, atomicZ, B, density, nu_omg[n], sigma = sig_c)
        eps_h = kj_epsilon_hot ( f, amu, atomicZ, B, density, harmonicNumber, $
                kPar, kPer, T_eV, nuOmg = nu_omg[n], sigma = sig_h )

        eps_zz_c[n] = eps_c[2,2]
        eps_zz_h[n] = eps_h[2,2]

        sig_zz_c[n] = sig_c[2,2]
        sig_zz_h[n] = sig_h[2,2]

        eps_xx_c[n] = eps_c[0,0]
        eps_xx_h[n] = eps_h[0,0]

        sig_xx_c[n] = sig_c[0,0]
        sig_xx_h[n] = sig_h[0,0]

    endfor

    p=plot(nu_omg,eps_zz_c,layout=[2,2,1])
    p=plot(nu_omg,imaginary(eps_zz_c),color='r',/over)

    p=plot(nu_omg,eps_zz_h,/over,thick=2)
    p=plot(nu_omg,imaginary(eps_zz_h),color='r',/over,thick=2)

    p=plot(nu_omg,sig_zz_c,layout=[2,2,2], /current)
    p=plot(nu_omg,imaginary(sig_zz_c),color='r',/over)

    p=plot(nu_omg,sig_zz_h,/over,thick=2)
    p=plot(nu_omg,imaginary(sig_zz_h),color='r',/over,thick=2)

    p=plot(nu_omg,eps_xx_c,layout=[2,2,3], /current)
    p=plot(nu_omg,imaginary(eps_xx_c),color='r',/over)

    p=plot(nu_omg,eps_xx_h,/over,thick=2)
    p=plot(nu_omg,imaginary(eps_xx_h),color='r',/over,thick=2)

    p=plot(nu_omg,sig_xx_c,layout=[2,2,4],/current)
    p=plot(nu_omg,imaginary(sig_xx_c),color='r',/over)

    p=plot(nu_omg,sig_xx_h,/over,thick=2)
    p=plot(nu_omg,imaginary(sig_xx_h),color='r',/over,thick=2)

    stop

end
