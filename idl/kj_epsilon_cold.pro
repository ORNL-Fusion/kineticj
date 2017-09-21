function kj_epsilon_cold, f, amu, atomicZ, B, density, nu_omg, sigma = sigma

    @constants

    if B lt 0 then stop ; B is magnitude of B
   
    w = 2 * !Pi * f
    m = amu * _amu * complex( 1, nu_omg) 
    q = atomicZ * _e

    wp = sqrt( density * q^2 / (m * _e0) ) 

    ; Stix pg 7
    ; ---------

    wc = q * B / m

    ;P = 1-wp^2/( w^2 )
    ;R = 1-wp^2 /( w*(w+wc) )
    ;L = 1-wp^2 /( w*(w-wc) )
    ;S = 0.5*(R+L)
    ;D = 0.5*(R-L)

    ;epsilon_stix = dComplexArr(3,3)
    ;
    ;epsilon_stix[0,0] = +S
    ;epsilon_stix[0,1] = -_II * D
    ;epsilon_stix[0,2] = 0
    ;
    ;epsilon_stix[1,0] = +_II*D 
    ;epsilon_stix[1,1] = +S
    ;epsilon_stix[1,2] = 0
    ;
    ;epsilon_stix[2,0] = 0 
    ;epsilon_stix[2,1] = 0
    ;epsilon_stix[2,2] = +P

    ;; Brambilla pg 165
    ;; ---------------- 

    ;S_bram = 1-wp^2/(w^2-wc^2)
    ;D_bram = wc/w * wp^2/(w^2-wc^2)
    ;P_bram = 1-wp^2/w^2

    ;epsilon_bram = dComplexArr(3,3)
    ;
    ;epsilon_bram[0,0] = +S_bram 
    ;epsilon_bram[0,1] = -_ii*D_bram
    ;epsilon_bram[0,2] = 0 
    ;
    ;epsilon_bram[1,0] = +_ii*D_bram 
    ;epsilon_bram[1,1] = +S_bram 
    ;epsilon_bram[1,2] = 0 
    ;
    ;epsilon_bram[2,0] = 0 
    ;epsilon_bram[2,1] = 0 
    ;epsilon_bram[2,2] = +P_bram 


    ; Swanson pg 23 - 24
    ; ------------------

    eps_swan = q/abs(q)
    wc_swan = abs(q) * B / m

    K1 = 1d0 - wp^2 / (w^2 - wc_swan^2)
    K2 = (1/_ii) * eps_swan * wc_swan * wp^2 / (w * (w^2 - wc_swan^2) )
    K3 = 1d0 - wp^2 / w^2

    epsilon_swan = dComplexArr(3,3)
    
    epsilon_swan[0,0] = +K1
    epsilon_swan[0,1] = +K2 
    epsilon_swan[0,2] = 0

    epsilon_swan[1,0] = -K2
    epsilon_swan[1,1] = +K1 
    epsilon_swan[1,2] = 0

    epsilon_swan[2,0] = 0 
    epsilon_swan[2,1] = 0 
    epsilon_swan[2,2] = +K3 

    if arg_present(sigma) then begin

        ;sigma_stix =  ( epsilon_stix - identity(3) ) * w * _e0 / _ii
        ;sigma_bram =  ( epsilon_bram - identity(3) ) * w * _e0 / _ii
        sigma_swan =  ( epsilon_swan - identity(3) ) * w * _e0 / _ii

        sigma = sigma_swan

    endif

    return, epsilon_swan

stop
end


