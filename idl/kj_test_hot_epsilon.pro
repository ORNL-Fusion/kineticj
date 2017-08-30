pro kj_test_hot_epsilon

    @constants

    ; This input deck is an output of AORSA2D
    ; for the purposes of verifying the sigmas.


 f=    42000000.000000000
  amu=    1.00000000
   Z=    1.00000000
    B=    4.01730919
     density=    5.0953423543494253E+017
      harm=            3
       kper=    1562.9423828125000
        kpar=    13.642897605895996
         T_eV=    245.05668640136719

    w = 2*!pi*f

    epsilon_bram = kj_hot_epsilon( f, amu, Z, B, density, harm, kPar, kPer, T_eV, $
        epsilon_cold = epsilon_cold, epsilon_swan_WD = epsilon_swan, epsilon_swan_ND = epsilon_swan_ND, $
        kx = 0, nuOmg = 0 )

    sigma_bram =  ( epsilon_bram - identity(3) ) * w * _e0 / _ii
    sigma_swan =  ( epsilon_swan_ND - identity(3) ) * w * _e0 / _ii
    sigma_cold =  ( epsilon_cold - identity(3) ) * w * _e0 / _ii

    print, 'Bram: '
    print, sigma_bram

    print, 'Swan: '
    print, sigma_swan


end
