pro kj_test_hot_epsilon

    @dlg_constants

    ; This input deck is an output of AORSA2D
    ; for the purposes of verifying the sigmas.

     f=    42000000.000000000
      amu=    2.00000000
       Z=    1.00000000
        B=    1.00000000
         density=    3.6771202663201116E+019
          harm=           12
           kper=    1767.1459960937500
            kpar=    9.5483818054199219
             T_eV=    1995.3190917968750
              nuOmg=    0.00000000

    w = 2*!pi*f

    epsilon_bram = kj_hot_epsilon( f, amu, Z, B, density, harm, kPar, kPer, T_eV, $
        epsilon_cold = epsilon_cold, epsilon_swan_WD = epsilon_swan, epsilon_swan_ND = epsilon_swan_ND, $
        kx = 0, nuOmg = nuOmg )

    sigma_bram =  ( epsilon_bram - identity(3) ) * w * _e0 / _ii
    sigma_swan =  ( epsilon_swan_ND - identity(3) ) * w * _e0 / _ii
    sigma_cold =  ( epsilon_cold - identity(3) ) * w * _e0 / _ii

    print, 'Bram: '
    print, sigma_bram

    print, 'Swan: '
    print, sigma_swan


end
