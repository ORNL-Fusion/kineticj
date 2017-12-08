pro kj_single_k_j

@dlg_constants

B = 1.0
amu = [_me_amu, _mi/_amu]
f = 42e6
Z = [-1.0,1.0]
density = [1e19,1e19]
T_eV = [2e3,2e3]
harmonicNumber = 1

lambdaPer = 0.05
lambdaPar = 0.57

kPar = 2*!pi/lambdaPar 
kPer = 2*!pi/lambdaPer 

kj_hot_epsilon, f, amu, Z, B, density, harmonicNumber, kPar, kPer, T_eV

stop

end
