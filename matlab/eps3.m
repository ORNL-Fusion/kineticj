function [eps] = eps3(x,B)

global f

% Cold plasma dielectric for a two species plasma

phys = dlg_constants();
me_amu = phys.('me_amu');

amu=[me_amu,2];
Z=[-1,1];
dens=[1,1]*4e18;
nu_omg=0;

eps = zeros(3,3);
sig = zeros(3,3);

for s=1:numel(amu)
    
    [this_eps,this_sig] = kj_epsilon_cold(f, amu(s), Z(s), B, dens(s), nu_omg);
    
    eps = eps + this_eps;
    sig = sig + this_sig;
    
end

end
