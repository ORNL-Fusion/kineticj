function [eps,sigma] = kj_epsilon_cold (f, amu, Z, B, n, nu_omg)

phys = dlg_constants();

c = phys.('c');
u0 = phys.('u0');
eps0 = phys.('eps0');
e = phys.('e');
amu0 = phys.('amu');

validateattributes(B,{'numeric'},{'nonnegative'});

w0 = 2 * pi * f;
w = w0 * complex( 1, nu_omg);
m = amu * amu0;
q = Z * e;

wp = sqrt( n * q^2 / (m * eps0) );

% Swanson pg 23 - 24
% ------------------

eps_swan = q/abs(q);
wc_swan = abs(q) * B / m;

K1 = 1 - wp^2 / (w^2 - wc_swan^2);
K2 = (1/1i) * eps_swan * wc_swan * wp^2 / (w * (w^2 - wc_swan^2) );
K3 = 1 - wp^2 / w^2;

epsilon_swan = complex(zeros(3,3));

epsilon_swan(1,1) = +K1;
epsilon_swan(1,2) = +K2;
epsilon_swan(1,3) = 0;

epsilon_swan(2,1) = -K2;
epsilon_swan(2,2) = +K1;
epsilon_swan(2,3) = 0;

epsilon_swan(3,1) = 0;
epsilon_swan(3,2) = 0;
epsilon_swan(3,3) = +K3;


sigma_swan =  ( epsilon_swan - eye(3) ) * w0 * eps0 / 1i;

sigma = sigma_swan;
eps = epsilon_swan;

end
