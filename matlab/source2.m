function [Sx,Sy,Sz] = source2(x)

% Cold plasma
% Z = 1; Z = {-1, 1}; amu = {me/amu0, 2}; n = 4 10^19; B = 0.8;
% kx = 40; ky = 0; kz = 0;

ExpVar = exp(1i*40*x);

Sx = complex(+499.56,+1058.6) * ExpVar;
Sy = complex(+2099.56,-1058.6) * ExpVar;
Sz = complex(+14184.4,0) * ExpVar;

end
