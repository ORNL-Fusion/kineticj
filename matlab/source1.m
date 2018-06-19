function [Sx,Sy,Sz] = source1(x)

% Vaccum ky=0.1;kz=0.05;

ExpVar = exp(5*1i*pi*x);

Sx = -0.0742344 * ExpVar;
Sy = +246.666   * ExpVar;
Sz = +2.46666   * ExpVar;

end
