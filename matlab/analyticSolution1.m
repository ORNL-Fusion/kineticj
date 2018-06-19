function [Ex,Ey,Ez] = analyticSolution1(x)

ExpVar = exp(5*1i*pi*x);

Ex = 1.00*ExpVar;
Ey = 1.00*ExpVar;
Ez = 0.01*ExpVar;

end
