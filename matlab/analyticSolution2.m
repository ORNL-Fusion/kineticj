function [Ex,Ey,Ez] = analyticSolution2(x)

ExpVar = exp(40*1i*x);

Ex = 1.00*ExpVar;
Ey = 1.00*ExpVar;
Ez = 0.01*ExpVar;

end
