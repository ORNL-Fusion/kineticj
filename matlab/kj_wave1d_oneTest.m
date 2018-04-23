function [] = kj_wave1d_oneTest()

n1 = 4;
n2 = 12;

n=fix(2.^linspace(n1,n2,n2-n1+1));

N=numel(n);

all_err = zeros(N,1);

for ii=1:N
    
    % Vacuum test
    % [E,err_l2] = kj_cold1d(n(ii));
    
    % Cold plasma test
    f = 13e6;xMin = 0;xMax = 2*pi/3;ky=0;kz=0;B=0.5;
    [E,err_l2] = kj_cold1d(f,xMin,xMax,n(ii),ky,kz,B);
    
    all_err(ii) = err_l2;
    
end

ideal_err = (1./n).^2;
fac = all_err(end)/ideal_err(end);

figure()
loglog(n,all_err)
hold on
loglog(n,ideal_err*fac)


end

function [Ex,Ey,Ez] = analyticSolution(x)

% ExpVar = exp(5*1i*pi*x);

ExpVar = exp(1i*15*x);

Ex = 1.00*ExpVar;
Ey = 1.00*ExpVar;
Ez = 0.01*ExpVar;

end

function [jA_x,jA_y,jA_z] = jA(x)

jA_x = x*0;
jA_y = x*0;
jA_z = x*0;

% loc = (x(end)-x(1))/2+x(1);
% sig = (x(end)-x(1))/10;
% 
% jA_x = 1*exp(-(x-loc).^2/sig^2)*(1+1i);
% jA_y = 1*exp(-(x-loc).^2/sig^2)*(1+1i);
% jA_z = 1*exp(-(x-loc).^2/sig^2)*(1+1i);

end

function [Sx,Sy,Sz] = source(x)

Sx = 0;
Sy = 0;
Sz = 0;

% Vaccum
% ky=0.1;kz=0.05;

% Sx = -1.64038467985 * exp(5*1i*pi*x);
% Sy = +245.097529329 * exp(5*1i*pi*x);
% Sz = +1.67636059316 * exp(5*1i*pi*x);

% Cold plasma
% ky=0;kz=0;

        Sx = complex(+21.0723,+72.0662) * exp(1i*15*x);
        Sy = complex(+246.072,-72.0662) * exp(1i*15*x);
        Sz = complex(+710.668,0) * exp(1i*15*x);

end