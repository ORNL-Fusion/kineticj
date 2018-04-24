function tests = kj_wave1d_test

tests = functiontests(localfunctions);

end

function [] = kj_wave1d_oneTest(testCase)

global f

n1 = 4;
n2 = 12;

n=fix(2.^linspace(n1,n2,n2-n1+1));

N=numel(n);

all_err = zeros(N,1);

for ii=1:N
    
    % Vacuum test
    
    f = 13e6;
    xMin = -1;
    xMax = +1;
    ky=0.0;
    kz=0.0;
    
    S = @source1;
    EA = @analyticSolution1;
    eps = @eps1;
    [ExA,EyA,EzA] = EA(xMin);
    lBC = {'periodic',[ExA,EyA,EzA]};
    rBC = {'periodic',[0,0,0]}; % not used
    
    [E,err_l2] = kj_wave1d(f,xMin,xMax,n(ii),lBC,rBC, ky,kz,'',eps,S,EA);
    
    %     % Cold plasma test
    %
    %     f = 13e6; xMin = 0; xMax = 2*pi/3; ky=0; kz=0;
    
    %
    %     S = @source2 EA = @analyticSolution2 eps = @eps2
    %
    %     [E,err_l2] = kj_wave1d(f,xMin,xMax,n(ii),ky,kz,B);
    
    all_err(ii) = err_l2;
    
end

ideal_err = (1./n).^2;
fac = all_err(end)/ideal_err(end);

doPlots = 0;

if doPlots
    figure()
    loglog(n,all_err)
    hold on
    loglog(n,ideal_err*fac)
end

% Extract the error vs 1/h slope and compare to expected

p = polyfit(log10(n.'),log10(all_err),1);

actSlope = p(1);
expSlope = -2.0;

verifyEqual(testCase,actSlope,expSlope,'RelTol',1e-1)

end

function [Ex,Ey,Ez] = analyticSolution1(x)

ExpVar = exp(5*1i*pi*x);

Ex = 1.00*ExpVar;
Ey = 1.00*ExpVar;
Ez = 0.01*ExpVar;

end

function [jA_x,jA_y,jA_z] = jA(x)

jA_x = x*0;
jA_y = x*0;
jA_z = x*0;

end

function [Sx,Sy,Sz] = source1(x)

% Vaccum ky=0.1;kz=0.05;

expVar = exp(5*1i*pi*x);

Sx = -0.0742344 * expVar;
Sy = +246.666   * expVar;
Sz = +2.46666   * expVar;

end

function [Sx,Sy,Sz] = source(x)

Sx = 0;
Sy = 0;
Sz = 0;

% Vaccum ky=0.1;kz=0.05;

% Sx = -1.64038467985 * exp(5*1i*pi*x); Sy = +245.097529329 *
% exp(5*1i*pi*x); Sz = +1.67636059316 * exp(5*1i*pi*x);

% Cold plasma ky=0;kz=0;

Sx = complex(+21.0723,+72.0662) * exp(1i*15*x);
Sy = complex(+246.072,-72.0662) * exp(1i*15*x);
Sz = complex(+710.668,0) * exp(1i*15*x);

end


function [eps] = eps1(x)

% Vacuum dielectric

eps = eye(3);

end

function [eps] = eps2(x)

global f

% Cold plasma dielectric for a two species plasma

amu=[me_amu,2];
Z=[-1,1];
dens=[1,1]*2e18;
nu_omg=0;
B=0.5;

eps = zeros(3,3);
sig = zeros(3,3);

for s=1:numel(amu)
    
    [this_eps,this_sig] = kj_epsilon_cold(f, amu(s), Z(s), B, dens(s), nu_omg);
    
    eps = eps + this_eps;
    sig = sig + this_sig;
    
end

end


