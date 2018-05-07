function tests = kj_wave1d_test

tests = functiontests(localfunctions);

end

%% Test 1: Error vs 1/h for vacuum eps. 

function [] = kj_wave1d_vacuum_test(testCase)

disp('Test 1: Error vs 1/h for vacuum eps');

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
    rBC = {'periodic',[0,0,0]}; % not used since only the left is used.
    
    [E,err_l2] = kj_wave1d(f,xMin,xMax,n(ii),lBC,rBC, ky,kz,'',eps,S,EA);
    
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


%% Test 2: Error vs 1/h for cold plasma eps

function [] = kj_wave1d_cold_plasma_test(testCase)

disp('Test 2: Error vs 1/h for cold plasma eps');

global f

n1 = 6;
n2 = 12;

n=fix(2.^linspace(n1,n2,n2-n1+1));

N=numel(n);

all_err = zeros(N,1);

for ii=1:N
        
    f = 13e6;
    xMin = 0;
    xMax = pi/4;
    ky=0.0;
    kz=0.0;
    
    S = @source2;
    EA = @analyticSolution2;
    eps = @eps2;
    [ExA,EyA,EzA] = EA(xMin);
    lBC = {'periodic',[ExA,EyA,EzA]};
    rBC = {'periodic',[0,0,0]}; % not used
    
    [E,err_l2] = kj_wave1d(f,xMin,xMax,n(ii),lBC,rBC,ky,kz,'',eps,S,EA);
    
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

verifyEqual(testCase,actSlope,expSlope,'RelTol',1.5e-1)

end


%% Test 3: Compare kx actual with dispersion

function [] = kj_wave1d_cold_plasma_dispersion_test(testCase)

disp('Test 3: Compare kx actual with dispersion');

global f

phys = dlg_constants();
me_amu = phys('me_amu');

f = 13e6;
xMin = -1;
xMax = +1;
nPts = 512;
ky=0;
kz=0;
amu=[me_amu,2];
Z=[-1,1];
dens=[1,1]*4e18;
nu_omg=0;

lBC = {'dirichlet',[0,0,0]};
rBC = {'dirichlet',[0,0,0]};

jA = @jA3;

B = linspace(0.1,1.0,20);

power = zeros(nPts,numel(B));

dx = (xMax-xMin)/(nPts-1);
kxaxis = (linspace(1,nPts,nPts)-1) / ( dx * nPts ) * 2*pi;

for b=1:numel(B)
    
    eps = @(x) eps3(x,B(b));
    
    E = kj_wave1d(f,xMin,xMax,nPts,lBC,rBC,ky,kz,jA,eps);
    Ex = E(1:nPts);
    power(:,b) = abs(fft(Ex)).^2;
    power(:,b) = power(:,b)/sum(power(:,b));
    
end

% figure()
% contourf(B,kxaxis,power);
% ylim([0,50])

load('tests/kj_wave1d.mat','actPower');

verifyEqual(testCase,actPower,power,'RelTol',1e-3)

end


function [Ex,Ey,Ez] = analyticSolution1(x)

ExpVar = exp(5*1i*pi*x);

Ex = 1.00*ExpVar;
Ey = 1.00*ExpVar;
Ez = 0.01*ExpVar;

end

function [Ex,Ey,Ez] = analyticSolution2(x)

ExpVar = exp(40*1i*x);

Ex = 1.00*ExpVar;
Ey = 1.00*ExpVar;
Ez = 0.01*ExpVar;

end

function [jA_x,jA_y,jA_z] = jA(x)

jA_x = x*0;
jA_y = x*0;
jA_z = x*0;

end

function [jA_x,jA_y,jA_z] = jA3(x)

xCenter = pi/2;
c = pi/10;
gauss = 1*exp(- (x-xCenter)^2 / (2*c^2) );

jA_x = gauss * complex(1,1);
jA_y = gauss * complex(1,1);
jA_z = gauss * complex(1,1);

end

function [Sx,Sy,Sz] = source1(x)

% Vaccum ky=0.1;kz=0.05;

ExpVar = exp(5*1i*pi*x);

Sx = -0.0742344 * ExpVar;
Sy = +246.666   * ExpVar;
Sz = +2.46666   * ExpVar;

end

function [Sx,Sy,Sz] = source2(x)

% Cold plasma
% Z = 1; Z = {-1, 1}; amu = {me/amu0, 2}; n = 4 10^19; B = 0.8;
% kx = 40; ky = 0; kz = 0;

ExpVar = exp(1i*40*x);

Sx = complex(+499.56,+1058.6) * ExpVar;
Sy = complex(+2099.56,-1058.6) * ExpVar;
Sz = complex(+14184.4,0) * ExpVar;

end


function [eps] = eps1(x)

% Vacuum dielectric

eps = eye(3);

end

function [eps] = eps2(x)

global f

% Cold plasma dielectric for a two species plasma

phys = dlg_constants();
me_amu = phys('me_amu');

amu=[me_amu,2];
Z=[-1,1];
dens=[1,1]*4e19;
nu_omg=0;
B=0.8;

eps = zeros(3,3);
sig = zeros(3,3);

for s=1:numel(amu)
    
    [this_eps,this_sig] = kj_epsilon_cold(f, amu(s), Z(s), B, dens(s), nu_omg);
    
    eps = eps + this_eps;
    sig = sig + this_sig;
    
end

end


function [eps] = eps3(x,B)

global f

% Cold plasma dielectric for a two species plasma

phys = dlg_constants();
me_amu = phys('me_amu');

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

