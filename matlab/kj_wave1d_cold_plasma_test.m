%% Test 2: Error vs 1/h for cold plasma eps

function [Amat] = kj_wave1d_cold_plasma_test(testCase)

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
    
    [E,err_l2,A,b] = kj_wave1d(f,xMin,xMax,n(ii),lBC,rBC,ky,kz,'',eps,S,EA);

    Amat{ii}.A = A;
    Amat{ii}.b = b;
    
    all_err(ii) = err_l2;
    
    ii
    
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


