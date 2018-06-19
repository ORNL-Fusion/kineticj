%% Test 3: Compare kx actual with dispersion

function [Amat] = kj_wave1d_cold_plasma_dispersion_test(testCase)

disp('Test 3: Compare kx actual with dispersion');

global f

phys = dlg_constants();
me_amu = phys.('me_amu');

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
    
    [E,err_L2,A,bvec] = kj_wave1d(f,xMin,xMax,nPts,lBC,rBC,ky,kz,jA,eps);
    Amat{b}.A = A;
    Amat{b}.b = bvec;

    Ex = E(1:nPts);
    power(:,b) = abs(fft(Ex)).^2;
    power(:,b) = power(:,b)/sum(power(:,b));
    
end

% figure()
% contourf(B,kxaxis,power);
% ylim([0,50])

% ------------------------------------
% check whether the file can be opened
% ------------------------------------
fid = fopen('tests/kj_wave1d.mat');
isok = (fid > 0);
if (isok),
  fclose(fid);
  load('tests/kj_wave1d.mat','actPower');
  verifyEqual(testCase,actPower,power,'RelTol',1e-3)
end;

end


