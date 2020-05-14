%% Example for calling kj_wave1d()

function [stat] = kj_wave1d_example()

global f

% Number of grid points
n = 1024;

% Frequency [Hz]
f = 3e6;

% Domain Range [m]
xMin = -4;
xMax = +4;

% Ignorable direction k values [1/m]
ky = 0.0;
kz = @kz2;

% Set the function that returns the current source (see below)
S = @source2;

% Set the function that returns the cold plasma dielectric (see below)
eps = @eps2;

% Set the spatial damping profile
damping = @damping2;

% Set the boundary conditions to 'periodic' or 'dirichlet'
lBC = {'periodic',[0,0,0]};
rBC = {'periodic',[0,0,0]}; % not used

% Call kj_wave1d()
[E,err,x,eps] = kj_wave1d(f,xMin,xMax,n,lBC,rBC,ky,kz,'',eps,S,'',damping);

M = numel(E);

n = M/3;

ex = E(0*n+1:1*n);
ey = E(1*n+1:2*n);
ez = E(2*n+1:3*n);


figure()
subplot(4,1,1)
plot(x,real(ex));
hold on 
plot(x,imag(ex),'color','r');
plot(x,abs(ex),'color','black');
subplot(4,1,2)
plot(x,real(ey));
hold on 
plot(x,imag(ey),'color','r');
plot(x,abs(ey),'color','black');
subplot(4,1,3)
plot(x,real(ez));
hold on 
plot(x,imag(ez),'color','r');
plot(x,abs(ez),'color','black');
subplot(4,1,4)
plot(x,damping(x));


% Plot solution

kj_plot_cmplx_3vec(x,E)

% Look at the real/imag parts of the dielectric to understand the damping

figure()
subplot(3,3,1)
plot(x,real(eps.exx));
hold on 
plot(x,imag(eps.exx),'color','r');
subplot(3,3,2)
plot(x,real(eps.exy));
hold on 
plot(x,imag(eps.exy),'color','r');
subplot(3,3,3)
plot(x,real(eps.exz));
hold on 
plot(x,imag(eps.exz),'color','r');

subplot(3,3,4)
plot(x,real(eps.eyx));
hold on 
plot(x,imag(eps.eyx),'color','r');
subplot(3,3,5)
plot(x,real(eps.eyy));
hold on 
plot(x,imag(eps.eyy),'color','r');
subplot(3,3,6)
plot(x,real(eps.eyz));
hold on 
plot(x,imag(eps.eyz),'color','r');

subplot(3,3,7)
plot(x,real(eps.ezx));
hold on 
plot(x,imag(eps.ezx),'color','r');
subplot(3,3,8)
plot(x,real(eps.ezy));
hold on 
plot(x,imag(eps.ezy),'color','r');
subplot(3,3,9)
plot(x,real(eps.ezz));
hold on 
plot(x,imag(eps.eyz),'color','r');

figure
plot(x,real(eps.kz))
hold on
plot(x,imag(eps.kz),'color','r')

figure
plot(x,real(eps.w))
hold on
plot(x,imag(eps.w),'color','r')

end

%% Define the kz profile

function [kz1] = kz2(x)

kz0 = 20i;

xMax = +4.;
xMin = -4.;
damping_width = (xMax-xMin)/10;
lOffSet = xMin + damping_width;
rOffSet = xMax - damping_width; 
lSide = exp(-(x-lOffSet).^2./(damping_width/2).^2);
rSide = exp(-(x-rOffSet).^2./(damping_width/2).^2);

lii = x <= lOffSet;
rii = x >= rOffSet;

lSide(lii) = 1;
rSide(rii) = 1;

kz1 = (1-(lSide + rSide)) .* kz0;
kz1 = kz1.*0 + kz0; % just turn this off for now

end

%% Setup the cold plasma dielectric

function [eps] = eps2(x)

global f

% Cold plasma dielectric for a two species (e,D) plasma

phys = dlg_constants();
me_amu = phys('me_amu');

amu=[me_amu,2];
Z=[-1,1];
dens=[1,1]*1e17;
nu_omg=damping2(x);
% nu_omg=0;
B=0.1;

eps = zeros(3,3);
sig = zeros(3,3);

for s=1:numel(amu)
    
    % Get the eps and sig
    [this_eps,this_sig] = kj_epsilon_cold(f, amu(s), Z(s), B, dens(s), nu_omg);
    
    eps = eps + this_eps;
    sig = sig + this_sig;
    
end

end

%% Setup the source (RHS)

function [Sx,Sy,Sz] = source2(x)

% Returns the x, y and z components of the volume source current as a function of position. 

% Cold plasma
% Z = 1; Z = {-1, 1}; amu = {me/amu0, 2}; n = 4 10^19; B = 0.8;
% kx = 40; ky = 0; kz = 0;

offset = 0;
width = 0.06;

ExpVar = exp(-(x-offset).^2./(width).^2);

Sx = 0.0;
Sy = 1 * ExpVar;
Sz = 0.0;

end

%% Setup the damping profile (nu_omg)

function [result] = damping2(x)

amplitude = 4e3;

xMax = +4.;
xMin = -4.;
damping_width = (xMax-xMin)/10;
lSide = amplitude * exp(-(x-xMax).^2./(damping_width).^2);
rSide = amplitude * exp(-(x-xMin).^2./(damping_width).^2);

result = lSide + rSide;

end