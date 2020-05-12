%% Example for calling kj_wave1d()

function [stat] = kj_wave1d_example()

global f

% Number of grid points
n = 1024;

% Frequency [Hz]
f = 2.52e6;

% Domain Range [m]
xMin = -4;
xMax = +4;

% Ignorable direction k values [1/m]
ky=0.0;
kz=20.0i;

% Set the function that returns the current source (see below)
S = @source2;

% Set the function that returns the cold plasma dielectric (see below)
eps = @eps2;

% Set the boundary conditions to 'periodic' or 'dirichlet'
lBC = {'periodic',[0,0,0]};
rBC = {'periodic',[0,0,0]}; % not used

% Call kj_wave1d()
[E,err,x] = kj_wave1d(f,xMin,xMax,n,lBC,rBC,ky,kz,'',eps,S);

% Plot solution

kj_plot_cmplx_3vec(x,E)

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

amplitude = 2;

xMax = +4.;
xMin = -4.;
damping_width = (xMax-xMin)/5;
lSide = amplitude * exp(-(x-xMax).^2./(damping_width).^2);
rSide = amplitude * exp(-(x-xMin).^2./(damping_width).^2);

result = lSide + rSide;

end