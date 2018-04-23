function [E,err_L2] = kj_wave1d(f,xMin,xMax,N,ky,kz,jA,eps,S,EA)
% KJ_WAVE1D  1D cold plasma wave solver.
%   [E] = KJ_WAVE1D(f,xMin,xMax,nPts) takes a frequency in Hz (f), domain
%   extents (xMin,xMax), number of points (nPts) and returns a 1D vector of
%   length nPtsx3 with the Ex,Ey,Ez components of the solution. This
%   assumes ky=kz=0, vacuum, and zero source current.
%
%   [E] = KJ_WAVE1D(xMin,xMax,nPts,ky,kz) specifies the wavenumber in the y
%   and z directions.
%
%   [E] = KJ_WAVE1D(xMin,xMax,nPts,ky,kz,jA,eps) specifies the driving
%   antenna current (jA) and the background magnetic field strength. Here
%   [jA_x,jA_y,jA_z]=jA(x) and
%   [[exx,exy,exz],[eyx,eyy,eyz],[ezx,ezy,ezz]]=eps(x) are function handles
%   which accept the location and returns the 3 components of jA and the
%   dielectric tensor.
%
%   [E,err] = KJ_WAVE1D(xMin,xMax,nPts,ky,kz,jA,S,EA) returns the L2 error
%   (err) between the solution and an analytic solution, used for testing
%   with the Method of Manufactured Solutions with both [S_x,S_y,S_z]=S(x)
%   and [ExA,EyA,EzA]=EA(x) being function handles which accept position x
%   and return the 3 components of the source and analytic solution
%   respectively.


if ~exist('xMin','var')     || isempty(xMin)
    xMin = -1;
end
if ~exist('xMax','var')     || isempty(xMax)
    xMax = +1;
end
if ~exist('f','var')        || isempty(f)
    f = 13e6;
end
if ~exist('nPts','var')     || isempty(N)
    N = 512;
end
if ~exist('ky','var')       || isempty(ky)
    ky = 0.0;
end
if ~exist('kz','var')       || isempty(kz)
    kz = 0.0;
end
if ~exist('eps','var')      || isempty(eps)
    [eps] = @(x) eye(3);
end
if ~exist('jA','var')       || isempty(jA)
    [jA_x,jA_y,jA_z] = @(x) 0;
end
if ~exist('S','var')        || isempty(S)
    [S_x,S_y,S_z] = @(x) 0;
end
if ~exist('S','var')        || isempty(S)
    [S_x,S_y,S_z] = @(x) 0;
end

% Boundary conditions

lbc = 'periodic';
rbc = 'periodic';
%lbc = 'dirichlet'; rbc = 'dirichlet';

ExL = 0;
EyL = 0;
EzL = 0;

ExR = 0;
EyR = 0;
EzR = 0;

periodic = 0;
ldirichlet = 0;
rdirichlet = 0;

if strcmp(lbc,'periodic') || strcmp(rbc,'periodic')
    n = N-1;
    periodic = 1;
else
    n = N;
end

if strcmp(lbc,'dirichlet')
    ldirichlet = 1;
end
if strcmp(rbc,'dirichlet')
    rdirichlet = 1;
end

% Physical constants

phys = dlg_constants();

c = phys('c');
u0 = phys('u0');
eps0 = phys('eps0');
me_amu = phys('me_amu');
amu0 = phys('amu');


% Parameters

h = (xMax-xMin) / (N-1);

x = linspace(xMin,xMax-h,n);

w = 2*pi*f;
k0 = w/c;

N = n*3;

A = sparse(N,N);
b = complex(zeros(N,1));


% Dielectric tensor

epsc = eye(3); % Vacuum

amu=[me_amu,2];
Z=[-1,1];
dens=[1,1]*2e18;
nu_omg=0;

epsc = zeros(3,3);
sigc = zeros(3,3);

for s=1:numel(amu)
    
    [this_epsc,this_sigc] = kj_epsilon_cold(f, amu(s), Z(s), B, dens(s), nu_omg); % Cold plasma
    
    epsc = epsc + this_epsc;
    sigc = sigc + this_sigc;
    
end


exx = epsc(1,1);
exy = epsc(1,2);
exz = epsc(1,3);

eyx = epsc(2,1);
eyy = epsc(2,2);
eyz = epsc(2,3);

ezx = epsc(3,1);
ezy = epsc(3,2);
ezz = epsc(3,3);


% Current source

[jA_x,jA_y,jA_z] = jA(x);



% Get analytic solution

[ExA,EyA,EzA] = analyticSolution(x);


% Fill matrix with grouping by component

for jj=1:N
    
    if jj == 1
        
        if periodic
            
            A(jj +0*n,jj +0*n) = 1;
            A(jj +1*n,jj +1*n) = 1;
            A(jj +2*n,jj +2*n) = 1;
            
            b(jj +0*n) = ExA(jj);
            b(jj +1*n) = EyA(jj);
            b(jj +2*n) = EzA(jj);
            
        elseif ldirichlet
            
            A(jj +0*n,jj +0*n) = 1;
            A(jj +1*n,jj +1*n) = 1;
            A(jj +2*n,jj +2*n) = 1;
            
            b(jj +0*n) = ExL;
            b(jj +1*n) = EyL;
            b(jj +2*n) = EzL;
            
        end
        
    elseif jj == N
        
        
        if periodic
            
        elseif rdirichlet
            
            A(jj +0*n,jj +0*n) = 1;
            A(jj +1*n,jj +1*n) = 1;
            A(jj +2*n,jj +2*n) = 1;
            
            b(jj +0*n) = ExR;
            b(jj +1*n) = EyR;
            b(jj +2*n) = EzR;
            
        end
        
    else
        
        jm = jj-1;
        jp = jj+1;
        
        if periodic && jp == n+1
            jp = 1;
        end
        
        % Ex
        
        A(jj +0*n,jm +0*n) = 0;
        A(jj +0*n,jj +0*n) = -exx * k0^2 + ky^2 + kz^2;
        A(jj +0*n,jp +0*n) = 0;
        
        A(jj +0*n,jm +1*n) = -1i * ky / (2*h);
        A(jj +0*n,jj +1*n) = -exy * k0^2;
        A(jj +0*n,jp +1*n) = +1i * ky / (2*h);
        
        A(jj +0*n,jm +2*n) = -1i * kz / (2*h);
        A(jj +0*n,jj +2*n) = -exz * k0^2;
        A(jj +0*n,jp +2*n) = +1i * kz / (2*h);
        
        % Ey
        
        A(jj +1*n,jm +0*n) = -1i * ky / (2*h);
        A(jj +1*n,jj +0*n) = -eyx * k0^2;
        A(jj +1*n,jp +0*n) = +1i * ky / (2*h);
        
        A(jj +1*n,jm +1*n) = -1/h^2;
        A(jj +1*n,jj +1*n) = +2/h^2 - eyy * k0^2 + kz^2;
        A(jj +1*n,jp +1*n) = -1/h^2;
        
        A(jj +1*n,jm +2*n) = 0;
        A(jj +1*n,jj +2*n) = -eyz * k0^2 - ky * kz;
        A(jj +1*n,jp +2*n) = 0;
        
        % Ez
        
        A(jj +2*n,jm +0*n) = -1i * kz / (2*h);
        A(jj +2*n,jj +0*n) = -ezx * k0^2;
        A(jj +2*n,jp +0*n) = +1i * kz / (2*h);
        
        A(jj +2*n,jm +1*n) = 0;
        A(jj +2*n,jj +1*n) = -ezy * k0^2 - ky * kz;
        A(jj +2*n,jp +1*n) = 0;
        
        A(jj +2*n,jm +2*n) = -1/h^2;
        A(jj +2*n,jj +2*n) = +2/h^2 - ezz * k0^2 + ky^2;
        A(jj +2*n,jp +2*n) = -1/h^2;
        
        % RHS
        
        [Sx,Sy,Sz] = source(x(jj));
        
        b(jj +0*n) = 1i*w*u0*jA_x(jj) + Sx;
        b(jj +1*n) = 1i*w*u0*jA_y(jj) + Sy;
        b(jj +2*n) = 1i*w*u0*jA_z(jj) + Sz;
        
    end
    
end

% Solve

E = A\b;

Ex = E(0*n+1:1*n);
Ey = E(1*n+1:2*n);
Ez = E(2*n+1:3*n);

if
    
    EA = [ExA,EyA,EzA].';
    
    err_L2 = abs(sqrt(mean((E-EA).^2)));
    
    doPlots = 0;
    if doPlots
        
        kj_plot_cmplx_3vec(E,EA)
        
        figure()
        ax1 = subplot(3,1,1);
        plot(ax1,x,real(Ex))
        hold on
        plot(ax1,x,imag(Ey))
        
        ax2 = subplot(3,1,2);
        plot(ax2,x,real(Ey))
        hold on
        plot(ax2,x,imag(Ex))
        
        ax3 = subplot(3,1,3);
        plot(ax3,x,real(Ez))
        hold on
        plot(ax3,x,imag(Ez))
        
    end
    
    if periodic
        Ex = [Ex.',Ex(1)];
        Ey = [Ey.',Ey(1)];
        Ez = [Ez.',Ez(1)];
    end
    
end


    function [Ex,Ey,Ez] = analyticSolution(x)
        
        Ex = 0;
        Ey = 0;
        Ez = 0;
        
    end

    function [jA_x,jA_y,jA_z] = jA(x)
        
        jA_x = x*0;
        jA_y = x*0;
        jA_z = x*0;
        
        
    end

    function [Sx,Sy,Sz] = source(x)
        
        Sx = 0;
        Sy = 0;
        Sz = 0;
        
    end
