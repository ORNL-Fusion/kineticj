function [E,err_L2, A,b] = kj_wave1d(f,xMin,xMax,N,lBC,rBC,ky,kz,jA,eps,S,EAnalytic,terms)
% KJ_WAVE1D  1D cold plasma wave solver.
%   [E] = KJ_WAVE1D(f,xMin,xMax,nPts,lbc,rbc) takes a frequency in Hz (f),
%   domain extents (xMin,xMax), number of points (nPts) and returns a 1D
%   vector of length nPtsx3 with the Ex,Ey,Ez components of the solution.
%   This assumes ky=kz=0, vacuum, and zero source current. The boundary
%   conditions are specified on the left (lbc) and right (rbc) using key
%   value pairs of the type and value, e.g.,
%       lBC = {'periodic',[Ex,Ey,Ez]} rBC = {'dirichlet',[Ex,Ey,Ez]}
%
%   [E] = KJ_WAVE1D(xMin,xMax,nPts,lbc,rbc,ky,kz) specifies the wavenumber
%   in the y and z directions.
%
%   [E] = KJ_WAVE1D(xMin,xMax,nPts,lbc,rbc,ky,kz,jA,eps) specifies the
%   driving antenna current (jA) and the background magnetic field
%   strength. Here [jA_x,jA_y,jA_z]=jA(x) and
%   [[exx,exy,exz],[eyx,eyy,eyz],[ezx,ezy,ezz]]=eps(x) are function handles
%   which accept the location and returns the 3 components of jA and the
%   dielectric tensor.
%
%   [E,err] = KJ_WAVE1D(xMin,xMax,nPts,lbc,rbc,ky,kz,jA,S,EAnalytic)
%   returns the L2 error (err) between the solution and an analytic
%   solution, used for testing with the Method of Manufactured Solutions
%   with both [S_x,S_y,S_z]=S(x) and [ExA,EyA,EzA]=EAnalytic(x) being
%   function handles which accept position x and return the 3 components of
%   the source and analytic solution respectively.
%
%   [E,err] = KJ_WAVE1D(xMin,xMax,nPts,lbc,rbc,ky,kz,jA,S,EAnalytic,terms)
%   allows specification of which of the 4 terms in the wave equation to
%   include (for preconditioner purposes). Specified as multipliers of
%   those terms = [1,1,1,1];

err_L2 = 0;
A = [];
b = [];
E = [];

compareWithAnalytic = 0;

if ~exist('xMin','var')     || isempty(xMin)
    xMin = -1;
end
if ~exist('xMax','var')     || isempty(xMax)
    xMax = +1;
end
if ~exist('f','var')        || isempty(f)
    f = 13e6;
end
if ~exist('N','var')        || isempty(N)
    N = 512;
end
if ~exist('ky','var')       || isempty(ky)
    ky = 0.0;
end
if ~exist('kz','var')       || isempty(kz)
    kz = 0.0;
end
if ~exist('eps','var')      || isempty(eps)
    eps = @(x) eye(3);
end
if ~exist('jA','var')       || isempty(jA)
    jA = @(x) deal(0,0,0);
end
if ~exist('S','var')        || isempty(S)
    S = @(x) deal(0,0,0);
end
if ~exist('EAnalytic','var')        || isempty(EAnalytic)
    EAnalytic = @(x) deal(0,0,0);
else
    compareWithAnalytic = 1;
end
if ~exist('lBC','var')        || isempty(lBC)
    lBC = {'dirichlet',[0,0,0]};
end
if ~exist('rBC','var')        || isempty(rBC)
    rBC = {'dirichlet',[0,0,0]};
end
if ~exist('terms','var')      || isempty(terms)
    terms = [1,1,1,1];
end


% Setup BC logic

periodic = NaN;
ldirichlet = NaN;
rdirichlet = NaN;

if strcmp(lBC{1},'periodic') || strcmp(rBC{1},'periodic')
    n = N-1;
    periodic = 1;
    ldirichlet = 0;
    rdirichlet = 0;
else
    n = N;
end

if strcmp(lBC{1},'dirichlet')
    ldirichlet = 1;
    periodic = 0;
end
if strcmp(rBC{1},'dirichlet')
    rdirichlet = 1;
    periodic = 0;
end

ExL = lBC{2}(1);
EyL = lBC{2}(2);
EzL = lBC{2}(3);

if rdirichlet
    
    ExR = rBC{2}(1);
    EyR = rBC{2}(2);
    EzR = rBC{2}(3);
    
end

assert(periodic==periodic,'Error setting BC (periodic)');
assert(ldirichlet==ldirichlet,'Error setting BC (ldirichlet)');
assert(rdirichlet==rdirichlet,'Error setting BC (rdirichlet)');

% disp('Boundary Conditions:'); lBCStr = sprintf("lBC: %s,
% [%f,%f,%f]",lBC{1},lBC{2}(1),lBC{2}(2),lBC{2}(3)); rBCStr = sprintf("rBC:
% %s, [%f,%f,%f]",rBC{1},rBC{2}(1),rBC{2}(2),rBC{2}(3)); disp(lBCStr)
% disp(rBCStr)


% Physical constants

phys = dlg_constants();

c = phys.('c');
u0 = phys.('u0');
eps0 = phys.('eps0');
me_amu = phys.('me_amu');
amu0 = phys.('amu');


% Parameters

h = (xMax-xMin) / (N-1);

x = linspace(xMin,xMax,N);

w = 2*pi*f;
k0 = w/c;

NDOF = n*3;

A = sparse(NDOF,NDOF);
b = complex(zeros(NDOF,1));

A_T1A = sparse(NDOF,NDOF); % Term T1A: -Laplacian(E)
A_T1B = sparse(NDOF,NDOF); % Term T1B: +Grad(Div(E))

A_T2A = sparse(NDOF,NDOF); % Term T2A: -k0^2 (I.E)
A_T2B = sparse(NDOF,NDOF); % Term T2B: -k0^2 ((eps-I).E)

A_BC = sparse(NDOF,NDOF);

% Fill matrix with grouping by component

% Construct BC and RHS

for jj=1:N
    
    if jj == 1
        
        if periodic
            
            A_BC(jj +0*n,jj +0*n) = 1;
            A_BC(jj +1*n,jj +1*n) = 1;
            A_BC(jj +2*n,jj +2*n) = 1;
            
            b(jj +0*n) = ExL;
            b(jj +1*n) = EyL;
            b(jj +2*n) = EzL;
            
        elseif ldirichlet
            
            A_BC(jj +0*n,jj +0*n) = 1;
            A_BC(jj +1*n,jj +1*n) = 1;
            A_BC(jj +2*n,jj +2*n) = 1;
            
            b(jj +0*n) = ExL;
            b(jj +1*n) = EyL;
            b(jj +2*n) = EzL;
            
        end
        
    elseif jj == N
        
        
        if periodic
            
        elseif rdirichlet
            
            A_BC(jj +0*n,jj +0*n) = 1;
            A_BC(jj +1*n,jj +1*n) = 1;
            A_BC(jj +2*n,jj +2*n) = 1;
            
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
        
        % RHS
        
        [jA_x,jA_y,jA_z] = jA(x(jj));
        
        [Sx,Sy,Sz] = S(x(jj));
        
        b(jj +0*n) = 1i*w*u0*jA_x + Sx;
        b(jj +1*n) = 1i*w*u0*jA_y + Sy;
        b(jj +2*n) = 1i*w*u0*jA_z + Sz;
        
    end
    
end

% Construct inner rows

for jj=1:N
    
    this_eps = eps(x(jj));
    
    exx = this_eps(1,1);
    exy = this_eps(1,2);
    exz = this_eps(1,3);
    
    eyx = this_eps(2,1);
    eyy = this_eps(2,2);
    eyz = this_eps(2,3);
    
    ezx = this_eps(3,1);
    ezy = this_eps(3,2);
    ezz = this_eps(3,3);
    
    if jj == 1
        % BC dealt with above
    elseif jj == N
        % BC dealt with above
    else
        
        jm = jj-1;
        jp = jj+1;
        
        if periodic && jp == n+1
            jp = 1;
        end
        
        %% Term T1A: -Laplacian(E)
        
        % Ex
        
        A_T1A(jj +0*n,jm +0*n) = -1/h^2;
        A_T1A(jj +0*n,jj +0*n) = 2.0/h^2 + ky^2 + kz^2;
        A_T1A(jj +0*n,jp +0*n) = -1/h^2;
        
        A_T1A(jj +0*n,jm +1*n) = 0;
        A_T1A(jj +0*n,jj +1*n) = 0;
        A_T1A(jj +0*n,jp +1*n) = 0;
        
        A_T1A(jj +0*n,jm +2*n) = 0;
        A_T1A(jj +0*n,jj +2*n) = 0;
        A_T1A(jj +0*n,jp +2*n) = 0;
        
        % Ey
        
        A_T1A(jj +1*n,jm +0*n) = 0;
        A_T1A(jj +1*n,jj +0*n) = 0;
        A_T1A(jj +1*n,jp +0*n) = 0;
        
        A_T1A(jj +1*n,jm +1*n) = -1/h^2;
        A_T1A(jj +1*n,jj +1*n) = +2/h^2 + ky^2 + kz^2;
        A_T1A(jj +1*n,jp +1*n) = -1/h^2;
        
        A_T1A(jj +1*n,jm +2*n) = 0;
        A_T1A(jj +1*n,jj +2*n) = 0;
        A_T1A(jj +1*n,jp +2*n) = 0;
        
        % Ez
        
        A_T1A(jj +2*n,jm +0*n) = 0;
        A_T1A(jj +2*n,jj +0*n) = 0;
        A_T1A(jj +2*n,jp +0*n) = 0;
        
        A_T1A(jj +2*n,jm +1*n) = 0;
        A_T1A(jj +2*n,jj +1*n) = 0;
        A_T1A(jj +2*n,jp +1*n) = 0;
        
        A_T1A(jj +2*n,jm +2*n) = -1/h^2;
        A_T1A(jj +2*n,jj +2*n) = +2/h^2 + ky^2 + kz^2;
        A_T1A(jj +2*n,jp +2*n) = -1/h^2;
        
        
        %% Term T1B: +Grad(Div(E))
        
        % Ex
        
        A_T1B(jj +0*n,jm +0*n) = 1/h^2;
        A_T1B(jj +0*n,jj +0*n) = -2.0/h^2;
        A_T1B(jj +0*n,jp +0*n) = 1/h^2;
        
        A_T1B(jj +0*n,jm +1*n) = -i*ky/(2*h);
        A_T1B(jj +0*n,jj +1*n) = 0;
        A_T1B(jj +0*n,jp +1*n) = +i*ky/(2*h);
        
        A_T1B(jj +0*n,jm +2*n) = -i*kz/(2*h);
        A_T1B(jj +0*n,jj +2*n) = 0;
        A_T1B(jj +0*n,jp +2*n) = +i*kz/(2*h);
        
        % Ey
        
        A_T1B(jj +1*n,jm +0*n) = -i*ky/(2*h);
        A_T1B(jj +1*n,jj +0*n) = 0;
        A_T1B(jj +1*n,jp +0*n) = +i*ky/(2*h);
        
        A_T1B(jj +1*n,jm +1*n) = 0;
        A_T1B(jj +1*n,jj +1*n) = -ky^2;
        A_T1B(jj +1*n,jp +1*n) = 0;
        
        A_T1B(jj +1*n,jm +2*n) = 0;
        A_T1B(jj +1*n,jj +2*n) = -ky * kz;
        A_T1B(jj +1*n,jp +2*n) = 0;
        
        % Ez
        
        A_T1B(jj +2*n,jm +0*n) = -i*kz/(2*h);
        A_T1B(jj +2*n,jj +0*n) = 0;
        A_T1B(jj +2*n,jp +0*n) = +i*kz/(2*h);
        
        A_T1B(jj +2*n,jm +1*n) = 0;
        A_T1B(jj +2*n,jj +1*n) = -ky*kz;
        A_T1B(jj +2*n,jp +1*n) = 0;
        
        A_T1B(jj +2*n,jm +2*n) = 0;
        A_T1B(jj +2*n,jj +2*n) = -kz^2;
        A_T1B(jj +2*n,jp +2*n) = 0;
        
        
        %% Term T2A: -k0^2 (I.E)
        
        % Ex
        
        A_T2A(jj +0*n,jm +0*n) = 0;
        A_T2A(jj +0*n,jj +0*n) = -k0^2;
        A_T2A(jj +0*n,jp +0*n) = 0;
        
        A_T2A(jj +0*n,jm +1*n) = 0;
        A_T2A(jj +0*n,jj +1*n) = 0;
        A_T2A(jj +0*n,jp +1*n) = 0;
        
        A_T2A(jj +0*n,jm +2*n) = 0;
        A_T2A(jj +0*n,jj +2*n) = 0;
        A_T2A(jj +0*n,jp +2*n) = 0;
        
        % Ey
        
        A_T2A(jj +1*n,jm +0*n) = 0;
        A_T2A(jj +1*n,jj +0*n) = 0;
        A_T2A(jj +1*n,jp +0*n) = 0;
        
        A_T2A(jj +1*n,jm +1*n) = 0;
        A_T2A(jj +1*n,jj +1*n) = -k0^2;
        A_T2A(jj +1*n,jp +1*n) = 0;
        
        A_T2A(jj +1*n,jm +2*n) = 0;
        A_T2A(jj +1*n,jj +2*n) = 0;
        A_T2A(jj +1*n,jp +2*n) = 0;
        
        % Ez
        
        A_T2A(jj +2*n,jm +0*n) = 0;
        A_T2A(jj +2*n,jj +0*n) = 0;
        A_T2A(jj +2*n,jp +0*n) = 0;
        
        A_T2A(jj +2*n,jm +1*n) = 0;
        A_T2A(jj +2*n,jj +1*n) = 0;
        A_T2A(jj +2*n,jp +1*n) = 0;
        
        A_T2A(jj +2*n,jm +2*n) = 0;
        A_T2A(jj +2*n,jj +2*n) = -k0^2;
        A_T2A(jj +2*n,jp +2*n) = 0;
        
        
        %% Term T2B: -k0^2 ((eps-I).E)
        
        % Ex
        
        A_T2B(jj +0*n,jm +0*n) = 0;
        A_T2B(jj +0*n,jj +0*n) = -(-1+exx)*k0^2;
        A_T2B(jj +0*n,jp +0*n) = 0;
        
        A_T2B(jj +0*n,jm +1*n) = 0;
        A_T2B(jj +0*n,jj +1*n) = -exy*k0^2;
        A_T2B(jj +0*n,jp +1*n) = 0;
        
        A_T2B(jj +0*n,jm +2*n) = 0;
        A_T2B(jj +0*n,jj +2*n) = -exz*k0^2;
        A_T2B(jj +0*n,jp +2*n) = 0;
        
        % Ey
        
        A_T2B(jj +1*n,jm +0*n) = 0;
        A_T2B(jj +1*n,jj +0*n) = -eyx*k0^2;
        A_T2B(jj +1*n,jp +0*n) = 0;
        
        A_T2B(jj +1*n,jm +1*n) = 0;
        A_T2B(jj +1*n,jj +1*n) = -(-1+eyy)*k0^2;
        A_T2B(jj +1*n,jp +1*n) = 0;
        
        A_T2B(jj +1*n,jm +2*n) = 0;
        A_T2B(jj +1*n,jj +2*n) = -eyz*k0^2;
        A_T2B(jj +1*n,jp +2*n) = 0;
        
        % Ez
        
        A_T2B(jj +2*n,jm +0*n) = 0;
        A_T2B(jj +2*n,jj +0*n) = -ezx*k0^2;
        A_T2B(jj +2*n,jp +0*n) = 0;
        
        A_T2B(jj +2*n,jm +1*n) = 0;
        A_T2B(jj +2*n,jj +1*n) = -ezy*k0^2;
        A_T2B(jj +2*n,jp +1*n) = 0;
        
        A_T2B(jj +2*n,jm +2*n) = 0;
        A_T2B(jj +2*n,jj +2*n) = -(-1+ezz)*k0^2;
        A_T2B(jj +2*n,jp +2*n) = 0;
        
        
        %         %% Original Curl x Curl x version
        %
        %         % Ex
        %
        %         A(jj +0*n,jm +0*n) = 0; A(jj +0*n,jj +0*n) = -exx * k0^2
        %         + ky^2 + kz^2; A(jj +0*n,jp +0*n) = 0;
        %
        %         A(jj +0*n,jm +1*n) = -1i * ky / (2*h); A(jj +0*n,jj +1*n)
        %         = -exy * k0^2; A(jj +0*n,jp +1*n) = +1i * ky / (2*h);
        %
        %         A(jj +0*n,jm +2*n) = -1i * kz / (2*h); A(jj +0*n,jj +2*n)
        %         = -exz * k0^2; A(jj +0*n,jp +2*n) = +1i * kz / (2*h);
        %
        %         % Ey
        %
        %         A(jj +1*n,jm +0*n) = -1i * ky / (2*h); A(jj +1*n,jj +0*n)
        %         = -eyx * k0^2; A(jj +1*n,jp +0*n) = +1i * ky / (2*h);
        %
        %         A(jj +1*n,jm +1*n) = -1/h^2; A(jj +1*n,jj +1*n) = +2/h^2
        %         - eyy * k0^2 + kz^2; A(jj +1*n,jp +1*n) = -1/h^2;
        %
        %         A(jj +1*n,jm +2*n) = 0; A(jj +1*n,jj +2*n) = -eyz * k0^2
        %         - ky * kz; A(jj +1*n,jp +2*n) = 0;
        %
        %         % Ez
        %
        %         A(jj +2*n,jm +0*n) = -1i * kz / (2*h); A(jj +2*n,jj +0*n)
        %         = -ezx * k0^2; A(jj +2*n,jp +0*n) = +1i * kz / (2*h);
        %
        %         A(jj +2*n,jm +1*n) = 0; A(jj +2*n,jj +1*n) = -ezy * k0^2
        %         - ky * kz; A(jj +2*n,jp +1*n) = 0;
        %
        %         A(jj +2*n,jm +2*n) = -1/h^2; A(jj +2*n,jj +2*n) = +2/h^2
        %         - ezz * k0^2 + ky^2; A(jj +2*n,jp +2*n) = -1/h^2;
        
    end
    
end

% Assemble the components of the equation

A = A_BC + A_T1A*terms(1) + A_T1B*terms(2) + A_T2A*terms(3) + A_T2B*terms(4);

% Solve

E = A\b;

Ex = E(0*n+1:1*n);
Ey = E(1*n+1:2*n);
Ez = E(2*n+1:3*n);

if periodic
    Ex = [Ex.',Ex(1)];
    Ey = [Ey.',Ey(1)];
    Ez = [Ez.',Ez(1)];
    E  = [Ex,Ey,Ez].';
end

if compareWithAnalytic
    
    [ExA,EyA,EzA] = EAnalytic(x);
    
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
    
end


end
