function [E,err_L2,err_max] = kj_cold1d(f,xMin,xMax,nPts,ky,kz,MMSSourceFunction,MMSAnalyticFunction)

if ~exist('xMin','var')     || isempty(xMin)
    xMin = -1;
end
if ~exist('xMax','var')     || isempty(xMax)
    xMax = +1;
end
if ~exist('f','var')        || isempty(f)
    f = 13e6;
end
if ~exist('nPts','var')     || isempty(nPts)
    nPts = 512;
end
if ~exist('ky','var')       || isempty(ky)
    ky = 0;
end
if ~exist('kz','var')       || isempty(kz)
    kz = 0;
end

% Boundary conditions

% lbc = 'periodic';
% rbc = 'periodic';
lbc = 'dirichlet';
rbc = 'dirichlet';

ExL = 0;
EyL = 0;
EzL = 0;

ExR = 0;
EyR = 0;
EzR = 0;

periodic = 0;
ldirichlet = 1;
rdirichlet = 1;

if strcmp(lbc,'periodic') || strcmp(rbc,'periodic')
    n = nPts-1;
    periodic = 1;
else
    n = nPts;
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

h = (xMax-xMin) / (nPts-1);

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
dens=[1,1]*4e19;
nu_omg=0;
B=0.5;

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

jA_x = complex(zeros(n,1));
jA_y = complex(zeros(n,1));
jA_z = complex(zeros(n,1));

loc = (xMax-xMin)/2+xMin;
sig = (xMax-xMin)/10;

jA_x = 1*exp(-(x-loc).^2/sig^2)*(1+1i);
jA_y = 1*exp(-(x-loc).^2/sig^2)*(1+1i);
jA_z = 1*exp(-(x-loc).^2/sig^2)*(1+1i);

% Get analytic solution

[ExA,EyA,EzA] = analyticSolutionCP(x);

ExA = ExA*0;
EyA = EyA*0;
EzA = EzA*0;

% Fill matrix with grouping by component

for jj=1:nPts
    
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
        
    elseif jj == nPts
        
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
        
        %         Sx = -1.64038467985 * exp(5*1i*pi*x(jj));
        %         Sy = +245.097529329 * exp(5*1i*pi*x(jj));
        %         Sz = +1.67636059316 * exp(5*1i*pi*x(jj));
        
        %         Sx = complex(-302.689,-446.826) * exp(1i*10*x(jj));
        %         Sy = complex(-202.689,+446.826) * exp(1i*10*x(jj));
        %         Sz = 4.57071 * exp(1i*10*x(jj));
        
        Sx=0;Sy=0;Sz=0;
        
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

EA = [ExA,EyA,EzA].';

err_L2 = abs(sqrt(mean((E-EA).^2)));
err_max = max(abs(E-EA));


kj_plot_cmplx_3vec(E,EA)


end

function [Ex,Ey,Ez] = analyticSolution(x)

ExpVar = exp(5*1i*pi*x);

Ex = 1.00*ExpVar;
Ey = 1.00*ExpVar;
Ez = 0.01*ExpVar;

end

function [Ex,Ey,Ez] = analyticSolutionCP(x)

ExpVar = exp(1i*10*x);

Ex = 1.00*ExpVar;
Ey = 1.00*ExpVar;
Ez = 0.01*ExpVar;

end