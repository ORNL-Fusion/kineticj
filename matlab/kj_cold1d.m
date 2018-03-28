function [E,err_L2,err_max] = kj_cold1d(nPts)

phys = dlg_constants();

c = phys('c');
u0 = phys('u0');
eps0 = phys('eps0');

%%% Inputs

f = 13e6;
%nPts = 256;

n = nPts-1; % for periodic

ky = 0;
kz = 0;

xMin = -1;
xMax = +1;

h = (xMax-xMin) / (n-1);

x = linspace(xMin,xMax-h,n);

w = 2*pi*f;
k0 = w^2/c^2;

N = n*3;

A = complex(zeros(N,N));
b = complex(zeros(N,1));

% Dielectric tensor

exx = 1;
exy = 0;
exz = 0;

eyx = 0;
eyy = 1;
eyz = 0;

ezx = 0;
ezy = 0;
ezz = 1;


% Current source

jA_x = complex(zeros(n,1));
jA_y = complex(zeros(n,1));
jA_z = complex(zeros(n,1));


% Get analytic solution

[ExA,EyA,EzA] = analyticSolution(x);

% Fill matrix with grouping by component

for jj=1:n
    
    if jj == 1
        
        % With periodic BC's we have to set an offset value
        % Here we set it to be the analytic solution for 
        % comparison purposes. 
        
        A(jj +0*n,jj +0*n) = 1;
        A(jj +1*n,jj +1*n) = 1;
        A(jj +2*n,jj +2*n) = 1;
        
        b(jj +0*n) = ExA(1);
        b(jj +1*n) = EyA(1);
        b(jj +2*n) = EzA(1);
        
    else
        
        % Ex
        
        jm = jj-1;
        jp = jj+1;
        
        if jm == 0
            jm = n;
        end
        
        if jp == n+1
            jp = 1;
        end
        
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
        
%             Sx = +44.8839 * exp(5*1i*pi*x(jj));
%             Sy = +267.695 * exp(5*1i*pi*x(jj));
%             Sz = -204.362 * exp(5*1i*pi*x(jj));
        
        Sx = -0.00551074 * exp(5*1i*pi*x(jj));
        Sy = +246.735 * exp(5*1i*pi*x(jj));
        Sz = +2.46735 * exp(5*1i*pi*x(jj));
        
        b(jj +0*n) = 1i*w*u0*jA_x(jj) + Sx;
        b(jj +1*n) = 1i*w*u0*jA_y(jj) + Sy;
        b(jj +2*n) = 1i*w*u0*jA_z(jj) + Sz;
        
    end
    
end

E = linsolve(A,b);

Ex = E(0*n+1:1*n);
Ey = E(1*n+1:2*n);
Ez = E(2*n+1:3*n);

EA = [ExA,EyA,EzA]';

err_L2 = norm(E-EA)/N;
err_max = max(abs(E-EA));

% kj_plot_cmplx_3vec(E,EA)

end

function [Ex,Ey,Ez] = analyticSolution(x)

Ex = 1.00*exp(5*1i*pi*x);
Ey = 1.00*exp(5*1i*pi*x);
Ez = 0.01*exp(5*1i*pi*x);

end