function [stat] = kj_iterate_gmes()

phys = dlg_constants();

c = phys('c');
u0 = phys('u0');
eps0 = phys('eps0');

% Load initial guess for x (stacked E=[Er,Et,Ez] field)

initialSolutionDir = 'template-ar';

arS = ar2_read_solution(initialSolutionDir);

rIn = arS('r');

E_r_init = arS('E_r')';
E_t_init = arS('E_t')';
E_z_init = arS('E_z')';

% % Perturb the correct initial guess by smoothing it to remove some of the
% % IBW. 
% 
% width = 50;
% smooth = hanning(width)/sum(hanning(width));
% 
% E_r_init = conv(E_r_init,smooth,'same');
% E_t_init = conv(E_t_init,smooth,'same');
% E_z_init = conv(E_z_init,smooth,'same');

E_init = [E_r_init,E_t_init,E_z_init]';

[M,N] = size(E_init);
n = M/3;

f1=figure();
f1.Name = 'E_init';
ax1 = subplot(3,1,1);
plot(ax1,rIn,real(E_r_init))
hold on
plot(ax1,rIn,imag(E_r_init))
ax2 = subplot(3,1,2);
plot(ax2,rIn,real(E_t_init))
hold on
plot(ax2,rIn,imag(E_t_init))
ax3 = subplot(3,1,3);
plot(ax3,rIn,real(E_z_init))
hold on
plot(ax3,rIn,imag(E_z_init))


% Load b (RHS for that guess)

arR = ar2_read_rundata(initialSolutionDir);

f2=figure();
f2.Name = 'RHS';
ax1 = subplot(3,1,1);
plot(ax1,arR('r'),real(arR('jA_r')))
hold on
plot(ax1,arR('r'),imag(arR('jA_r')))
ax2 = subplot(3,1,2);
plot(ax2,arR('r'),real(arR('jA_t')))
hold on
plot(ax2,arR('r'),imag(arR('jA_t')))
ax3 = subplot(3,1,3);
plot(ax3,arR('r'),real(arR('jA_z')))
hold on
plot(ax3,arR('r'),imag(arR('jA_z')))

f = arR('freq');
nPhi = cast(arR('nPhi'),'single');
kz = arR('kz_1d');
w = 2 * pi * f;


jA = [arR('jA_r')',arR('jA_t')',arR('jA_z')']';

RHS = -i * w * u0 * jA;


% Test my LHS function by applying it to the AORSA solution and then
% comparing LHS with RHS.

[myLHS,LHS_t1,LHS_t2] = kj_LHS(E_init);

res = myLHS;% - RHS;

res_r = res(0*n+1:1*n);
res_t = res(1*n+1:2*n);
res_z = res(2*n+1:3*n);

f5=figure();
f5.Name = 'Residual Test';
ax1 = subplot(3,1,1);
plot(ax1,rIn,real(res_r))
hold on
plot(ax1,rIn,imag(res_r))
ax2 = subplot(3,1,2);
plot(ax2,rIn,real(res_t))
hold on
plot(ax2,rIn,imag(res_t))
ax3 = subplot(3,1,3);
plot(ax3,rIn,real(res_z))
hold on
plot(ax3,rIn,imag(res_z))


% Call GMRES (using the nested function handle defined below)

b = RHS;
restart = [];
tol = [];
maxit = 5;
M1 = [];
M2 = [];
x0 = E_init;

[x,flag,relres,ite,resvec] = gmres(@kj_LHS,b,restart,tol,maxit,M1,M2,x0);

f3=figure();
f3.Name = 'Residual';
semilogy(0:maxit,resvec/norm(b),'-o');

E_r_final = x(0*n+1:1*n);
E_t_final = x(1*n+1:2*n);
E_z_final = x(2*n+1:3*n);

f4=figure();
f4.Name = 'E_final';
ax1 = subplot(3,1,1);
plot(ax1,arS('r'),real(E_r_final))
hold on
plot(ax1,arS('r'),imag(E_r_final))
ax2 = subplot(3,1,2);
plot(ax2,arS('r'),real(E_t_final))
hold on
plot(ax2,arS('r'),imag(E_t_final))
ax3 = subplot(3,1,3);
plot(ax3,arS('r'),real(E_z_final))
hold on
plot(ax3,arS('r'),imag(E_z_final))

stat = 0;

% Setup A*x = LHS evaluation function handle as nested function

    function [LHS,LHS_t1,LHS_t2] = kj_LHS (E)
        
        Er = E(0*n+1:1*n);
        Et = E(1*n+1:2*n);
        Ez = E(2*n+1:3*n);
                        
        [M,N] = size(rIn);
        
        LHS = zeros(size(E));
        
        LHS_r = zeros(1,n);
        LHS_t = zeros(1,n);
        LHS_z = zeros(1,n);
        
        LHS_t1_r = zeros(1,n);
        LHS_t1_t = zeros(1,n);
        LHS_t1_z = zeros(1,n);
        
        LHS_t2_r = zeros(1,n);
        LHS_t2_t = zeros(1,n);
        LHS_t2_z = zeros(1,n);
        
        % Call to kineticj for this E to get jP
            
        [jP_r,jP_t,jP_z] = kj_runkj(Er,Et,Ez);
        
        h = rIn(2)-rIn(1);
        
        for j=1+1:M-1
            
            r = rIn(j);
            
            % - curl curl E
            %
            % See the $KINETICJ/mathematica/kj_curl_curl_finite_difference file for the
            % derivation of this code chunk.
            
            curlcurl_r=i.*r.^(-2).*((-1).*i.*(nPhi.^2+kz.^2.*r.^2).*Er(j)+nPhi.*Et( ...
                j)+(1/2).*h.^(-1).*r.*((-1).*nPhi.*Et((-1)+j)+nPhi.*Et(1+j)+ ...
                kz.*r.*((-1).*Ez((-1)+j)+Ez(1+j))));
            
            
            curlcurl_t=r.^(-2).*((-1).*i.*nPhi.*Er(j)+Et(j)+(-1/2).*h.^(-2).*r.*( ...
                2.*h.^2.*i.^2.*kz.^2.*r.*Et(j)+h.*((-1).*Et((-1)+j)+Et(1+j)) ...
                +2.*r.*(Et((-1)+j)+(-2).*Et(j)+Et(1+j))+h.*i.*nPhi.*(Er((-1) ...
                +j)+(-1).*Er(1+j)+(-2).*h.*i.*kz.*Ez(j))));
            
            
            curlcurl_z=(1/2).*h.^(-2).*(h.*r.^(-2).*(i.*(kz.*r.*((-1).*r.*Er((-1)+ ...
                j)+r.*Er(1+j)+2.*h.*(Er(j)+i.*nPhi.*Et(j)))+(-2).*h.*i.* ...
                nPhi.^2.*Ez(j))+r.*(Ez((-1)+j)+(-1).*Ez(1+j)))+(-2).*(Ez(( ...
                -1)+j)+(-2).*Ez(j)+Ez(1+j)));
            
            % - curl curl E (again to check - from my old RS notes)
            
            ...
            
            
            % + w^2/c^2 * ( E + i/(w*eps_0)*Jp )
            
            term2_r = w^2/c^2 * ( Er(j) + i/(w*eps0) * jP_r(j) );
            term2_t = w^2/c^2 * ( Et(j) + i/(w*eps0) * jP_t(j) );
            term2_z = w^2/c^2 * ( Ez(j) + i/(w*eps0) * jP_z(j) );
            
            % LHS = - curl curl E + w^2/c^2 * ( E + i/(w*eps_0)*Jp )
            
            LHS_r(j) = curlcurl_r + term2_r;
            LHS_t(j) = curlcurl_t + term2_t;
            LHS_z(j) = curlcurl_z + term2_z;
            
            LHS_t1_r(j) = curlcurl_r;
            LHS_t1_t(j) = curlcurl_t;
            LHS_t1_z(j) = curlcurl_z;
            
            LHS_t2_r(j) = term2_r;
            LHS_t2_t(j) = term2_t;
            LHS_t2_z(j) = term2_z;
            
        end
        
        % Return A*x vector for GMRES
        
        LHS = [LHS_r,LHS_t,LHS_z]';
        LHS_t1 = [LHS_t1_r,LHS_t1_t,LHS_t1_z]';
        LHS_t2 = [LHS_t2_r,LHS_t2_t,LHS_t2_z]';
        
    end

    function [w] = hanning(N)
        
        alpha = 0.5;
        k = linspace(0,N-1,N);
        w = alpha - (1-alpha)*cos(2*pi*k/N);
        
    end

end