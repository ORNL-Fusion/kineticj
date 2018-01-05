function [stat] = kj_iterate_gmres()

phys = dlg_constants();

c = phys('c');
u0 = phys('u0');
eps = phys('eps0');

% Load initial guess for x (stacked E=[Er,Et,Ez] field)

initialSolutionDir = 'template-ar';

arS = ar2_read_solution(initialSolutionDir);

figure
ax1 = subplot(3,1,1);
plot(ax1,arS('r'),real(arS('E_r')))
hold on
plot(ax1,arS('r'),imag(arS('E_r')))
ax2 = subplot(3,1,2);
plot(ax2,arS('r'),real(arS('E_t')))
hold on
plot(ax2,arS('r'),imag(arS('E_t')))
ax3 = subplot(3,1,3);
plot(ax3,arS('r'),real(arS('E_z')))
hold on
plot(ax3,arS('r'),imag(arS('E_z')))

Einit = [arS('E_r')',arS('E_t')',arS('E_z')']';


% Load b (RHS for that guess)

arR = ar2_read_rundata(initialSolutionDir);

figure
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
rIn = arS('r');

jA = [arR('jA_r')',arR('jA_t')',arR('jA_z')']';

RHS = -i * w * u0 * jA;


% Call GMRES (using the nested function handle defined below)

b = RHS;
restart = [];
tol = [];
maxit = [];
M1 = [];
M2 = [];
x0 = Einit;

[x,flag,relres,ite,resvec] = gmres(@kj_LHS,b,restart,tol,maxit,M1,M2,x0);


% Setup A*x = LHS evaluation function handle as nested function

    function [LHS] = kj_LHS (E)
        
        [M,N] = size(E);
        
        n = M/3;
        
        Er = E(0*n+1:1*n);
        Et = E(1*n+1:2*n);
        Ez = E(2*n+1:3*n);
                        
        [M,N] = size(rIn);
        
        LHS = zeros(size(E));
        
        LHS_r = zeros(n);
        LHS_t = zeros(n);
        LHS_z = zeros(n);
        
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
            
            
            % + w^2/c^2 * ( E + i/(w*eps_0)*Jp )
            
            term2_r = w^2/c^2 * ( Er(j) + i/(w*eps0) * jP_r(j) );
            term2_t = w^2/c^2 * ( Et(j) + i/(w*eps0) * jP_t(j) );
            term2_z = w^2/c^2 * ( Ez(j) + i/(w*eps0) * jP_z(j) );
            
            % LHS = - curl curl E + w^2/c^2 * ( E + i/(w*eps_0)*Jp )
            
            LHS_r(j) = curlcurl_r + term2_r;
            LHS_t(j) = curlcurl_t + term2_t;
            LHS_z(j) = curlcurl_z + term2_z;
            
        end
        
        % Return A*x vector for GMRES
        
        %LHS = [LHS_r,LHS_t,LHS_z]';
        
    end


end