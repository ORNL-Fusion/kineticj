function [stat] = kj_iterate_gmes()

readonly = 0;

phys = dlg_constants();

c = phys('c');
u0 = phys('u0');
eps0 = phys('eps0');

solutionFile = 'gmres-solution.mat';

useAorsa = 0;

% Load initial guess for x (stacked E=[Er,Et,Ez] field)

initialSolutionDir = 'template';

if (useAorsa)
    sol = ar2_read_solution(initialSolutionDir);
    run = ar2_read_rundata(initialSolutionDir);
else
    sol = rs_read_solution(initialSolutionDir);
    run = rs_read_rundata(initialSolutionDir);
end

rIn = sol('r');

f = run('freq');
nPhi = cast(run('nPhi'),'single');
kz = run('kz');

E0_r = sol('E_r')';
E0_t = sol('E_t')';
E0_z = sol('E_z')';

J0_r = sol('jP_r')';
J0_t = sol('jP_t')';
J0_z = sol('jP_z')';


% % Perturb the correct initial guess by smoothing it to remove some of the
% % IBW.
%
% width = 50;
% smooth = hanning(width)/sum(hanning(width));
%
% E_r_init = conv(E_r_init,smooth,'same');
% E_t_init = conv(E_t_init,smooth,'same');
% E_z_init = conv(E_z_init,smooth,'same');

E0 = [E0_r,E0_t,E0_z]';
J0 = [J0_r,J0_t,J0_z]';

[M,N] = size(E0);
n = M/3;


% RHS

w = 2 * pi * f;

jA = [run('jA_r')',run('jA_t')',run('jA_z')']';

RHS = -i * w * u0 * jA;


% Evaluate Jacobian (J = dResidual / dE)

J = complex(zeros(M,M));
dE = complex(zeros(M,1));

f_stageResidual = @kj_runResidual;
f_readResidual = @kj_readResidual;

readonlyFileName = 'kj-rs-readonly.mat';

if readonly
    
    load(readonlyFileName,'loc0','loc');
    
else
    
    disp('staging residuals ...');
    
    loc = [];
    
    for ii=1:M
        
        E_r = E0_r;
        E_t = E0_t;
        E_z = E0_z;
        
        ix = mod(ii-1,n)+1;
        
        dE_r = 0.0; dE_t = 0.0; dE_z = 0.0;
        
        if (ii > 0*n) && (ii <= 1*n)
            dE_r = E_r(ix) * 0.1 + complex(1,1)*1e-7;
            dE(ii) = dE_r;
        end
        
        if (ii > 1*n) && (ii <= 2*n)
            dE_t = E_t(ix) * 0.1 + complex(1,1)*1e-7;
            dE(ii) = dE_t;
        end
        
        if (ii > 2*n) && (ii <= 3*n)
            dE_z = E_z(ix) * 0.1 + complex(1,1)*1e-7;
            dE(ii) = dE_z;
        end
         
        E_r(ix) = E_r(ix) + dE_r;
        E_t(ix) = E_t(ix) + dE_t;
        E_z(ix) = E_z(ix) + dE_z;
        
        E1 = [E_r,E_t,E_z]';
        
        [thisDir] = f_stageResidual(E1); % comes out as [resx1,resx2,...resxn,resy1,resy2,...,resyn,etc]
        
        loc = [loc;thisDir];
        
    end
    
    disp('done');
    
    % Run all the residual calculations in IDL
    
    disp('computing residuals ...');
    
    [loc0] = f_stageResidual(E0);
    
    [stat] = kj_runResiduals();
    
    save(readonlyFileName,'loc0','loc','dE');
    
    disp('done');
    
end

[res0] = f_readResidual(loc0);

disp('reading residuals ...');

for ii=1:M
        
    res1 = f_readResidual(loc(ii,:)); % comes out as [resx1,resx2,...resxn,resy1,resy2,...,resyn,etc]
        
    dres = res0 - res1;
    
    % Jacobian row
    % J could be evaluated with a higher order diferencing scheme.
    
    dRes_dE = dres ./ dE(ii);
    
    J(ii,:) = dRes_dE; % for collected components ordering indicated above
    
end

disp('done');

save('kj-rs-jacobian.mat','J');


% Solve using Newtons method

E1 = E0 - inv(J) * res0';

% Check residual of final solution

[loc1] = f_stageResidual(E1);

% *** Need the ability to run a single case, or only run if not run ***

loadPreviousSolution = 0;

if loadPreviousSolution
    
    load(solutionFile);
    E0 = E_final;
    
end



f1=figure();
f1.Name = 'E_init';
ax1 = subplot(3,1,1);
plot(ax1,rIn,real(E0_r))
hold on
plot(ax1,rIn,imag(E0_r))
ax2 = subplot(3,1,2);
plot(ax2,rIn,real(E0_t))
hold on
plot(ax2,rIn,imag(E0_t))
ax3 = subplot(3,1,3);
plot(ax3,rIn,real(E0_z))
hold on
plot(ax3,rIn,imag(E0_z))



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




% Test my LHS function by applying it to the AORSA solution and then
% comparing LHS with RHS.

[myLHS,LHS_t1,LHS_t2] = kj_LHS(E0);

jP_ar = ([sum(arS('jP_r'),3)',sum(arS('jP_t'),3)',sum(arS('jP_z'),3)']');

LHS_t2_ar = w^2/c^2 .* ( E0 + i/(w*eps0).*jP_ar );

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

f6=figure();
f6.Name = 'Residual terms : LHS_t1, LHS_t2, RHS';
ax1 = subplot(3,1,1);
hold on
plot(ax1,rIn,real(LHS_t1(1:n)))
plot(ax1,rIn,real(LHS_t2(1:n)))
plot(ax1,rIn,real(RHS(1:n)))
plot(ax1,rIn,real(LHS_t2_ar(1:n)))

ax2 = subplot(3,1,2);
hold on
plot(ax2,rIn,real(LHS_t1(1+n:2*n)))
plot(ax2,rIn,real(LHS_t2(1+n:2*n)))
plot(ax2,rIn,real(RHS(1+n:2*n)))

ax3 = subplot(3,1,3);
hold on
plot(ax3,rIn,real(LHS_t1(1+2*n:3*n)))
plot(ax3,rIn,real(LHS_t2(1+2*n:3*n)))
plot(ax3,rIn,real(RHS(1+2*n:3*n)))

% Read the AORSA matrix and try GMRES on that also

aorsaMatrixFile = '/Users/dg6/scratch/aorsa2d/colestock-kashuba-reference/matrix.nc';
arA = dlg_read_netcdf(aorsaMatrixFile);

ar_A = complex(arA('reA'),arA('imA'));

[T,B] = balance(ar_A);

ar_x = linsolve(ar_A, RHS);

ar_Ek_alp = ar_x(0*n+1:1*n);
ar_Ek_bet = ar_x(1*n+1:2*n);
ar_Ek_prl = ar_x(2*n+1:3*n);

ar_Ek_alp_0 = arS('Ek_alp');
ar_Ek_bet_0 = arS('Ek_bet');
ar_Ek_prl_0 = arS('Ek_prl');

f7=figure();
ax1 = subplot(3,1,1);
hold on
plot(ar_Ek_alp,ar_Ek_alp_0);

% Call GMRES (using the nested function handle defined below)

b = RHS;
restart = [];
tol = [];
maxit = 250;
M1 = [];
M2 = [];
x0 = E0;

[x,flag,relres,ite,resvec] = gmres(@kj_LHS,b,restart,tol,maxit,M1,M2,x0);

E_final = x;

save(solutionFile, 'E_final');

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

    function [LHS] = kj_LHS2 (E)
        
        Er = E(0*n+1:1*n);
        Et = E(1*n+1:2*n);
        Ez = E(2*n+1:3*n);
                        
        [M,N] = size(rIn);
        
        % Update jP(E) via call to kinetic-j
            
        [jP_r,jP_t,jP_z] = kj_runkj(Er,Et,Ez);
        
        % Update E(jP) via call to full-wave solve
        
        
        h = rIn(2)-rIn(1);
        
        % term1 = -curlxcurl(E)
        
        t1_r = -(i*nPhi./rIn.^2.*gradient(rIn.*Et,h) + nPhi^2./rIn.^2.*Er + kz^2*Er + i*kz.*gradient(Ez,h));
        t1_t = -(-kz*nPhi./rIn.*Ez + kz^2.*Et - gradient(gradient(rIn.*Et,h)./rIn,h) + i*nPhi.*gradient(Er./rIn,h));
        t1_z = -(i*kz./rIn.*gradient(rIn.*Er,h) - 1./rIn.*gradient(rIn.*gradient(Ez,h),h) + nPhi^2./rIn.^2.*Ez - nPhi*kz./rIn.*Et);
        
        % term2 = +w^2/c^2 * ( E + i/(w*e0) * jP )
        
        t2_r = w^2/c^2 .* ( Er + i/(w*eps0) .* jP_r );
        t2_t = w^2/c^2 .* ( Et + i/(w*eps0) .* jP_t );
        t2_z = w^2/c^2 .* ( Ez + i/(w*eps0) .* jP_z );
        
        LHS_r = t1_r' + t2_r';
        LHS_t = t1_t' + t2_t';
        LHS_z = t1_z' + t2_z';
        
        % Return A*x vector for GMRES
        
        LHS = [LHS_r,LHS_t,LHS_z]';
        
        LHS_t1 = [t1_r',t1_t',t1_z']';
        LHS_t2 = [t2_r',t2_t',t2_z']';
        
    end

    function [LHS,LHS_t1,LHS_t2] = kj_LHS (E)
        
        Er = E(0*n+1:1*n);
        Et = E(1*n+1:2*n);
        Ez = E(2*n+1:3*n);
                        
        [M,N] = size(rIn);
        
        % Call to kineticj for this E to get jP
            
        [jP_r,jP_t,jP_z] = kj_runkj(Er,Et,Ez);
        
        h = rIn(2)-rIn(1);
        
        % term1 = -curlxcurl(E)
        
        t1_r = -(i*nPhi./rIn.^2.*gradient(rIn.*Et,h) + nPhi^2./rIn.^2.*Er + kz^2*Er + i*kz.*gradient(Ez,h));
        t1_t = -(-kz*nPhi./rIn.*Ez + kz^2.*Et - gradient(gradient(rIn.*Et,h)./rIn,h) + i*nPhi.*gradient(Er./rIn,h));
        t1_z = -(i*kz./rIn.*gradient(rIn.*Er,h) - 1./rIn.*gradient(rIn.*gradient(Ez,h),h) + nPhi^2./rIn.^2.*Ez - nPhi*kz./rIn.*Et);
        
        % term2 = +w^2/c^2 * ( E + i/(w*e0) * jP )
        
        t2_r = w^2/c^2 .* ( Er + i/(w*eps0) .* jP_r );
        t2_t = w^2/c^2 .* ( Et + i/(w*eps0) .* jP_t );
        t2_z = w^2/c^2 .* ( Ez + i/(w*eps0) .* jP_z );
        
        LHS_r = t1_r' + t2_r';
        LHS_t = t1_t' + t2_t';
        LHS_z = t1_z' + t2_z';
        
        % Return A*x vector for GMRES
        
        LHS = [LHS_r,LHS_t,LHS_z]';
        
        LHS_t1 = [t1_r',t1_t',t1_z']';
        LHS_t2 = [t2_r',t2_t',t2_z']';
        
    end

    function [J] = kj_update(E)
        
        Er = E(0*n+1:1*n);
        Et = E(1*n+1:2*n);
        Ez = E(2*n+1:3*n);
        
        [jP_r,jP_t,jP_z] = kj_runkj(Er,Et,Ez);
        
        J = [jP_r,jP_t,jP_z]';
        
    end

    function [res] = kj_residual(E,J,RHS)
        
        Er = E(0*n+1:1*n);
        Et = E(1*n+1:2*n);
        Ez = E(2*n+1:3*n);
        
        jP_r = J(0*n+1:1*n);
        jP_t = J(1*n+1:2*n);
        jP_z = J(2*n+1:3*n);
        
        RHS_r = RHS(0*n+1:1*n);
        RHS_t = RHS(1*n+1:2*n);
        RHS_z = RHS(2*n+1:3*n);
        
        h = rIn(2)-rIn(1);
        
        % term1 = -curlxcurl(E)
        
        t1_r = -(i*nPhi./rIn.^2.*gradient(rIn.*Et,h) + nPhi^2./rIn.^2.*Er + kz^2*Er + i*kz.*gradient(Ez,h));
        t1_t = -(-kz*nPhi./rIn.*Ez + kz^2.*Et - gradient(gradient(rIn.*Et,h)./rIn,h) + i*nPhi.*gradient(Er./rIn,h));
        t1_z = -(i*kz./rIn.*gradient(rIn.*Er,h) - 1./rIn.*gradient(rIn.*gradient(Ez,h),h) + nPhi^2./rIn.^2.*Ez - nPhi*kz./rIn.*Et);
        
        % term2 = +w^2/c^2 * ( E + i/(w*e0) * jP )
        
        t2_r = w^2/c^2 .* ( Er + i/(w*eps0) .* jP_r );
        t2_t = w^2/c^2 .* ( Et + i/(w*eps0) .* jP_t );
        t2_z = w^2/c^2 .* ( Ez + i/(w*eps0) .* jP_z );
        
        LHS_r = t1_r' + t2_r';
        LHS_t = t1_t' + t2_t';
        LHS_z = t1_z' + t2_z';
        
        % Return A*x vector for GMRES
        
        LHS = [LHS_r,LHS_t,LHS_z]';
        
        res = RHS - LHS;
        
        f6=figure();
        f6.Name = 'Residual terms : LHS_t1, LHS_t2, RHS';
        ax1 = subplot(3,1,1);
        hold on
        plot(ax1,rIn,real(t1_r))
        plot(ax1,rIn,real(t2_r))
        plot(ax1,rIn,real(RHS_r))
        plot(ax1,rIn,real(RHS_r - (t1_r+t2_r)),'LineWidth',3)
        
        ax2 = subplot(3,1,2);
        hold on
        plot(ax2,rIn,real(t1_t))
        plot(ax2,rIn,real(t2_t))
        plot(ax2,rIn,real(RHS_t))
        plot(ax2,rIn,real(RHS_t - (t1_t+t2_t)),'LineWidth',3)
        
        ax3 = subplot(3,1,3);
        hold on
        plot(ax3,rIn,real(t1_z))
        plot(ax3,rIn,real(t2_z))
        plot(ax3,rIn,real(RHS_z))
        plot(ax3,rIn,real(RHS_z - (t1_z+t2_z)),'LineWidth',3)
        
    end

    function [w] = hanning(N)
        
        alpha = 0.5;
        k = linspace(0,N-1,N);
        w = alpha - (1-alpha)*cos(2*pi*k/N);
        
    end

end