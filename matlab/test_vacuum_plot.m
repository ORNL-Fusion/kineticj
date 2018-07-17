% ------------------------------------------------- 
% script to test  preconditioners such as laplacian
% ------------------------------------------------- 
testCase = 1;
phys = dlg_constants();

Amat = kj_wave1d_vacuum_test(testCase);

mval = 30;
maxit = 30;
rtol = 1e-6;

use_diagonal = 0;
clf;

ncase = numel(Amat);

icase = ncase;
  A = Amat{icase}.A;
  b = Amat{icase}.b;
  n = size(A,1);
  m = round(n/3);
  D = spdiag( diag(A,0),0);

  xMin = -1;
  xMax = 1;
  h = (xMax-xMin)/(n-1);

  c = phys.('c');
  u0 = phys.('u0');
  eps0 = phys.('eps0');
  me_amu = phys.('me_amu');
  amu0 = phys.('amu');

 f = 13e6;

  w = 2*pi*f;
  k0 = w/c;


  disp(sprintf('icase=%d, n=%d', ... 
                icase,    n ));
  % ---------------------
  % 2nd derivative matrix
  % ---------------------
  C = spdiag( 2*ones(m,1),0) + ...
      spdiag(-1*ones(m-1,1),1) + ...
      spdiag(-1*ones(m-1,1),-1);
  C = C/(h*h);
  [ia,ja,aa] = find(C);
  C = sparse( [ia,ia+m,ia+2*m],  [ja,ja+m,ja+2*m], [aa,aa,aa]);

  % [L,U] = lu( C + 1i*k0*k0*speye(n,n) );
  U = chol(C);
  L = U';

  [x,flag,relres,iter,resvec] = ...
         gmres(A,b, mval, rtol, maxit,  L,U);

  x_lap = x;
  resvec_lap = resvec;

  % --------------
  % symmetric part
  % --------------
  [L,U] = lu( 0.5*(A+A') );
  [x,flag,relres,iter,resvec] = ...
         gmres(A,b, mval, rtol, maxit,  L,U);
  
  x_sym = x;
  resvec_sym = resvec;

  % -----------------------
  % diagonal preconditioner
  % -----------------------
  if (use_diagonal),
    [x, flag, relres,iter,resvec] = ...
         gmres(A,b, mval, rtol, maxit );
    x_diag = x;
    resvec_diag = resvec;
    res_diag = norm(b-A*x_diag);
  end;


  res_lap = norm(b-A*x_lap);
  res_sym = norm(b-A*x_sym);
  norm_b = norm(b);
 
  if (use_diagonal),
    disp(sprintf('res_lap=%e, res_sym=%e, res_diag=%e, norm_b=%e', ...
                  res_lap,    res_sym,    res_diag,    norm_b));

    semilogy( 1:length(resvec_lap), resvec_lap ./resvec_lap(1), 'r.-', ...
            1:length(resvec_sym), resvec_sym ./resvec_sym(1), 'b.-', ...
            1:length(resvec_diag), resvec_diag ./resvec_diag(1), 'g.-' );
    legend('laplacian','Hermitian','diagonal');
   else
    disp(sprintf('res_lap=%e, res_sym=%e,  norm_b=%e', ...
                  res_lap,    res_sym,     norm_b));

    semilogy( 1:length(resvec_lap), resvec_lap ./resvec_lap(1), 'r.-', ...
            1:length(resvec_sym), resvec_sym ./resvec_sym(1), 'b.-' );
    legend('laplacian','Hermitian');

   end;
  title(sprintf('Vacuum problem, n=%d',n));
  % pause(1);

