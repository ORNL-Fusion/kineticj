% ------------------------------------------------- 
% script to test  preconditioners such as laplacian
% ------------------------------------------------- 
testCase = 1;
Amat = kj_wave1d_cold_plasma_test(testCase);

mval = 30;
maxit = 30;
rtol = 1e-6;

use_diagonal = 0;
clf;

ncase = numel(Amat);
nrow = ceil( sqrt(ncase) );
ncol = nrow;


for icase=1:ncase,
  A = Amat{icase}.A;
  b = Amat{icase}.b;
  n = size(A,1);
  D = spdiag( diag(A,0),0);


  disp(sprintf('icase=%d, n=%d', ... 
                icase,    n ));
  % ---------------------
  % 2nd derivative matrix
  % ---------------------
  C = spdiag( 2*ones(n,1),0) + ...
      spdiag(-1*ones(n-1,1),1) + ...
      spdiag(-1*ones(n-1,1),-1);

  L = chol(C)';

  [x,flag,relres,iter,resvec] = ...
         gmres(A,b, mval, rtol, maxit,  L,L');

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
 
  subplot(3,3,icase);
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
    % legend('laplacian','Hermitian');

   end;
  title(sprintf('icase=%d,n=%d',icase,n));
  % pause(1);
end;

