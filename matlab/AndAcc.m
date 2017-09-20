function [x,iter,res_hist] = AndAcc(g,x,mMax,itmax,atol,rtol,droptol,beta,AAstart)
% This performs fixed-point iteration with or without Anderson
% acceleration for a given fixed-point map g and initial
% approximate solution x.
%
% Required inputs:
%   g = fixed-point map (function handle); form gval = g(x).
%   x = initial approximate solution (column vector).
%
% Optional inputs:
%   mMax = maximum number of stored residuals (non-negative integer).
%           NOTE: mMax = 0 => no acceleration.
%   itmax = maximum allowable number of iterations.
%   atol = absolute error tolerance.
%   rtol = relative error tolerance.
%   droptol = tolerance for dropping stored residual vectors to improve
%           conditioning: If droptol > 0, drop residuals if the
%           condition number exceeds droptol; if droptol <= 0,
%           do not drop residuals.
%   beta = damping factor: If beta > 0 (and beta ~= 1), then the step is
%           damped by beta; otherwise, the step is not damped.
%           NOTE: beta can be a function handle; form beta(iter), where iter is
%           the iteration number and 0 < beta(iter) <= 1.
%   AAstart = acceleration delay factor: If AAstart > 0, start acceleration
%           when iter = AAstart.
%
% Output:
%   x = final approximate solution.
%   iter = final iteration number.
%   res_hist = residual history matrix (iteration numbers and residual norms).
%
% Homer Walker (walker@wpi.edu), 10/14/2011.
% Set the method parameters.

if nargin < 2, error('AndAcc requires at least two arguments.'); end
if nargin < 3, mMax = min{10, size(x,1)}; end
if nargin < 4, itmax = 100; end
if nargin < 5, atol = 1.e-10; end
if nargin < 6, rtol = 1.e-10; end
if nargin < 7, droptol = 1.e10; end
if nargin < 8, beta = 1; end
if nargin < 9, AAstart = 0; end
% Initialize the storage arrays.
res_hist = []; % Storage of residual history.
DG = []; % Storage of g-value differences.
% Initialize printing.
if mMax == 0
    fprintf('\n No acceleration.');
elseif mMax > 0
    fprintf('\n Anderson acceleration, mMax = %d \n',mMax);
else
    error('AndAcc.m: mMax must be non-negative.');
end
fprintf('\n iter res_norm \n');
% Initialize the number of stored residuals.
mAA = 0;
% Top of the iteration loop.

for iter = 0:itmax
    
    % Plot the delta
    
    [jr,jt,jz] = kj_x_to_vec(x);
    idx = 1:size(jr,1);
    ax1=subplot(2,2,1);
    hold(ax1,'on');
    plot(idx',real(jr),'black',idx',imag(jr),'r');
    ax2=subplot(2,2,2);
    hold(ax2,'on');
    plot(idx',real(jt),'black',idx',imag(jt),'r');
    ax3=subplot(2,2,3);
    hold(ax3,'on');
    plot(idx',real(jz),'black',idx',imag(jz),'r');
    
    % Apply g and compute the current residual norm.
    gval = g(x,iter);
    fval = gval - x;
    res_norm = norm(fval);
    fprintf(' %d %e \n', iter, res_norm);
    res_hist = [res_hist;[iter,res_norm]];
    
    ax4=subplot(2,2,4);
    hold(ax4,'on');
    plot(res_hist);
    
    % Set the residual tolerance on the initial iteration.
    if iter == 0, tol = max(atol,rtol*res_norm); end
    % Test for stopping.
    if res_norm <= tol,
        fprintf('Terminate with residual norm = %e \n\n', res_norm);
        break;
    end
    if mMax == 0 || iter < AAstart,
        % Without acceleration, update x <- g(x) to obtain the next
        % approximate solution.
        x = gval;
    else
        % Apply Anderson acceleration.
        % Update the df vector and the DG array.
        if iter > AAstart,
            df = fval-f_old;
            if mAA < mMax,
                DG = [DG gval-g_old];
            else
                DG = [DG(:,2:mAA) gval-g_old];
            end
            mAA = mAA + 1;
        end
        f_old = fval;
        g_old = gval;
        if mAA == 0
            % If mAA == 0, update x <- g(x) to obtain the next approximate solution.
            x = gval;
        else
            % If mAA > 0, solve the least-squares problem and update the
            % solution.
            if mAA == 1
                % If mAA == 1, form the initial QR decomposition.
                R(1,1) = norm(df);
                Q = R(1,1)\df;
            else
                % If mAA > 1, update the QR decomposition.
                if mAA > mMax
                    % If the column dimension of Q is mMax, delete the first column and
                    % update the decomposition.
                    [Q,R] = qrdelete(Q,R,1);
                    mAA = mAA - 1;
                    % The following treats the qrdelete quirk described below.
                    if size(R,1) ~= size(R,2),
                        Q = Q(:,1:mAA-1); R = R(1:mAA-1,:);
                    end
                    % Explanation: If Q is not square, then qrdelete(Q,R,1) reduces the
                    % column dimension of Q by 1 and the column and row
                    % dimensions of R by 1. But if Q *is* square, then the
                    % column dimension of Q is not reduced and only the column
                    % dimension of R is reduced by one. This is to allow for
                    % MATLAB?s default "thick" QR decomposition, which always
                    % produces a square Q.
                end
                % Now update the QR decomposition to incorporate the new
                % column.
                for j = 1:mAA - 1
                    R(j,mAA) = Q(:,j)'*df;
                    df = df - R(j,mAA)*Q(:,j);
                end
                R(mAA,mAA) = norm(df);
                Q = [Q,R(mAA,mAA)\df];
            end
            if droptol > 0
                % Drop residuals to improve conditioning if necessary.
                condDF = cond(R);
                while condDF > droptol && mAA > 1
                    fprintf(' cond(D) = %e, reducing mAA to %d \n', condDF, mAA-1);
                    [Q,R] = qrdelete(Q,R,1);
                    DG = DG(:,2:mAA);
                    mAA = mAA - 1;
                    % The following treats the qrdelete quirk described above.
                    if size(R,1) ~= size(R,2),
                        Q = Q(:,1:mAA); R = R(1:mAA,:);
                    end
                    condDF = cond(R);
                end
            end
            % Solve the least-squares problem.
            gamma = R\(Q'*fval);
            % Update the approximate solution.
            x = gval - DG*gamma;
            % Apply damping if beta is a function handle or if beta > 0
            % (and beta ~= 1).
            if isa(beta,'function_handle'),
                x = x - (1-beta(iter))*(fval - Q*R*gamma);
            else
                if beta > 0 && beta ~= 1,
                    x = x - (1-beta)*(fval - Q*R*gamma);
                end
            end
        end
    end
end

% Bottom of the iteration loop.
if res_norm > tol && iter == itmax,
    fprintf('\n Terminate after itmax = %d iterations. \n', itmax);
    fprintf(' Residual norm = %e \n\n', res_norm);
end

end
