function [fac] = kj_damping(it,N)

    th = it / N * pi;
    th(th>pi) = 0;
    fac = 1-cos(th);
    
    fprintf('Iteration: %i, N: %i, beta: %f\n',it,N,fac);

    return, fac

end