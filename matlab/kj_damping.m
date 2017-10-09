function [fac] = kj_damping(it,N)

    maxFac = 0.1;
    
    th = it / N * pi;
    th(th>pi) = pi;
    fac = (1-cos(th))/2.0 * maxFac;
    
    fprintf('Iteration: %i, N: %i, beta: %f, th: %f\n',it,N,fac,th*180/pi);

    return, fac

end