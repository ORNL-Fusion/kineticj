function kj_iterate_aa ()

f1 = figure;

useAR = 0;

% Get the initial guess at x, either zero (startIteration=0)
% or resume from startIteration = N. 

startIteration = 0;

if startIteration == 0
    if useAR
        startRunDir = 'template-ar';
    else
        startRunDir = 'template-rs';
    end
else
    startRunDir = char(strcat('run',sprintf('%4.4i',startIteration)));
end

deltaFile = strcat(startRunDir,'/kj-delta-in.nc');

jr_re = ncread(deltaFile,'jP_r_re');
jr_im = ncread(deltaFile,'jP_r_im');
jt_re = ncread(deltaFile,'jP_t_re');
jt_im = ncread(deltaFile,'jP_t_im');
jz_re = ncread(deltaFile,'jP_z_re');
jz_im = ncread(deltaFile,'jP_z_im');

if startIteration == 0
    jr = complex(jr_re,jr_re)*0;
    jt = complex(jt_re,jt_im)*0;
    jz = complex(jz_re,jz_im)*0;
else
    jr = complex(jr_re,jr_re);
    jt = complex(jt_re,jt_im);
    jz = complex(jz_re,jz_im);
end


x = kj_vec_to_x(jr,jt,jz);

% Iterate x

g = @kj_g;
beta = @kj_damping;
beta_N = 5;
maxIterations = 1000;
atol = 1.0e-10;
rtol = 1.0e-10;
AAstart = 0;
mMax = 20;

[x,iter,res_hist] = AndAcc(g,x,mMax,maxIterations,...
    atol,rtol,1.0e10,beta,AAstart,f1,beta_N);

[jr2,jt2,jz2] = kj_x_to_vec(x);

end