function kj_iterate_aa ()

f1 = figure;

% Get the initial guess at x (i.e., zero)

templateDir = 'template-rs';

deltaFile = strcat(templateDir,'/kj-delta.nc');

jr_re = ncread(deltaFile,'jP_r_re');
jr_im = ncread(deltaFile,'jP_r_im');
jt_re = ncread(deltaFile,'jP_t_re');
jt_im = ncread(deltaFile,'jP_t_im');
jz_re = ncread(deltaFile,'jP_z_re');
jz_im = ncread(deltaFile,'jP_z_im');

jr = complex(jr_re,jr_re)*0;
jt = complex(jt_re,jt_im)*0;
jz = complex(jz_re,jz_im)*0;

x = kj_vec_to_x(jr,jt,jz);

% Iterate x

g = @kj_g;
beta = @kj_damping;
beta_N = 10;

[x,iter,res_hist] = AndAcc(g,x,10,100,1.0e-10,1.0e-10,1.0e10,beta,0,f1,beta_N);

% for it=1:3
% 
%     x = kj_g(x,it);
% 
% end

[jr2,jt2,jz2] = kj_x_to_vec(x);

end