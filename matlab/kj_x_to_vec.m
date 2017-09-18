function [jr,jt,jz] = kj_x_to_vec(x)

N = size(x,1)/6;

jr_re = x(0*N+1:0*N+N);
jr_im = x(1*N+1:1*N+N);
jt_re = x(2*N+1:2*N+N);
jt_im = x(3*N+1:3*N+N);
jz_re = x(4*N+1:4*N+N);
jz_im = x(5*N+1:5*N+N);

jr = complex(jr_re,jr_im);
jt = complex(jt_re,jt_im);
jz = complex(jz_re,jz_im);

end