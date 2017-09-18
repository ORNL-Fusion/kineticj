function [x] = kj_vec_to_x(jr,jt,jz)

jr_re = real(jr);
jr_im = imag(jr);
jt_re = real(jt);
jt_im = imag(jt);
jz_re = real(jz);
jz_im = imag(jz);

x = [jr_re',jr_im',jt_re',jt_im',jz_re',jz_im']';

end