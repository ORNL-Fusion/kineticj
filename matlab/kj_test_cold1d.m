function [] = kj_test_cold1d()

n1 = 3;
n2 = 14;

n=fix(2.^linspace(n1,n2,n2-n1+1));

N=numel(n);

all_err = zeros(N,1);
all_err_max = zeros(N,1);

for ii=1:N
    
    [E,err_l2,err_max] = kj_cold1d(n(ii));
    
    all_err(ii) = err_l2;
    all_err_max(ii) = err_max;
    
end


loglog(n,all_err)
hold on
% loglog(n,all_err_max)
loglog(n,(1./n).^2)

end