function [] = kj_cold1d_twoTest()

f = 13e6;
xMin = -1;
xMax = +1;
nPts = 512;
ky=0;
kz=0;

B = linspace(0.1,1.0,100);

power = zeros(nPts,numel(B));

dx = (xMax-xMin)/(nPts-1);
kxaxis = linspace(1,nPts,nPts) / ( dx * nPts ) * 2*pi;

for b=1:numel(B)
    
    E = kj_cold1d(f,xMin,xMax,nPts,ky,kz,B(b));
    Ex = E(1:nPts);
    power(:,b) = abs(fft(Ex)).^2;
    power(:,b) = power(:,b)/sum(power(:,b));

end

figure()
nLev=20;
contourf(B,kxaxis,power);
ylim([0,50])

figure()
semilogy(kxaxis,power(:,50))
xlim([0 80])
end