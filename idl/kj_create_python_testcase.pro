pro kj_create_python_testcase

@dlg_constants

f = 13e6
q = 1
amu = 1
kPar = 0  
T_eV = 5e3

xMin=0.1
xMax=10.0
b0 = 1.5
nPts = 1000

xRange = xMax-xMin 
dx = xRange / (nPts-1)
x = fIndGen(nPts)*dx+xMin

b = b0/x

p=plot(x,b)

m = (amu*_amu)
wc = q*_e*b/m
w = 2*!pi*f
vth = sqrt(2*T_eV*_e/m)

rhoL = m * vTh / ( abs(q*_e) * b ) 

res2 = w-kPar*vth-2*wc
loc = where(res2 eq max(res2))

kx = (!pi/(2*rhoL[loc]))[0]
p=plot(x,w/wc)

res = w-kPar*vth-1*wc
xrange=[0,6]
p=plot(x,1/abs(res),/ylog,xrange=xrange)

for n=-5,5 do begin
    res = w-kPar*vth-n*wc
    p=plot(x,1/abs(res),/over)
endfor

kj_create_single_k_input, b0=b, kx=kx, f_Hz=f, x=x, nPts = nPts, /writeOutput, $
        E1Multiplier=1, $
        E2Multiplier=-_ii, $
        E3Multiplier=0

stop

end
