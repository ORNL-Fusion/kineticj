pro kj_sigma_vs_t, runKJ=runKJ

@constants

n = 100
tMin = 0.1
tMax = 10e3
T_eV = 10d0^(findGen(n)/(n-1)*(alog10(tMax)-alog10(tMin))+alog10(tMin)) 

f = 13d6
Z = 1d0
amu =  2.0;_me_amu
B = 1d0
density = 2d19
harmonicNumber = 1
kPar = 100  
kPer = 10 
kj_nStepsPerCycle = 30.0
kj_nRFCycles = 10.0 
 
eps = ComplexArr(3,3,n)
eps_cold = ComplexArr(3,3,n)

sig = ComplexArr(3,3,n)
sig_cold = ComplexArr(3,3,n)

w = 2*!pi*f
for i=0,n-1 do begin

    thisTeV = T_eV[i]
    this_eps = kj_hot_epsilon(f, amu, Z, B, density, $
            harmonicNumber, kPar, kPer, thisTeV, $
            epsilon_cold = this_eps_cold )
    
    eps[*,*,i] = this_eps 
    eps_cold[*,*,i] = this_eps_cold 

    sig[*,*,i] = (eps[*,*,i] - identity(3)) * w * _e0 / _ii
    sig_cold[*,*,i] = (eps_cold[*,*,i] - identity(3)) * w * _e0 / _ii

endfor

; Stage and run kj over this range of temperatures 

n_kj = 100
tMin = 0.1 
tMax = 10e3
T_eV_kj = 10d0^(findGen(n_kj)/(n_kj-1)*(alog10(tMax)-alog10(tMin))+alog10(tMin)) 
;T_eV_kj = [0.1,1.0,10.0,100.0,1000.0, 10000.0]

nT = n_elements(T_eV_kj)

sig1 = ComplexArr(nT)
sig2 = ComplexArr(nT)
sig3 = ComplexArr(nT)

TemplateRunDir = 'benchmark-perp'
RootDir = '/home/dg6/scratch/kineticJ/temp_scan'
cd, RootDir

for t=0,nT-1 do begin
    
    ThisRunDir = string(t,format='(i5.5)')

    if keyword_set(runKJ) then begin

        file_delete, ThisRunDir, /recursive, /allow_nonexistent

        file_copy, TemplateRunDir, ThisRunDir, /recursive

    endif

    cd, ThisRunDir

    kj_create_single_k_input, b0=B, kPar=kPar, kPer=kPer, f_Hz=f, n_m3=density, $
            Er=Er, Et=Et, Ez=Ez, x=x, writeOutput=runKJ

    if keyword_set(runKJ) then begin

        ; Stage the input wave fields

        ; Adjust the kj.cfg config file parameters

        kj = kj_read_cfg('./')
        kj['T_keV'] = T_eV_kj[t]*1e-3 
        kj['species_amu'] = amu 
        kj['species_Z'] = Z
        kj['ky'] = kPar
        kj['nStepsPerCycle'] = kj_nStepsPerCycle 
        kj['nRFCycles'] = kj_nRFCycles
        kj['nPx'] = 10
        kj['nPy'] = 65
        kj['nPz'] = 10

        kj_write_kj_cfg, kj, './'

        ; Run kj

        RunCommand = '~/code/kineticj/bin/kineticj'
        spawn, RunCommand, StdOut, StdErr
        print, StdOut
        print, StdErr

    endif

    ; Read in kj results

    kj_read_jp_old, x=kj_x, j1x=kj_j1x, j1y=kj_j1y, j1z=kj_j1z
    
    kj_Er = interpol(Er,x,kj_x)
    kj_Et = interpol(Et,x,kj_x)
    kj_Ez = interpol(Ez,x,kj_x)
    
    sig1[t] = (kj_j1x/kj_Er)[0]
    sig2[t] = (kj_j1y/kj_Et)[0]
    sig3[t] = (kj_j1z/kj_Ez)[0]

    cd, RootDir

endfor

layout=[3,3]
pos = 1
thick = 2 
style = '--'
transparency = 50

plotThis = sig
plotThis_cold = sig_cold

p=plot(T_eV,plotThis[0,0,*],layout=[[layout],pos],title='plotThis[0,0]',/xlog,yRange=[-1,1]*max(abs(plotThis[0,0,*])))
p=plot(T_eV,imaginary(plotThis[0,0,*]),color='r',/over)
p=plot(T_eV,plotThis_cold[0,0,*],/over,thick=thick,transparency=transparency,LineStyle=style)
p=plot(T_eV,imaginary(plotThis_cold[0,0,*]),color='r',/over,thick=thick,transparency=transparency,LineStyle=style)

p=plot(T_eV_kj, sig1, /over, thick=2)
p=plot(T_eV_kj, imaginary(sig1), color='r', /over, thick=2)

++pos 
p=plot(T_eV,plotThis[0,1,*],layout=[[layout],pos],title='plotThis[0,1]',/current,/xlog)
p=plot(T_eV,imaginary(plotThis[0,1,*]),color='r',/over)
p=plot(T_eV,plotThis_cold[0,1,*],/over,thick=thick,transparency=transparency,LineStyle=style)
p=plot(T_eV,imaginary(plotThis_cold[0,1,*]),color='r',/over,thick=thick,transparency=transparency,LineStyle=style)

p=plot(T_eV_kj, sig3, /over, thick=2)
p=plot(T_eV_kj, imaginary(sig3), color='r', /over, thick=2)

++pos 
p=plot(T_eV,plotThis[0,2,*],layout=[[layout],pos],title='plotThis[0,2]',/current,/xlog)
p=plot(T_eV,imaginary(plotThis[0,2,*]),color='r',/over)
p=plot(T_eV,plotThis_cold[0,2,*],/over,thick=thick,transparency=transparency,LineStyle=style)
p=plot(T_eV,imaginary(plotThis_cold[0,2,*]),color='r',/over,thick=thick,transparency=transparency,LineStyle=style)

p=plot(T_eV_kj, sig2, /over, thick=2)
p=plot(T_eV_kj, imaginary(sig2), color='r', /over, thick=2)

++pos
p=plot(T_eV,plotThis[1,0,*],layout=[[layout],pos],title='plotThis[1,0]',/current,/xlog)
p=plot(T_eV,imaginary(plotThis[1,0,*]),color='r',/over)
p=plot(T_eV,plotThis_cold[1,0,*],/over,thick=thick,transparency=transparency,LineStyle=style)
p=plot(T_eV,imaginary(plotThis_cold[1,0,*]),color='r',/over,thick=thick,transparency=transparency,LineStyle=style)

++pos 
p=plot(T_eV,plotThis[1,1,*],layout=[[layout],pos],title='plotThis[1,1]',/current,/xlog)
p=plot(T_eV,imaginary(plotThis[1,1,*]),color='r',/over)
p=plot(T_eV,plotThis_cold[1,1,*],/over,thick=thick,transparency=transparency,LineStyle=style)
p=plot(T_eV,imaginary(plotThis_cold[1,1,*]),color='r',/over,thick=thick,transparency=transparency,LineStyle=style)

++pos 
p=plot(T_eV,plotThis[1,2,*],layout=[[layout],pos],title='plotThis[1,2]',/current,/xlog)
p=plot(T_eV,imaginary(plotThis[1,2,*]),color='r',/over)
p=plot(T_eV,plotThis_cold[1,2,*],/over,thick=thick,transparency=transparency,LineStyle=style)
p=plot(T_eV,imaginary(plotThis_cold[1,2,*]),color='r',/over,thick=thick,transparency=transparency,LineStyle=style)

++pos
p=plot(T_eV,plotThis[2,0,*],layout=[[layout],pos],title='plotThis[2,0]',/current,/xlog)
p=plot(T_eV,imaginary(plotThis[2,0,*]),color='r',/over)
p=plot(T_eV,plotThis_cold[2,0,*],/over,thick=thick,transparency=transparency,LineStyle=style)
p=plot(T_eV,imaginary(plotThis_cold[2,0,*]),color='r',/over,thick=thick,transparency=transparency,LineStyle=style)

++pos 
p=plot(T_eV,plotThis[2,1,*],layout=[[layout],pos],title='plotThis[2,1]',/current,/xlog)
p=plot(T_eV,imaginary(plotThis[2,1,*]),color='r',/over)
p=plot(T_eV,plotThis_cold[2,1,*],/over,thick=thick,transparency=transparency,LineStyle=style)
p=plot(T_eV,imaginary(plotThis_cold[2,1,*]),color='r',/over,thick=thick,transparency=transparency,LineStyle=style)

++pos 
p=plot(T_eV,plotThis[2,2,*],layout=[[layout],pos],title='plotThis[2,2]',/current,/xlog)
p=plot(T_eV,imaginary(plotThis[2,2,*]),color='r',/over)
p=plot(T_eV,plotThis_cold[2,2,*],/over,thick=thick,transparency=transparency,LineStyle=style)
p=plot(T_eV,imaginary(plotThis_cold[2,2,*]),color='r',/over,thick=thick,transparency=transparency,LineStyle=style)


stop

end
