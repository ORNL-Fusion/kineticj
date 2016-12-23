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
BUnit = [0,0,1]
density = 2d19
harmonicNumber = 1
kPar = 100.0  
kPer = 10.0 

kx = kPer 
ky = 0
kz = kPar

kj_nStepsPerCycle = 100.0
kj_nRFCycles = 10.0 

eps = ComplexArr(3,3,n)
eps_swan = ComplexArr(3,3,n)
eps_cold = ComplexArr(3,3,n)

sig = ComplexArr(3,3,n)
sig_swan = ComplexArr(3,3,n)
sig_cold = ComplexArr(3,3,n)

w = 2*!pi*f
for i=0,n-1 do begin

    thisTeV = T_eV[i]
    this_eps = kj_hot_epsilon(f, amu, Z, B, density, $
            harmonicNumber, kPar, kPer, thisTeV, $
            epsilon_cold = this_eps_cold, kx = kx, epsilon_swan = this_eps_swan )
    
    eps[*,*,i] = this_eps 
    eps_swan[*,*,i] = this_eps_swan 
    eps_cold[*,*,i] = this_eps_cold 

    sig[*,*,i] = (eps[*,*,i] - identity(3)) * w * _e0 / _ii
    sig_swan[*,*,i] = (eps_swan[*,*,i] - identity(3)) * w * _e0 / _ii
    sig_cold[*,*,i] = (eps_cold[*,*,i] - identity(3)) * w * _e0 / _ii

endfor

; Stage and run kj over this range of temperatures 

n_kj = 10
tMin = 100 
tMax = 10e3
T_eV_kj = 10d0^(findGen(n_kj)/(n_kj-1)*(alog10(tMax)-alog10(tMin))+alog10(tMin)) 
;T_eV_kj = [0.1,1.0,10.0,100.0,1000.0, 10000.0]

nT_kj = n_elements(T_eV_kj)

sig_kj = ComplexArr(3,3,nT_kj)

TemplateRunDir = 'benchmark-perp'
RootDir = '/home/dg6/scratch/kineticJ/temp_scan'
cd, RootDir

for t=0,nT_kj-1 do begin
    
    ThisRunDir = string(t,format='(i5.5)')

    if keyword_set(runKJ) then begin

        file_delete, ThisRunDir, /recursive, /allow_nonexistent

        file_copy, TemplateRunDir, ThisRunDir, /recursive

    endif

    cd, ThisRunDir

    for row=0,2 do begin

        RowString = string(row,format='(i1.1)')
        This_E_FileName = 'data/kj_single_k_' + RowString
        This_jP2_FileName = 'jP2_' + RowString + '.nc'

        if row eq 0 then begin
            E1=1
            E2=0
            E3=0
        endif

        if row eq 1 then begin
            E1=0
            E2=1
            E3=0
        endif

        if row eq 2 then begin
            E1=0
            E2=0
            E3=1
        endif

        kj_create_single_k_input, b0=B, bUnit=bUnit, kx=kx, f_Hz=f, n_m3=density, $
                Er=Er, Et=Et, Ez=Ez, x=x, writeOutput=runKJ, $
                E1Multiplier=E1, E2Multiplier=E2, E3Multiplier=E3, fileName=This_E_FileName

        if keyword_set(runKJ) then begin

            ; Stage the input wave fields

            ; Adjust the kj.cfg config file parameters

            kj = kj_read_cfg('./')
            kj['eField_fName'] = this_E_FileName
            kj['xGridMin'] = x[0]+(x[-1]-x[0])/2 - (x[1]-x[0])
            kj['xGridMax'] = x[0]+(x[-1]-x[0])/2 + (x[1]-x[0])
            kj['T_keV'] = T_eV_kj[t]*1e-3 
            kj['species_amu'] = float(amu)
            kj['species_Z'] = float(Z)
            kj['ky'] = float(ky)
            kj['kz'] = float(kz) 
            kj['nStepsPerCycle'] = float(kj_nStepsPerCycle) 
            kj['nRFCycles'] = float(kj_nRFCycles)
            kj['nPx'] = 21 
            kj['nPy'] = 21 
            kj['nPz'] = 65 

            ; Set dt (kj_nStepsPerCycle) such that we sample the shortest wavelength
            ; at the highest velocity with adequote sampling

            nvTh = 3
            vThMax = nvTh * sqrt( 2.0 * T_eV_kj[t] * _e / ( amu * _amu ) )
            parSamples = ( 2 * !pi / kPar ) / ( 1 / f / kj_nStepsPerCycle * vThMax )
            print, 'T_eV: ', T_eV_kj[t], '  parSamples: ', parSamples
            perSamples = ( 2 * !pi / kPer ) / ( 1 / f / kj_nStepsPerCycle * vThMax )
            print, 'T_eV: ', T_eV_kj[t], '  perSamples: ', perSamples

            nMinSamples = 5
            par_dt = ( 2 * !pi / kPar ) / nMinSamples / vThMax 
            nHarmonic = 3
            wc = ( Z * _e ) * B / ( amu * _amu )
            per_dt = 1 / ( wc / 2 * !pi ) / ( nHarmonic * 2 )

            par_nStepsPerCycle = 1 / ( wc / 2 * !pi ) / par_dt 
            per_nStepsPerCycle = 1 / ( wc / 2 * !pi ) / per_dt 

            print, 'par_nStepsPerCycle: ', par_nStepsPerCycle
            print, 'per_nStepsPerCycle: ', per_nStepsPerCycle

            vPhsPer = w / kPer
            vPhsPar = w / kPar

            kj_write_kj_cfg, kj, './'

            ; Run kj

            RunCommand = '~/code/kineticj/bin/kineticj'
            spawn, RunCommand, StdOut, StdErr
            print, StdOut
            print, StdErr

            file_move, 'jP2.nc', This_jP2_FileName 

        endif

        ; Read in kj results

        ;kj_read_jp_old, x=kj_x, j1x=kj_j1x, j1y=kj_j1y, j1z=kj_j1z, /oldFormat
        kj_read_jp_old, x=kj_x, j1x=kj_j1x, j1y=kj_j1y, j1z=kj_j1z, fileName=This_jP2_FileName
   
        kj_Er = interpol(Er,x,kj_x)
        kj_Et = interpol(Et,x,kj_x)
        kj_Ez = interpol(Ez,x,kj_x)
   
        if row eq 0 then begin 
            sig_kj[row,0,t] = (kj_j1x/kj_Er)[0]
            sig_kj[row,1,t] = (kj_j1y/kj_Er)[0]
            sig_kj[row,2,t] = (kj_j1z/kj_Er)[0]
        endif

        if row eq 1 then begin 
            sig_kj[row,0,t] = (kj_j1x/kj_Et)[0]
            sig_kj[row,1,t] = (kj_j1y/kj_Et)[0]
            sig_kj[row,2,t] = (kj_j1z/kj_Et)[0]
        endif

        if row eq 2 then begin 
            sig_kj[row,0,t] = (kj_j1x/kj_Ez)[0]
            sig_kj[row,1,t] = (kj_j1y/kj_Ez)[0]
            sig_kj[row,2,t] = (kj_j1z/kj_Ez)[0]
        endif

    endfor

    cd, RootDir

    sig_kj[*,*,t] = transpose(sig_kj[*,*,t])

endfor


layout=[3,3]
pos = 1
thick = 2 
style = '--'
transparency = 50
xFS = 6
yFS = 6
margin = [0.15,0.15,0.1,0.15]

plotThis = sig
plotThis_cold = sig_cold

p=plot(alog10(T_eV),plotThis[0,0,*],layout=[[layout],pos],$
        title='$\sigma_{xx}$',yRange=[-1,1]*max(abs(plotThis[0,0,*])),/buffer,$
        font_size=12, xTitle='log10( T [eV] )', xTickFont_size=xFS, yTickFont_size=yFS, $
        xMinor = 0, axis_style=1, yTitle='$\sigma_{xx} [S/m]$', margin=margin )
p=plot(alog10(T_eV),imaginary(plotThis[0,0,*]),color='r',/over)
p=plot(alog10(T_eV),plotThis_cold[0,0,*],/over,thick=thick,transparency=transparency,LineStyle=style)
p=plot(alog10(T_eV),imaginary(plotThis_cold[0,0,*]),color='r',/over,thick=thick,transparency=transparency,LineStyle=style)

p=plot(alog10(T_eV_kj), sig_kj[0,0,*], /over, thick=2)
p=plot(alog10(T_eV_kj), imaginary(sig_kj[0,0,*]), color='r', /over, thick=2)

p=plot(alog10(T_eV), sig_swan[0,0,*], /over, thick=1, color='m', lineStyle='--')
p=plot(alog10(T_eV), imaginary(sig_swan[0,0,*]), color='m', /over, thick=1, lineStyle='--')

++pos 
p=plot(alog10(T_eV),plotThis[0,1,*],layout=[[layout],pos],/current, $
        title='$\sigma_{xy}$',yRange=[-1,1]*max(abs(plotThis[0,1,*])),/buffer,$
        font_size=12, xTitle='log10( T [eV] )', xTickFont_size=xFS, yTickFont_size=yFS, $
        xMinor = 0, axis_style=1, yTitle='$\sigma_{xy} [S/m]$', margin=margin )
p=plot(alog10(T_eV),imaginary(plotThis[0,1,*]),color='r',/over)
p=plot(alog10(T_eV),plotThis_cold[0,1,*],/over,thick=thick,transparency=transparency,LineStyle=style)
p=plot(alog10(T_eV),imaginary(plotThis_cold[0,1,*]),color='r',/over,thick=thick,transparency=transparency,LineStyle=style)

p=plot(alog10(T_eV_kj), sig_kj[0,1,*], /over, thick=2)
p=plot(alog10(T_eV_kj), imaginary(sig_kj[0,1,*]), color='r', /over, thick=2)

p=plot(alog10(T_eV), sig_swan[0,1,*], /over, thick=1, color='m', lineStyle='--')
p=plot(alog10(T_eV), imaginary(sig_swan[0,1,*]), color='m', /over, thick=1, lineStyle='--')

++pos 
p=plot(alog10(T_eV),plotThis[0,2,*],layout=[[layout],pos],/current, $
        title='$\sigma_{xz}$',yRange=[-1,1]*max(abs(plotThis[0,2,*])),/buffer,$
        font_size=12, xTitle='log10( T [eV] )', xTickFont_size=xFS, yTickFont_size=yFS, $
        xMinor = 0, axis_style=1, yTitle='$\sigma_{xz} [S/m]$', margin=margin )

p=plot(alog10(T_eV),imaginary(plotThis[0,2,*]),color='r',/over)
p=plot(alog10(T_eV),plotThis_cold[0,2,*],/over,thick=thick,transparency=transparency,LineStyle=style)
p=plot(alog10(T_eV),imaginary(plotThis_cold[0,2,*]),color='r',/over,thick=thick,transparency=transparency,LineStyle=style)

p=plot(alog10(T_eV_kj), sig_kj[0,2,*], /over, thick=2)
p=plot(alog10(T_eV_kj), imaginary(sig_kj[0,2,*]), color='r', /over, thick=2)

p=plot(alog10(T_eV), sig_swan[0,2,*], /over, thick=1, color='m', lineStyle='--')
p=plot(alog10(T_eV), imaginary(sig_swan[0,2,*]), color='m', /over, thick=1, lineStyle='--')

++pos
p=plot(alog10(T_eV),plotThis[1,0,*],layout=[[layout],pos],/current, $
        title='$\sigma_{yx}$',yRange=[-1,1]*max(abs(plotThis[1,0,*])),/buffer,$
        font_size=12, xTitle='log10( T [eV] )', xTickFont_size=xFS, yTickFont_size=yFS, $
        xMinor = 0, axis_style=1, yTitle='$\sigma_{yx} [S/m]$', margin=margin )

p=plot(alog10(T_eV),imaginary(plotThis[1,0,*]),color='r',/over)
p=plot(alog10(T_eV),plotThis_cold[1,0,*],/over,thick=thick,transparency=transparency,LineStyle=style)
p=plot(alog10(T_eV),imaginary(plotThis_cold[1,0,*]),color='r',/over,thick=thick,transparency=transparency,LineStyle=style)

p=plot(alog10(T_eV_kj), sig_kj[1,0,*], /over, thick=2)
p=plot(alog10(T_eV_kj), imaginary(sig_kj[1,0,*]), color='r', /over, thick=2)

p=plot(alog10(T_eV), sig_swan[1,0,*], /over, thick=1, color='m', lineStyle='--')
p=plot(alog10(T_eV), imaginary(sig_swan[1,0,*]), color='m', /over, thick=1, lineStyle='--')

++pos 
p=plot(alog10(T_eV),plotThis[1,1,*],layout=[[layout],pos],/current, $
        title='$\sigma_{yy}$',yRange=[-1,1]*max(abs(plotThis[1,1,*])),/buffer,$
        font_size=12, xTitle='log10( T [eV] )', xTickFont_size=xFS, yTickFont_size=yFS, $
        xMinor = 0, axis_style=1, yTitle='$\sigma_{yy} [S/m]$', margin=margin )

p=plot(alog10(T_eV),imaginary(plotThis[1,1,*]),color='r',/over)
p=plot(alog10(T_eV),plotThis_cold[1,1,*],/over,thick=thick,transparency=transparency,LineStyle=style)
p=plot(alog10(T_eV),imaginary(plotThis_cold[1,1,*]),color='r',/over,thick=thick,transparency=transparency,LineStyle=style)

p=plot(alog10(T_eV_kj), sig_kj[1,1,*], /over, thick=2)
p=plot(alog10(T_eV_kj), imaginary(sig_kj[1,1,*]), color='r', /over, thick=2)

p=plot(alog10(T_eV), sig_swan[1,1,*], /over, thick=1, color='m', lineStyle='--')
p=plot(alog10(T_eV), imaginary(sig_swan[1,1,*]), color='m', /over, thick=1, lineStyle='--')

++pos 
p=plot(alog10(T_eV),plotThis[1,2,*],layout=[[layout],pos],/current, $
        title='$\sigma_{yz}$',yRange=[-1,1]*max(abs(plotThis[1,2,*])),/buffer,$
        font_size=12, xTitle='log10( T [eV] )', xTickFont_size=xFS, yTickFont_size=yFS, $
        xMinor = 0, axis_style=1, yTitle='$\sigma_{yz} [S/m]$', margin=margin )

p=plot(alog10(T_eV),imaginary(plotThis[1,2,*]),color='r',/over)
p=plot(alog10(T_eV),plotThis_cold[1,2,*],/over,thick=thick,transparency=transparency,LineStyle=style)
p=plot(alog10(T_eV),imaginary(plotThis_cold[1,2,*]),color='r',/over,thick=thick,transparency=transparency,LineStyle=style)

p=plot(alog10(T_eV_kj), sig_kj[1,2,*], /over, thick=2)
p=plot(alog10(T_eV_kj), imaginary(sig_kj[1,2,*]), color='r', /over, thick=2)

p=plot(alog10(T_eV), sig_swan[1,2,*], /over, thick=1, color='m', lineStyle='--')
p=plot(alog10(T_eV), imaginary(sig_swan[1,2,*]), color='m', /over, thick=1, lineStyle='--')

++pos
p=plot(alog10(T_eV),plotThis[2,0,*],layout=[[layout],pos],/current, $
        title='$\sigma_{zx}$',yRange=[-1,1]*max(abs(plotThis[2,0,*])),/buffer,$
        font_size=12, xTitle='log10( T [eV] )', xTickFont_size=xFS, yTickFont_size=yFS, $
        xMinor = 0, axis_style=1, yTitle='$\sigma_{zx} [S/m]$', margin=margin )

p=plot(alog10(T_eV),imaginary(plotThis[2,0,*]),color='r',/over)
p=plot(alog10(T_eV),plotThis_cold[2,0,*],/over,thick=thick,transparency=transparency,LineStyle=style)
p=plot(alog10(T_eV),imaginary(plotThis_cold[2,0,*]),color='r',/over,thick=thick,transparency=transparency,LineStyle=style)

p=plot(alog10(T_eV_kj), sig_kj[2,0,*], /over, thick=2)
p=plot(alog10(T_eV_kj), imaginary(sig_kj[2,0,*]), color='r', /over, thick=2)

p=plot(alog10(T_eV), sig_swan[2,0,*], /over, thick=1, color='m', lineStyle='--')
p=plot(alog10(T_eV), imaginary(sig_swan[2,0,*]), color='m', /over, thick=1, lineStyle='--')

++pos 
p=plot(alog10(T_eV),plotThis[2,1,*],layout=[[layout],pos],/current, $
        title='$\sigma_{zy}$',yRange=[-1,1]*max(abs(plotThis[2,1,*])),/buffer,$
        font_size=12, xTitle='log10( T [eV] )', xTickFont_size=xFS, yTickFont_size=yFS, $
        xMinor = 0, axis_style=1, yTitle='$\sigma_{zy} [S/m]$', margin=margin )

p=plot(alog10(T_eV),imaginary(plotThis[2,1,*]),color='r',/over)
p=plot(alog10(T_eV),plotThis_cold[2,1,*],/over,thick=thick,transparency=transparency,LineStyle=style)
p=plot(alog10(T_eV),imaginary(plotThis_cold[2,1,*]),color='r',/over,thick=thick,transparency=transparency,LineStyle=style)

p=plot(alog10(T_eV_kj), sig_kj[2,1,*], /over, thick=2)
p=plot(alog10(T_eV_kj), imaginary(sig_kj[2,1,*]), color='r', /over, thick=2)

p=plot(alog10(T_eV), sig_swan[2,1,*], /over, thick=1, color='m', lineStyle='--')
p=plot(alog10(T_eV), imaginary(sig_swan[2,1,*]), color='m', /over, thick=1, lineStyle='--')

++pos 
p=plot(alog10(T_eV),plotThis[2,2,*],layout=[[layout],pos],/current, $
        title='$\sigma_{zz}$',yRange=[-1,1]*max(abs(plotThis[2,2,*])),/buffer,$
        font_size=12, xTitle='log10( T [eV] )', xTickFont_size=xFS, yTickFont_size=yFS, $
        xMinor = 0, axis_style=1, yTitle='$\sigma_{zz} [S/m]$', margin=margin )

p=plot(alog10(T_eV),imaginary(plotThis[2,2,*]),color='r',/over)
p=plot(alog10(T_eV),plotThis_cold[2,2,*],/over,thick=thick,transparency=transparency,LineStyle=style)
p=plot(alog10(T_eV),imaginary(plotThis_cold[2,2,*]),color='r',/over,thick=thick,transparency=transparency,LineStyle=style)

p=plot(alog10(T_eV_kj), sig_kj[2,2,*], /over, thick=2)
p=plot(alog10(T_eV_kj), imaginary(sig_kj[2,2,*]), color='r', /over, thick=2)

p=plot(alog10(T_eV), sig_swan[2,2,*], /over, thick=1, color='m', lineStyle='--')
p=plot(alog10(T_eV), imaginary(sig_swan[2,2,*]), color='m', /over, thick=1, lineStyle='--')

p.save, 'kj_sigma_vs_t.png', resolution=300, /transparent
p.save, 'kj_sigma_vs_t.pdf'

stop

end
