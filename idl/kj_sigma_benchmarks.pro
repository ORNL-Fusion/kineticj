pro kj_sigma_benchmarks, runKJ=runKJ, $
        benchmark = _benchmark

if keyword_set(_benchmark) then benchmark = _benchmark else benchmark = 1

@constants

n = 300
n_kj = 10 

kj_nPts_grid = 301
kj_nPts_eval = 1

if benchmark eq 1 then begin

    ; Benchmark 1
    ; -----------
    ; T Scan, D, ky = 0

    f = 13d6
    Z = 1d0
    amu =  2.0;_me_amu
    BUnit = [0,0,1]
    density = 2d19
    harmonicNumber = 4
  
    B_T = [1d0]
    B_T_kj = B_T

    ; Analytic calculation range

    tMin = 0.1
    tMax = 10e3
    T_eV = 10d0^(findGen(n)/(n-1)*(alog10(tMax)-alog10(tMin))+alog10(tMin)) 

    ; KJ calculation range

    tMin_kj = 100 
    tMax_kj = 10e3
    T_eV_kj = 10d0^(findGen(n_kj)/(n_kj-1)*(alog10(tMax_kj)-alog10(tMin_kj))+alog10(tMin_kj)) 
   
    kx = 10.0
    ky = 0.0 
    kz = 100.0

    ; KJ config parameters

    kj_nPx = 21
    kj_nPy = 21
    kj_nPz = 65
    kj_nStepsPerCyclotronPeriod = 100.0
    kj_nRFCycles = 10.0 

endif else if benchmark eq 2 then begin

    ; Benchmark 2
    ; -----------
    ; Te Scan, D, ky != 0

    f = 13d6
    Z = 1d0
    amu =  2.0;_me_amu
    BUnit = [0,0,1]
    density = 2d19
    harmonicNumber = 4

    B_T = [1d0]
    B_T_kj = B_T
   
    ; Analytic calculation range

    tMin = 0.1
    tMax = 10e3
    T_eV = 10d0^(findGen(n)/(n-1)*(alog10(tMax)-alog10(tMin))+alog10(tMin)) 
   
    ; KJ calculation range

    tMin_kj = 100 
    tMax_kj = 10e3
    T_eV_kj = 10d0^(findGen(n_kj)/(n_kj-1)*(alog10(tMax_kj)-alog10(tMin_kj))+alog10(tMin_kj)) 
 
    kx = 10.0
    ky = 23.0 
    kz = 100.0

    ; KJ config parameters

    kj_nPx = 21
    kj_nPy = 21
    kj_nPz = 65
    kj_nStepsPerCyclotronPeriod = 100.0
    kj_nRFCycles = 10.0 

endif else if benchmark eq 3 then begin

    ; Benchmark 3
    ; -----------
    ; B Scan over a few electron-cyclotron resonances
    
    f = 32d9
    Z = -1d0
    amu =  _me_amu
    B = 1d0
    BUnit = [0,0,1]
    density = 5d19
    harmonicNumber = 6
    
    T_eV = [1e3] 
    T_eV_kj = T_eV

    ; Analytic calculation range

    bMin = 0.9
    bMax = 1.5
    B_T = fIndGen(n)/(n-1)*(bMax-bMin)+bMin 

    ; KJ calculation range

    bMin_kj = 1.0
    bMax_kj = 1.3
    B_T_kj = fIndGen(n_kj)/(n_kj-1)*(bMax_kj-bMin_kj)+bMin_kj 

    kx = 10.0
    ky = 0.0 
    kz = 100.0

    ; KJ config parameters

    kj_nPx = 11
    kj_nPy = 11
    kj_nPz = 15
    kj_nStepsPerCyclotronPeriod = 100.0
    kj_nRFCycles = 100.0 

    ; Diagnose the scenario

    beta_ = density * _kB * T_eV / ( B_T^2 / ( 2 *_u0 ) )
    vTh = sqrt( 2.0 * T_eV * _e / ( amu * _amu ) )
    wc = ( Z * _e ) * B_T / ( amu * _amu )
    wp = sqrt ( density * _e^2 / ( amu * _amu * _e0 ) )
    w_wc = 2*!pi*f / wc

endif else if benchmark eq 4 then begin

    ; Benchmark 4
    ; -----------
    ; B Scan over a few ion-cyclotron resonances
    ; Within a single spatial domain.
    
    f = 32d9
    Z = -1d0
    amu =  _me_amu
    B = 1d0
    BUnit = [0,0,1]
    density = 5d19
    harmonicNumber = 6
    
    T_eV = [1.0e3] 
    T_eV_kj = T_eV

    kj_nPts_eval = n_kj 

    ; Analytic calculation range

    bMin = 0.9
    bMax = 1.5
    B_T = fIndGen(n)/(n-1)*(bMax-bMin)+bMin 

    ; KJ calculation range

    bMin_kj = 0.9
    bMax_kj = 1.5
    B_T_kj = fIndGen(kj_nPts_grid)/(kj_nPts_grid-1)*(bMax_kj-bMin_kj)+bMin_kj 

    kx = 10.0
    ky = 0.0 
    kz = 100.0

    ; KJ config parameters

    kj_nPx = 11
    kj_nPy = 11
    kj_nPz = 11
    kj_nStepsPerCyclotronPeriod = 100.0
    kj_nRFCycles = 100.0 

    ; Diagnose the scenario

    beta_ = density * _kB * T_eV / ( B_T^2 / ( 2 *_u0 ) )
    vTh = sqrt( 2.0 * T_eV * _e / ( amu * _amu ) )
    wc = ( Z * _e ) * B_T / ( amu * _amu )
    wc_kj = ( Z * _e ) * B_T_kj / ( amu * _amu )
    wp = sqrt ( density * _e^2 / ( amu * _amu * _e0 ) )
    w_wc = 2*!pi*f / wc
    w_wc_kj = 2*!pi*f / wc_kj

endif else if benchmark eq 5 then begin

    ; Benchmark 3
    ; -----------
    ; B Scan over a few electron-cyclotron resonances
    
    f = 32d9
    Z = -1d0
    amu =  _me_amu
    B = 1d0
    BUnit = [0,0,1]
    density = 5d19
    harmonicNumber = 6
    
    T_eV = [1e3] 
    T_eV_kj = T_eV

    ; Analytic calculation range

    bMin = 0.3
    bMax = 1.2
    B_T = fIndGen(n)/(n-1)*(bMax-bMin)+bMin 

    ; KJ calculation range

    bMin_kj = 0.3
    bMax_kj = 1.2
    B_T_kj = fIndGen(n_kj)/(n_kj-1)*(bMax_kj-bMin_kj)+bMin_kj 

    kx = 10.0
    ky = 0.0 
    kz = 100.0

    ; KJ config parameters

    kj_nPx = 11
    kj_nPy = 11
    kj_nPz = 15
    kj_nStepsPerCyclotronPeriod = 100.0
    kj_nRFCycles = 1.0 

    ; Diagnose the scenario

    beta_ = density * _kB * T_eV / ( B_T^2 / ( 2 *_u0 ) )
    vTh = sqrt( 2.0 * T_eV * _e / ( amu * _amu ) )
    wc = ( Z * _e ) * B_T / ( amu * _amu )
    wp = sqrt ( density * _e^2 / ( amu * _amu * _e0 ) )
    w_wc = 2*!pi*f / wc


endif

nT = n_elements(T_eV)
nB = n_elements(B_T)

nT_kj = n_elements(T_eV_kj)
nB_kj = n_elements(B_T_kj)

if benchmark eq 4 then nB_kj = 1

kPar = kz  
kPer = sqrt( kx^2 + ky^2 ) 

eps = ComplexArr(3,3,n)
eps_swan = ComplexArr(3,3,n)
eps_cold = ComplexArr(3,3,n)

sig = ComplexArr(3,3,n)
sig_swan = ComplexArr(3,3,n)
sig_cold = ComplexArr(3,3,n)

w = 2*!pi*f
cnt = 0
print, 'Calculating analytic sigma ...'
for i=0,nT-1 do begin
for j=0,nB-1 do begin

    thisTeV = T_eV[i]
    thisB = B_T[j]
    this_eps = kj_hot_epsilon(f, amu, Z, thisB, density, $
            harmonicNumber, kPar, kPer, thisTeV, $
            epsilon_cold = this_eps_cold, kx = kx, $
            epsilon_swan_WD = this_eps_swan, $
            epsilon_swan_ND = this_eps_swan_ND  )
    
    eps[*,*,cnt] = this_eps 
    eps_swan[*,*,cnt] = this_eps_swan 
    eps_cold[*,*,cnt] = this_eps_cold 

    sig[*,*,cnt] = (eps[*,*,cnt] - identity(3)) * w * _e0 / _ii
    sig_swan[*,*,cnt] = (eps_swan[*,*,cnt] - identity(3)) * w * _e0 / _ii
    sig_cold[*,*,cnt] = (eps_cold[*,*,cnt] - identity(3)) * w * _e0 / _ii
    ++cnt

endfor
endfor
print, 'DONE'

; Stage and run kj over this range of temperatures 

if benchmark eq 4 then begin
    sig_kj = ComplexArr(3,3,kj_nPts_eval)
endif else begin
    sig_kj = ComplexArr(3,3,n_kj)
endelse

TemplateRunDir = 'template'
cd, current = RootDir
print, 'RootDir: ', RootDir

cnt = 0
for t=0,nT_kj-1 do begin
for b=0,nB_kj-1 do begin
    
    ThisRunDir = string(cnt,format='(i5.5)')

    if keyword_set(runKJ) then begin

        file_delete, ThisRunDir, /recursive, /allow_nonexistent

        file_copy, TemplateRunDir, ThisRunDir, /recursive

    endif

    print, ThisRunDir
    cd, ThisRunDir

    for row=0,2 do begin

        RowString = string(row,format='(i1.1)')
        this_inputFileName = 'input/input-data_' + RowString + '.nc'
        This_jP2_FileName = 'output/jP2_' + RowString + '.nc'

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

        if benchmark eq 4 then begin
            this_b0 = B_T_kj
        endif else begin
            this_b0 = B_T_kj[b]
        endelse

        kj_create_single_k_input, b0=this_b0, bUnit=bUnit, kx=kx, f_Hz=f, n_m3=density, $
                Er=Er, Et=Et, Ez=Ez, x=x_kjGrid, writeOutput=runKJ, $
                E1Multiplier=E1, E2Multiplier=E2, E3Multiplier=E3, fileName=this_inputFileName, $
                nPts = kj_nPts_grid

        if keyword_set(runKJ) then begin

            ; Stage the input wave fields

            ; Adjust the kj.cfg config file parameters

            kj = kj_read_cfg('./')
            kj['input_fName'] = this_inputFileName
            if benchmark eq 4 then begin
                kj['xGridMin'] = 0.5 
                kj['xGridMax'] = 2.5
            endif else begin
                kj['xGridMin'] = x_kjGrid[0]+(x_kjGrid[-1]-x_kjGrid[0])/2 - (x_kjGrid[1]-x_kjGrid[0])
                kj['xGridMax'] = x_kjGrid[0]+(x_kjGrid[-1]-x_kjGrid[0])/2 + (x_kjGrid[1]-x_kjGrid[0])
            endelse
            kj['T_keV'] = T_eV_kj[t]*1e-3 
            kj['species_amu'] = double(amu)
            kj['species_Z'] = double(Z)
            kj['ky'] = float(ky)
            kj['kz'] = float(kz) 
            kj['nStepsPerCyclotronPeriod'] = float(kj_nStepsPerCyclotronPeriod) 
            kj['nRFCycles'] = float(kj_nRFCycles)
            kj['nP_Vx'] = fix(kj_nPx) 
            kj['nP_Vy'] = fix(kj_nPy) 
            kj['nP_Vz'] = fix(kj_nPz) 
            kj['nXGrid'] = fix(kj_nPts_eval) 

            ; Set dt (kj_nStepsPerCyclotronPeriod) such that we sample the shortest wavelength
            ; at the highest velocity with adequote sampling

            nvTh = 3
            vThMax = nvTh * sqrt( 2.0 * T_eV_kj[t] * _e / ( amu * _amu ) )
            parSamples = ( 2 * !pi / kPar ) / ( 1 / f / kj_nStepsPerCyclotronPeriod * vThMax )
            print, 'T_eV: ', T_eV_kj[t], '  parSamples: ', parSamples
            perSamples = ( 2 * !pi / kPer ) / ( 1 / f / kj_nStepsPerCyclotronPeriod * vThMax )
            print, 'T_eV: ', T_eV_kj[t], '  perSamples: ', perSamples

            nMinSamples = 5
            par_dt = ( 2 * !pi / kPar ) / nMinSamples / vThMax 
            nHarmonic = 3
            wc = ( Z * _e ) * B_T_kj[b] / ( amu * _amu )
            per_dt = 1 / ( wc / 2 * !pi ) / ( nHarmonic * 2 )

            par_nStepsPerCyclotronPeriod = 1 / ( wc / 2 * !pi ) / par_dt 
            per_nStepsPerCyclotronPeriod = 1 / ( wc / 2 * !pi ) / per_dt 

            print, 'par_nStepsPerCyclotronPeriod: ', par_nStepsPerCyclotronPeriod
            print, 'per_nStepsPerCyclotronPeriod: ', per_nStepsPerCyclotronPeriod

            vPhsPer = w / kPer
            vPhsPar = w / kPar

            kj_write_kj_cfg, kj, './'

            ; Run kj

            RunCommand = '~/code/kineticj/bin/kineticj'
            spawn, RunCommand, StdOut, StdErr
            print, StdOut
            print, StdErr

            file_move, 'output/jP2.nc', This_jP2_FileName 

        endif

        ; Read in kj results

        kj_read_jp_old, x=kj_x, j1x=kj_j1x, j1y=kj_j1y, j1z=kj_j1z, fileName=This_jP2_FileName
   
        kj_Er = interpol(Er,x_kjGrid,kj_x)
        kj_Et = interpol(Et,x_kjGrid,kj_x)
        kj_Ez = interpol(Ez,x_kjGrid,kj_x)

        if benchmark eq 4 then begin

            if row eq 0 then begin 
                sig_kj[row,0,*] = (kj_j1x/kj_Er)
                sig_kj[row,1,*] = (kj_j1y/kj_Er)
                sig_kj[row,2,*] = (kj_j1z/kj_Er)
            endif

            if row eq 1 then begin 
                sig_kj[row,0,*] = (kj_j1x/kj_Et)
                sig_kj[row,1,*] = (kj_j1y/kj_Et)
                sig_kj[row,2,*] = (kj_j1z/kj_Et)
            endif

            if row eq 2 then begin 
                sig_kj[row,0,*] = (kj_j1x/kj_Ez)
                sig_kj[row,1,*] = (kj_j1y/kj_Ez)
                sig_kj[row,2,*] = (kj_j1z/kj_Ez)
            endif

        endif else begin   

            if row eq 0 then begin 
                sig_kj[row,0,cnt] = (kj_j1x/kj_Er)[0]
                sig_kj[row,1,cnt] = (kj_j1y/kj_Er)[0]
                sig_kj[row,2,cnt] = (kj_j1z/kj_Er)[0]
            endif

            if row eq 1 then begin 
                sig_kj[row,0,cnt] = (kj_j1x/kj_Et)[0]
                sig_kj[row,1,cnt] = (kj_j1y/kj_Et)[0]
                sig_kj[row,2,cnt] = (kj_j1z/kj_Et)[0]
            endif

            if row eq 2 then begin 
                sig_kj[row,0,cnt] = (kj_j1x/kj_Ez)[0]
                sig_kj[row,1,cnt] = (kj_j1y/kj_Ez)[0]
                sig_kj[row,2,cnt] = (kj_j1z/kj_Ez)[0]
            endif

        endelse

    endfor

    cd, RootDir

    ++cnt

endfor
endfor

for i=0,n_elements(sig_kj[0,0,*])-1 do begin

    sig_kj[*,*,i] = transpose(sig_kj[*,*,i])
  
endfor

; Plot results

layout=[3,3]
pos = 1
thick = 1 
style = '--'
transparency = 30
transparencyC = 0
xFS = 6
yFS = 6
margin = [0.18,0.15,0.1,0.15]
xRange = 0
kj_color1 = 'teal'
kj_color2 = 'tomato'

plotThis = sig
plotThis_cold = sig_cold

if benchmark eq 1 or benchmark eq 2 then begin
    xTitle ='log10( T [eV] )'
    x = alog10(T_eV)
    x_kj = alog10(T_eV_kj)
endif

if benchmark eq 3 then begin
    xTitle ='$\omega/\omega_{c}$'

    wc = abs(( Z * _e ) * B_T / ( amu * _amu ))
    wc_kj = abs(( Z * _e ) * B_T_kj / ( amu * _amu ))

    x = w / wc
    x_kj = w / wc_kj

    xRange = [0.8,1.2]

endif

if benchmark eq 4 then begin
    xTitle ='$\omega/\omega_{c}$'

    wc = abs(( Z * _e ) * B_T / ( amu * _amu ))

    kj_B = interpol(B_T_kj,x_kjGrid,kj_x)

    wc_kj = abs(( Z * _e ) * kj_B / ( amu * _amu ))

    x = w / wc
    x_kj = w / wc_kj

endif

if benchmark eq 5 then begin
    xTitle ='$\omega/\omega_{c}$'

    wc = abs(( Z * _e ) * B_T / ( amu * _amu ))
    wc_kj = abs(( Z * _e ) * B_T_kj / ( amu * _amu ))

    x = w / wc
    x_kj = w / wc_kj

    xRange = 0 

endif

yRange = [-1,1]*max(abs(sig[0,0,*]))*1.1

p=plot(x,plotThis[0,0,*],layout=[[layout],pos],$
        title='$\sigma_{xx}$',/buffer,$
        font_size=12, xTitle=xTitle, xTickFont_size=xFS, yTickFont_size=yFS, $
        xMinor = 0, axis_style=1, yTitle='$\sigma_{xx} [S/m]$', margin=margin, yRange=yRange, xRange=xRange )
p=plot(x,imaginary(plotThis[0,0,*]),color='r',/over)
p=plot(x,plotThis_cold[0,0,*],/over,thick=thick,transparency=transparencyC,LineStyle=style)
p=plot(x,imaginary(plotThis_cold[0,0,*]),color='r',/over,thick=thick,transparency=transparencyC,LineStyle=style)

p=plot(x_kj, sig_kj[0,0,*], /over, thick=3, color=kj_color1, transparency=transparency)
p=plot(x_kj, imaginary(sig_kj[0,0,*]), color=kj_color2, /over, thick=3, transparency=transparency)

p=plot(x, sig_swan[0,0,*], /over, thick=1, lineStyle='--')
p=plot(x, imaginary(sig_swan[0,0,*]), color='r', /over, thick=1, lineStyle='--')

++pos 
p=plot(x,plotThis[0,1,*],layout=[[layout],pos],/current, $
        title='$\sigma_{xy}$',yRange=[-1,1]*max(abs(plotThis[0,1,*])),/buffer,$
        font_size=12, xTitle=xTitle, xTickFont_size=xFS, yTickFont_size=yFS, $
        xMinor = 0, axis_style=1, yTitle='$\sigma_{xy} [S/m]$', margin=margin, xRange=xRange )
p=plot(x,imaginary(plotThis[0,1,*]),color='r',/over)
p=plot(x,plotThis_cold[0,1,*],/over,thick=thick,transparency=transparencyC,LineStyle=style)
p=plot(x,imaginary(plotThis_cold[0,1,*]),color='r',/over,thick=thick,transparency=transparencyC,LineStyle=style)

p=plot(x_kj, sig_kj[0,1,*], /over, thick=3, color=kj_color1, transparency=transparency)
p=plot(x_kj, imaginary(sig_kj[0,1,*]), color=kj_color2, /over, thick=3,transparency=transparency)

p=plot(x, sig_swan[0,1,*], /over, thick=1, lineStyle='--')
p=plot(x, imaginary(sig_swan[0,1,*]), color='r', /over, thick=1, lineStyle='--')

++pos 
p=plot(x,plotThis[0,2,*],layout=[[layout],pos],/current, $
        title='$\sigma_{xz}$',yRange=[-1,1]*max(abs(plotThis[0,2,*])),/buffer,$
        font_size=12, xTitle=xTitle, xTickFont_size=xFS, yTickFont_size=yFS, $
        xMinor = 0, axis_style=1, yTitle='$\sigma_{xz} [S/m]$', margin=margin, xRange=xRange )

p=plot(x,imaginary(plotThis[0,2,*]),color='r',/over)
p=plot(x,plotThis_cold[0,2,*],/over,thick=thick,transparency=transparencyC,LineStyle=style)
p=plot(x,imaginary(plotThis_cold[0,2,*]),color='r',/over,thick=thick,transparency=transparencyC,LineStyle=style)

p=plot(x_kj, sig_kj[0,2,*], /over, thick=3, color=kj_color1, transparency=transparency)
p=plot(x_kj, imaginary(sig_kj[0,2,*]), color=kj_color2, /over, thick=3,transparency=transparency)

p=plot(x, sig_swan[0,2,*], /over, thick=1, lineStyle='--')
p=plot(x, imaginary(sig_swan[0,2,*]), color='r', /over, thick=1, lineStyle='--')

++pos
p=plot(x,plotThis[1,0,*],layout=[[layout],pos],/current, $
        title='$\sigma_{yx}$',yRange=[-1,1]*max(abs(plotThis[1,0,*])),/buffer,$
        font_size=12, xTitle=xTitle, xTickFont_size=xFS, yTickFont_size=yFS, $
        xMinor = 0, axis_style=1, yTitle='$\sigma_{yx} [S/m]$', margin=margin, xRange=xRange )

p=plot(x,imaginary(plotThis[1,0,*]),color='r',/over)
p=plot(x,plotThis_cold[1,0,*],/over,thick=thick,transparency=transparencyC,LineStyle=style)
p=plot(x,imaginary(plotThis_cold[1,0,*]),color='r',/over,thick=thick,transparency=transparencyC,LineStyle=style)

p=plot(x_kj, sig_kj[1,0,*], /over, thick=3, color=kj_color1, transparency=transparency)
p=plot(x_kj, imaginary(sig_kj[1,0,*]), color=kj_color2, /over, thick=3,transparency=transparency)

p=plot(x, sig_swan[1,0,*], /over, thick=1, lineStyle='--')
p=plot(x, imaginary(sig_swan[1,0,*]), color='r', /over, thick=1, lineStyle='--')

++pos 
p=plot(x,plotThis[1,1,*],layout=[[layout],pos],/current, $
        title='$\sigma_{yy}$',yRange=[-1,1]*max(abs(plotThis[1,1,*])),/buffer,$
        font_size=12, xTitle=xTitle, xTickFont_size=xFS, yTickFont_size=yFS, $
        xMinor = 0, axis_style=1, yTitle='$\sigma_{yy} [S/m]$', margin=margin, xRange=xRange )

p=plot(x,imaginary(plotThis[1,1,*]),color='r',/over)
p=plot(x,plotThis_cold[1,1,*],/over,thick=thick,transparency=transparencyC,LineStyle=style)
p=plot(x,imaginary(plotThis_cold[1,1,*]),color='r',/over,thick=thick,transparency=transparencyC,LineStyle=style)

p=plot(x_kj, sig_kj[1,1,*], /over, thick=3, color=kj_color1, transparency=transparency)
p=plot(x_kj, imaginary(sig_kj[1,1,*]), color=kj_color2, /over, thick=3,transparency=transparency)

p=plot(x, sig_swan[1,1,*], /over, thick=1, lineStyle='--')
p=plot(x, imaginary(sig_swan[1,1,*]), color='r', /over, thick=1, lineStyle='--')

++pos 
p=plot(x,plotThis[1,2,*],layout=[[layout],pos],/current, $
        title='$\sigma_{yz}$',yRange=[-1,1]*max(abs(plotThis[1,2,*])),/buffer,$
        font_size=12, xTitle=xTitle, xTickFont_size=xFS, yTickFont_size=yFS, $
        xMinor = 0, axis_style=1, yTitle='$\sigma_{yz} [S/m]$', margin=margin, xRange=xRange )

p=plot(x,imaginary(plotThis[1,2,*]),color='r',/over)
p=plot(x,plotThis_cold[1,2,*],/over,thick=thick,transparency=transparencyC,LineStyle=style)
p=plot(x,imaginary(plotThis_cold[1,2,*]),color='r',/over,thick=thick,transparency=transparencyC,LineStyle=style)

p=plot(x_kj, sig_kj[1,2,*], /over, thick=3, color=kj_color1, transparency=transparency)
p=plot(x_kj, imaginary(sig_kj[1,2,*]), color=kj_color2, /over, thick=3,transparency=transparency)

p=plot(x, sig_swan[1,2,*], /over, thick=1, lineStyle='--')
p=plot(x, imaginary(sig_swan[1,2,*]), color='r', /over, thick=1, lineStyle='--')

++pos
p=plot(x,plotThis[2,0,*],layout=[[layout],pos],/current, $
        title='$\sigma_{zx}$',yRange=[-1,1]*max(abs(plotThis[2,0,*])),/buffer,$
        font_size=12, xTitle=xTitle, xTickFont_size=xFS, yTickFont_size=yFS, $
        xMinor = 0, axis_style=1, yTitle='$\sigma_{zx} [S/m]$', margin=margin, xRange=xRange )

p=plot(x,imaginary(plotThis[2,0,*]),color='r',/over)
p=plot(x,plotThis_cold[2,0,*],/over,thick=thick,transparency=transparencyC,LineStyle=style)
p=plot(x,imaginary(plotThis_cold[2,0,*]),color='r',/over,thick=thick,transparency=transparencyC,LineStyle=style)

p=plot(x_kj, sig_kj[2,0,*], /over, thick=3, color=kj_color1, transparency=transparency)
p=plot(x_kj, imaginary(sig_kj[2,0,*]), color=kj_color2, /over, thick=3,transparency=transparency)

p=plot(x, sig_swan[2,0,*], /over, thick=1, lineStyle='--')
p=plot(x, imaginary(sig_swan[2,0,*]), color='r', /over, thick=1, lineStyle='--')

++pos 
p=plot(x,plotThis[2,1,*],layout=[[layout],pos],/current, $
        title='$\sigma_{zy}$',yRange=[-1,1]*max(abs(plotThis[2,1,*])),/buffer,$
        font_size=12, xTitle=xTitle, xTickFont_size=xFS, yTickFont_size=yFS, $
        xMinor = 0, axis_style=1, yTitle='$\sigma_{zy} [S/m]$', margin=margin, xRange=xRange )

p=plot(x,imaginary(plotThis[2,1,*]),color='r',/over)
p=plot(x,plotThis_cold[2,1,*],/over,thick=thick,transparency=transparencyC,LineStyle=style)
p=plot(x,imaginary(plotThis_cold[2,1,*]),color='r',/over,thick=thick,transparency=transparencyC,LineStyle=style)

p=plot(x_kj, sig_kj[2,1,*], /over, thick=3, color=kj_color1, transparency=transparency)
p=plot(x_kj, imaginary(sig_kj[2,1,*]), color=kj_color2, /over, thick=3,transparency=transparency)

p=plot(x, sig_swan[2,1,*], /over, thick=1, lineStyle='--')
p=plot(x, imaginary(sig_swan[2,1,*]), color='r', /over, thick=1, lineStyle='--')

++pos 

yRange = [-1,1]*max(abs(sig[2,2,*]))*1.1

p=plot(x,plotThis[2,2,*],layout=[[layout],pos],/current, $
        title='$\sigma_{zz}$',yRange=yRange,/buffer,$
        font_size=12, xTitle=xTitle, xTickFont_size=xFS, yTickFont_size=yFS, $
        xMinor = 0, axis_style=1, yTitle='$\sigma_{zz} [S/m]$', margin=margin, xRange=xRange )

p=plot(x,imaginary(plotThis[2,2,*]),color='r',/over)
p=plot(x,plotThis_cold[2,2,*],/over,thick=thick,transparency=transparencyC,LineStyle=style)
p=plot(x,imaginary(plotThis_cold[2,2,*]),color='r',/over,thick=thick,transparency=transparencyC,LineStyle=style)

p=plot(x_kj, sig_kj[2,2,*], /over, thick=3, color=kj_color1, transparency=transparency)
p=plot(x_kj, imaginary(sig_kj[2,2,*]), color=kj_color2, /over, thick=3,transparency=transparency)

p=plot(x, sig_swan[2,2,*], /over, thick=1, lineStyle='--')
p=plot(x, imaginary(sig_swan[2,2,*]), color='r', /over, thick=1, lineStyle='--')

p.save, 'kj_sigma_vs_t.png', resolution=300, /transparent
p.save, 'kj_sigma_vs_t.pdf'

stop

end
