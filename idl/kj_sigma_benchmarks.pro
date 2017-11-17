pro kj_sigma_benchmarks, runKJ=runKJ, $
        benchmark = _benchmark

if keyword_set(_benchmark) then benchmark = _benchmark else benchmark = 1

cd, current = pwd

benchmarkDirString = 'benchmark'+StrTrim(string(benchmark),2)

if not file_test(benchmarkDirString,/directory) then $
    file_mkdir, benchmarkDirString 

cd, benchmarkDirString

@constants

n = 300
n_kj = 15

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

    kj_nPx = 11
    kj_nPy = 11
    kj_nPz = 85
    kj_nStepsPerCyclotronPeriod = 600.0
    kj_nRFCycles = 10.0 

    ; Diagnose the scenario

    beta_ = density * _kB * T_eV / ( B_T^2 / ( 2 *_u0 ) )
    vTh = sqrt( 2.0 * T_eV * _e / ( amu * _amu ) )
    wc = ( Z * _e ) * B_T / ( amu * _amu )
    wp = sqrt ( density * _e^2 / ( amu * _amu * _e0 ) )
    w_wc = 2*!pi*f / wc

endif else if benchmark eq 2 then begin

    ; Benchmark 2
    ; -----------
    ; Te Scan, D, ky != 0

    f = 13d6
    Z = 1d0
    amu =  2.0
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

    kj_nPx = 11
    kj_nPy = 11
    kj_nPz = 85
    kj_nStepsPerCyclotronPeriod = 500.0
    kj_nRFCycles = 10.0 

endif else if benchmark eq 3 then begin

    ; Benchmark 3
    ; -----------
    ; B Scan over a few electron-cyclotron resonances

    n_kj = 50 

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
    ; Same as benchmark 3 but
    ; within a single spatial domain, 
    ; although there is a bug in the z 
    ; off-diagnoal elements not present in b3.
    
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

    ; Benchmark 5
    ; -----------
    ; B Scan over a few electron-cyclotron resonances
    ; but the 2nd and 3rd harmonic responses are at least 
    ; 5 orders of magnitude lower, so resolving it is 
    ; not practical, nor likely ever required. This is
    ; a touch confusing though ... i.e., when is this ever
    ; going to be important?
    
    f = 32d9
    Z = -1d0
    amu = _me_amu 
    B = 1d0
    BUnit = [0,0,1]
    density = 5d19
    harmonicNumber = 6
    
    T_eV = [20e3] 
    T_eV_kj = T_eV

    ; Analytic calculation range

    bMin = 0.3
    bMax = 1.2
    B_T = fIndGen(n)/(n-1)*(bMax-bMin)+bMin 

    ; KJ calculation range

    bMin_kj = 0.5
    bMax_kj = 0.6
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
    this_eps = kj_epsilon_hot(f, amu, Z, thisB, density, $
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

TemplateRunDir = '../template'
cd, current = RootDir
print, 'RootDir: ', RootDir

tic
cnt = 0
for t=0,nT_kj-1 do begin
for b=0,nB_kj-1 do begin
    
    ThisRunDir = string(cnt,format='(i5.5)')

    if keyword_set(runKJ) then begin

        file_delete, ThisRunDir, /recursive, /allow_nonexistent

        file_copy, TemplateRunDir, ThisRunDir, /recursive

    endif

    print, ThisRunDir, ' of ', strTrim(string((nT_kj*nB_kj)),2)
    cd, ThisRunDir

    for row=0,2 do begin

        RowString = string(row,format='(i1.1)')
        this_inputFileName = 'input/input-data_' + RowString + '.nc'
        This_jP2_FileName = 'output/jP2_' + RowString + '.nc'

        if row eq 0 then begin
            print, 'Ex!=0'
            E1=1
            E2=0
            E3=0
        endif

        if row eq 1 then begin
            print, 'Ey!=0'
            E1=0
            E2=1
            E3=0
        endif

        if row eq 2 then begin
            print, 'Ez!=0'
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
            ;print, 'T_eV: ', T_eV_kj[t], '  parSamples: ', parSamples
            perSamples = ( 2 * !pi / kPer ) / ( 1 / f / kj_nStepsPerCyclotronPeriod * vThMax )
            ;print, 'T_eV: ', T_eV_kj[t], '  perSamples: ', perSamples

            nMinSamples = 5
            par_dt = ( 2 * !pi / kPar ) / nMinSamples / vThMax 
            nHarmonic = 3
            wc = ( Z * _e ) * B_T_kj[b] / ( amu * _amu )
            per_dt = 1 / ( wc / 2 * !pi ) / ( nHarmonic * 2 )

            par_nStepsPerCyclotronPeriod = 1 / ( wc / 2 * !pi ) / par_dt 
            per_nStepsPerCyclotronPeriod = 1 / ( wc / 2 * !pi ) / per_dt 

            ;print, 'par_nStepsPerCyclotronPeriod: ', par_nStepsPerCyclotronPeriod
            ;print, 'per_nStepsPerCyclotronPeriod: ', per_nStepsPerCyclotronPeriod

            vPhsPer = w / kPer
            vPhsPar = w / kPar

            kj_write_kj_cfg, kj, './'

            ; Run kj

            RunCommand = 'DYLD_LIBRARY_PATH=/Users/dg6/code/netcdf-cxx4/lib ~/code/kineticj/bin/kineticj'
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

toc

for i=0,n_elements(sig_kj[0,0,*])-1 do begin

    sig_kj[*,*,i] = transpose(sig_kj[*,*,i])
  
endfor

; Plot results

layout=[3,3]
pos = 1
thick = 1 
style = '--'
transparency = 40
transparencyC = 0
xFS = 8
yFS = 8
margin = [0.24,0.15,0.1,0.15]
xRange = 0
kj_color1 = 'teal'
kj_color2 = 'tomato'
yLog = 0

plotThis = sig
plotThis_cold = sig_cold

if benchmark eq 1 or benchmark eq 2 then begin
    xTitle ='log10( T [eV] )'
    x = alog10(T_eV)
    x_kj = alog10(T_eV_kj)
endif

if benchmark eq 2 then begin

    plotThis = sig_swan

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

    C = 1e-10

    plotThis = complex( alog10( C + abs( real_part(sig) )), alog10( C + abs( imaginary(sig) )) )
    plotThis_cold = complex( alog10( C + abs( real_part(sig_cold)) ), alog10( C + abs( imaginary(sig_cold) )) )
    sig_kj = complex( alog10( C + abs( real_part(sig_kj)) ), alog10( C + abs( imaginary(sig_kj) )) )
    sig_swan = complex( alog10( C + abs( real_part(sig_swan)) ), alog10( C + abs( imaginary(sig_swan) )) )

endif

label=['x','y','z']
pos = 1
for i=0,2 do begin
for j=0,2 do begin

    subscript = label[i]+label[j]
    p=plot(x,plotThis[i,j,*],layout=[[layout],pos],$
            title='$\sigma_{'+subscript+'}$',yRange=[-1,1]*max(abs(plotThis[i,j,*])),/buffer,$
            font_size=12, xTitle=xTitle, xTickFont_size=xFS, yTickFont_size=yFS, $
            xMinor = 0, axis_style=1, yTitle='$\sigma_{'+subscript+'} [S/m]$', margin=margin, $
            xRange=xRange, /current )
    p=plot(x,imaginary(plotThis[i,j,*]),color='r',/over)

    p=plot(x,plotThis_cold[i,j,*],/over,thick=thick,transparency=transparencyC,LineStyle=style)
    p=plot(x,imaginary(plotThis_cold[i,j,*]),color='r',/over,thick=thick,transparency=transparencyC,LineStyle=style)
    
    p=plot(x_kj, sig_kj[i,j,*], /over, thick=3, color=kj_color1, transparency=transparency)
    p=plot(x_kj, imaginary(sig_kj[i,j,*]), color=kj_color2, /over, thick=3, transparency=transparency)
    
    pos++

endfor
endfor

p.save, benchmarkDirString+'.png', resolution=300, /transparent
p.save, benchmarkDirString+'.pdf'

cd, pwd

stop

end
