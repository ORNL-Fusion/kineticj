pro kj_stix_current, overPlotAR2 = _overPlotAR2, useRS=_useRS, hot=_hot

if keyword_set(_overPlotAR2) then overPlotAR2 = 1 else overPlotAR2 = 0
if keyword_set(_useRS) then useRS = _useRS else useRS = 0
if keyword_set(_hot) then hot = _hot else hot = 0

useAR = 1
if useRS then useAR = 0 

@constants

; Get the plasma parameters

arP = ar2_read_runData('./',1)

f = arP.freq
w = 2 * !pi * f
B = sqrt( arP.br^2 + arP.bt^2 + arP.bz^2 )
r = arP.r
nPhi = arP.nPhi
nX = n_elements(arP.br)
density = arP.densitySpec
temp = arP.tempSpec
nuOmg = arP.nuOmg

ar2 = ar2_read_ar2input('./')

amu = ar2.amu
atomicZ = ar2.atomicZ

nS = n_elements(amu)

; Get the E field

if useAR then begin
    print, 'Reading AR Solution'
    arS = ar2_read_solution('./',1)
    solution = arS
endif

if useRS then begin
    print, 'Reading RS Solution'
    rsS = rsfwc_read_solution('./')
    solution = rsS
endif

if hot then begin
    print, 'Using HOT dielectric'
endif else begin
    print, 'Using COLD dielectric'
endelse

; For each species

kx = 0
;T_eV = 2e3
kPar = nPhi / r 
harmonicNumber = 3

windowWidth = 100

jr = complexArr(nX,nS)
jt = complexArr(nX,nS)
jz = complexArr(nX,nS)

run = 1

if run then begin

for s=0,nS-1 do begin
for i=windowWidth/2,nX-windowWidth/2-1 do begin

    ;print, 'Spec: ', s
    ;print, 'iX: ', i

    iL = i-windowWidth/2
    iR = i+windowWidth/2 

    ; Extract and window the E field
    
    N = n_elements(solution.E_r[iL:iR])

    Er = solution.E_r[iL:iR] * hanning(N)
    Et = solution.E_t[iL:iR] * hanning(N)
    Ez = solution.E_z[iL:iR] * hanning(N)
    
    ; Forward FFT
    
    Ekr = fft(Er,/center)
    Ekt = fft(Et,/center)
    Ekz = fft(Ez,/center)
    
    dX = r[1]-r[0]
    kxAxis = fIndGen(N) / ( dX * N ) * 2*!Pi 
    dk = kxAxis[1]-kxAxis[0]
    kxAxis = kxAxis - kxAxis[-1]/2 - dk/2
    
    xRange = [-1,1]*600
    
    ;p=plot(kxAxis,abs(Ekr)^2,layout=[1,3,1],xRange=xRange)
    ;p=plot(kxAxis,abs(Ekt)^2,layout=[1,3,2],/current,xRange=xRange)
    ;p=plot(kxAxis,abs(Ekz)^2,layout=[1,3,3],/current,xRange=xRange)

    ; Get sigma for each k 
  
    jkr = complexArr(N)
    jkt = complexArr(N)
    jkz = complexArr(N)

    ;for k=0,N-1 do begin

        print, 'TEMP: ',temp[i,0,s]
        epsilon = kj_hot_epsilon( f, amu[s], atomicZ[s], B[i], $
                density[i,0,s], harmonicNumber, kPar[i], kxAxis, temp[i,0,s], $
                kx = kx, nuOmg = nuOmg[i,0,s], epsilon_cold = epsilon_cold );

        epsilon_cold = complex(rebin(real_part(epsilon_cold),3,3,N),rebin(imaginary(epsilon_cold),3,3,N))

        sigma =  ( epsilon - rebin(identity(3),3,3,N) ) * w * _e0 / _ii
        sigma_cold =  ( epsilon_cold - rebin(identity(3),3,3,N) ) * w * _e0 / _ii

        ; Choose hot or cold sigma 

        _sigma = sigma_cold
        if hot then begin
            _sigma = sigma
        endif

        ; Calculate k-space plasma current
        ; This would have to be generalize for non magnetically aligned coordinates.

        ; Try flipping the indexing order also.

        jkr = (_sigma[0,0,*] * Ekr + _sigma[1,0,*] * Ekz + _sigma[2,0,*] * Ekt)[*]
        jkz = (_sigma[0,1,*] * Ekr + _sigma[1,1,*] * Ekz + _sigma[2,1,*] * Ekt)[*]
        jkt = (_sigma[0,2,*] * Ekr + _sigma[1,2,*] * Ekz + _sigma[2,2,*] * Ekt)[*]

    ;endfor
    
    ; Inverse FFT for configuration-space plasma current
    
    thisjr = fft(jkr,/center,/inverse)
    thisjt = fft(jkt,/center,/inverse)
    thisjz = fft(jkz,/center,/inverse)

    jr[i,s] = thisjr[N/2]
    jt[i,s] = thisjt[N/2]
    jz[i,s] = thisjz[N/2]

endfor
endfor

    save, jr, jt, jz, fileName='kj_sc.sav', /variables

endif else begin

    restore, 'kj_sc.sav'

endelse

Er = solution.E_r
Et = solution.E_t
Ez = solution.E_z

p=plot(r,Er,layout=[1,3,1])
p=plot(r,imaginary(Er),color='r',/over)

p=plot(r,Et,layout=[1,3,2],/current)
p=plot(r,imaginary(Et),color='r',/over)

p=plot(r,Ez,layout=[1,3,3],/current)
p=plot(r,imaginary(Ez),color='r',/over)

current=0
for s=0,nS-1 do begin

    p=plot(r,jr[*,s],layout=[nS,3,1+s],current=current)
    p=plot(r,imaginary(jr[*,s]),color='r',/over)
    if overPlotAR2 then begin
        p=plot(solution.r,solution.JP_r[*,0,s],thick=4,transparency=80,/over)            
        p=plot(solution.r,imaginary(solution.JP_r[*,0,s]),color='r',thick=4,transparency=80,/over)            
    endif

    current = (current + 1)<1
    p=plot(r,jt[*,s],layout=[nS,3,1+1*nS+s],current=current)
    p=plot(r,imaginary(jt[*,s]),color='r',/over)
    if overPlotAR2 then begin
        p=plot(solution.r,solution.JP_t[*,0,s],thick=4,transparency=80,/over)            
        p=plot(solution.r,imaginary(solution.JP_t[*,0,s]),color='r',thick=4,transparency=80,/over)            
    endif

    p=plot(r,jz[*,s],layout=[nS,3,1+2*nS+s],current=current)
    p=plot(r,imaginary(jz[*,s]),color='r',/over)
    if overPlotAR2 then begin
        p=plot(solution.r,solution.JP_z[*,0,s],thick=4,transparency=80,/over)            
        p=plot(solution.r,imaginary(solution.JP_z[*,0,s]),color='r',thick=4,transparency=80,/over)            
    endif

endfor

stop
end
