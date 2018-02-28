pro kj_stix_current, overPlotRefSolution = _overPlotRefSolution, useRS=_useRS, hot=_hot, $
        jr = jr, jt = jt, jz = jz, $
        kjDeltaFileName_in = _kjDeltaFileName_in, $
        kjDeltaFileName_out = _kjDeltaFileName_out, $
        referenceSolutionDir = _referenceSolutionDir, rgrid = rgrid

if keyword_set(_overPlotRefSolution) then overPlotRefSolution = 1 else overPlotRefSolution = 0
if keyword_set(_useRS) then useRS = _useRS else useRS = 0
if keyword_set(_hot) then hot = _hot else hot = 0
if keyword_set(_kjDeltaFileName_in) then kjDeltaFileName_in = _kjDeltaFileName_in else kjDeltaFileName_in = 'kj-delta-in.nc'
if keyword_set(_kjDeltaFileName_out) then kjDeltaFileName_out = _kjDeltaFileName_out else kjDeltaFileName_out = 'kj-delta-out.nc'
if keyword_set(_referenceSolutionDir) then referenceSolutionDir = _referenceSolutionDir else referenceSolutionDir = './'

useAR = 1
if useRS then useAR = 0 

@dlg_constants

; Get the plasma parameters

ar2 = ar2_read_ar2input('./')
amu = ar2.amu
atomicZ = ar2.atomicZ
nS = n_elements(amu)

; Get the E field

if useAR then begin

    arR = ar2_read_runData('./',1)

    f = arR['freq']

    br = arR['brU']*arR['bmod']
    bt = arR['btU']*arR['bmod']
    bz = arR['bzU']*arR['bmod']

    B = sqrt( br^2 + bt^2 + bz^2 )
    r = arR['r']
    nPhi = arR['nPhi']
    kz = rs['kz']
    density = reform(arR['densitySpec'])
    temp = reform(arR['tempSpec']) 
    nu_omg = reform(arR['nuOmg'])

    solution = ar2_read_solution('./',1)
    solution_ref = ar2_read_solution(referenceSolutionDir,1)

endif

if useRS then begin

    rs = rs_read_runData('./')

    f = rs['freq']
    br = rs['br']
    bt = rs['bt']
    bz = rs['bz']
    B = sqrt( br^2 + bt^2 + bz^2 )
    r = rs['r']
    nPhi = rs['nPhi']
    kz = rs['kz']
    density = rs['densitySpec']
    temp = density*0 

    ; Since RS doesn't have temperature, we import from ar2Input.nc
    ; ----------
    nS = n_elements(density[0,*])
    z = r*0
    for s=0,nS-1 do begin
        temp[*,s] = reform(ar2.temp_eV[*,*,s])
    endfor 
    ; ----------

    ;temp = temp * 0.0001

    nu_omg = rs['nuOmg']

    solution = rsfwc_read_solution('./')
    solution_ref = rsfwc_read_solution(referenceSolutionDir)
    
endif

w = 2 * !pi * f
nX = n_elements(r)

if hot then begin
    print, 'Using HOT dielectric'
endif else begin
    print, 'Using COLD dielectric'
endelse

; For each species

kPar = nPhi / r
harmonicNumber = 3

windowWidth = 512

jr = complexArr(nX,nS)
jt = complexArr(nX,nS)
jz = complexArr(nX,nS)

jrAll = complexArr(nX,nX,nS)
jtAll = complexArr(nX,nX,nS)
jzAll = complexArr(nX,nX,nS)

run = 1
savFileName = 'kj-stix.sav'

if run then begin

kxArr = fltArr(nX,windowWidth+1)
NArr = fltArr(nX)
ekrArr = complexArr(nX,windowWidth+1)
ektArr = complexArr(nX,windowWidth+1)
ekzArr = complexArr(nX,windowWidth+1)
sig2 = complexArr(3,3,nX,windowWidth+1,nS)
sigc = complexArr(3,3,nX,windowWidth+1,nS)
sig2_abp = complexArr(3,3,nX,windowWidth+1,nS)
sigc_abp = complexArr(3,3,nX,windowWidth+1,nS)

; Move the rotation matrix outside of the loop 

R_abp_to_rtz = fltArr(3,3,nX) 
for i=1,nX-2 do begin
    thisBUnitVec = [br[i],bt[i],bz[i]]/B[i]
    R_abp_to_rtz[*,*,i] = get_rotmat_abp_to_rtz(thisBUnitVec)
endfor

; Main loop

for s=0,nS-1 do begin

    for i=1,nX-2 do begin
    
        iL = (i-windowWidth/2)>0
        iR = (i+windowWidth/2)<(nX-1)
    
        thisWindowWidth = min([i-iL,iR-i])*2+1
    
        if thisWindowWidth gt windowWidth+1 then stop
    
        iL = i - (thisWindowWidth-1)/2
        iR = i + (thisWindowWidth-1)/2
    
        _iL = 0 
        _iR = thisWindowWidth - 1
    
        ; Extract and window the E field
        
        N = n_elements((solution['E_r'])[iL:iR])
    
        er = +(solution['E_r'])[iL:iR] * hanning(n)
        et = +(solution['E_t'])[iL:iR] * hanning(n)
        ez = +(solution['E_z'])[iL:iR] * hanning(n)

        ;; Test for FFT
        ;kk = 200.0  
        ;er = exp(_ii * kk * r[iL:iR]) 
       
        ; forward fft
        
        ekr = fft(er,/center)
        ekt = fft(et,/center)
        ekz = fft(ez,/center)
       
        ekrArr[i,_iL:_iR] = ekr
        ektArr[i,_iL:_iR] = ekt
        ekzArr[i,_iL:_iR] = ekz
  
        dx = r[1]-r[0]
        kxaxis = findgen(n) / ( dx * n ) * 2*!pi 
        dk = kxaxis[1]-kxaxis[0]
        kxaxis = kxaxis - kxaxis[-1]/2 
        if (n_elements(kxaxis) mod 2) eq 0 then begin
            kxaxis = kxaxis - dk/2
        endif
      
        kxArr[i,_iL:_iR] = kxAxis
        nArr[i] = n_elements(kxAxis)
    
        kPer = sqrt(kxAxis^2+kz^2)
    
        ; Get sigma for each k 
      
        jkr = complexArr(N)
        jkt = complexArr(N)
        jkz = complexArr(N)
    
        ;for k=0,N-1 do begin
    
            epsilon_bram = kj_epsilon_hot( f, amu[s], atomicZ[s], B[i], $
                    density[i,s], harmonicNumber, kPar[i], kPer, temp[i,s], $
                    kx = 0, nuOmg = nu_omg[i,s]);, epsilon_cold = epsilon_cold) ;, $
                    ;epsilon_swan_ND = epsilon_swan );
   
            epsilon_cold = kj_epsilon_cold(f, amu[s], atomicZ[s], B[i], $
                    density[i,s], nu_omg[i,s])

            epsilon_cold = complex(rebin(real_part(epsilon_cold),3,3,N),rebin(imaginary(epsilon_cold),3,3,N))
    
            sigma_bram =  ( epsilon_bram - rebin(identity(3),3,3,N) ) * w * _e0 / _ii
            ;sigma_swan =  ( epsilon_swan - rebin(identity(3),3,3,N) ) * w * _e0 / _ii
            sigma_cold =  ( epsilon_cold - rebin(identity(3),3,3,N) ) * w * _e0 / _ii
    
            ; Rotate sigma from ABP to RTZ

            sigma_abp_bram = sigma_bram
            sigma_abp_cold = sigma_cold
   
            ; Choose hot or cold sigma 
    
            if hot then begin

                for k=0,N-1 do begin
    	            sigma_bram[*,*,k] = rotateEpsilon ( sigma_bram[*,*,k], thisBUnitVec, R = R_abp_to_rtz[*,*,i] )
    	            ;sigma_swan[*,*,k] = rotateEpsilon ( sigma_swan[*,*,k], thisBUnitVec, R = R_abp_to_rtz[*,*,i] )
                endfor
 
                ;_sigma = sigma_swan
                _sigma = sigma_bram

            endif else begin

                for k=0,N-1 do begin
    	            sigma_cold[*,*,k] = rotateEpsilon ( sigma_cold[*,*,k], thisBUnitVec, R = R_abp_to_rtz[*,*,i] )
                endfor

                _sigma = sigma_cold

            endelse

            for k=0,N-1 do begin
    	            sigma_cold[*,*,k] = rotateEpsilon ( sigma_cold[*,*,k], thisBUnitVec, R = R_abp_to_rtz[*,*,i] )
            endfor


            sig2[*,*,i,_iL:_iR,s] = _sigma
            sigc[*,*,i,_iL:_iR,s] = sigma_cold
            sig2_abp[*,*,i,_iL:_iR,s] = sigma_abp_bram
            sigc_abp[*,*,i,_iL:_iR,s] = sigma_abp_cold
    
            ; Calculate k-space plasma current
    
            jkr = +(_sigma[0,0,*] * Ekr + _sigma[1,0,*] * Ekt + _sigma[2,0,*] * Ekz)[*]
            jkt = +(_sigma[0,1,*] * Ekr + _sigma[1,1,*] * Ekt + _sigma[2,1,*] * Ekz)[*]
            jkz = +(_sigma[0,2,*] * Ekr + _sigma[1,2,*] * Ekt + _sigma[2,2,*] * Ekz)[*]

        ;endfor
        
        ; Inverse FFT for configuration-space plasma current
        
        thisjr = fft(jkr,/center,/inverse)
        thisjt = fft(jkt,/center,/inverse)
        thisjz = fft(jkz,/center,/inverse)

        ; Extract the central point from the FFT

        jr[i,s] = thisjr[N/2]
        jt[i,s] = thisjt[N/2]
        jz[i,s] = thisjz[N/2]
   
        ; Store all the j(x) such that we can average later. 

        jrAll[i,iL:iR,s] = thisjr
        jtAll[i,iL:iR,s] = thisjt
        jzAll[i,iL:iR,s] = thisjz

    endfor

    ;; Now do the windowed average over spatial points
   
    ;windowWidth2 = 10 
    ;for i=1,nX-2 do begin

    ;    iL = (i-windowWidth2/2)>0
    ;    iR = (i+windowWidth2/2)<(nX-1)
    ;
    ;    thisWindowWidth = min([i-iL,iR-i])*2+1
    ;
    ;    if thisWindowWidth gt windowWidth+1 then stop
    ;
    ;    iL = i - (thisWindowWidth-1)/2
    ;    iR = i + (thisWindowWidth-1)/2
 
    ;    N = n_elements(jrAll[iL:iR,0,0])
    ;    win = kj_hanning(N,/sym)
    ;    
    ;    win = win/total(win)
    ;    win2D = rebin(win,N,n_elements(jrAll[0,*,0]))

    ;    thisJr = total(reform(jrAll[iL:iR,*,s]) * win2D,1)
    ;    thisJt = total(reform(jtAll[iL:iR,*,s]) * win2D,1)
    ;    thisJz = total(reform(jzAll[iL:iR,*,s]) * win2D,1)

    ;    jr[i,s] = thisjr[N/2]
    ;    jt[i,s] = thisjt[N/2]
    ;    jz[i,s] = thisjz[N/2]

    ;    ;if i eq 350 then begin
    ;    ;for j=0,N-1 do begin
    ;    ;    p=plot(jtAll[iL+j,*,s],/over)
    ;    ;endfor
    ;    ;p=plot(jtAll[iL+N/2,*,s],/over,color='b')
    ;    ;p=plot(thisJt,/over,color='r')
    ;    ;;stop
    ;    ;endif
    ;endfor

endfor

    save, jr, jt, jz, sig2, fileName=savFileName, /variables

endif else begin

    restore, savFileName

endelse

Er = solution['E_r']
Et = solution['E_t']
Ez = solution['E_z']

sign = -1

delta_r = sign * (reform(solution_ref['jP_r_spec']) - jr)
delta_t = sign * (reform(solution_ref['jP_t_spec']) - jt)
delta_z = sign * (reform(solution_ref['jP_z_spec']) - jz)

doPlots = 1

if doPlots then begin

p=plot(r,Er,layout=[1,3,1],title='Er')
p=plot(r,imaginary(Er),color='r',/over)

p=plot(r,Et,layout=[1,3,2],/current,title='Et')
p=plot(r,imaginary(Et),color='r',/over)

p=plot(r,Ez,layout=[1,3,3],/current,title='Ez')
p=plot(r,imaginary(Ez),color='r',/over)

current=0
for s=0,nS-1 do begin

    p=plot(r,jr[*,s],layout=[nS,3,1+s],current=current)
    p=plot(r,imaginary(jr[*,s]),color='r',/over)
    if overPlotRefSolution then begin
        p=plot(solution_ref['r'],(solution['jP_r_spec'])[*,0,s],thick=4,transparency=80,/over)            
        p=plot(solution_ref['r'],imaginary((solution['jP_r_spec'])[*,0,s]),color='r',thick=4,transparency=80,/over)            
        p=plot(solution_ref['r'],delta_r[*,s],thick=2,/over,lineStyle='--')            
        p=plot(solution_ref['r'],imaginary(delta_r[*,s]),color='r',/over,lineStyle='--',thick=2)            
    endif

    current = (current + 1)<1
    p=plot(r,jt[*,s],layout=[nS,3,1+1*nS+s],current=current)
    p=plot(r,imaginary(jt[*,s]),color='r',/over)
    if overPlotRefSolution then begin
        p=plot(solution_ref['r'],(solution['jP_t_spec'])[*,0,s],thick=4,transparency=80,/over)            
        p=plot(solution_ref['r'],imaginary((solution['jP_t_spec'])[*,0,s]),color='r',thick=4,transparency=80,/over)            
        p=plot(solution_ref['r'],delta_t[*,s],thick=2,/over,lineStyle='--')            
        p=plot(solution_ref['r'],imaginary(delta_t[*,s]),color='r',/over,lineStyle='--',thick=2)            
    endif

    p=plot(r,jz[*,s],layout=[nS,3,1+2*nS+s],current=current)
    p=plot(r,imaginary(jz[*,s]),color='r',/over)
    if overPlotRefSolution then begin
        p=plot(solution_ref['r'],(solution['jP_z_spec'])[*,0,s],thick=4,transparency=80,/over)            
        p=plot(solution_ref['r'],imaginary((solution['jP_z_spec'])[*,0,s]),color='r',thick=4,transparency=80,/over)            
        p=plot(solution_ref['r'],delta_z[*,s],thick=2,/over,lineStyle='--')            
        p=plot(solution_ref['r'],imaginary(delta_z[*,s]),color='r',/over,lineStyle='--',thick=2)            
    endif

endfor

endif

plotSpectra = 0

if plotSpectra then begin
if not useRS then begin

    ; Compare the global spectrum and sigmas for that spectrum
    
    xrange = [-200,800]
    
    eka = fft(solution.ealp,/center)
    ekb = fft(solution.ebet,/center)
    ekp = fft(solution.eprl,/center)
    N = n_elements(solution.eprl)
    
    dx = r[1]-r[0]
    kxaxis = findgen(n) / ( dx * N ) * 2*!pi 
    dk = kxaxis[1]-kxaxis[0]
    kxaxis = kxaxis - kxaxis[-1]/2
    if (n_elements(kxaxis) mod 2) eq 0 then begin
        kxaxis = kxaxis - dk/2
    endif

    s = 1
    iX = nX/2 ; This should be about the same as the AR window. 
    
    thiskr = ekrArr[iX,*]
    thiskt = ektArr[iX,*]
    thiskz = ekzArr[iX,*]

    thiska = thiskr*0
    thiskb = thiskt*0
    thiskp = thiskz*0

    for i=0,n_elements(thiskr)-1 do begin

        thisk_rtz = [thiskr[i],thiskt[i],thiskz[i]]
        thisk_abp = transpose(R_abp_to_rtz[*,*,iX]) # thisk_rtz
        
        thiska[i] = thisk_abp[0]
        thiskb[i] = thisk_abp[1]
        thiskp[i] = thisk_abp[2]

    endfor

    p=plot(solution.kr,solution.ealpk,layout=[1,3,1],xrange=xrange)
    p=plot(solution.kr,imaginary(solution.ealpk),color='r',/over)
    p=plot(kxaxis,eka,/over,thick=3,transparency=50)
    p=plot(kxaxis,imaginary(eka),color='r',/over,thick=3,transparency=50)
    p=plot(kxArr[iX,*],thiska,/over,color='b',thick=4,trans=50)
    p=plot(kxArr[iX,*],imaginary(thiska),/over,color='magenta',thick=4,trans=50)
    
    p=plot(solution.kr,solution.ebetk,layout=[1,3,2],/current,xrange=xrange)
    p=plot(solution.kr,imaginary(solution.ebetk),color='r',/over)
    p=plot(kxaxis,ekb,/over,thick=3,transparency=50)
    p=plot(kxaxis,imaginary(ekb),color='r',/over,thick=3,transparency=50)
    p=plot(kxArr[iX,*],thiskb,/over,color='b',thick=4,trans=50)
    p=plot(kxArr[iX,*],imaginary(thiskb),/over,color='magenta',thick=4,trans=50)
    
    p=plot(solution.kr,solution.eprlk,layout=[1,3,3],/current,xrange=xrange)
    p=plot(solution.kr,imaginary(solution.eprlk),color='r',/over)
    p=plot(kxaxis,ekp,thick=3,transparency=50,/over)
    p=plot(kxaxis,imaginary(ekp),color='r',/over,thick=3,transparency=50)
    p=plot(kxArr[iX,*],thiskp,/over,color='b',thick=4,trans=50)
    p=plot(kxArr[iX,*],imaginary(thiskp),/over,color='magenta',thick=4,trans=50)
    
    current = 0 
    for ii=0,2 do begin
    for jj=0,2 do begin
        p=plot(solution.kr,solution.sig[ii,jj,ix,*,s],layout=[3,3,ii*3+jj+1],current=current)
        p=plot(solution.kr,imaginary(solution.sig[ii,jj,ix,*,s]),color='r',/over)
        p=plot(kxArr[iX,0:nArr[iX]-1],sig2_abp[ii,jj,ix,0:nArr[iX]-1,s],/over,thick=3,trans=50,linestyle='--')
        p=plot(kxArr[iX,0:nArr[iX]-1],imaginary(sig2_abp[ii,jj,ix,0:nArr[iX]-1,s]),color='r',/over,thick=3,trans=50,linestyle='--')
        p=plot(kxArr[iX,0:nArr[iX]-1],sigc_abp[ii,jj,ix,0:nArr[iX]-1,s],/over,linestyle=':',thick=2)
        p=plot(kxArr[iX,0:nArr[iX]-1],imaginary(sigc_abp[ii,jj,ix,0:nArr[iX]-1,s]),color='r',/over,linestyle=':',thick=2)
        current = (current + 1)<1
    endfor
    endfor 

   stop 
endif
endif

updateDeltaFile = 0

if updateDeltaFile then begin

    ; Update the kj-delta-out file by adding the new delta to the 
    ; kj-delta-in file
    
    kj_in = ncdf_parse(kjDeltaFileName_in,/read)
    
    jP_r_re_in = kj_in['jP_r_re','_DATA']
    jP_r_im_in = kj_in['jP_r_im','_DATA']
    jP_t_re_in = kj_in['jP_t_re','_DATA']
    jP_t_im_in = kj_in['jP_t_im','_DATA']
    jP_z_re_in = kj_in['jP_z_re','_DATA']
    jP_z_im_in = kj_in['jP_z_im','_DATA']
    
    jP_r_in = complex(jP_r_re_in,jP_r_im_in)
    jP_t_in = complex(jP_t_re_in,jP_t_im_in)
    jP_z_in = complex(jP_z_re_in,jP_z_im_in)
    
    nc_id = nCdf_create (kjDeltaFileName_out, /clobber )
    
    	nCdf_control, nc_id, /fill
    	
    	nr_id = nCdf_dimDef ( nc_id, 'nR', n_elements(r) )
    	scalar_id = nCdf_dimDef ( nc_id, 'scalar', 1 )
    
    	freq_id = nCdf_varDef ( nc_id, 'freq', scalar_id, /float )
    	r_id = nCdf_varDef ( nc_id, 'r', nr_id, /float )
    
    	jr_re_id = nCdf_varDef ( nc_id, 'jP_r_re', nr_id, /float )
    	jr_im_id = nCdf_varDef ( nc_id, 'jP_r_im', nr_id, /float )
    	jt_re_id = nCdf_varDef ( nc_id, 'jP_t_re', nr_id, /float )
    	jt_im_id = nCdf_varDef ( nc_id, 'jP_t_im', nr_id, /float )
    	jz_re_id = nCdf_varDef ( nc_id, 'jP_z_re', nr_id, /float )
    	jz_im_id = nCdf_varDef ( nc_id, 'jP_z_im', nr_id, /float )
    
        nCdf_control, nc_id, /enDef
    
    	nCdf_varPut, nc_id, freq_id, f
    
    	nCdf_varPut, nc_id, r_id, r 
       
    	nCdf_varPut, nc_id, jr_re_id, real_part( jP_r_in + total(delta_r,2) )
    	nCdf_varPut, nc_id, jr_im_id, imaginary( jP_r_in + total(delta_r,2) )
    	nCdf_varPut, nc_id, jt_re_id, real_part( jP_t_in + total(delta_t,2) )
    	nCdf_varPut, nc_id, jt_im_id, imaginary( jP_t_in + total(delta_t,2) )
    	nCdf_varPut, nc_id, jz_re_id, real_part( jP_z_in + total(delta_z,2) )
    	nCdf_varPut, nc_id, jz_im_id, imaginary( jP_z_in + total(delta_z,2) )
    
    nCdf_close, nc_id

end

; Also write an actual Jp file

kjJpFileName_out = 'kj-jp.nc'

nc_id = nCdf_create (kjJpFileName_out, /clobber )

	nCdf_control, nc_id, /fill
	
	nr_id = nCdf_dimDef ( nc_id, 'nR', n_elements(r) )
	scalar_id = nCdf_dimDef ( nc_id, 'scalar', 1 )

	freq_id = nCdf_varDef ( nc_id, 'freq', scalar_id, /float )
	r_id = nCdf_varDef ( nc_id, 'r', nr_id, /float )

	jr_re_id = nCdf_varDef ( nc_id, 'jP_r_re', nr_id, /float )
	jr_im_id = nCdf_varDef ( nc_id, 'jP_r_im', nr_id, /float )
	jt_re_id = nCdf_varDef ( nc_id, 'jP_t_re', nr_id, /float )
	jt_im_id = nCdf_varDef ( nc_id, 'jP_t_im', nr_id, /float )
	jz_re_id = nCdf_varDef ( nc_id, 'jP_z_re', nr_id, /float )
	jz_im_id = nCdf_varDef ( nc_id, 'jP_z_im', nr_id, /float )

    nCdf_control, nc_id, /enDef

	nCdf_varPut, nc_id, freq_id, f

	nCdf_varPut, nc_id, r_id, r 
   
	nCdf_varPut, nc_id, jr_re_id, real_part( total(jr,2) )
	nCdf_varPut, nc_id, jr_im_id, imaginary( total(jr,2) )
	nCdf_varPut, nc_id, jt_re_id, real_part( total(jt,2) )
	nCdf_varPut, nc_id, jt_im_id, imaginary( total(jt,2) )
	nCdf_varPut, nc_id, jz_re_id, real_part( total(jz,2) )
	nCdf_varPut, nc_id, jz_im_id, imaginary( total(jz,2) )

nCdf_close, nc_id


rgrid = r

end
