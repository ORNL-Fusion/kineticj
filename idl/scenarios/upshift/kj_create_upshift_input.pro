@xyz_cyl 
@interpb
@kj_zfunction

pro kj_create_upshift_input

    @constants    

    FieldsOutFileName = 'data/kj_upshift.nc'
    f_Hz = 30e6
    n_e = 1.1d14
    bPolFactor = 1.0
    EqdskFile = 'g130608.00355.EFIT02.mds.corrected.qscale_1.00000'
    E_keV = 2.0
    Np = 500
    AtomicZ = -1
    Theta_Toroidal_Sign = -1 ; Switch this depending on the magnetic field tordoial direction.

    ParticlesOutFileName = 'data/f_E_keV_'+string(E_keV,format='(f3.1)')+$
            '_n_e_'+string(n_e,format='(e7.1)')+'_Np_'+string(Np,format='(i4.4)')+'.nc'

    m = 15
    n = 81
    nPhi = -12 

    c0_CYL = [1.2,0.0,0.0]

    c0_XYZ = Coords_CYL_to_XYZ(c0_CYL)
    dS = 0.01
    TraceNPts = 2000

    g = ReadGEqdsk ( EqdskFile, $
            bPolFactor = bPolFactor, $
            fieldLineIn = c0_CYL, $
            FieldLine_XYZ = FieldLine_XYZ, $
            FieldLine_CYL = FieldLine_CYL, $
            B_AlongFieldLine_XYZ = B_AlongFieldLine_XYZ, $
            B_AlongFieldLine_CYL = B_AlongFieldLine_CYL, $
            SafetyFactor = SafetyFactor, $
            FieldLineTraceDir = 1, $
            FieldLineTraceDS = dS, $
            FieldLineTraceNSteps = TraceNPts, $
            BInterpS = BInterpS ) 

    g = ReadGEqdsk ( EqdskFile, $
            bPolFactor = bPolFactor, $
            fieldLineIn = c0_CYL, $
            FieldLine_XYZ = FieldLine_XYZ_, $
            FieldLine_CYL = FieldLine_CYL_, $
            B_AlongFieldLine_XYZ = B_AlongFieldLine_XYZ_, $
            B_AlongFieldLine_CYL = B_AlongFieldLine_CYL_, $
            SafetyFactor = SafetyFactor_, $
            FieldLineTraceDir = -1, $
            FieldLineTraceDS = dS, $
            FieldLineTraceNSteps = TraceNPts ) 

    B_XYZ = [ [reverse(B_AlongFieldLine_XYZ_[*,0:TraceNPts-1],2)], [B_AlongFieldLine_XYZ[*,1:TraceNPts-1]] ]
    B_CYL = [ [reverse(B_AlongFieldLine_CYL_[*,0:TraceNPts-1],2)], [B_AlongFieldLine_CYL[*,1:TraceNPts-1]] ]

    c_XYZ = [ [reverse(FieldLine_XYZ_[*,0:TraceNPts-1],2)], [FieldLine_XYZ[*,1:TraceNPts-1]] ]
    c_CYL = [ [reverse(FieldLine_CYL_[*,0:TraceNPts-1],2)], [FieldLine_CYL[*,1:TraceNPts-1]] ]

    s_Coord = (fIndGen(TraceNpts*2-1)-TraceNPts+1)*dS 


    nS = n_elements(s_Coord)
	SPointsMin = -3.0
	SPointsMax = +5.0
	SPointsN = 40
	SPoints = fIndGen(SPointsN)/(SPointsN-1)*(SPointsMax-SPointsMin)+SPointsMin

	SPointsIndex = IntArr(SPointsN)
	for i=0,SPointsN-1 do begin
    	SPointsIndex[i] = where(abs(s_Coord-SPoints[i]) eq min(abs(s_Coord-SPoints[i])))
	endfor

    BMag = sqrt(total(B_XYZ^2,1))
    BMag_CYL = sqrt(total(B_CYL^2,1))

    print, 'Safety Factor: ', SafetyFactor, SafetyFactor_

    vColors = BytScl(BMag, top=253)+1
    p = plot3d((c_XYZ[0,*])[*],(c_XYZ[1,*])[*],(c_XYZ[2,*])[*],$
            thick = 3.0, $
            aspect_ratio = 1.0, aspect_z = 1.0, $
            rgb_table = 7, vert_colors = vColors,$
            depth_cue = [0,2], /perspective )

    nFrames = 48 
    for _n = 0, nFrames*2/3-1 do begin
            thisR = g.rLim
            thisZ = g.zLim
            thisT = 360./(nFrames-1)*_n+fltArr(n_elements(thisR))
            Frame_Cyl = transpose([[thisT[*]],[thisR[*]],[thisZ[*]]])
            Frame_Rec = cv_coord(from_cylin=Frame_Cyl,/degrees,/to_rect)
            pF = plot3d((Frame_Rec[0,*])[*],(Frame_Rec[1,*])[*],(Frame_Rec[2,*])[*],/over,thick=5,transparency=70)
    endfor


    ; Create a Fourier basis function of given m,n,nPhi

    nR = 129
    nZ = 257

    rMin = 0.18
    rMax = 1.60
    
    zMin = +1.65
    zMax = -1.65

    x = fIndGen(nR)/(nR-1)*2*!pi
    y = fIndGen(nZ)/(nZ-1)*2*!pi

    r = fIndGen(nR)/(nR-1)*(rMax-rMin)+rMin
    z = fIndGen(nZ)/(nZ-1)*(zMax-zMin)+zMin

    x2d = rebin(x,nR,nZ)
    y2d = transpose(rebin(y,nZ,nR))

    ii = complex(0,1)
    E = exp( ii * (n * x2d + m * y2d) )   

    bru_2d = fltArr(nR,nZ)
    btu_2d = fltArr(nR,nZ)
    bzu_2d = fltArr(nR,nZ)

    for i=0,nR-1 do begin
            for j=0,nZ-1 do begin

                this_c_CYL = [r[i],0,z[j]]
                this_b_CYL = bHere_CYL(bInterpS,this_c_CYL, bMag=this_bMag)
                this_bu_CYL = this_b_CYL / this_bMag
                
                bru_2d[i,j] = this_bu_CYL[0]
                btu_2d[i,j] = this_bu_CYL[1]
                bzu_2d[i,j] = this_bu_CYL[2]

            endfor
    endfor

    Er_2d = E * bru_2d
    Et_2d = E * btu_2d
    Ez_2d = E * bzu_2d

    c = contour ( Er_2d, r, z, rgb_table = 17, /fill, aspect_ratio = 1.0, layout=[3,1,1], current=0 )
    p = plot ( g.rlim, g.zLim, /over, thick = 3, transparency=70)
    p = plot ( c_CYL[0,*],c_CYL[2,*], /over, thick = 3, rgb_table =  7, vert_colors = vColors )  

    c = contour ( Et_2d, r, z, rgb_table = 17, /fill, aspect_ratio = 1.0, layout=[3,1,2], current=1 )
    p = plot ( g.rlim, g.zLim, /over, thick = 3, transparency=70)
    p = plot ( c_CYL[0,*],c_CYL[2,*], /over, thick = 3, rgb_table =  7, vert_colors = vColors )  

    c = contour ( Ez_2d, r, z, rgb_table = 17, /fill, aspect_ratio = 1.0, layout=[3,1,3], current=1 )
    p = plot ( g.rlim, g.zLim, /over, thick = 3, transparency=70)
    p = plot ( c_CYL[0,*],c_CYL[2,*], /over, thick = 3, rgb_table =  7, vert_colors = vColors )  


    c = contour ( E, r, z, rgb_table = 17, /fill, aspect_ratio = 1.0 )
    p = plot ( g.rlim, g.zLim, /over, thick = 3, transparency=70)
    p = plot ( c_CYL[0,*],c_CYL[2,*], /over, thick = 3, rgb_table =  7, vert_colors = vColors )  

    kr = 2*!pi*n/(rMax-rMin)
    kz = 2*!pi*m/(zMax-zMin)

    kDotR_AlongS  = fltArr(nS)
    kb = fltArr(nS)
    Eb = complexArr(nS)
    Eb_check = complexArr(nS)
    for s=0,nS-1 do begin

        this_r = c_CYL[0,s]
        this_t = c_CYL[1,s]
        this_z = c_CYL[2,s]

        this_kr = kr
        ;this_kt = Theta_Toroidal_Sign * nPhi/this_r
        this_kt = nPhi/this_r
        this_kz = kz

        this_k_CYL = [this_kr,this_kt,this_kz]
        this_k_XYZ = vector_CYL_to_XYZ(c_CYL[*,s],this_k_CYL)

        this_x = (this_r-rMin)/(rMax-rMin)*2*!pi
        this_y = (this_z-zMin)/(zMax-zMin)*2*!pi

        kb[s] = this_k_XYZ[0]*b_XYZ[0,s]/bMag[s] + $
                this_k_XYZ[1]*b_XYZ[1,s]/bMag[s] + $
                this_k_XYZ[2]*b_XYZ[2,s]/bMag[s]

        kDotR_AlongS[s] =  $
                 c_CYL[0,s]*this_k_CYL[0] $
                +c_CYL[1,s]*this_k_CYL[1]*c_CYL[0,s]$
                +c_CYL[2,s]*this_k_CYL[2]

        ;this_E = exp(ii*n*this_x)*exp(ii*m*this_y)*exp(ii*nPhi*this_t)
        ;this_E = exp(ii*(n*this_x+m*this_y+nPhi*this_t))
        ;this_E = exp(ii*(c_CYL[0,s]*this_kr+c_CYL[1,s]*nPhi+c_CYL[2,s]*this_kz))
        ;this_E = exp(ii*( $
        ;         c_CYL[0,s]*this_k_CYL[0] $
        ;        +c_CYL[1,s]*this_k_CYL[1]*c_CYL[0,s]$
        ;        +c_CYL[2,s]*this_k_CYL[2] $
        ;        ))
        this_E = exp(ii*kDotR_AlongS[s])

        this_bu_CYL = b_CYL[*,s]/bMag_CYL[s]

        this_Er = this_bu_CYL[0] * this_E 
        this_Et = this_bu_CYL[1] * this_E 
        this_Ez = this_bu_CYL[2] * this_E 

        this_E_CYL = [this_Er,this_Et,this_Ez]

        Eb_check[s] = this_E
        Eb[s] = this_bu_CYL[0]*this_E_CYL[0] $
           + this_bu_CYL[1]*this_E_CYL[1] $
           + this_bu_CYL[2]*this_E_CYL[2]  

    endfor


	SPoints_sig33 = ComplexArr(SPointsN)
	SPoints_sig33_cold = ComplexArr(SPointsN)

	for Pt=0,SPointsN-1 do begin

    	w    = 2*!pi*f_Hz
    	E_eV = E_keV * 1e3
    	b0   = bMag_CYL[SPointsIndex[Pt]]
    	vTh  = sqrt(2d0*E_eV*_e/me)
    	wpe  = sqrt(n_e*_e^2/(me*e0))
    	wce  = abs(AtomicZ)*_e*b0/me
    	kPar = kb[SPointsIndex[Pt]]
    	vPhs = w/kPar
    	lambdaPar = 2*!Pi/kPar
		print, 'vTh distance: ', vTh * (1/f_Hz)

    	kPer = 0
    	lambda=0
    	In = 1
    	sum = !null
    	for l=0,0 do begin
    		zeta_n=(w-l*wce)/(kPar*vTh)

    	    print, 'Z function argument for Mathematica: ', zeta_n

			Z_n = kj_zfunction(zeta_n,Zp=Zp_n)

    		if sum then begin
    			sum = sum + In*zeta_n*Zp_n
    		endif else begin
    			sum = In*zeta_n*Zp_n
    		endelse
    	endfor

    	K3 = 1d0 - wpe^2 * exp(-lambda) / (w*kPar*vTh) * sum 
    	II = complex(0,1)
    	SPoints_sig33[Pt] = -(K3 - 1d0)*II*w*e0
    	
    	stixP = 1-wpe^2/w^2
    	SPoints_sig33_cold[Pt] = -(stixP-1d0)*II*w*e0

	endfor

	save, SPoints, SPoints_sig33, SPoints_sig33_cold, FileName='AnalyticSig33.sav'

   	print, 'Analytic Sig33(k): ',SPoints_sig33
   	print, 'Analytic Sig33(cold)', SPoints_sig33_cold


    p = plot(s_Coord, BMag, layout=[1,4,1], thick = 2, color='b', title='bMag Along S')
    p = plot(s_Coord, real_part(Eb), layout=[1,4,2], /current, thick=2, title='ePar Along S',xRange=[-3,5])
    p = plot(s_Coord, imaginary(Eb), layout=[1,4,2], /over, thick=2, color='r', transparency=50)

    p = plot(s_Coord, real_part(Eb_check), layout=[1,4,2], /over, thick=4, title='ePar Along S', transparency=70, lineStyle='--')
    p = plot(s_Coord, imaginary(Eb_check), layout=[1,4,2], /over, thick=4, color='r', transparency=70, lineStyle='--')

    p = plot(s_Coord, kDotR_AlongS, layout=[1,4,3], /current, thick=2, color='purple',title='kDotR Along S')
    p = plot(s_Coord, kb, layout=[1,4,4], /current, thick=2, color='g', title='kPar Along S',xrange=[-3,5])

    ; Write the input file

    nc_id = nCdf_create ( FieldsOutFileName, /clobber )

	nCdf_control, nc_id, /fill
	
	nr_id = nCdf_dimDef ( nc_id, 'nR', nS )
	scalar_id = nCdf_dimDef ( nc_id, 'scalar', 1 )

	freq_id = nCdf_varDef ( nc_id, 'freq', scalar_id, /float )
	r_id = nCdf_varDef ( nc_id, 'r', nr_id, /float )

	B0_r_id = nCdf_varDef ( nc_id, 'B0_r', nr_id, /float )
	B0_p_id = nCdf_varDef ( nc_id, 'B0_p', nr_id, /float )
	B0_z_id = nCdf_varDef ( nc_id, 'B0_z', nr_id, /float )

	e_r_re_id = nCdf_varDef ( nc_id, 'e_r_re', nr_id, /float )
	e_r_im_id = nCdf_varDef ( nc_id, 'e_r_im', nr_id, /float )
	e_p_re_id = nCdf_varDef ( nc_id, 'e_p_re', nr_id, /float )
	e_p_im_id = nCdf_varDef ( nc_id, 'e_p_im', nr_id, /float )
	e_z_re_id = nCdf_varDef ( nc_id, 'e_z_re', nr_id, /float )
	e_z_im_id = nCdf_varDef ( nc_id, 'e_z_im', nr_id, /float )

	jP_r_re_id = nCdf_varDef ( nc_id, 'jP_r_re', nr_id, /float )
	jP_r_im_id = nCdf_varDef ( nc_id, 'jP_r_im', nr_id, /float )
	jP_p_re_id = nCdf_varDef ( nc_id, 'jP_p_re', nr_id, /float )
	jP_p_im_id = nCdf_varDef ( nc_id, 'jP_p_im', nr_id, /float )
	jP_z_re_id = nCdf_varDef ( nc_id, 'jP_z_re', nr_id, /float )
	jP_z_im_id = nCdf_varDef ( nc_id, 'jP_z_im', nr_id, /float )

	nCdf_control, nc_id, /enDef

	nCdf_varPut, nc_id, freq_id, f_Hz

	nCdf_varPut, nc_id, r_id, s_Coord 

	nCdf_varPut, nc_id, B0_r_id, bMag
	nCdf_varPut, nc_id, B0_p_id, bMag*0
	nCdf_varPut, nc_id, B0_z_id, bMag*0

	nCdf_varPut, nc_id, e_r_re_id,real_part(Eb)
	nCdf_varPut, nc_id, e_r_im_id,imaginary(Eb)
	nCdf_varPut, nc_id, e_p_re_id,real_part(Eb)*0
	nCdf_varPut, nc_id, e_p_im_id,imaginary(Eb)*0
	nCdf_varPut, nc_id, e_z_re_id,real_part(Eb)*0
	nCdf_varPut, nc_id, e_z_im_id,imaginary(Eb)*0

	nCdf_varPut, nc_id, jP_r_re_id,real_part(Eb)*0
	nCdf_varPut, nc_id, jP_r_im_id,imaginary(Eb)*0 
	nCdf_varPut, nc_id, jP_p_re_id,real_part(Eb)*0 
	nCdf_varPut, nc_id, jP_p_im_id,imaginary(Eb)*0 
	nCdf_varPut, nc_id, jP_z_re_id,real_part(Eb)*0 
	nCdf_varPut, nc_id, jP_z_im_id,imaginary(Eb)*0 

nCdf_close, nc_id


    create_test_particle_f, $
            /weighted_maxwellian_XYZ, $
            energy_keV = E_keV, $
            density_m3 = n_e, $
            rsfwc_1d = FieldsOutFileName, $
            OutputFileName = ParticlesOutFileName, $
            n_particles = Np


    stop

end
