pro kj_create_upshift_input

    f_Hz = 30e6
    bPolFactor = 1.0
    EqdskFile = 'g130608.00355.EFIT02.mds.corrected.qscale_1.00000'

    rtz = [1.5,0.0,0.0]
    dS = 0.01
    dir = 1
    TraceNPts = 1000

    g = ReadGEqdsk ( EqdskFile, $
            bPolFactor = bPolFactor, $
            fieldLineIn = rtz, $
            FieldLine_XYZ = FieldLine_XYZ, $
            B_AlongFieldLine_XYZ = B_AlongFieldLine_XYZ, $
            B_AlongFieldLine_CYL = B_AlongFieldLine_CYL, $
            SafetyFactor = SafetyFactor, $
            FieldLineTraceDir = 1, $
            FieldLineTraceDS = dS, $
            FieldLineTraceNSteps = TraceNPts ) 

    g = ReadGEqdsk ( EqdskFile, $
            bPolFactor = bPolFactor, $
            fieldLineIn = rtz, $
            FieldLine_XYZ = FieldLine_XYZ_, $
            B_AlongFieldLine_XYZ = B_AlongFieldLine_XYZ_, $
            B_AlongFieldLine_CYL = B_AlongFieldLine_CYL_, $
            SafetyFactor = SafetyFactor_, $
            FieldLineTraceDir = -1, $
            FieldLineTraceDS = dS, $
            FieldLineTraceNSteps = TraceNPts ) 

    B_XYZ = [ [reverse(B_AlongFieldLine_XYZ[*,0:TraceNPts-1],2)], [B_AlongFieldLine_XYZ_[*,1:TraceNPts-1]] ]
    c_XYZ = [ [reverse(FieldLine_XYZ[*,0:TraceNPts-1],2)], [FieldLine_XYZ_[*,1:TraceNPts-1]] ]

    BMag = sqrt(total(B_XYZ^2,1))

    print, 'Safety Factor: ', SafetyFactor, SafetyFactor_

    p = plot3d((c_XYZ[0,*])[*],(c_XYZ[1,*])[*],(c_XYZ[2,*])[*],$
            thick = 3.0, $
            aspect_ratio = 1.0, aspect_z = 1.0, $
            rgb_table = 7, vert_colors = BytScl(BMag, top=253)+1,$
            depth_cue = [0,2], /perspective )

    nFrames = 48 
    for n = 0, nFrames*2/3-1 do begin
            thisR = g.rLim
            thisZ = g.zLim
            thisT = 360./(nFrames-1)*n+fltArr(n_elements(thisR))
            Frame_Cyl = transpose([[thisT[*]],[thisR[*]],[thisZ[*]]])
            Frame_Rec = cv_coord(from_cylin=Frame_Cyl,/degrees,/to_rect)
            pF = plot3d((Frame_Rec[0,*])[*],(Frame_Rec[1,*])[*],(Frame_Rec[2,*])[*],/over,thick=5,transparency=70)
    endfor

    p = plot(BMag)

stop
    ; Write the input file

    nc_id = nCdf_create ('kj_upshift.nc', /clobber )

	nCdf_control, nc_id, /fill
	
	nr_id = nCdf_dimDef ( nc_id, 'nR', n_elements(x) )
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

	nCdf_varPut, nc_id, r_id, x

	nCdf_varPut, nc_id, B0_r_id, br
	nCdf_varPut, nc_id, B0_p_id, bt 
	nCdf_varPut, nc_id, B0_z_id, bz

	nCdf_varPut, nc_id, e_r_re_id,real_part(Er) 
	nCdf_varPut, nc_id, e_r_im_id,imaginary(Er) 
	nCdf_varPut, nc_id, e_p_re_id,real_part(Et)
	nCdf_varPut, nc_id, e_p_im_id,imaginary(Et)
	nCdf_varPut, nc_id, e_z_re_id,real_part(Ez)
	nCdf_varPut, nc_id, e_z_im_id,imaginary(Ez)

	nCdf_varPut, nc_id, jP_r_re_id,real_part(Jpr) 
	nCdf_varPut, nc_id, jP_r_im_id,imaginary(Jpr) 
	nCdf_varPut, nc_id, jP_p_re_id,real_part(Jpt) 
	nCdf_varPut, nc_id, jP_p_im_id,imaginary(Jpt) 
	nCdf_varPut, nc_id, jP_z_re_id,real_part(Jpz) 
	nCdf_varPut, nc_id, jP_z_im_id,imaginary(Jpz) 

nCdf_close, nc_id



    stop

end
