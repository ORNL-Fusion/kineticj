pro kj_rs_run_residuals

    compile_opt idl2, logical_predicate

    ; Get list of runs to do

    runList = file_search('./','2018*',/fully_qualify_path)
    
    N = n_elements(runList)

    cd, current=root

    n_proc = 24 

    br = objarr(n_proc) 

    print, 'Setting up IDL parallel pool ...'

    for i=0,n_proc-1 do begin

        br[i] = IDL_IDLbridge()
        br[i]->SetVar, '!path', !path

    endfor

    print, 'DONE'

    print, 'Doing work ...'

    pos = 0 
    while pos lt N do begin  

        for i=0,n_proc-1 do begin

            if (br[i]->status() eq 0) then begin

                br[i]->SetVar, 'dir', runList[pos]
                br[i]->Execute, 'cd, dir'
                br[i]->Execute, 'res=kj_rs_residual()',/nowait

                ++pos

                print, pos, N

            endif

        endfor

    endwhile

end
