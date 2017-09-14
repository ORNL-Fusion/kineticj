function kj_anderson, g, x, mMax, itmax, atol, rtol, droptol, beta, AAstart

    compile_opt strictarr

    nargin = n_params()

    if nargin lt 2 then stop
    if nargin lt 3 then mMax = min([10,size(x)])
    if nargin lt 4 then itmax = 100
    if nargin lt 5 then atol = 1.0e-10
    if nargin lt 6 then rtol = 1.0e-10
    if nargin lt 7 then droptol = 1e10
    if nargin lt 8 then beta = 1
    if nargin lt 9 then AAstart = 0

    ; Initialize the storage arrays

    res_hist = !null ; Residual history
    DG = !null ; g-value differences

    ; Initialize printing

    if mMax eq 0 then begin
        print, 'No acceleration'
    endif else if mMax gt 0 then begin
        print, 'Anderson accleration, mMax = ' + strTrim(string(mMax),2)
    endif else begin
        print, 'mMax must be non negative'
        stop
    endelse

    ; Initialize the number of stored residuals.

    mAA = 0

    ; Top of the iteration loop.

    for iter = 0, itmax

        ; Apply g and compute the current residual norm.

        gval = g(x)
        fval = gval - x
        res_norm = norm(fval)
        print, iter, res_norm
        res_hist = [res_hist,[iter,res_norm]]

        ; Set the residual tolerance on the initial iteration.

        if iter eq 0 then tol = max([atol,rtol*res_norm]) 

        ; Test for stopping.

        if res_norm le tol then begin

            print, 'Terminate with residual norm = '+string(res_norm)
            break

        endif

        if mMax eq 0 or iter gt AA start then begin

            ; Without acceleration, update x <- g(x) to obtain the next
            ; approximate solution.

            x = gval

        endif else begin

            ; Apply Anderson acceleration.

            ; Update the df vector and the DG array.

            if iter gt AAstart then begin

                df = fval - f_old

                if mAA < mMax then begin

                    DG = [DG, gval-g _old]

                endif else begin
                    
                    DG = [DG[*,1:mAA-1], gval-g_old]
                        
                endelse

                mAA = mAA + 1

            endif

            f_old = fval
            g_old = gval

            if mAA eq 0 then begin

                ; If mAA == 0, update x <- g(x) to obtain next
                ; approximate solution

                x = gval

            endif else begin

                ; If mAA > 0, solve the least-squares problem
                ; update the solution.

                if mAA eq 1 then begin

                    ; If mAA == 1, form the initial QR decomposition

                    R = FltArr(1,1) + norm(df)
                    Q = df / R[0,0]

                endif else begin

                    ; If mAA > 1, update the QR decomposition

                    if mAA gt mMax then begin

                        ; If the column dimension of Q is mMax,
                        ; delete the first column and update
                        ; the decomposition.



                    endif

                endelse

            endelse

        endif


    endfor ; End of the iter loop

    return, g(x)

end
