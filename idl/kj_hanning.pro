function kj_hanning, N, symmetric=symmetric

    k = fIndGen(N)

    if keyword_set(symmetric) then begin
        f = 0.5-(1-0.5)*cos(2*!pi*k/(N-1))
    endif else begin
        f = 0.5-(1-0.5)*cos(2*!pi*k/(N))
    endelse

    return, f

end
