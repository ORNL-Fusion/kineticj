pro kj_iterate_anderson

    x = fIndgen(10)
    g = kj_update()

    x2 = kj_anderson(g,x)

    stop
    
end
