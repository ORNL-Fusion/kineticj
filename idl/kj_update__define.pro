function kj_update::_overloadFunction, x
    compile_opt strictarr

    g = x*2

    return, g

end


pro kj_update__define
    compile_opt strictarr

    !null = { kj_update, inherits IDL_Object }

end
