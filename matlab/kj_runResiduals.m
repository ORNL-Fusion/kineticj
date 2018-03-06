function [s] = kj_runResiduals()
    
    s=system('IDL_STARTUP="/Users/dg6/idlStartup.pro" /usr/local/bin/idl -quiet run_kj_rs_run_residuals');

    return, s
end