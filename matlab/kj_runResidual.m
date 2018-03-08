function [s] = kj_runResidual(loc)
    
    rootDir = pwd();
    
    cd(char(loc));

    s=system('IDL_STARTUP="/Users/dg6/idlStartup.pro" /usr/local/bin/idl -quiet run_kj_rs_residual');

    cd(rootDir);
    
    return, s
end