function [res] = kj_readResidual(loc)

rootDir = pwd();

cd(char(loc));

resFile = 'output/kj-rs-res.nc';

while exist(resFile,'file')==0
    
    % Wait for file to be created
    
    disp('Waiting for residual file ...');
    disp(loc);
    pause(1);
    
end


% Read LHS from file

res = dlg_read_netcdf(resFile);

res_r = res('res_r');
res_t = res('res_t');
res_z = res('res_z');

res = [res_r',res_t',res_z'];


cd(rootDir);

return, res

end