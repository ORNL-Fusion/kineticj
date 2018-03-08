function [thisDir] = kj_stageResidual(E,uid)

if ~exist('uid','var') || isempty(uid)
  uidStr = string(java.util.UUID.randomUUID);   
end


%uidStr = num2str(uid,'%3.3i');

% Get the run directory

rootDir = pwd();

% Set the template directory

templateDir = 'template';

% Stage iteration

%thisDir = char(strcat('run',sprintf('%4.4i',it)));

now = datetime('now');
thisDir = strcat( datestr(now,'yyyy-mm-dd-HH-MM-SS'),'-',uidStr);

copyfile( templateDir, char(thisDir) );

% Change to run directory

cd(char(thisDir));

% Update the input E field to KJ

[M,N] = size(E);
n = M/3;

Er = E(0*n+1:1*n);
Et = E(1*n+1:2*n);
Ez = E(2*n+1:3*n);

eFieldFile = 'output/rs-solution.nc';

ncwrite(eFieldFile,'E_r_re',real(Er));
ncwrite(eFieldFile,'E_r_im',imag(Er));
ncwrite(eFieldFile,'E_t_re',real(Et));
ncwrite(eFieldFile,'E_t_im',imag(Et));
ncwrite(eFieldFile,'E_z_re',real(Ez));
ncwrite(eFieldFile,'E_z_im',imag(Ez));

% Run IDL routine rs_lhs 

% s=system('IDL_STARTUP="/Users/dg6/idlStartup.pro" /usr/local/bin/idl -quiet run_kj_rs_residual &');

% %Read LHS from file
% 
% resFile = 'output/kj-rs-res.nc';
% res = dlg_read_netcdf(resFile);
% 
% res_r = res('res_r');
% res_t = res('res_t');
% res_z = res('res_z');
% 
% res = [res_r',res_t',res_z'];

cd(rootDir);


end