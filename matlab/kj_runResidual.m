function [res] = kj_runResidual(Er,Et,Ez,uid)

if ~exist('uid','var') || isempty(uid)
  uid = 0;   
end
uidStr = num2str(uid,'%3.3i'); 

% Get the run directory

rootDir = pwd();

% Set the template directory

templateDir = 'template';

% Stage iteration

%thisDir = char(strcat('run',sprintf('%4.4i',it)));

now = datetime('now');
thisDir = strcat( datestr(now,'yyyy-mm-dd-HH-MM-SS'),'-',uidStr);

copyfile( templateDir, thisDir );

% Change to run directory

cd(thisDir);

% Update the input E field to KJ

file = dir('output/rs-solution.nc');

eFieldFile = strcat('output/',file.name());

ncwrite(eFieldFile,'E_r_re',real(Er));
ncwrite(eFieldFile,'E_r_im',imag(Er));
ncwrite(eFieldFile,'E_t_re',real(Et));
ncwrite(eFieldFile,'E_t_im',imag(Et));
ncwrite(eFieldFile,'E_z_re',real(Ez));
ncwrite(eFieldFile,'E_z_im',imag(Ez));

% Run IDL routine rs_lhs 

!IDL_STARTUP="/Users/dg6/idlStartup.pro" /usr/local/bin/idl run_kj_rs_residual

% Read LHS from file

resFile = 'output/kj-rs-res.nc';
res = dlg_read_netcdf(resFile);

jP_r = kj('jP_r');
jP_t = kj('jP_t');
jP_z = kj('jP_z');

cd(rootDir);


end