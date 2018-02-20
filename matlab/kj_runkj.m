function [jP_r,jP_t,jP_z] = kj_runkj(Er,Et,Ez,uid)

useAR = 1;

if ~exist('uid','var') || isempty(uid)
  uid = 0;   
end
uidStr = num2str(uid,'%3.3i'); 

% Assuming an existing RS run with an existing zeroed delta file
% and that all the RS runs have a KJ input, just with the first
% one have zero value. Use the delta as x in x=g(x). 

% Get the run directory

rootDir = pwd();

% Set the template directory

if useAR
    templateDir = 'template-ar';
else
    templateDir = 'template-rs';
end

% Stage iteration

%thisDir = char(strcat('run',sprintf('%4.4i',it)));

now = datetime('now');
thisDir = strcat( datestr(now,'yyyy-mm-dd-HH-MM-SS'),'-',uidStr);

copyfile( templateDir, thisDir );

% Change to run directory

cd(thisDir);

% Update the input E field to KJ


files = dir('output/solution*.nc');

[M,N] = size(files);

if (M > 1)
    print, 'Too many solution files, please cleanup';
    exit
end

eFieldFile = strcat('output/',files(1).name());

ncwrite(eFieldFile,'E_r_re',real(Er));
ncwrite(eFieldFile,'E_r_im',imag(Er));
ncwrite(eFieldFile,'E_t_re',real(Et));
ncwrite(eFieldFile,'E_t_im',imag(Et));
ncwrite(eFieldFile,'E_z_re',real(Ez));
ncwrite(eFieldFile,'E_z_im',imag(Ez));

% Run KJ 

%!IDL_STARTUP="/Users/dg6/idlStartup.pro" /usr/local/bin/idl run_kj &> kj.idl.log
!IDL_STARTUP="/Users/dg6/idlStartup.pro" /usr/local/bin/idl run_kj

% Read Jp from file

kjFile = 'kj-jp.nc';
kj = dlg_read_netcdf(kjFile);

jP_r = kj('jP_r');
jP_t = kj('jP_t');
jP_z = kj('jP_z');

cd(rootDir);

end