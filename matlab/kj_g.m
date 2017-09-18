function [x] = kj_g (x,it)

if ~exist('it','var') || isempty(it)
  it=1;
end

% Assuming an existing RS run with an existing zeroed delta file
% and that all the RS runs have a KJ input, just with the first
% one have zero value. Use the delta as x in x=g(x). 

% Get the run directory

rootDir = pwd();

% Set the template directory

templateDir = 'template-rs';

% Stage iteration

thisDir = char(strcat('run',string(it)))

copyfile( templateDir, thisDir );

% Change to run directory

cd(thisDir);

% Update the delta file with the new x values

deltaFile = 'kj-delta.nc';

[jr,jt,jz] = kj_x_to_vec(x);

ncwrite(deltaFile,'jP_r_re',real(jr));
ncwrite(deltaFile,'jP_r_im',imag(jr));
ncwrite(deltaFile,'jP_t_re',real(jt));
ncwrite(deltaFile,'jP_t_im',imag(jt));
ncwrite(deltaFile,'jP_z_re',real(jz));
ncwrite(deltaFile,'jP_z_im',imag(jz));

% Run RS with IDL

!IDL_STARTUP="/Users/dg6/idlStartup.pro" /usr/local/bin/idl run_rs

% Run KJ to get the new delta

!IDL_STARTUP="/Users/dg6/idlStartup.pro" /usr/local/bin/idl run_kj

% Read new delta from file

deltaFile = 'kj-stix.nc';

jr_re = ncread(deltaFile,'jP_r_re');
jr_im = ncread(deltaFile,'jP_r_im');
jt_re = ncread(deltaFile,'jP_t_re');
jt_im = ncread(deltaFile,'jP_t_im');
jz_re = ncread(deltaFile,'jP_z_re');
jz_im = ncread(deltaFile,'jP_z_im');

jr = complex(jr_re,jr_im);
jt = complex(jt_re,jt_im);
jz = complex(jz_re,jz_im);

x = kj_vec_to_x(jr,jt,jz);

cd(rootDir);

end