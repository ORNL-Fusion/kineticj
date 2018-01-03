function [jr,jt,jz] = kj_runkj(Er,Et,Ez)

useAR = 0;

if ~exist('it','var') || isempty(it)
  it=1;
end

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

thisDir = char(strcat('run',sprintf('%4.4i',it)));

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

eFieldFile = strcat('output/',files(1));

ncwrite(eFieldFile,'er_re',real(Er));
ncwrite(eFieldFile,'er_im',imag(Er));
ncwrite(eFieldFile,'et_re',real(Et));
ncwrite(eFieldFile,'et_im',imag(Et));
ncwrite(eFieldFile,'ez_re',real(Ez));
ncwrite(eFieldFile,'ez_im',imag(Ez));

% Run KJ 

%!IDL_STARTUP="/Users/dg6/idlStartup.pro" /usr/local/bin/idl run_kj &> kj.idl.log
!IDL_STARTUP="/Users/dg6/idlStartup.pro" /usr/local/bin/idl run_kj

% Read Jp from file

% Read new delta from file

kj_JpFile = 'kj-jp.nc';

jr_re = ncread(deltaFile,'jP_r_re');
jr_im = ncread(deltaFile,'jP_r_im');
jt_re = ncread(deltaFile,'jP_t_re');
jt_im = ncread(deltaFile,'jP_t_im');
jz_re = ncread(deltaFile,'jP_z_re');
jz_im = ncread(deltaFile,'jP_z_im');

jr = complex(jr_re,jr_im);
jt = complex(jt_re,jt_im);
jz = complex(jz_re,jz_im);

cd(rootDir);

end

end