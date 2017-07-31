# Kinetic-j

## Description
The Kinetic-j program takes a time harmonic electric wave field, and calculates the kinetic plasma current. This calculation is done in configuration-space, as opposed to the more traditional Fourier-space evaluation forms provided in the standard plasma wave physics texts. 

## Program Files
'''
.
├── bin
├── idl
├── include
├── machine-makefiles
├── obj
└── src
'''

## Dependencies 

### Libconfig (https://github.com/hyperrealm/libconfig)
The human readable input files are of the format defined by the libconfig API. Libconfig can be found at the above location on github.

### NetCDF-CXX4 (https://github.com/Unidata/netcdf-cxx4)
We utilize the netcdf file format for both inputs and outputs. We further utilize the more recent CXX4 components of the NetCDF API, which meand a netcdf installation built with the --enable-cxx4 flag. The above official netcdf-cxx4 github provides further information.  

### CUDA (https://developer.nvidia.com/cuda-downloads)
To enable the GPU capability the CUDA API is required. This will likely be already available on appropriate machines, if not, see the link above. 

## Installation
We utilize a simple machine specific makefile to specify the locations of the To build Kinetic-j on a machine wnew machineThe "machine-makefiles" directory contains machine specific 


### Build on gpufusion.ornl.gov
```
source env-gpufusion.sh
make clean
make
```

## Other Information

### Calculate the Guiding Center terms file

From the https://github.com/dlg0/orbit_tracer repo, run the `gc_terms.pro` routine on an `ar2Input.nc` file ...
```
cd ~/scratch/aorsa2d/colestock-kashuba-reference
idl
IDL>@constants
IDL>gc_terms, _me_amu, -1, ar2='input/ar2Input.nc', /oneD
```
