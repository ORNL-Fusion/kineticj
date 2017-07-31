# Kinetic-j

## Description
The Kinetic-j program takes a time harmonic electric wave field, and calculates the kinetic plasma current. This calculation is done in configuration-space, as opposed to the more traditional Fourier-space evaluation forms provided in the standard plasma wave physics texts. 

## Program Files
The program is organized into C++ and IDL files. The C++ are the buildable program (located in src and include directories), with the IDL scripts (located in the idl directory) being used to generate inputs, run benchmarks, and post process outputs. The bin and obj directories are empty, being only placeholders for build locations. The mathematica directory containts some relevant analysis and derivations of the equations upon which the program is based. 
```
. (Makefile etc)
├── bin (empty)
├── idl (IDL scripts)
├── include (C++ header files)
├── machine-makefiles (machine specific makefile setting files)
├── mathematica (mathematica scripts)
├── obj (empty)
└── src (C++ source files)
```

## Dependencies 

### Libconfig (https://github.com/hyperrealm/libconfig)
The human readable input files are of the format defined by the libconfig API. Libconfig can be found at the above location on github.

### NetCDF-CXX4 (https://github.com/Unidata/netcdf-cxx4)
We utilize the netcdf file format for both inputs and outputs. We further utilize the more recent CXX4 components of the NetCDF API, which meand a netcdf installation built with the --enable-cxx4 flag. The above official netcdf-cxx4 github provides further information.  

### CUDA (https://developer.nvidia.com/cuda-downloads)
To enable the GPU capability the CUDA API is required. This will likely be already available on appropriate machines, if not, see the link above. Note that the code is designed to be compiled using the THRUST API wihtin the CUDA SDK, so even if you do not have a GPU, you will still need the CUDA API, just using the CPU or OPENMP THRUST targets. 

## Installation
We utilize a simple machine specific makefile to specify the locations of the above dependencies. To build Kinetic-j on a new machine you will need to create a copy of an existing file (e.g., `machine-makefiles/Makefile.dlg-macpro`) with the name of your machine as the extension. This can be done via 

`cp machine-makefiles/Makefile.dlg-macpro machine-makefiles/$(uname -n)`

and then editing the resulting file appropriately. Then just 

```
make clean
make
```
### Build Options
Note: The code is designed to be built using CUDA, even if you do not have a GPU, i.e., the THRUST API within the CUDA SDK provides for single thread CPU (CPP), multi thread OpenMP (OMP), and CUDA (CUDA) targets via the appropriate choice of THRUST policy.

Within the top level `Makefile` the THRUST policy (which selects either CPU, OPENMP, or GPU) is selected by uncommenting the desired choice as 

```
THRUST_POLICY:=THRUST_DEVICE_SYSTEM_CUDA
#THRUST_POLICY:=THRUST_DEVICE_SYSTEM_CPP
#THRUST_POLICY:=THRUST_DEVICE_SYSTEM_OMP
 ```

## Specific Machine Build Notes

### gpufusion.ornl.gov

```
source env-gpufusion.sh
make clean
make
```

## Other Information

### Generate the Guiding Center terms file
To speed up the guiding center orbit calculation (if selected), some of the terms in the ODE are precalculated and tabulated in a file. This is done using code from the https://github.com/dlg0/orbit_tracer repo. After cloning that repo, run the `gc_terms.pro` routine on an `ar2Input.nc` file, e.g., 

```
cd ~/scratch/aorsa2d/colestock-kashuba-reference
idl
IDL>@constants
IDL>gc_terms, _me_amu, -1, ar2='input/ar2Input.nc', /oneD
```
