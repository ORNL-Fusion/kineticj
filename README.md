# Kinetic-j
NOTE: This file is best viewed at https://github.com/ORNL-Fusion/kineticj/

## Description
The Kinetic-j program takes a time harmonic electric wave field, and calculates the kinetic plasma current. This calculation is done in configuration-space, as opposed to the more traditional Fourier-space evaluation forms provided in the standard plasma wave physics texts. 

## Program Files
The program is organized into C++ and IDL files. The C++ are the buildable program (located in src and include directories), with the IDL scripts (located in the idl directory) being used to generate inputs, run benchmarks, and post process outputs. The bin and obj directories are empty, being only placeholders for build locations. The mathematica directory contains some relevant analysis and derivations of the equations upon which the program is based. 
```
. (Makefile etc)
├── benchmarks (template and output for the sigma benchmarking)
├── bin (empty)
├── idl (IDL utility scripts)
├── include (C++ header files)
├── machine-makefiles (machine specific makefile setting files)
├── mathematica (mathematica scripts)
├── matlab (matlab utility scripts)
├── obj (empty)
├── python (python scripts)
├── src (C++ source files)
├── template (single test case)
└── tests (regression test files)
```

## Dependencies 

### Libconfig (https://github.com/hyperrealm/libconfig)
The human readable input files are of the format defined by the libconfig API. Libconfig can be found at the above location on github.

### NetCDF-C (https://github.com/Unidata/netcdf-c)
We utilize the netcdf file format for both inputs and outputs.

### NetCDF-CXX4 (https://github.com/Unidata/netcdf-cxx4)
We further utilize the more recent CXX-4 components of the NetCDF API (note that the --enable-cxx4 flag is not required as suggested in the netcdf-cxx4 github readme file).  

### THRUST (https://github.com/thrust/thrust) or CUDA (https://developer.nvidia.com/cuda-downloads)
The THRUST headers are required and are available either standalone from the THRUST repo (for building on systems without GPUs), or from within the CUDA SDK (for building on systems with GPUs). So if you do not have an nvidia compute GPU, you will still need the THRUST library (unbuilt).

## Installation
We utilize a simple machine specific makefile to specify the locations of the above dependencies. To build Kinetic-j on a new machine you will need to create a copy of an existing file (e.g., `machine-makefiles/Makefile.dlg-macpro`) with the name of your machine as the extension. This can be done via 

`cp machine-makefiles/Makefile.dlg-macpro machine-makefiles/$(uname -n)`

and then editing the resulting file appropriately. Then just 

```
make clean
make
```
### Build Options
The THRUST API provides for single thread CPU (CPP), multi thread OpenMP (OMP), and CUDA (CUDA) targets via the appropriate choice of THRUST policy. Within the top level `Makefile` the THRUST policy (which selects either CPU, OPENMP, or GPU) is selected by uncommenting the desired choice as 

```
# Select only one.

USE_SERIAL:= 0
USE_CUDA  := 0
USE_OPENMP:= 1
 ```
Also, various code features and debugging levels can be enabled at build time by enabling any of the following within the top level `Makefile`. 

```
# Features
CPPFLAGS += -DCYLINDRICAL_INPUT_FIELDS=1 # Are the input (NetCDF) fields in Cylindrical or not. 
CPPFLAGS += -D_PARTICLE_BOUNDARY=1 # 1 = particle absorbing walls, 2 = periodic, 3 = reflective
CPPFLAGS += -DGC_ORBITS=0 # Use Guiding Center evalulation (Not yet functional)
CPPFLAGS += -DCOMPLEX_WRF=0 # Use the complex omega as the time integral decay time, as opposed to the Hanning window.  

# Non THRUST Implementations
CPPFLAGS += -DDO_CPU_ITERATOR_APPROACH=0
CPPFLAGS += -DDO_CPU_APPROACH=0

# Timing Info
CPPFLAGS += -DUSEPAPI=0
CPPFLAGS += -DLOWMEM_USEPAPI=0
CPPFLAGS += -DCLOCK=1
CPPFLAGS += -DPRINT_INFO=1

# Debug Info Level
CPPFLAGS += -DDEBUG_GC=0
CPPFLAGS += -DDEBUG_EVAL_VGC=0
CPPFLAGS += -DDEBUG_EVAL_APAR=0
CPPFLAGS += -DDEBUGLEVEL=0
CPPFLAGS += -DDEBUG_LINES=0
CPPFLAGS += -DDEBUG_INTERP=0
CPPFLAGS += -DDEBUG_MAXWELLIAN=0
CPPFLAGS += -DDEBUG_FORCE_TERM=0
CPPFLAGS += -DDEBUG_MOVE=0
CPPFLAGS += -DDEBUG_ROTATION=0
CPPFLAGS += -DDEBUG_INTVECARRAY=0
CPPFLAGS += -DDEBUG_READ_E_FIELD=0

# File Write Options
CPPFLAGS += -DLOWMEM_ORBIT_WRITE=0
CPPFLAGS += -DF1_WRITE=0
```

## Specific Machine Build Notes

### edison.nersc.go
```
source env-edison.sh
make clean
make
```

### gpufusion.ornl.gov

```
source env-gpufusion.sh
make clean
make
```

## Running Kinetic-J
Set a `KINETICJ` environment variable to be the location of the cloned source ...
```
cd ~/code
git clone https://github.com/ORNL-Fusion/kineticj.git
export KINETICJ=~/code/kineticj
```

### Run the regression tests
To aid development testing we include some simple regression testing. This is run as follows ...
```
cd $KINETICJ/tests
python ../python/kj_test.py
```
or at NERSC ...
```
cd $KINETICJ
source env-edison.sh
cd $SCRATCH
mkdir kineticj
cp -r $KINETICJ/tests .
cd tests
salloc -N 1 -p debug
python $KINETICJ/python/kj_plot.py
exit
```
with the expected output being ...
```
dlg-macbookpro2:tests dg6$ python ../python/kj_test.py
benchmark1-00007    PASS
benchmark2-00013    PASS
benchmark3-00004    PASS
test4               PASS
```

### Run the test case
A standalone test case is also availble where we demonstrate the current response for a variable magnetic field (1/r) where the fundamental and 2nd harmonic ion cyclotron resonances are in the domain (at x=1.75 and x=3.5 respectively). We have provided an input file `template/input/input-data.nc` with an electric field with `kx=pi/(2*rho_L)` where `rho_L` is the Lamor radius at the location of the 2nd harmonic resonance, i.e., the finite Lamor radius interaction is captured per pg.270-271 of Stix. This case is run via the following commands ...

```
cd $KINETICJ
cd template
./bin/kineticj
python ../python/kj_plot.py
```

or for NERSC (Edison with a single node using OpenMP)

```
cd $KINETICJ
source env-edison.sh
cd $SCRATCH
mkdir kineticj
cp -r $KINETICJ/template .
cd template
salloc -N 1 -p debug
$KINETICJ/bin/kineticj
exit
python $KINETICJ/python/kj_plot.py
```
<img src="/template/template.png" width="300">

Changing the variables in `template/kj.cfg` will allow experimenting with running the code. 

```
  1 xGridMin = 1.0;
  2 xGridMax = 5.0;
  3 nXGrid = 200
  4 nRFCycles = 50.0;
  5 species_number = 0;
  6 species_amu = 1.0;
  7 species_Z = 1.0;
  8 runIdent = "template";
  9 nP_Vx = 5;
 10 nP_Vy = 5;
 11 nP_Vz = 5;
 12 nThermal = 3;
 13 nPhi = 0;
 14 ky = 0.0;
 15 T_keV = 5.0;
 16 input_fName = "input/input-data.nc";
 17 kz = 0.0;
 18 nStepsPerCyclotronPeriod = 30.0;
```

### Using IDL to run the sigma benchmarks

1. Clone https://github.com/dlg0/lib_dlg
2. Add `kineticj` and `lib_dlg` idl folders to the IDL path, i.e., edit the `export $IDL_STARTUP=~/idlStartup.pro` file, e.g., 
    ```
    !PATH=EXPAND_PATH('+~/code/kineticj/idl:+~/code/lib_dlg/idl:<IDL_DEFAULT>')
    ```
3. Change to the benchmarks folder, start IDL, run one of the 3 benchmarks.
    ```
    cd $KINETICJ
    cd benchmarks
    ./cleanBenchmarks.sh
    idl
    IDL>kj_sigma_benchmarks, runKJ=1, benchmark=1
    ```
    or at NERSC (Edison with a single node using OpenMP)
    ```
    cd $KINETICJ
    source env-edison.sh
    cd $SCRATCH
    mkdir kineticj
    cp -r $KINETICJ/benchmarks .
    cd benchmarks
    ./cleanBenchmarks.sh
    salloc -N 1 -p debug
    idl
    IDL>kj_sigma_benchmarks, runKJ=1, benchmark=1
    ```
4. Examine the `kineticj/benchmarks/benchmark1/benchmark.png` file for output. 

## Other Information

### Generate the Guiding Center terms file
To speed up the guiding center orbit calculation (if selected), some of the terms in the ODE are precalculated and tabulated in a file. This is done using code from the https://github.com/dlg0/orbit_tracer repo. After cloning that repo, run the `gc_terms.pro` routine on an `ar2Input.nc` file, e.g., 

```
cd ~/scratch/aorsa2d/colestock-kashuba-reference
idl
IDL>@constants
IDL>gc_terms, _me_amu, -1, ar2='input/ar2Input.nc', /oneD
```
