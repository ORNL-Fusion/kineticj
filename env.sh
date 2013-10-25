module switch PrgEnv-pgi PrgEnv-gnu
module switch gcc gcc/4.6.3
module load cudatoolkit
module load cray-netcdf/4.2.1.1
export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
