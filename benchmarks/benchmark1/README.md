This benchmark is a scan of Te for a D plasma, with ky = 0. Run on gpufusion with ...
```
qsub -I
source ~/code/kineticj/env-gpufusion.sh
cd ~/scratch/kineticj
cp -r ~/code/kineticj/benchmarks .
cd benchmarks/benchmark1
idl
IDL>kj_sigma_benchmarks, benchmark = 1, runKJ=1
```
