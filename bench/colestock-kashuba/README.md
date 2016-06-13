#Colestock-Kashuba Test Case

##Overview
This is a 1-D benchmark with AORSA to test the iterative generation of an Ion Bernstein Wave (IBW), and its subsequent damping. This tests the following aspects of Kinetic-J:

1. Mulitple Species
2. Perpendicualr Kinetics
3. Parallel Kinetics (the IBW is damped via electron Landau damping)
4. The guiding center implementation for electrons to avoid the electron gyro-frequency.

##To Run
```
cd bench/colestock-kashuba
../../bin/kineticj
idl
IDL>kj_plot_current
```
