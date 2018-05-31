# Notes on `kj_wave1d.m`
This is a 1D, finite-difference, frequency-domain (FD-FD) solver for both vacuum and cold plasma problems. 

## Run the tests
Several test cases are implemented and run via ...
```
runtests('kj_wave1d_test.m')
Running kj_wave1d_test
Test 1: Error vs 1/h for vacuum eps
Test 2: Error vs 1/h for cold plasma eps
Test 3: Compare kx actual with dispersion
```
### Test 1
This test examines the error between an analytic solution and the numeric solution for vacuum EM with the check being on the slope of that curve, i.e., ~2 for the second order spatial discretization scheme. The analytic solution is constructed via the method of manufactured solutions. 

### Test 2
Same as Test 1 except using a cold plasma dielectric tensor instead of the vacuum (identity) tensor).

### Test 3
Does a scan over the B field (an input to the cold plasma dispersion solver) and compares the resulting k spectra scan with a matlab saved file of the same. 

## Run the example
An example on how to run the code is shown in `kj_wave1d_test.m` and can be run using ...
```
kj_wave1d_example
```
