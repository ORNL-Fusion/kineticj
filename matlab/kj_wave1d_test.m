function tests = kj_wave1d_test

% tests = functiontests(localfunctions);

testCase = 1;
tests{1} = kj_wave1d_vaccum_test(testCase);

testCase = 2;
tests{2} = kj_wave1d_cold_plasma_test(testCase);

testCase = 3;
tests{3} = kj_wave1d_cold_plasma_dispersion_test(testCase);

end












