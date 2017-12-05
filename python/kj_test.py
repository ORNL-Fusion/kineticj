from netCDF4 import Dataset
import numpy as np
from subprocess import call
import os

tests = ["benchmark1-00007","benchmark2-00013","benchmark3-00004","test4"]

rtol = 1e-5

KJPATH = os.environ['KINETICJ']
cmd = os.path.expanduser(KJPATH+"/bin/kineticj")
print cmd
args = ""
cwd = os.getcwd()

for test in tests:
   
    os.chdir(test)
    call([cmd,args]) 

    ncCorrect = Dataset('output/jP2-correct.nc')
    ncThisrun = Dataset('output/jP2.nc')

    j1xre1 = ncThisrun.variables['j1xc_re'][:]
    j1xim1 = ncThisrun.variables['j1xc_im'][:]
    j1yre1 = ncThisrun.variables['j1yc_re'][:]
    j1yim1 = ncThisrun.variables['j1yc_im'][:]
    j1zre1 = ncThisrun.variables['j1zc_re'][:]
    j1zim1 = ncThisrun.variables['j1zc_im'][:]

    j1xre2 = ncCorrect.variables['j1xc_re'][:]
    j1xim2 = ncCorrect.variables['j1xc_im'][:]
    j1yre2 = ncCorrect.variables['j1yc_re'][:]
    j1yim2 = ncCorrect.variables['j1yc_im'][:]
    j1zre2 = ncCorrect.variables['j1zc_re'][:]
    j1zim2 = ncCorrect.variables['j1zc_im'][:]

    result = []
    
    result.append( np.allclose(j1xre1,j1xre2,rtol=rtol) )
    result.append( np.allclose(j1xim1,j1xim2,rtol=rtol) )
    result.append( np.allclose(j1yre1,j1yre2,rtol=rtol) )
    result.append( np.allclose(j1yim1,j1yim2,rtol=rtol) )
    result.append( np.allclose(j1zre1,j1zre2,rtol=rtol) )
    result.append( np.allclose(j1zim1,j1zim2,rtol=rtol) )

    if not any(result):
        output = test.ljust(20)+ 'FAIL'
        #print test+": FAIL"
    else:
        output = test.ljust(20)+ 'PASS'
        print output
        #print test+": PASS"

    os.chdir(cwd)

