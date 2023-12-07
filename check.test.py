
import os, time
import numpy as np
import xarray as xr
from testNameList import nameList
#testIDList   = ['tpv8', 'tpv104','tpv1053d','tpv1053d.6c','meng2023a','meng2023cb']
fileNameList = ['fault.dyna.r.nc','frt.txt0','frt.txt1', 'frt.txt2','frt.txt3']
refRoot  = 'test.reference.results'
testRoot = 'test'

def compare_nc_files(fn1,fn2,threshold=1e-3):
    isTheSame = 'SUCCESS '+fn1+' '+fn2
    f1 = xr.open_dataset(fn1)
    f2 = xr.open_dataset(fn2)
    metadata_equal = f1.identical(f2)
    data_equal = (f1==f2).all().items()
    
    # Compare variables
    for var in f1.variables:
        var1 = f1[var]
        var2 = f2[var]
        if var1.dims != var2.dims:
            isTheSame = 'FAIL var dim '+fn1+' '+fn2
        if not np.allclose(var1,var2,rtol=threshold, atol=threshold):
            isTheSame = 'FAIL var numbers '+fn1+' '+fn2
        
    if not metadata_equal:# and data_equal:
        isTheSame = 'FAIL metadata '+fn1+' '+fn2
        
    print(isTheSame)
    
    return isTheSame

def compare_txt_files(fn1,fn2,threshold=1e-3):
    isTheSame = 'SUCCESS '+fn1+' '+fn2
    with open(fn1,'r') as f1, open(fn2,'r') as f2:
        result1 = f1.read().split()
        result2 = f2.read().split()
    if len(result1) != len(result2):
        isTheSame = 'FAIL '+fn1+' '+fn2
    
    for num1,num2 in zip(result1,result2):
        fnum1, fnum2 = float(num1),float(num2)
        if abs(fnum1-fnum2) > threshold:
            isTheSame = 'FAIL '+fn1+' '+fn2
    print(isTheSame)

    

for testid in nameList:
    for filename in fileNameList:
        refPath  = refRoot+'/'+testid+'/'+filename
        testPath = testRoot+'/'+testid+'/'+filename
        if os.path.exists(refPath):
            if 'nc' in filename:
                compare_nc_files(refPath, testPath, 1e-3)
            elif 'frt' in filename:
                compare_txt_files(refPath, testPath, 1e-3)

