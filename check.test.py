
import os
import xarray as xr

testIDList   = ['tpv104','tpv1053d','tpv1053d.6c','meng2023a','meng2023cb']
fileNameList = ['fault.dyna.r.nc','frt.txt0','frt.txt2']
refRoot  = 'test.reference.results'
testRoot = 'test'

def compare(fn1,fn2):
    f1 = xr.open_dataset(fn1)
    f2 = xr.open_dataset(fn2)
    metadata_equal = f1.identical(f2)
    data_equal = (f1==f2).all().items()

    if metadata_equal and data_equal:
        print('SUCCESS')
    else:
        print('FAIL')


for testid in testIDList:
    for filename in fileNameList:
        refPath  = refRoot+'/'+testid+'/'+filename
        testPath = testRoot+'/'+testid+'/'+filename
        os.system('diff -s '+refPath+' '+testPath)
        #compare(refPath, testPath)
        #os.system('ncdiff '+refPath+' '+testPath+' tmp.nc')

