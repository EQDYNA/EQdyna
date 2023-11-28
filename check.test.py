
import os

testIDList   = ['tpv104','tpv1053d','meng2023a','meng2023cb']
fileNameList = ['fault.dyna.r.nc']
refRoot  = 'test.reference.results'
testRoot = 'test'

for testid in testIDList:
    for filename in fileNameList:
        refPath  = refRoot+'/'+testid+'/'+filename
        testPath = testRoot+'/'+testid+'/'+filename
        os.system('diff -s '+refPath+' '+testPath)
        os.system('cmp -bl '+refPath+' '+testPath)

