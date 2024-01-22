#! /bin/bash 
import os, time
from testNameList import nameList, coreNumList
# This script will perform tests on default test cases.
# Currently, it includes tpv104 and tpv1053d with smaller domains and coarse element sizes.
#
#testList = ['tpv8', 'tpv104','tpv1053d','tpv1053d.6c','meng2023a','meng2023cb']
#cpuNumList  = [4, 4, 4, 6, 4, 4]
MPIRUN='mpirun'

os.system('rm -rf test')
os.system('rm -rf bin/eqdyna')
os.system('mkdir test')

os.system('./install-eqdyna.sh -m ubuntu')
os.chdir('test')

startTime = time.time()
def runTest(testDir, compSet, coreNum):
    cmd = 'create.newcase '+testDir+' '+compSet
    os.system(cmd)
    os.chdir(testDir)
    os.system('./case.setup')
    os.system('mpirun -np '+str(coreNum)+' eqdyna')
    os.system('python3 plotRuptureDynamics')
    os.chdir('..')
    
for testName, coreNum in zip(nameList, coreNumList):
    runTest(testName, testName, coreNum)

os.chdir('..')
os.system('python3 check.test.py')

endTime = time.time()

print('Time consumed for all the tests are ', endTime-startTime, ' s')
