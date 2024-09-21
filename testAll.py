#! /bin/bash 
import os, time
from testNameList import nameList, coreNumList
# This script will perform tests on test cases defined in testNameList.

MPIRUN='mpirun' # please modify MPIRUN to fit your system accordingly. 
print('testAll: MPIRUN is ', MPIRUN)
print('testAll: please modify MPIRUN to fit your system accordingly.')

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
    os.system(MPIRUN+' -np '+str(coreNum)+' eqdyna')
    os.system('python3 plotRuptureDynamics')
    os.chdir('..')
    
for testName, coreNum in zip(nameList, coreNumList):
    runTest(testName, testName, coreNum)

os.chdir('..')
os.system('python3 check.test.py')

endTime = time.time()

print('Time consumed for all the tests are ', endTime-startTime, ' s')
