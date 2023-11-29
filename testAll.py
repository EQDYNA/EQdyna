#! /bin/bash 
import os
# This script will perform tests on default test cases.
# Currently, it includes tpv104 and tpv1053d with smaller domains and coarse element sizes.
#
testList = ['tpv104','tpv1053d','tpv1053d.6c','meng2023a','meng2023cb']
cpuNumList  = [4, 4, 6, 4, 4]
MPIRUN='mpirun.mpich'

os.system('rm -rf test')
os.system('rm -rf bin/eqdyna')
os.system('mkdir test')

os.system('./install-eqdyna.sh -m ubuntu')
os.chdir('test')

def runTest(testDir, compSet, coreNum):
    cmd = 'create.newcase '+testDir+' '+compSet
    os.system(cmd)
    os.chdir(testDir)
    os.system('./case.setup')
    os.system('mpirun.mpich -np '+str(coreNum)+' eqdyna')
    os.system('python3 plotRuptureDynamics')
    os.chdir('..')
    
for testName, coreNum in zip(testList, cpuNumList):
    runTest(testName, 'test.'+testName, coreNum)

os.chdir('..')
os.system('python3 check.test.py')

