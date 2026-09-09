#!/usr/bin/env python3
"""
Compares a completed test/ run against the golden test.reference.results/
tree (PROJECT_RULES.md rules 2, 3, 4, 5, 7).

"Pass" means every comparison below prints SUCCESS and the process exits 0.
A printed FAIL, a missing test file, or zero comparisons run is a hard
failure (rule 2) -- nothing here is allowed to fail silently or leave a
false-green exit code on screen.

threshold=1e-3 is the one calibrated tolerance (rule 5), used by both
compare_nc_files and compare_txt_files. scripts/compareTwoNc.py is a
separate ad hoc diff utility and is not held to this threshold.

compare_nc_files/compare_txt_files are plain functions with no import-time
side effects; the comparison loop itself only runs under __main__ so this
file can be imported by tests (see testsys/unit/test_check_comparisons.py)
without touching the real test/ or test.reference.results/ trees.
"""
import os, sys
import numpy as np
import xarray as xr

fileNameList = ['fault.dyna.r.nc', 'frt.txt0', 'frt.txt1', 'frt.txt2', 'frt.txt3']
refRoot = 'test.reference.results'
testRoot = 'test'
THRESHOLD = 1e-3


def compare_nc_files(fn1, fn2, threshold=THRESHOLD):
    isTheSame = 'SUCCESS ' + fn1 + ' ' + fn2
    f1 = xr.open_dataset(fn1)
    f2 = xr.open_dataset(fn2)
    try:
        metadata_equal = f1.identical(f2)
        for var in f1.variables:
            var1 = f1[var]
            var2 = f2[var]
            if var1.dims != var2.dims:
                isTheSame = 'FAIL var dim ' + fn1 + ' ' + fn2
            elif not np.allclose(var1, var2, rtol=threshold, atol=threshold):
                isTheSame = 'FAIL var numbers ' + fn1 + ' ' + fn2
        if not metadata_equal and isTheSame.startswith('SUCCESS'):
            isTheSame = 'FAIL metadata ' + fn1 + ' ' + fn2
    finally:
        f1.close()
        f2.close()
    print(isTheSame)
    return isTheSame


def compare_txt_files(fn1, fn2, threshold=THRESHOLD):
    with open(fn1, 'r') as f1, open(fn2, 'r') as f2:
        result1 = f1.read().split()
        result2 = f2.read().split()

    if len(result1) != len(result2):
        isTheSame = f'FAIL length mismatch ({len(result1)} vs {len(result2)}) {fn1} {fn2}'
        print(isTheSame)
        return isTheSame

    isTheSame = 'SUCCESS ' + fn1 + ' ' + fn2
    for num1, num2 in zip(result1, result2):
        fnum1, fnum2 = float(num1), float(num2)
        if abs(fnum1 - fnum2) > threshold:
            isTheSame = 'FAIL ' + fn1 + ' ' + fn2
            break
    print(isTheSame)
    return isTheSame


def main():
    from testNameList import nameList

    results = []
    for testid in nameList:
        print(' ')
        for filename in fileNameList:
            refPath = refRoot + '/' + testid + '/' + filename
            testPath = testRoot + '/' + testid + '/' + filename
            if not os.path.exists(refPath):
                continue
            if not os.path.exists(testPath):
                msg = 'FAIL missing ' + testPath
                print(msg)
                results.append(msg)
            elif 'nc' in filename:
                results.append(compare_nc_files(refPath, testPath, THRESHOLD))
            elif 'frt' in filename:
                results.append(compare_txt_files(refPath, testPath, THRESHOLD))

    failures = [r for r in results if r.startswith('FAIL')]
    print(' ')
    if not results:
        print('check.test.py: FAIL - no comparisons were run '
              '(missing test.reference.results/ or test/ trees?)')
        return 1
    print(f'check.test.py: {len(results) - len(failures)}/{len(results)} comparisons SUCCESS')
    if failures:
        print('check.test.py: FAIL -', len(failures), 'comparison(s) failed')
        return 1
    print('check.test.py: SUCCESS')
    return 0


if __name__ == '__main__':
    sys.exit(main())
