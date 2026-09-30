"""Shared helpers for the testsys tiers (see PROJECT_RULES.md rule 3).

Kept deliberately tiny: printing a result is not a gate, so every tier
funnels its outcome through here into a single non-zero-on-failure exit,
per rule 3 ("the script ran" and "the script passed" are different
questions).
"""
import os
import shutil
import subprocess

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# The machine whose src/fortran/makefile branch builds and probes use: what install-eqdyna.sh
# exported, else ubuntu. Hardcoding ubuntu built src/ with `-I/usr/include` on ls6 and failed.
MACHINE = os.environ.get('EQDYNA_TEST_MACHINE') or os.environ.get('MACHINE', 'ubuntu')


def banner(msg):
    print(f'\n==== testsys: {msg} ====')


FSRC = os.path.join(REPO_ROOT, 'src', 'fortran')


def make_var(name):
    """One src/fortran/makefile variable for this machine ($MACHINE, else ubuntu),
    e.g. NETCDF_LIB. Probes ask make rather than spell out ubuntu's paths, which
    on ls6 failed as `ld: cannot find -lnetcdf`. Raises if it is empty."""
    env = dict(os.environ, MACHINE=MACHINE)
    out = subprocess.run(['make', '-s', 'print-' + name], cwd=FSRC, env=env,
                         capture_output=True, text=True, check=True).stdout.split()
    if not out:
        raise RuntimeError('src/fortran/makefile defines no %s for MACHINE=%s'
                           % (name, MACHINE))
    return out


def stage(workdir, name):
    """Copy src/fortran/<name> into workdir; return the copy to compile.

    gfortran looks for modules in the SOURCE FILE'S directory before -I, so
    compiling src/fortran/x.f90 in place reads the ifort globalvar.mod the real
    ls6 build left beside it ("Reading module ... Unexpected EOF")."""
    dst = os.path.join(workdir, name)
    shutil.copyfile(os.path.join(FSRC, name), dst)
    return dst
