#! /usr/bin/env python3
"""
One shared helper (pathway item 95; not a test -- the regression tier runs
only `test_*.py`): launch a binary under mpirun with each rank's stdout and
stderr redirected to a rank-owned FILE, and return what the ranks wrote.

WHY: item 89 measured that Open MPI discards already-flushed bytes it has not
yet forwarded when MPI_Abort tears the job down, so a guard asserting on text
from mpirun's PIPE after an abortRun refusal can flake. A rank-owned file
cannot lose bytes the rank flushed. Technique: test_stop_exit_status.py
probe B. `exec` keeps the binary's own exit status; `$$` (the shell's pid,
which becomes the binary's) names each rank's file.
"""
import glob
import os
import shlex
import subprocess

PATTERN = 'eqdyna-stdout.*.txt'


def run_rank_files(mpirun, binary, cwd, np=1, env=None, timeout=120):
    """(returncode, text): the ranks' own output files (new since this call)
    followed by whatever mpirun itself printed. Raises RuntimeError when no
    rank file appears -- then the binary never ran, and a text assertion on
    mpirun's own chatter would prove nothing. subprocess.TimeoutExpired
    propagates to the caller."""
    before = set(glob.glob(os.path.join(glob.escape(cwd), PATTERN)))
    shell = 'exec %s > eqdyna-stdout.$$.txt 2>&1' % shlex.quote(binary)
    r = subprocess.run([mpirun, '-np', str(np), 'sh', '-c', shell], cwd=cwd, env=env,
                       capture_output=True, timeout=timeout)
    launcher = (r.stdout + r.stderr).decode(errors='replace')
    files = sorted(set(glob.glob(os.path.join(glob.escape(cwd), PATTERN))) - before)
    if not files:
        raise RuntimeError('no rank-owned stdout file under %s -- %s never ran '
                           '(mpirun rc=%d): %s' % (cwd, binary, r.returncode,
                                                   launcher[-500:]))
    text = ''.join(open(f, errors='replace').read() for f in files)
    return r.returncode, text + launcher
