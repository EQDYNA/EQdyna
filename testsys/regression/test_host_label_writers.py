#! /usr/bin/env python3
"""
Regression guard (board row 162, owner decision 2026-11-03): no writer may
record this box's real hostname.  Every perf/profile/evidence writer records
profile_emit.host_label() (EQDYNA_HOST_LABEL, default 'workstation'); the
Fortran profile writer does the same via readHostname.

Two halves, both able to go RED:

  STATIC   every tracked .py/.f90/.sh/.c file under src/, testsys/, scripts/
           (profile_emit.py aside, whose docstring names the banned calls) is
           scanned, comment lines skipped, for a call that reads the real
           hostname (python socket/os.uname/platform.node, the Fortran
           hostnm, the kernel hostname file).  A NEGATIVE CONTROL first feeds
           the scanner a planted line per pattern and requires a hit, so a
           pattern that rotted into matching nothing fails loudly.
  BEHAVE   host_label() in a clean subprocess returns 'workstation' with the
           variable unset, the override when set, and never the real
           hostname.

Complements test_no_hardcoded_paths.py, which catches a committed artifact
AFTER the leak; this catches the writer BEFORE it runs.  Never loosen by
adding a file to an allowlist -- route the writer through host_label().
"""
import os
import re
import socket
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SELF_RELPATH = 'testsys/regression/test_host_label_writers.py'
HELPER_RELPATH = 'src/python/eqdyna/profile_emit.py'
SCAN_DIRS = ('src/', 'testsys/', 'scripts/')
SCAN_EXT = ('.py', '.f90', '.sh', '.c')

# Fragments so this file's own source never matches its own patterns.
BANNED = [
    ('python socket', r'socket\.' + 'get' + r'hostname\s*\('),
    ('python fqdn', r'socket\.' + 'get' + r'fqdn\s*\('),
    ('python uname', r'\.' + 'uname' + r'\(\)\.' + 'node' + 'name'),
    ('python platform.node', r'platform\.' + 'node' + r'\s*\('),
    ('fortran hostnm', r'\b' + 'hostnm' + r'\b'),
    ('kernel hostname file', 'kernel/' + 'host' + 'name'),
    ('shell hostname', r'\$\(\s*host' + r'name\b'),
]
PLANTED = {
    'python socket': 'h = socket.' + 'get' + 'hostname()',
    'python fqdn': 'h = socket.' + 'get' + 'fqdn()',
    'python uname': 'h = os.' + 'uname' + '().' + 'node' + 'name',
    'python platform.node': 'h = platform.' + 'node' + '()',
    'fortran hostnm': 'call ' + 'hostnm' + '(s)',
    'kernel hostname file': "open(file='/proc/sys/kernel/" + 'host' + "name')",
    'shell hostname': 'H=$( host' + 'name -s)',
}


def scan_text(text):
    """[(label, lineno)] for every non-comment line matching a banned call."""
    hits = []
    for n, line in enumerate(text.splitlines(), 1):
        st = line.lstrip()
        if st.startswith('#') or st.startswith('!'):
            continue
        for label, pat in BANNED:
            if re.search(pat, line):
                hits.append((label, n))
    return hits


def tracked_files():
    r = subprocess.run(['git', 'ls-files', '--'] + list(SCAN_DIRS), cwd=ROOT,
                       capture_output=True, text=True, timeout=60)
    if r.returncode != 0:
        raise RuntimeError('git ls-files failed: %s' % r.stderr)
    return [p for p in r.stdout.splitlines() if p.endswith(SCAN_EXT)]


def main():
    problems = []

    # negative control: every pattern must fire on its planted line
    for label, _ in BANNED:
        if not scan_text(PLANTED[label]):
            problems.append('scanner pattern %r matches nothing (rotted)' % label)
        if scan_text('# ' + PLANTED[label]):
            problems.append('scanner pattern %r fires on a comment line' % label)

    files = tracked_files()
    if len(files) < 100:
        problems.append('only %d tracked source files scanned -- git ls-files '
                        'or SCAN_DIRS broke' % len(files))
    leaks = []
    for p in files:
        if p in (SELF_RELPATH, HELPER_RELPATH):
            continue
        try:
            text = open(os.path.join(ROOT, p), errors='replace').read()
        except FileNotFoundError:
            continue  # deleted in the working tree, not scanned
        for label, n in scan_text(text):
            leaks.append('%s:%d (%s)' % (p, n, label))
    if leaks:
        problems.append('%d writer line(s) read the real hostname; route each '
                        'through profile_emit.host_label():\n    %s'
                        % (len(leaks), '\n    '.join(leaks)))

    # behavioural half, in a clean interpreter (-I: no cwd/env import tricks)
    real = socket.gethostname()
    code = ('import sys; sys.path.insert(0, %r); '
            'from eqdyna.profile_emit import host_label; print(host_label())'
            % os.path.join(ROOT, 'src', 'python'))
    env = {k: v for k, v in os.environ.items() if k != 'EQDYNA_HOST_LABEL'}
    got = subprocess.run([sys.executable, '-I', '-c', code], env=env,
                         capture_output=True, text=True, timeout=60)
    if got.returncode != 0 or got.stdout.strip() != 'workstation':
        problems.append('host_label() default is %r (rc %d), want '
                        "'workstation': %s" % (got.stdout.strip(), got.returncode,
                                               got.stderr.strip()[-200:]))
    if real and real in got.stdout:
        problems.append('host_label() returned the real hostname')
    env['EQDYNA_HOST_LABEL'] = 'nodeX'
    got2 = subprocess.run([sys.executable, '-I', '-c', code], env=env,
                          capture_output=True, text=True, timeout=60)
    if got2.stdout.strip() != 'nodeX':
        problems.append('EQDYNA_HOST_LABEL override ignored: %r' % got2.stdout.strip())

    if problems:
        for p in problems:
            print('FAIL test_host_label_writers: %s' % p)
        return 1
    print('PASS test_host_label_writers: %d tracked source files scanned, 0 '
          'read the real hostname; %d patterns each fire on a planted line; '
          "host_label() default 'workstation', override honoured"
          % (len(files), len(BANNED)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
