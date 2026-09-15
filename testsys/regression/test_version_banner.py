#! /usr/bin/env python3
"""
Regression guard: the runtime banner must match VERSION (rules 11, 14).

Guards a real drift: at the v5.6.0 release `src/eqdyna3d.f90` still printed
"Welcome to EQdyna 5.3.3" while VERSION said 5.5.0 -- two releases stale. It
is the first thing every run prints and the first thing a user quotes in a bug
report, so a stale banner silently misattributes results to the wrong version.

Nothing enforced it because the number lives in a Fortran string literal that
no build step touches. This test is the enforcement: it is a pure text
comparison, needs no build, and runs in milliseconds.

Exits non-zero on any failure (rule 2).
"""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
VERSION_FILE = os.path.join(ROOT, 'VERSION')
BANNER_FILE = os.path.join(ROOT, 'src', 'eqdyna3d.f90')
BANNER_RE = re.compile(r'Welcome to EQdyna\s+([0-9]+\.[0-9]+\.[0-9]+)')


def main():
    print('Regression guard: runtime banner must match VERSION')
    version = open(VERSION_FILE, errors='replace').read().strip()
    if not version:
        print('FAIL: VERSION is empty')
        return 1

    text = open(BANNER_FILE, errors='replace').read()
    found = BANNER_RE.findall(text)
    if not found:
        print('FAIL: no "Welcome to EQdyna <x.y.z>" banner in %s'
              % os.path.relpath(BANNER_FILE, ROOT))
        return 1
    if len(set(found)) != 1:
        print('FAIL: the banner names several versions: %s' % sorted(set(found)))
        return 1

    banner = found[0]
    if banner != version:
        print('FAIL: banner says %s, VERSION says %s.' % (banner, version))
        print('  %s is the first line of every run and the number users quote'
              % os.path.relpath(BANNER_FILE, ROOT))
        print('  in bug reports. Update the banner in the same commit as'
              ' VERSION (rule 11).')
        return 1

    print('  banner %s == VERSION %s' % (banner, version))
    print('\nPASS')
    return 0


if __name__ == '__main__':
    sys.exit(main())
