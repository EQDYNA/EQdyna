#! /usr/bin/env python3
"""
Regression guard: Docker.guide.md must not pin a specific version (rules 2, 11).

THE INCIDENT this guards against: `dunyuliu/eqdyna:v5.3.1` on Docker Hub sat
seven minor versions behind the repo before `publish.yml`/`ghcr.io` replaced
that hand-`docker commit` workflow (see the root `Dockerfile`'s own comment).
`Docker.guide.md` had the SAME failure mode waiting to happen the ordinary
way: it hardcoded `ghcr.io/eqdyna/eqdyna:v5.8.2` in four places (found
2026-09-16/17, while the repo was already at v5.9.0) -- a doc that quotes one
version number as its running example goes stale the next time anyone tags a
release, silently, because nothing reads it and nothing runs it.

WHAT THIS PINS: no `vX.Y.Z`-shaped version string appears anywhere in
Docker.guide.md. `:latest` and a generic `:vX.Y.Z` placeholder (for the "pin a
specific release yourself" instruction) are both fine -- neither can drift,
because neither names an actual version. A real tag string
(`v` + digits + `.` + digits + `.` + digits) is the only thing this refuses.

Cheap (rule 9): one regex over one file, no subprocess, no network.
"""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
GUIDE = os.path.join(ROOT, 'Docker.guide.md')

# A real, concrete version tag: v + digits.digits.digits. Does NOT match the
# placeholder spelling "vX.Y.Z" used in the guide's own "pin a version"
# instruction, since X/Y/Z are letters, not digits.
PINNED_VERSION = re.compile(r'\bv\d+\.\d+\.\d+\b')


def main():
    if not os.path.exists(GUIDE):
        print('FAIL test_docker_guide_no_pinned_version: %s does not exist' % GUIDE)
        return 1

    text = open(GUIDE, errors='replace').read()
    hits = []
    for lineno, line in enumerate(text.splitlines(), 1):
        for m in PINNED_VERSION.finditer(line):
            hits.append((lineno, m.group(0), line.strip()))

    if hits:
        print('FAIL test_docker_guide_no_pinned_version: %d pinned version '
              'string(s) found -- Docker.guide.md must reference `:latest` or '
              'the literal placeholder `:vX.Y.Z`, never a real tag, or this '
              'doc goes stale the next time a release is cut (rule 2):' % len(hits))
        for lineno, hit, line in hits:
            print('  Docker.guide.md:%d: %r in: %s' % (lineno, hit, line))
        return 1

    print('  PASS  no pinned version string in Docker.guide.md '
          '(%d bytes checked)' % len(text))
    print('SUCCESS test_docker_guide_no_pinned_version')
    return 0


if __name__ == '__main__':
    sys.exit(main())
