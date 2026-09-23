#! /usr/bin/env python3
"""
Image gate: `.git/config` inside the just-built EQdyna image must carry NO
persisted credential (2026-09-23 incident).

THE INCIDENT. `actions/checkout@v4` in `.github/workflows/publish.yml`
defaults `persist-credentials` to `true`, which writes the live
`GITHUB_TOKEN` into the checkout's `.git/config` as an
`http.<origin>/.extraheader` entry. The cleanup for that is a POST step
(`action.yml`: `post: dist/index.js`), which runs AFTER every main step in
the job -- including `docker build`. `Dockerfile` does `COPY . /opt/eqdyna`
and `.dockerignore` deliberately keeps `.git` (two in-image regression
guards read real history from it), so at the default, the token would be
baked into a layer of every published image. Fixed by adding
`with: persist-credentials: false` to that checkout step, pinned by
`testsys/regression/test_publish_image_fetch_depth.py`. THIS script is the
second layer: it runs INSIDE the just-built image, before the push, and
checks the artifact itself rather than trusting the workflow source that
produced it -- catching a future regression on this same setting, or a
different mechanism (e.g. an `ssh-key` checkout input) that leaves a
credential in the same file by a different route.

Usage:
    python3 testsys/check_git_config_no_credential.py [path-to-git-config]

`path-to-git-config` defaults to `.git/config` relative to the current
working directory (the publish gate runs this with `cwd=/opt/eqdyna`, so the
default resolves to `/opt/eqdyna/.git/config`, the file `COPY . /opt/eqdyna`
actually shipped). A path is accepted as an argument so this script can be
pointed at a constructed fixture for its own regression/mutation test,
without needing a real GitHub Actions credential leak to reproduce one.

MARKERS CHECKED, case-insensitively, and why each is here:
  - `extraheader`     -- the config KEY `persist-credentials: true` writes
                         (`http.<url>.extraheader`), holding an
                         `AUTHORIZATION: basic ...`-style header.
  - `authorization`    -- the header name itself, in case a future encoding
                         puts it somewhere other than an extraheader value.
  - `x-access-token`   -- the basic-auth username checkout uses when
                         encoding the token (`x-access-token:<token>`,
                         base64'd, but the literal also appears unencoded in
                         some checkout code paths and in insteadOf rewrites).
  - `ghp_`, `ghs_`, `gho_`, `ghu_`, `ghr_`, `github_pat_`
                        -- GitHub token prefixes (classic PAT, GitHub App
                         installation token, OAuth, user-to-server, refresh,
                         fine-grained PAT respectively). A token can end up
                         in `.git/config` by a route OTHER than
                         `persist-credentials` (e.g. hand-written into an
                         `insteadOf` URL), so these are checked regardless of
                         which config key carries them.
  - `sshcommand`       -- `core.sshCommand`, what checkout's `ssh-key` input
                         (not used by this workflow today, but a plausible
                         future edit) persists instead of an http header;
                         its presence means a private key was wired in, even
                         though the key material itself lives in a separate
                         file this script does not have to find to flag the
                         risk.

Exits 2 if the target file does not exist -- a gate with nothing to scan is
not a passing gate, it is a gate that never ran (rule 2: fail loudly, never
skip). Exits 1 and names every marker found, with the offending line, if any
marker is present. Exits 0, and prints the byte/line count of what it
scanned, only if the file exists and none of the markers appear.
"""
import os
import sys

MARKERS = (
    'extraheader',
    'authorization',
    'x-access-token',
    'ghp_',
    'ghs_',
    'gho_',
    'ghu_',
    'ghr_',
    'github_pat_',
    'sshcommand',
)


def scan(path):
    """Returns (line_count, [(marker, line_no, line_text), ...])."""
    with open(path, errors='replace') as fh:
        lines = fh.readlines()
    hits = []
    for i, raw in enumerate(lines, start=1):
        lower = raw.lower()
        for marker in MARKERS:
            if marker in lower:
                hits.append((marker, i, raw.strip()))
    return len(lines), hits


def main(argv):
    path = argv[1] if len(argv) > 1 else os.path.join('.git', 'config')

    if not os.path.isfile(path):
        print('REFUSED check_git_config_no_credential: %r does not exist -- '
              'there is nothing to scan, which is not the same as nothing '
              'being found. If `.git/config` is genuinely absent from this '
              'image, re-derive this gate rather than let it report clean '
              'on a file it never opened.' % path)
        return 2

    line_count, hits = scan(path)

    if hits:
        print('FAIL check_git_config_no_credential: %d credential marker(s) '
              'found in %r (%d line(s) scanned):' % (len(hits), path, line_count))
        for marker, lineno, text in hits:
            print('  line %d: matched %r -- %s' % (lineno, marker, text))
        print('  A persisted git credential in this file means the image '
              'build baked a token (or SSH key reference) into a layer. '
              'See .github/workflows/publish.yml\'s checkout step '
              '(persist-credentials) and re-check for an accidental '
              'insteadOf/credential.helper rewrite.')
        return 1

    print('PASS check_git_config_no_credential: %r scanned (%d line(s)), no '
          'credential marker among %d checked (%s)'
          % (path, line_count, len(MARKERS), ', '.join(MARKERS)))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
