#! /usr/bin/env python3
"""
Regression guard: fatal exits must carry a non-zero status (rules 2, 9).

Guards a past defect with two faces, both of which silently exit ZERO:

  * `stop` with no code           -- 9 sites
  * `stop 'some message'`         -- 4 sites

Under gfortran both print and then exit 0 (the probe below verifies this on
the toolchain actually in use). Two of the offenders were the hard-stop
alignment gates -- checkPMLAlignment and checkFaultMPIAlignment -- whose
entire purpose is to refuse a mesh that would produce a wrong answer. A
batch script, testsys/run.py, or CI checking $? scored those refusals as
passes.

All fatal paths now go through `abortRun(code, reason)` in src/fortran/errorCodes.f90,
which prints a structured block and calls MPI_Abort so the whole job dies
with `code` as its exit status instead of one rank stopping and the rest
blocking forever in their next collective.

This test fails if either form reappears in src/. A `stop` that legitimately
ends a successful run is allowed, but must say so on the same line with the
marker NORMAL-EXIT -- that keeps "this run succeeded" visibly distinct from
"someone forgot the exit code", which is the confusion that caused the bug.

Cheap (rule 9): one ~20-byte compile, a source scan, and two `mpirun -np 2`
refusals of an empty case (see probe_real_binary for why two); ~1 s measured.
Exits non-zero on any failure (rule 2).
"""
import glob
import os
import re
import shlex
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = os.path.join(ROOT, 'src', 'fortran')
REGISTRY = os.path.join(SRC, 'errorCodes.f90')

ALLOW_MARKER = 'NORMAL-EXIT'

# `stop` as a statement: at the start of a statement, or after `then`.
# Captures what follows so bare vs quoted vs numeric can be told apart.
#
# The leading `\s*` is load-bearing: every real site is indented, so an
# anchor of bare `^` matched nothing and the scan reported a vacuous PASS.
# Caught by the self-check below, which is why that self-check exists.
STOP_STMT = re.compile(
    r"(?:^\s*|\bthen\b\s*|\)\s*)\bstop\b(?P<rest>[^\n]*)", re.IGNORECASE)


def probe_compiler():
    """Report the exit status of `stop 'string'` on this toolchain."""
    fc = os.environ.get('FC', 'gfortran')
    with tempfile.TemporaryDirectory() as d:
        f90, exe = os.path.join(d, 'p.f90'), os.path.join(d, 'p')
        with open(f90, 'w') as fh:
            fh.write("program p\n  stop 'refused'\nend program p\n")
        try:
            subprocess.run([fc, '-o', exe, f90], check=True,
                           capture_output=True, timeout=120)
        except (OSError, subprocess.SubprocessError) as exc:
            print('  probe skipped (%s unavailable: %s)' % (fc, exc))
            return None
        rc = subprocess.run([exe], capture_output=True).returncode
        print("  %s: `stop 'string'` exits %d%s"
              % (fc, rc, '  <- the defect this guards' if rc == 0 else ''))
        return rc


def scan_sources():
    """Return [(relpath, lineno, form, code_text)] for every offending stop."""
    offenders = []
    for path in sorted(glob.glob(os.path.join(SRC, '*.f90'))):
        rel = os.path.relpath(path, ROOT)
        with open(path, errors='replace') as fh:
            for lineno, raw in enumerate(fh, 1):
                if raw.lstrip().startswith('!'):
                    continue
                code, _, comment = raw.partition('!')
                m = STOP_STMT.search(code)
                if not m:
                    continue
                rest = m.group('rest').strip()
                if rest.startswith(("'", '"')):
                    form = "stop '<literal>'"
                elif rest == '':
                    if ALLOW_MARKER in comment:
                        continue           # deliberate successful termination
                    form = 'bare stop'
                else:
                    continue               # `stop <integer>` still exits non-zero
                offenders.append((rel, lineno, form, code.strip()))
    return offenders


RANGES = [
    (1, 9,  'Generic'),
    (11, 19, 'Configuration and parameters'),
    (21, 29, 'Input files'),
    (31, 39, 'Fault geometry'),
    (41, 49, 'Mesh generation and element quality'),
    (51, 59, 'MPI and domain decomposition'),
    (61, 69, 'Numerics and runtime state'),
    (71, 79, 'External libraries'),
]
TROUBLESHOOTING = os.path.join(ROOT, 'docs', 'user', 'troubleshooting.md')
BEGIN = '<!-- BEGIN EXIT CODES (generated from src/fortran/errorCodes.f90; do not edit by hand) -->'
END = '<!-- END EXIT CODES -->'

CODE_DECL = re.compile(
    r'^\s*integer,\s*parameter\s*::\s*(ERR_\w+)\s*=\s*(\d+)\s*(?:!\s*(.*?))?\s*$')


def parse_codes():
    """[(value, name, one-line meaning)] from the registry, sorted by value."""
    out = []
    for raw in open(REGISTRY, errors='replace'):
        m = CODE_DECL.match(raw)
        if m:
            out.append((int(m.group(2)), m.group(1), (m.group(3) or '').strip()))
    return sorted(out)


def render_table(codes):
    """The docs/user/troubleshooting.md section, derived entirely from the
    registry's own comments."""
    L = [BEGIN, '',
         'When a run is refused or fails, EQdyna prints a `FATAL` block naming the code',
         'and the reason, and exits with that code. Codes are kept in 1-125 so the number',
         'in the source is the number the shell reports (a status above 255 wraps).', '',
         'How faithfully the number survives depends on the launcher. `mpirun` (Open MPI,',
         'hydra) reports it directly. `srun` reports the **maximum** status across tasks,',
         'so a straggler killed while `MPI_Abort` tears the job down yields 137 or 143 and',
         'masks the code; `sbatch` reports the wrapper script\'s status unless the script',
         'ends with `exit $?`. On ls6 and grace the launcher is `ibrun` -> `srun`. So treat',
         'a **non-zero status** as the reliable signal and the printed `FATAL` block as',
         'authoritative; the specific number is advisory under `srun`.', '']
    for lo, hi, title in RANGES:
        rows = [c for c in codes if lo <= c[0] <= hi]
        if not rows:
            continue
        L += ['**%s** (%d-%d)' % (title, lo, hi), '',
              '| code | name | meaning |', '|-----:|------|---------|']
        L += ['| %d | `%s` | %s |' % (v, n, d) for v, n, d in rows]
        L += ['']
    L += ['Exit status 0 means the run completed. Every fatal path goes through',
          '`abortRun` in `src/fortran/errorCodes.f90`, which calls `MPI_Abort` so the whole job',
          'ends instead of one rank stopping while the others block in a collective.',
          '', END]
    return '\n'.join(L)


def check_or_update_docs(update=False):
    """The docs/user/troubleshooting.md table is generated, so it cannot
    drift from the registry."""
    s = open(TROUBLESHOOTING, errors='replace').read()
    if BEGIN not in s or END not in s:
        return ["docs/user/troubleshooting.md has no generated exit-code section"
                " (run this test with --update to insert one)"]
    head, rest = s.split(BEGIN, 1)
    _, tail = rest.split(END, 1)
    current = BEGIN + rest.split(END, 1)[0] + END
    wanted = render_table(parse_codes())
    if current == wanted:
        return []
    if update:
        open(TROUBLESHOOTING, 'w').write(head + wanted + tail)
        print('  docs/user/troubleshooting.md exit-code table regenerated')
        return []
    return ["docs/user/troubleshooting.md's exit-code table no longer matches "
            "src/fortran/errorCodes.f90 -- rerun this test with --update"]


def self_check():
    """Prove the scanner can still SEE each defect before trusting its PASS.

    A source scan that matches nothing reports success indistinguishably
    from a clean tree. That is not hypothetical: the first version of this
    file anchored on a bare `^`, matched none of the indented real sites,
    and printed PASS over thirteen live offenders. Every case below is
    written the way it actually appeared in src/.
    """
    must_flag = [
        ('                    stop',                     'bare stop'),
        ('        stop 1004',                            None),   # fine: non-zero
        ("        stop 'checkPMLAlignment failed'",       "stop '<literal>'"),
        ('        stop "Stopped"',                        "stop '<literal>'"),
        ('        if (iFault>1) stop',                    'bare stop'),
        ('    stop ! NORMAL-EXIT: success',               None),   # explicitly allowed
        ('    call abortRun(ERR_NETCDF, "x")',            None),   # not a stop at all
    ]
    problems = []
    for line, expected in must_flag:
        code, _, comment = line.partition('!')
        m = STOP_STMT.search(code)
        got = None
        if m:
            rest = m.group('rest').strip()
            if rest.startswith(("'", '"')):
                got = "stop '<literal>'"
            elif rest == '' and ALLOW_MARKER not in comment:
                got = 'bare stop'
        if got != expected:
            problems.append('scanner returned %r, expected %r, for: %s'
                            % (got, expected, line.strip()))
    return problems


def check_registry():
    """Every ERR_* code must be unique and in the shell-safe range 1-125."""
    problems = []
    if not os.path.exists(REGISTRY):
        return ['src/fortran/errorCodes.f90 is missing']
    seen = {}
    pat = re.compile(r'^\s*integer,\s*parameter\s*::\s*(ERR_\w+)\s*=\s*(\d+)',
                     re.IGNORECASE)
    for lineno, raw in enumerate(open(REGISTRY, errors='replace'), 1):
        m = pat.match(raw)
        if not m:
            continue
        name, value = m.group(1), int(m.group(2))
        # 126/127 mean "cannot execute" and 128+n means "killed by signal n";
        # anything over 255 wraps, so the source number stops matching what
        # the shell reports.
        if not 1 <= value <= 125:
            problems.append('%s = %d is outside the shell-safe range 1-125'
                            % (name, value))
        if value in seen:
            problems.append('%s = %d collides with %s' % (name, value, seen[value]))
        else:
            seen[value] = name
    if not seen:
        problems.append('no ERR_* codes found in src/fortran/errorCodes.f90')
    return problems, len(seen)


def probe_real_binary():
    """Run the REAL MPI binary on a refusable input, in TWO invocations.

    The compiler probe above only proves what a serial `stop 'string'` does.
    The claims this file is actually about -- that a refused EQdyna run exits
    non-zero through mpirun rather than hanging, AND that it says why -- need
    the real binary. An empty directory has no bGlobal.txt, so the run must
    refuse with ERR_INPUT_FILE_MISSING.

    Why two invocations (pathway item 89). Both claims used to be asserted
    against ONE `capture_output=True` run of `mpirun -np 2`, and that made the
    FATAL-text half INTERMITTENT: CI run 35833383579 got exit 21 and no hang,
    but an empty capture. The text is not lost in the Fortran -- abortRun
    flushes unit 6 (src/fortran/errorCodes.f90:162) BEFORE MPI_Abort
    (:165) -- and not in Python. It is lost in Open MPI's I/O forwarding,
    which discards bytes a rank has already written and flushed when
    MPI_Abort tears the job down. Measured here on Open MPI 4.1.1 with a
    rank that writes 20000 lines and then aborts: piped capture returned
    2264946 of 2879970 bytes when the abort was immediate, and 2879970 of
    2879970 on every run when a 1 s sleep preceded the abort. Same writer,
    same capture, different abort timing. How much survives is a scheduling
    race, so on a loaded runner "all of it" is a possible outcome and so is
    "none of it".

    So the two claims are asserted where each is deterministic:

      A. `mpirun -np 2 <binary>` with piped capture -- the real user path.
         Asserts: the refusal exits NON-ZERO, and does not hang. Says
         nothing about the forwarded text; a dropped forward is reported as
         a note, not a failure.
      B. the same refusal with each rank's stdout redirected to its OWN FILE
         by the launched shell, so no forwarding is involved and the bytes
         are in the filesystem the instant the rank's flush returns.
         Asserts: some rank's file carries the FATAL block, and the refusal
         exits NON-ZERO here too.

    No claim is weakened: a build that printed no FATAL block, or that
    forgot the flush (an unflushed buffer dies with the SIGKILLed rank and
    leaves the file empty), or that exited 0, or that hung, still fails --
    it just fails in B, or in both. The only thing no longer gated is
    whether Open MPI's forwarder happens to drain in time, which is a
    property of the launcher and not of EQdyna.

    Reported as SKIPPED, never as PASS, when the binary or mpirun is absent,
    so an unbuilt tree cannot look like a green check.
    """
    # bin/eqdyna FIRST: install-eqdyna.sh moves the binary there, so after a
    # normal build src/fortran/eqdyna does not exist and this probe was
    # skipping on every ordinary tree -- the strongest check in this file,
    # silently absent, on a repo that HAD a built binary two directories away.
    # src/fortran/eqdyna is still accepted for a bare `make` with no install.
    mpirun = os.environ.get('EQDYNA_MPIRUN', 'mpirun')
    binary = None
    for cand in (os.path.join(ROOT, 'bin', 'eqdyna'),
                 os.path.join(SRC, 'eqdyna')):
        if os.path.exists(cand):
            binary = cand
            break
    if binary is None:
        print('  real-binary probe SKIPPED (no bin/eqdyna or %s; build it to enable)'
              % os.path.relpath(os.path.join(SRC, 'eqdyna'), ROOT))
        return None
    want = None
    for raw in open(REGISTRY, errors='replace'):
        m = CODE_DECL.match(raw)
        if m and m.group(1) == 'ERR_INPUT_FILE_MISSING':
            want = int(m.group(2))
    if want is None:
        return ['ERR_INPUT_FILE_MISSING is not in the registry']
    problems = []

    # --- A: the real user path -- exit status and no-hang ---------------
    with tempfile.TemporaryDirectory() as d:
        try:
            p = subprocess.run([mpirun, '-np', '2', binary], cwd=d,
                               capture_output=True, timeout=120)
        except FileNotFoundError:
            print('  real-binary probe SKIPPED (%s not found)' % mpirun)
            return None
        except subprocess.TimeoutExpired:
            return ['the refused run HUNG under %s instead of aborting -- this'
                    ' is the exact failure MPI_Abort exists to prevent' % mpirun]
    out = (p.stdout + p.stderr).decode(errors='replace')
    if p.returncode == 0:
        problems.append('a refused run exited 0 under %s' % mpirun)
    elif p.returncode != want:
        # Non-zero but not the registry code: still a correct refusal, and
        # under srun/ibrun the number is expected to be unreliable.
        print('  real-binary probe A: exited %d (non-zero, but not %d -- '
              'acceptable, launchers may remap)' % (p.returncode, want))
    else:
        print('  real-binary probe A: mpirun -np 2 on an empty case exits %d '
              '(ERR_INPUT_FILE_MISSING), no hang' % p.returncode)
    if 'EQdyna: FATAL' not in out:
        # NOT a failure: see the docstring. Kept visible so the launcher's
        # drop rate stays observable instead of being silently absorbed.
        print('  note: %s forwarded %d byte(s) and no FATAL text on this run '
              '(Open MPI drops un-drained output at MPI_Abort teardown); the '
              'text is asserted in probe B' % (mpirun, len(out)))

    # --- B: same refusal, capture that cannot lose bytes ----------------
    # `exec` so the shell is replaced and mpirun still sees the binary's own
    # exit status; `$$` so each rank names its own file under any launcher.
    shell = 'exec %s > eqdyna-stdout.$$.txt 2>&1' % shlex.quote(binary)
    with tempfile.TemporaryDirectory() as d:
        try:
            q = subprocess.run([mpirun, '-np', '2', 'sh', '-c', shell], cwd=d,
                               capture_output=True, timeout=120)
        except subprocess.TimeoutExpired:
            problems.append('the refused run HUNG under %s with redirected'
                            ' output instead of aborting' % mpirun)
            return problems
        files = sorted(glob.glob(os.path.join(d, 'eqdyna-stdout.*.txt')))
        texts = [open(f, errors='replace').read() for f in files]
    nbytes = sum(len(t) for t in texts)
    if not files:
        problems.append('no rank wrote a stdout file under %s -- the binary'
                        ' never ran, so probe B proved nothing' % mpirun)
    elif not [t for t in texts if 'EQdyna: FATAL' in t]:
        problems.append('the refused run printed no FATAL block (%d rank'
                        ' file(s), %d byte(s), none containing it)'
                        % (len(files), nbytes))
    else:
        print('  real-binary probe B: FATAL block present in a rank-owned file'
              ' (%d file(s), %d byte(s)), exit %d'
              % (len(files), nbytes, q.returncode))
    if q.returncode == 0:
        problems.append('a refused run exited 0 under %s (redirected probe B)'
                        % mpirun)
    return problems


def main():
    print('Regression guard: fatal exits must carry a non-zero status')
    probe_compiler()
    real = probe_real_binary()

    rc = 0
    if real is None:
        # The docstring for probe_real_binary promises this is "SKIPPED,
        # never PASS" -- `if real:` alone treated None the same as an empty
        # (passing) list, which folded a missing prerequisite into a green
        # exit. bin/eqdyna (or src/fortran/eqdyna) and mpirun are both
        # supposed to exist in this tier's declared, built environment
        # (testsys/run.py and CI build the binary first), so their absence
        # here is itself a failure, not an optional environment.
        print('\nFAIL: real-binary probe SKIPPED -- bin/eqdyna (or '
              'src/fortran/eqdyna) and mpirun are required in this tier\'s '
              'build environment; build the binary (./install-eqdyna.sh) '
              'and ensure mpirun is on PATH to enable this probe.')
        rc = 1
    elif real:
        print('\nFAIL: the real binary does not refuse correctly:')
        for r in real:
            print('  ' + r)
        rc = 1

    sc = self_check()
    if sc:
        print('\nFAIL: the scanner itself is broken -- its PASS would be vacuous:')
        for p in sc:
            print('  ' + p)
        rc = 1
    else:
        print('  self-check: scanner classifies all known defect forms correctly')

    result = check_registry()
    if isinstance(result, list):
        problems, ncodes = result, 0
    else:
        problems, ncodes = result
    if problems:
        print('\nFAIL: error-code registry:')
        for p in problems:
            print('  ' + p)
        rc = 1
    else:
        print('  registry: %d codes, all unique and within 1-125' % ncodes)

    drift = check_or_update_docs(update='--update' in sys.argv)
    if drift:
        print('\nFAIL: documentation drift:')
        for d in drift:
            print('  ' + d)
        rc = 1
    else:
        print('  docs/user/troubleshooting.md exit-code table matches the registry')

    offenders = scan_sources()
    if offenders:
        print('\nFAIL: %d stop statement(s) that exit 0:' % len(offenders))
        for rel, lineno, form, code in offenders:
            print('  %s:%d  [%s]  %s' % (rel, lineno, form, code))
        print('\nUse `call abortRun(ERR_..., \'what was wrong and what to do\')`'
              '\n  from src/fortran/errorCodes.f90. If a `stop` really does end a'
              ' SUCCESSFUL run,\n  mark it with a trailing `! %s: ...` comment.'
              % ALLOW_MARKER)
        rc = 1

    if rc == 0:
        nfiles = len(glob.glob(os.path.join(SRC, '*.f90')))
        print('\nPASS: no exit-0 fatal paths in %d src/*.f90 files' % nfiles)
    return rc


if __name__ == '__main__':
    sys.exit(main())
