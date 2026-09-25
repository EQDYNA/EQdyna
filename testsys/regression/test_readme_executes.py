#! /usr/bin/env python3
"""
The stranger-clone gate, made mechanical (mission: README EXECUTES, not just
resolves). `test_readme_commands.py` (this directory) only checks that a
README command's script/tier/module/compset NAME exists -- it would pass a
README whose `python3 -m eqdyna ...` line raised ModuleNotFoundError, because
"import eqdyna resolves" is a different question from "import eqdyna
SUCCEEDS with the exact PYTHONPATH the README tells a user to export". This
file answers the second question by actually running the README, verbatim,
as a brand-new clone in a bare environment would.

TWO MODES, one file, one engine -- so the parser and the execution engine are
never two implementations to keep in sync.

  FAST (default; this is what `run.py regression` / ci_shard run, and what
  is registered in a shard): parses README.md's fenced blocks (cheap), and
  runs the mutation self-test on a tiny SYNTHETIC README (a handful of
  `echo`/`export` lines, sub-second, no clone, no network, no install). This
  is "the fast regression check that runs only the parser and the mutation
  self-test" the mission calls for. Under a second; safe for the seconds
  budget (rule 9).

  FULL (`EQDYNA_README_GATE=full` in the environment; this is what
  `python3 testsys/run.py readme` sets): additionally clones THIS commit
  into a temp dir and executes README.md's Requirements/Install/Quick
  start/Python-solver fenced blocks verbatim, in one bash session, in a
  bare env -i-equivalent environment (fresh HOME, PATH=/usr/bin:/bin plus
  only what the README itself exports). ~2-4 minutes -- this is why it is
  never the default; a fresh checkout's `run.py regression` must stay
  seconds-scale. The mutation self-test still runs FIRST in this mode too
  (mission: "built in and run first") -- if the harness itself cannot tell
  a broken README from a working one, there is no point spending the 2-4
  minutes on the real one.

WHY env -i-EQUIVALENT, NOT THE REAL SHELL. This developer's own login shell
carries PYTHONPATH, EQDYNAROOT and a project venv already exported (this
very box's ambient `sys.executable` is `.../venv_cotopaxi/bin/python3`, not
`/usr/bin/python3`) -- running the README through THAT environment would
silently paper over exactly the class of bug this gate exists to catch (a
README step that only works because the tester's shell was already
configured). Passing `env={'HOME':..., 'PATH': '/usr/bin:/bin'}` to
subprocess.run IS `env -i` plus those two vars -- no child inherits
anything else.

PIP LINES: RUN THEM, not marked. `pip install numpy netCDF4 matplotlib
xarray` (Requirements) and the install step are idempotent by pip's own
contract (rerunning a satisfied install is a fast no-op; rerunning an
unsatisfied one is exactly what a stranger clone needs), and this box has
network reachability for PyPI. There is nothing here for a `# needs root`
marker to protect against -- the ONLY line in the README that needs root is
the `apt-get` line, and it is already marked. Marking the pip lines too
would make the gate skip the exact kind of step whose omission caused the
real bug class in this file's motivating incident.

Anti-vacuous-green discipline (papercuts): every verdict below names the
block and line it refers to, never a bare pass/fail.
"""
import glob
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
README_PATH = os.path.join(ROOT, 'README.md')

RUN_SECTIONS = ('Requirements', 'Install', 'Quick start', 'Python solver')
NEEDS_ROOT_MARKER = '# needs root'

FAILURES = []


def check(ok, what):
    print('%s -- %s' % ('PASS' if ok else 'FAIL', what))
    if not ok:
        FAILURES.append(what)


# --------------------------------------------------------------- parser ----

def parse_fenced_blocks(text):
    """Every ``` ... ``` fenced block in document order, tagged with the
    nearest preceding `## `/`### ` heading. Returns a list of
    dict(section=str, start_line=int, lines=[(lineno, text), ...])."""
    lines = text.splitlines()
    section = None
    blocks = []
    in_block = False
    current = None
    for i, line in enumerate(lines, start=1):
        if line.strip().startswith('```'):
            if not in_block:
                in_block = True
                current = {'section': section, 'start_line': i, 'lines': []}
            else:
                in_block = False
                blocks.append(current)
                current = None
            continue
        if in_block:
            current['lines'].append((i, line))
        else:
            m = re.match(r'^#{2,3}\s+(.*)$', line)
            if m:
                section = m.group(1).strip()
    return blocks


def selected_blocks(blocks, sections=RUN_SECTIONS):
    return [b for b in blocks if b['section'] in sections]


# ----------------------------------------------------------- the engine ----

def build_session_script(blocks, clone_substitution=None):
    """Turn the selected fenced blocks into ONE bash script: state (cd,
    export) carries across blocks because it is one process. Returns
    (script_text, executed_steps, skipped_steps); each step is
    dict(section, line, text).

    `clone_substitution`, if given, is (needle, replacement, announce) --
    the ONE literal substring substitution the mission permits (the GitHub
    clone URL), applied to a matching line and printed via an explicit
    `echo` in the generated script so the substitution is visible in the
    run log, not silently baked in.
    """
    out = ['set -u', 'set -o pipefail']   # unbound-variable use (e.g. a removed export) is a
                       # hard failure, not a silent empty string -- this is
                       # what makes the mutation self-test's "export
                       # removed" case actually fail.
    executed, skipped = [], []
    for b in blocks:
        section = b['section']
        out.append('echo "==== README block: %s (source line %d) ===="'
                    % (section, b['start_line']))
        for lineno, raw in b['lines']:
            line = raw.rstrip()
            if not line.strip():
                continue
            if line.rstrip().endswith(NEEDS_ROOT_MARKER):
                skipped.append(dict(section=section, line=lineno, text=line))
                out.append('echo "SKIP (needs root) README.md:%d"' % lineno)
                continue
            cmd = line
            if clone_substitution is not None:
                needle, replacement, _ = clone_substitution
                if needle in cmd:
                    cmd = cmd.replace(needle, replacement)
                    out.append('echo "SUBSTITUTED README.md:%d: %s -> %s"'
                                % (lineno, shlex.quote(needle), shlex.quote(replacement)))
            executed.append(dict(section=section, line=lineno, text=line))
            tag = 'README.md:%d' % lineno
            out.append('echo "RUN %s"' % tag)
            out.append(cmd)
            out.append('__rc=$?')
            out.append(
                'if [ "$__rc" -ne 0 ]; then '
                'echo "GATE-FAIL section=%s line=%d rc=$__rc"; exit "$__rc"; fi'
                % (json_ish(section), lineno))
    return '\n'.join(out) + '\n', executed, skipped


def json_ish(s):
    """Shell-safe-enough single-token rendering of a section name for the
    GATE-FAIL tag (no shell metacharacters in any of our section names)."""
    return re.sub(r'[^A-Za-z0-9_.-]', '_', s)


def run_session(script_text, cwd, home):
    """Runs the ONE bash session in an env -i-equivalent environment (fresh
    HOME, PATH=/usr/bin:/bin, nothing else -- no PYTHONPATH, no EQDYNAROOT,
    no venv). Returns subprocess.CompletedProcess."""
    os.makedirs(cwd, exist_ok=True)
    os.makedirs(home, exist_ok=True)
    script_path = os.path.join(cwd, '.readme_gate_session.sh')
    with open(script_path, 'w') as fh:
        fh.write(script_text)
    env = {'HOME': home, 'PATH': '/usr/bin:/bin'}
    return subprocess.run(['bash', script_path], cwd=cwd, env=env,
                           capture_output=True, text=True)


def attribute_failure(stdout, stderr, returncode):
    """Names the section/line the generated bash session was executing when
    it failed -- however it failed.

    Two failure shapes need covering, found while building this file's own
    mutation self-test:
      1. A command exits nonzero: our own `__rc` check catches it and prints
         the GATE-FAIL tag -- the precise, authoritative case.
      2. A bash-level abort (`set -u` dereferencing an unbound variable, a
         syntax error) kills the script BEFORE our `__rc` check ever runs.
         This is not a hypothetical: the mission's own 'export removed'
         mutation hits it, because the very next line in the SAME block
         already dereferences the now-unset var, so bash aborts mid-block,
         never reaching our wrapper. Falls back to the last block header /
         RUN tag this file's own generated script had echoed before the
         abort -- still exact, just sourced from the trail instead of a tag.
    """
    combined = stdout + '\n' + stderr
    m = re.search(r'GATE-FAIL section=(\S+) line=(\d+) rc=(-?\d+)', combined)
    if m:
        return dict(section=m.group(1), line=int(m.group(2)),
                    rc=int(m.group(3)), via='GATE-FAIL')
    if returncode == 0:
        return None
    section, line = None, None
    for ln in stdout.splitlines():
        hm = re.match(r'==== README block: (.*) \(source line \d+\) ====', ln)
        if hm:
            section = json_ish(hm.group(1))
        rm = re.match(r'RUN README\.md:(\d+)', ln)
        if rm:
            line = int(rm.group(1))
    return dict(section=section, line=line, rc=returncode, via='last-RUN-fallback')


# ------------------------------------------------------- mutation self-test

SYNTHETIC_README = """# Synthetic

## Requirements

```
echo requirements-marker-skipped  # needs root
echo requirements-ran
```

## Install

```
export GREETING=hello
echo "install saw: $GREETING"
```

## Quick start

```
echo "quickstart saw: $GREETING"
```

### Python solver

```
echo "pysolver saw: $GREETING"
```
"""


def _mutate_drop_export(text):
    """Mirrors the mission's own example: 'the PYTHONPATH export removed'.
    A downstream block then references an unset var; `set -u` makes that a
    hard failure at exactly that line, not a silently-empty string."""
    return text.replace('export GREETING=hello\n', '')


def _mutate_typo_command(text):
    """Mirrors the mission's other example: 'a typo in create.newcase' --
    here, typo the synthetic README's own stand-in command so bash reports
    it as "command not found" (rc=127)."""
    return text.replace('echo "install saw: $GREETING"',
                         'echoo "install saw: $GREETING"')


def _run_synthetic(text, tmp_parent):
    blocks = selected_blocks(parse_fenced_blocks(text))
    script, executed, skipped = build_session_script(blocks)
    work = os.path.join(tmp_parent, 'work')
    home = os.path.join(tmp_parent, 'home')
    proc = run_session(script, work, home)
    return proc, executed, skipped


def check_mutation_self_test():
    """Built-in, run-first harness self-check (mission step 4): an intact
    synthetic README passes; each of two independently-broken copies fails
    AND the failure is attributed to the exact section/line that was broken
    -- proving find_gate_fail's attribution, not just bash's exit code."""
    tmp = tempfile.mkdtemp(prefix='readme_gate_selftest_')
    try:
        proc, executed, skipped = _run_synthetic(SYNTHETIC_README, os.path.join(tmp, 'intact'))
        check(proc.returncode == 0,
              'mutation self-test (intact synthetic README): exit 0 '
              '(got %d); stdout tail: %s' % (proc.returncode, proc.stdout.strip().splitlines()[-3:]))
        check(len(executed) == 5 and len(skipped) == 1,
              'mutation self-test (intact): 5 executed + 1 skipped-by-marker '
              'line (got executed=%d skipped=%d)' % (len(executed), len(skipped)))
        check('pysolver saw: hello' in proc.stdout,
              'mutation self-test (intact): state (GREETING) propagated '
              'across blocks in ONE session (stdout has "pysolver saw: hello")')

        broken_export = _mutate_drop_export(SYNTHETIC_README)
        proc2, _, _ = _run_synthetic(broken_export, os.path.join(tmp, 'broken_export'))
        fail2 = attribute_failure(proc2.stdout, proc2.stderr, proc2.returncode)
        check(proc2.returncode != 0,
              'mutation self-test (export removed): the mutated copy FAILS '
              '(got exit %d)' % proc2.returncode)
        # The removed export's own block (Install) dereferences $GREETING on
        # the very next line, so `set -u` aborts bash there -- before our
        # own GATE-FAIL wrapper runs. attribute_failure falls back to the
        # last RUN tag this script echoed, which is exactly that line.
        check(fail2 is not None and fail2['section'] == 'Install' and fail2['via'] == 'last-RUN-fallback',
              'mutation self-test (export removed): failure is attributed '
              'to the Install block via the last-RUN fallback (bash\'s own '
              '`set -u` abort, not our __rc wrapper), got %r; stderr tail: %s'
              % (fail2, proc2.stderr.strip().splitlines()[-3:]))

        broken_typo = _mutate_typo_command(SYNTHETIC_README)
        proc3, _, _ = _run_synthetic(broken_typo, os.path.join(tmp, 'broken_typo'))
        fail3 = attribute_failure(proc3.stdout, proc3.stderr, proc3.returncode)
        check(proc3.returncode != 0,
              'mutation self-test (typo\'d command): the mutated copy FAILS '
              '(got exit %d)' % proc3.returncode)
        check(fail3 is not None and fail3['section'] == 'Install' and fail3['rc'] == 127,
              'mutation self-test (typo\'d command): failure is attributed '
              'to the Install block with rc=127 (command not found), got %r'
              % fail3)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------- fast checks

def check_readme_has_the_four_sections():
    blocks = parse_fenced_blocks(open(README_PATH, errors='replace').read())
    found = {b['section'] for b in blocks} & set(RUN_SECTIONS)
    missing = sorted(set(RUN_SECTIONS) - found)
    check(not missing,
          'README.md: every gated section (%s) has at least one fenced '
          'block (missing: %r)' % (', '.join(RUN_SECTIONS), missing))


def check_e2e_full_never_fenced():
    """Hours-long steps are never in a fenced block a stranger-clone gate
    would run by default."""
    text = open(README_PATH, errors='replace').read()
    blocks = parse_fenced_blocks(text)
    in_fence = [b for b in blocks for _, l in b['lines'] if 'e2e-full' in l]
    check(not in_fence,
          'README.md: "e2e-full" never appears inside a fenced code block '
          '(found in %d block(s)) -- an hours-long step must stay prose/'
          'inline, never something a copy-paste-the-block user runs by '
          'accident' % len(in_fence))
    check('e2e-full' in text,
          'README.md: "e2e-full" is still documented somewhere (as prose), '
          'so this check is not vacuously true')


def check_only_marked_lines_are_skipped():
    """Self-consistency: the ONLY reason build_session_script ever skips a
    line is the literal `# needs root` suffix -- assert that invariant
    against the real README's selected blocks directly, not by trusting the
    function's own code to have no other branch."""
    blocks = selected_blocks(parse_fenced_blocks(open(README_PATH, errors='replace').read()))
    _, executed, skipped = build_session_script(blocks)
    bad = [s for s in skipped if not s['text'].rstrip().endswith(NEEDS_ROOT_MARKER)]
    check(not bad,
          'README.md: every skipped line ends in the literal marker %r '
          '(bad: %r)' % (NEEDS_ROOT_MARKER, bad))
    check(len(skipped) >= 1,
          'README.md: at least one line is marked %r and skipped (got %d) '
          '-- if this drops to 0 the check above is vacuous' % (NEEDS_ROOT_MARKER, len(skipped)))


FAST_CHECKS = [check_readme_has_the_four_sections, check_e2e_full_never_fenced,
               check_only_marked_lines_are_skipped, check_mutation_self_test]


# ------------------------------------------------------------ full-mode ----

def parse_documented_outputs(text):
    """The README's own bullet list of quick-start outputs
    (`* \\`pattern\\`, ...`) between the Quick start fenced block and the
    next heading. Returns a list of glob patterns, e.g.
    ['cRuptureDynamics.png', 'frt.txt*', ...]."""
    lines = text.splitlines()
    patterns = []
    in_list = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith('* `'):
            in_list = True
        elif in_list and stripped and not stripped.startswith('*'):
            break
        if not in_list:
            continue
        for m in re.finditer(r'`([^`]+)`', line):
            token = m.group(1)
            # A star can sit anywhere before the extension (`frt.txt*` as
            # well as `faultst*.txt`) -- one shared pattern for all three
            # extensions the Quick start bullet list documents.
            if re.match(r'^[\w*.-]+\.(png|txt|nc)\*?$', token):
                patterns.append(token)
    seen = set()
    out = []
    for p in patterns:
        if p not in seen:
            seen.add(p)
            out.append(p)
    return out


def run_full_readme_gate():
    """The heavy real gate: clone THIS commit, run README.md's fenced
    blocks verbatim in a bare environment, verify the documented outputs.
    On failure, the temp dir is kept and its path printed (Box rule)."""
    sha = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT,
                          capture_output=True, text=True).stdout.strip()
    check(bool(sha), 'git rev-parse HEAD (this commit under test) resolved (%r)' % sha)

    text = open(README_PATH, errors='replace').read()
    blocks = selected_blocks(parse_fenced_blocks(text))
    check(len(blocks) >= 4,
          'README.md: >=4 fenced blocks selected across %s (got %d)'
          % (RUN_SECTIONS, len(blocks)))

    tmp = tempfile.mkdtemp(prefix='eqdyna_readme_gate_')
    work = os.path.join(tmp, 'work')
    home = os.path.join(tmp, 'home')
    clone_sub = ('git clone https://github.com/EQDYNA/EQdyna.git',
                 'git clone %s EQdyna && git -C EQdyna checkout %s' % (shlex.quote(ROOT), sha),
                 True)
    script, executed, skipped = build_session_script(blocks, clone_substitution=clone_sub)

    print('\n---- README-executes FULL gate ----')
    print('clone substitution: %r -> %r' % (clone_sub[0], clone_sub[1]))
    print('work dir: %s' % work)
    t0 = time.time()
    proc = run_session(script, work, home)
    wall = time.time() - t0
    print(proc.stdout)
    if proc.stderr:
        print('---- stderr ----')
        print(proc.stderr)
    print('---- wall time: %.1fs, exit %d ----' % (wall, proc.returncode))

    fail = attribute_failure(proc.stdout, proc.stderr, proc.returncode)
    if proc.returncode != 0:
        FAILURES.append('README FULL gate: exit %d, %s (kept: %s)'
                         % (proc.returncode, fail or 'no GATE-FAIL tag found', work))
        print('FAIL -- README FULL gate: exit %d at %s (work dir kept: %s)'
              % (proc.returncode, fail, work))
        return
    print('PASS -- README FULL gate: all %d executed line(s) across %s exited 0 '
          '(%d skipped-by-marker, %.1fs)'
          % (len(executed), RUN_SECTIONS, len(skipped), wall))

    case_dir = os.path.join(home, 'runs', 'tpv8')
    outputs = parse_documented_outputs(text)
    check(len(outputs) >= 4,
          'README.md: quick-start output bullet list parsed >=4 file '
          'pattern(s) (got %r)' % outputs)
    for pattern in outputs:
        matches = glob.glob(os.path.join(case_dir, pattern))
        nonempty = [p for p in matches if os.path.getsize(p) > 0]
        check(bool(nonempty),
              'quick start (%s): documented output %r exists and is '
              'non-empty (matched %d file(s), %d non-empty)'
              % (case_dir, pattern, len(matches), len(nonempty)))

    jax_case_dir = os.path.join(home, 'runs', 'tpv8-jax')
    frt0 = os.path.join(jax_case_dir, 'frt.txt0')
    ok = os.path.exists(frt0) and os.path.getsize(frt0) > 0
    check(ok, 'python solver (%s): python3 -m eqdyna wrote a non-empty '
               'frt.txt0 (%s)' % (jax_case_dir, frt0))

    if not FAILURES:
        shutil.rmtree(tmp, ignore_errors=True)
        print('work dir removed (gate passed): %s' % work)
    else:
        print('work dir KEPT (gate failed somewhere): %s' % tmp)


# --------------------------------------------------------------------- main

def main():
    full = os.environ.get('EQDYNA_README_GATE') == 'full'
    print('Regression guard: README.md must EXECUTE, not just resolve '
          '(mode: %s)' % ('FULL (real clone + real run)' if full
                           else 'FAST (parser + synthetic mutation self-test)'))
    for c in FAST_CHECKS:
        c()
    if FAILURES:
        print('\nFAIL test_readme_executes (%d check(s), fast phase) -- '
              'refusing to spend the FULL gate\'s 2-4 minutes on a harness '
              'that already failed its own self-test' % len(FAILURES))
        return 1
    if full:
        run_full_readme_gate()
    else:
        print('\n(FAST mode: skipping the real clone+run. Set '
              'EQDYNA_README_GATE=full, or run `python3 testsys/run.py '
              'readme`, for the real ~2-4 minute gate.)')
    if FAILURES:
        print('\nFAIL test_readme_executes (%d check(s))' % len(FAILURES))
        return 1
    print('\nSUCCESS test_readme_executes')
    return 0


if __name__ == '__main__':
    sys.exit(main())
