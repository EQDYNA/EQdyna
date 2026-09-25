#! /usr/bin/env python3
"""
Generates docs/user/parameters.md's parameter-reference list from
scripts/defaultParameters.py's `parameters` class -- the same class
case.setup and user_defined_params.py read defaults from, so this reference
cannot describe a knob that does not exist, or drift silently on a default
that changed underneath it (board row 126).

WHAT IS EXTRACTED, and what is deliberately not.
  Every top-level `name = value` (or `name1, name2 = v1, v2`) assignment
  directly in the `parameters` class body -- including inside its one
  `if nmat == 1 / elif nmat > 1` branch, tagged with that branch's own
  condition -- becomes one reference entry: the parameter's name, its
  literal default, and a short NOTE built only from:
    - the comment on the assignment's own line, and
    - a contiguous comment block immediately above or below it,
  excluding a bounded separator/title/separator "section banner" block
  (used instead to add a heading) and excluding any single comment LINE
  that trips the same internal-reference markers
  testsys/regression/test_user_docs_style.py refuses in a user-facing
  document (rule 26): defaultParameters.py's comments are written for
  developers and occasionally cite an internal board row. A filtered line
  is PRINTED at generation time -- never silently dropped (rule 2).

  A handful of names are internal working state, not something a user
  sets, and are excluded outright: on_fault_vars, mat, fx, fz, nfx, nfz,
  grav (SKIP_NAMES below). The per-node fill loop
  (`for ix, xcoor in enumerate(fx): ...`) is not descended into at all --
  it computes derived per-node arrays FROM the parameters already listed,
  it is not itself configuration.

Usage:
  python3 docs/user/gen_params.py            # regenerate parameters.md in place
  python3 docs/user/gen_params.py --check     # exit 1 if the committed page
                                               # would change; write nothing
"""
import ast
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SOURCE = os.path.join(ROOT, 'scripts', 'defaultParameters.py')
TARGET = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'parameters.md')

BEGIN = ('<!-- BEGIN PARAMETER REFERENCE (generated from '
          'scripts/defaultParameters.py by docs/user/gen_params.py; do not '
          'edit by hand) -->')
END = '<!-- END PARAMETER REFERENCE -->'

SKIP_NAMES = {'on_fault_vars', 'mat', 'fx', 'fz', 'nfx', 'nfz', 'grav'}

# Same class of internal reference test_user_docs_style.py refuses, applied
# here to ONE comment line at a time so a single contaminated line can be
# filtered out of an otherwise user-safe block instead of discarding it whole.
INTERNAL_MARKERS = [
    re.compile(r'\bPR\s*#\d+|\(#\d+\)'),
    re.compile(r'\brules?\s+\d+[a-z]?\b', re.I),
    re.compile(r'pathway_forward|board row|\bitem\s+\d+', re.I),
    re.compile(r'\b(mira|iris|lars|kai|haruto|nadia|sophia|zofia|victor|'
               r'wei-lin|wei lin|dunyu-liu|anya|marta|priya)\b'),
    re.compile(r'(?<![\w/.-])(?=[0-9a-f]*[a-f])(?=[0-9a-f]*[0-9])'
               r'(?:[0-9a-f]{7,12}|[0-9a-f]{40})(?![\w/.-])'),
    re.compile(r'\.(py|f90|sh|md|yml|txt):\d+'),
    # Bare (no :line) internal-source references, stricter than rule 26's
    # own file:line pattern: a Fortran source filename is never something a
    # user of this page edits or runs, so any mention of one is a
    # developer-rationale citation, not user content. Unlike a Python
    # filename (e.g. user_defined_params.py, which a user DOES edit and
    # this page's own intro names), so .py is not blanket-filtered here --
    # instead the specific internal-module citations actually present in
    # scripts/defaultParameters.py's comments are named explicitly.
    re.compile(r'\b[\w/]+\.f90\b'),
    re.compile(r'\btestsys/matrix\.py\b|\blib\.resolveViscoplasticParams\b|'
               r'\bNSTRESS_CONVENTION\b'),
]

SEP_RE = re.compile(r'^#{5,}\s*$')
BANNER_RE = re.compile(r'^#{2,}\s+(.+?)\s+#{2,}$')


def _is_internal(line):
    return any(rx.search(line) for rx in INTERNAL_MARKERS)


def _strip_comment(line):
    return line.strip().lstrip('#').strip()


def _find_banners(lines):
    """{0-indexed lineno of the title line: title text} for every
    separator/title/separator banner, plus the set of every 0-indexed line
    a banner occupies (never also read as a parameter's note)."""
    banners = {}
    consumed = set()
    for i in range(1, len(lines) - 1):
        if SEP_RE.match(lines[i - 1].strip()) and SEP_RE.match(lines[i + 1].strip()):
            m = BANNER_RE.match(lines[i].strip())
            if m:
                banners[i] = m.group(1).strip()
                consumed |= {i - 1, i, i + 1}
    return banners, consumed


def _comment_block(lines, start, step, consumed, taken):
    """Contiguous comment-only lines from `start` moving by `step`
    (+1 or -1), stopping at a non-comment line or one already claimed by a
    banner or a neighbouring parameter. Returned in reading order."""
    out, used = [], set()
    i = start
    while 0 <= i < len(lines):
        raw = lines[i]
        if not raw.strip().startswith('#') or i in consumed or i in taken:
            break
        out.append(raw)
        used.add(i)
        i += step
    if step < 0:
        out.reverse()
    return out, used


def _clean_note(raw_lines, filtered_log, owner):
    """Join a comment block into one note, or drop it ENTIRELY (never just
    the one offending line) when any line trips an internal-reference
    marker. Splicing around a removed mid-sentence line produces a
    grammatically broken, actively confusing sentence -- worse for a user
    than no note at all -- so the unit that is dropped is the whole block,
    logged so the omission is visible rather than silent (rule 2)."""
    texts = [_strip_comment(raw) for raw in raw_lines]
    texts = [t for t in texts if t]
    if not texts:
        return ''
    hit = next((t for t in texts if _is_internal(t)), None)
    if hit is not None:
        filtered_log.append(
            '%s: dropped whole block (%d line(s); contains %r)'
            % (owner, len(texts), hit))
        return ''
    return ' '.join(texts)


def extract(source_text=None):
    """Returns (rows, banners, filtered_log).
    rows    -- [(lineno, name, default_repr, note, conditions)] in source order.
    banners -- [(lineno, title)] in source order.
    """
    text = source_text if source_text is not None else open(SOURCE).read()
    lines = text.splitlines()
    banners, consumed = _find_banners(lines)
    tree = ast.parse(text)
    cls = next(n for n in tree.body
               if isinstance(n, ast.ClassDef) and n.name == 'parameters')

    rows = []
    taken = set()
    filtered_log = []

    def walk(body, conditions):
        for node in body:
            if isinstance(node, ast.If):
                walk(node.body, conditions + [ast.unparse(node.test)])
                walk(node.orelse, conditions)
                continue
            if not isinstance(node, ast.Assign):
                continue                           # For/Expr/etc: not a knob
            target = node.targets[0]
            if isinstance(target, ast.Tuple):
                names = [t.id for t in target.elts if isinstance(t, ast.Name)]
            elif isinstance(target, ast.Name):
                names = [target.id]
            else:
                continue                           # Subscript target etc.
            if isinstance(node.value, ast.Tuple) and len(node.value.elts) == len(names):
                values = node.value.elts
            else:
                values = [node.value] * len(names)

            lineno0, endlineno0 = node.lineno - 1, node.end_lineno - 1
            inline = ''
            if '#' in lines[endlineno0]:
                inline = lines[endlineno0].split('#', 1)[1].strip()

            lead, lead_used = _comment_block(lines, lineno0 - 1, -1, consumed, taken)
            trail, trail_used = _comment_block(lines, endlineno0 + 1, +1, consumed, taken)
            taken.update(lead_used | trail_used)

            for name, value in zip(names, values):
                if name in SKIP_NAMES:
                    continue
                try:
                    default = repr(ast.literal_eval(value))
                except Exception:
                    default = ast.get_source_segment(text, value)
                owner = '%s (scripts/defaultParameters.py line %d)' % (name, node.lineno)

                note_parts = []
                lead_note = _clean_note(lead, filtered_log, owner)
                if lead_note:
                    note_parts.append(lead_note)
                if inline:
                    if _is_internal(inline):
                        filtered_log.append('%s: filtered %r' % (owner, inline))
                    else:
                        note_parts.append(inline)
                trail_note = _clean_note(trail, filtered_log, owner)
                if trail_note:
                    note_parts.append(trail_note)
                note = ' '.join(note_parts)
                if _is_internal(note):
                    raise RuntimeError(
                        'gen_params: %s note still trips an internal-reference '
                        'marker after per-line filtering -- the joiner has a '
                        'bug, refusing to emit a user-facing doc with it: %r'
                        % (owner, note))
                rows.append((node.lineno, name, default, note, list(conditions)))

    walk(cls.body, [])
    ordered_banners = sorted(banners.items())
    return rows, ordered_banners, filtered_log


def render(rows, banners, filtered_log):
    events = [(ln, 'banner', title) for ln, title in banners]
    events += [(ln, 'row', (name, default, note, conditions))
               for (ln, name, default, note, conditions) in rows]
    events.sort(key=lambda e: e[0])

    out = [BEGIN, '',
           'Every entry below is an attribute of the `parameters` class in '
           '`scripts/defaultParameters.py`, read by `case.setup` and '
           'overridable per case in that case\'s own `user_defined_params.py`. '
           'Defaults shown are this class\'s own; a compset may set '
           'different values for its own physics.', '']
    for _, kind, data in events:
        if kind == 'banner':
            out += ['## %s' % data, '']
            continue
        name, default, note, conditions = data
        cond = ' (applies when %s)' % ' and '.join(conditions) if conditions else ''
        header = '* **`%s`**' % name
        one_line = default is not None and '\n' not in default and len(default) <= 100
        if one_line:
            header += ' -- default `%s`%s' % (default, cond)
            out.append(header)
            if note:
                out += ['', '  ' + note]
        else:
            out += [header + cond, '', '  Default:', '', '  ```']
            for l in (default or '').splitlines():
                out.append('  ' + l)
            out.append('  ```')
            if note:
                out += ['', '  ' + note]
        out.append('')
    out.append(END)
    return '\n'.join(out)


def check_or_update(update):
    if not os.path.exists(TARGET):
        return ['docs/user/parameters.md does not exist (run without --check to create it)']
    s = open(TARGET, errors='replace').read()
    if BEGIN not in s or END not in s:
        return ['docs/user/parameters.md has no generated parameter-reference '
                'section (run without --check to insert one)']
    head, rest = s.split(BEGIN, 1)
    _, tail = rest.split(END, 1)
    current = BEGIN + rest.split(END, 1)[0] + END
    rows, banners, filtered_log = extract()
    for f in filtered_log:
        print('  filtered: %s' % f)
    wanted = render(rows, banners, filtered_log)
    if current == wanted:
        return []
    if update:
        open(TARGET, 'w').write(head + wanted + tail)
        print('  docs/user/parameters.md parameter reference regenerated '
              '(%d entries, %d section(s))' % (len(rows), len(banners)))
        return []
    return ["docs/user/parameters.md's parameter reference no longer matches "
            "scripts/defaultParameters.py -- rerun 'python3 docs/user/gen_params.py'"]


def main():
    check = '--check' in sys.argv
    problems = check_or_update(update=not check)
    if problems:
        print('FAIL gen_params:')
        for p in problems:
            print(' -', p)
        return 1
    print('SUCCESS gen_params' + (' (check only, no write)' if check else ''))
    return 0


if __name__ == '__main__':
    sys.exit(main())
