#! /usr/bin/env python3
"""
Regression guard: docs/user/parameters.md's generated parameter-reference
section must match what docs/user/gen_params.py currently produces from
scripts/defaultParameters.py (board row 126) -- the same generated-doc
contract test_stop_exit_status.py already holds for
docs/user/troubleshooting.md's exit-code table.

Cheap (rule 9): an AST parse of one ~330-line file and a string compare,
well under 1 s. Exits non-zero on drift, a missing generated section, or a
gen_params.py bug that lets an internal-reference marker survive into the
rendered note (gen_params.py raises on that; this test lets that exception
surface as a failure rather than swallowing it).
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
GEN_PARAMS = os.path.join(ROOT, 'docs', 'user', 'gen_params.py')


def _load_gen_params():
    import importlib.util
    spec = importlib.util.spec_from_file_location('gen_params', GEN_PARAMS)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    print('Regression guard: docs/user/parameters.md must match gen_params.py')
    if not os.path.exists(GEN_PARAMS):
        print('FAIL: %s does not exist' % os.path.relpath(GEN_PARAMS, ROOT))
        return 1
    mod = _load_gen_params()
    problems = mod.check_or_update(update=False)
    if problems:
        print('FAIL test_params_reference_freshness:')
        for p in problems:
            print(' -', p)
        print('\nRun: python3 docs/user/gen_params.py')
        return 1
    rows, banners, filtered = mod.extract()
    print('  docs/user/parameters.md matches scripts/defaultParameters.py '
          '(%d parameter entries, %d section(s), %d internal-reference '
          'block/line(s) filtered at generation time)'
          % (len(rows), len(banners), len(filtered)))
    print('\nSUCCESS test_params_reference_freshness')
    return 0


if __name__ == '__main__':
    sys.exit(main())
