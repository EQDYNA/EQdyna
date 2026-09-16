"""
Milestone 7: standalone port of src/library_output.f90's `output_frt` --
the frt.txt writer, EXACT Fortran '(1x,22e18.7e4)' formatting.

Column layout (22 per fault node, ntotft==1, matching this port's scope),
read directly off output_frt's write statement, in order:
  1-3   meshCoor(1:3, nsmp(1,i,1))      -- slave-node x,y,z coordinates
  4     fnft(i,1)                       -- rupture time
  5-10  fric(71:76,i,1)                 -- slip S/D/N, sliprate S/D/max
  11    fric(47,i,1)                    -- peak/final slip rate
  12-14 fric(78:80,i,1)                 -- Tn, Ts, Td
  15-20 fric(31:36,i,1)                 -- vxm,vym,vzm (master), vxs,vys,vzs (slave)
  21    fric(20,i,1)                    -- RSF state variable (unused, friclaw=1)
  22    fric(23,i,1)                    -- theta_pc (normal-stress state)

Fortran's E18.7E4 edit descriptor: each value is printed normalized-
mantissa style, "0." followed by 7 significant digits (NOT "d.ddddddd" --
Fortran's E format always shows a leading "0", unlike Python/C's "%e"
which normalizes to one nonzero digit before the decimal point), then
"E", an explicit exponent sign, and a 4-digit exponent, right-justified
in an 18-character field with a single leading blank ('1x') before the
first field on the line. This module derives the Fortran-style 7-digit
mantissa from Python's `%.6e` (which gives 1+6=7 significant digits, just
with the decimal point one place to the left of Fortran's convention) by
shifting the exponent by +1 -- not a re-implementation of Fortran's
internal rounding, just an exact re-derivation of the SAME digits Python's
own correctly-rounded `%e` already produces, redistributed into Fortran's
"0.ddddddd" layout.

Verified via ROUND-TRIP against `testsys/parity/fixtures/test_tpv8_serial/
frt.txt0` (a real Fortran-written file, not synthetic): parse every row's
22 floats back out, re-format them with THIS module's writer, and assert
the reproduced text is byte-identical to the original file, via
the removed parity tier (test_standalone_meshgen.py, deleted 2026-09-15). This tests the FORMATTER only
(bit-exact reproduction of Fortran's own quantization) -- the actual
fnft/fric VALUES come from the time-stepping solver (Milestone 8), not
this milestone; M8's acceptance test is where the writer is exercised
end-to-end against fresh solver output.
"""
import numpy as np


def format_fortran_e(x, width=18, decimals=7, exp_digits=4):
    """Port of Fortran's Ew.dEe edit descriptor for a single value.

    x: a Python float (any finite double).
    Returns an 18-character string (default width), right-justified,
    e.g. '  -0.1500000E+0005' or '   0.0000000E+0000'.
    """
    if x == 0.0:
        digits = '0' * decimals
        exp = 0
        sign = '-' if np.signbit(x) else ''
    else:
        sign = '-' if x < 0.0 else ''
        ax = abs(x)
        s = '%.*e' % (decimals - 1, ax)  # 1 + (decimals-1) = `decimals` sig figs
        mantissa, exp_part = s.split('e')
        lead_digit = mantissa[0]
        rest = mantissa[2:] if len(mantissa) > 1 else ''
        digits = (lead_digit + rest).ljust(decimals, '0')[:decimals]
        exp = int(exp_part) + 1
    exp_sign = '+' if exp >= 0 else '-'
    field = '%s0.%sE%s%0*d' % (sign, digits, exp_sign, exp_digits, abs(exp))
    return field.rjust(width)


def format_frt_row(values, width=18, decimals=7, exp_digits=4):
    """Port of `write(unit,'(1x,22e18.7e4)') (values...)` for one row.
    values: iterable of floats (22 for frt.txt's column layout, but not
    hardcoded to 22 -- the '1x' + N-field structure is general).
    Returns the row's text WITHOUT a trailing newline.
    """
    return ' ' + ''.join(format_fortran_e(v, width, decimals, exp_digits) for v in values)


def build_frt_rows(meshCoor, nsmp, fnft, fric):
    """Assembles the 22-column value list for every fault node, in
    output_frt's exact column order (see module docstring). Does NOT
    write to disk -- see write_frt.

    meshCoor: (N+1,3) 1-indexed (row 0 unused).
    nsmp: (nftnd,2) int64 [slave_id, master_id], 1-indexed.
    fnft: (nftnd+1,) rupture time per fault node, 1-indexed (row 0 unused).
    fric: (nftnd+1,101) friction-state array, 1-indexed rows AND columns
        (row/col 0 unused), matching readInputFiles.read_on_fault_vars's
        convention.

    Returns a list of nftnd lists, each 22 floats.
    """
    nftnd = nsmp.shape[0]
    rows = []
    for i in range(1, nftnd + 1):
        slave = int(nsmp[i - 1, 0])
        row = [meshCoor[slave, 0], meshCoor[slave, 1], meshCoor[slave, 2], fnft[i]]
        row += [fric[i, s] for s in (71, 72, 73, 74, 75, 76)]
        row.append(fric[i, 47])
        row += [fric[i, 78], fric[i, 79], fric[i, 80]]
        row += [fric[i, s] for s in (31, 32, 33, 34, 35, 36)]
        row.append(fric[i, 20])
        row.append(fric[i, 23])
        rows.append(row)
    return rows


def write_frt(path, meshCoor, nsmp, fnft, fric):
    """Port of output_frt: writes `path` with one '(1x,22e18.7e4)'-
    formatted row per fault node, terminated by '\\n' (matching Fortran's
    list-directed record terminator for a formatted sequential file).
    Fortran's `if (nftnd(1) > 0)` guard: raises if nsmp is empty rather
    than silently writing an empty file (no fault, no frt.txt in Fortran
    either -- this port makes that an explicit error instead of a
    silent no-op, since a standalone run with zero fault nodes indicates
    a setup problem this port should not paper over)."""
    if nsmp.shape[0] == 0:
        raise ValueError('write_frt: nftnd==0 -- Fortran would skip writing '
                          'frt.txt entirely in this case; refusing to write '
                          'an empty/misleading file')
    rows = build_frt_rows(meshCoor, nsmp, fnft, fric)
    with open(path, 'w') as f:
        for row in rows:
            f.write(format_frt_row(row) + '\n')


def read_frt(path, n_cols=22):
    """Parses a Fortran-written frt.txt file back into a (nftnd, n_cols)
    float array -- used by the round-trip parity test, and generally
    useful for any downstream code that wants frt.txt's values without
    re-deriving output_frt's column layout."""
    rows = []
    with open(path) as f:
        for line in f:
            vals = [float(v) for v in line.split()]
            if len(vals) != n_cols:
                raise ValueError('read_frt: expected %d columns, got %d in line %r'
                                  % (n_cols, len(vals), line))
            rows.append(vals)
    return np.array(rows)
