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
Row 114 (board): output_onfault_st / output_offfault_st, the `faultst*.txt`
/ `body*.txt` SCEC-format station time series (library_output.f90:19-247).

Column derivation for `write_onfault_stations` (output_onfault_st): the
Fortran calls storeOnFaultStationQuantSCEC (faulting.f90:508) with
nsdSlipVector/nsdSliprateVector/nsdTractionVector AFTER solveSWTW/solveRSF
have run for this step, so:
  - fric[SLIP_STRIKE]/[SLIP_DIP]/[SLIPRATE_STRIKE]/[SLIPRATE_DIP] are written
    by getNsdSlipSliprateTraction (faulting.f90:108-112 / faulting.py's
    getNsdSlipSliprateTraction) BEFORE solveRSF's background-slip-rate add
    (faulting.f90:219-220) -- so for friclaw>=3 the creep term
    (fric[VINI_X]/[VINI_Z]) must be added back in AT RECORD TIME to match
    what solveRSF's own local nsdSlipVector/nsdSliprateVector carry into
    storeOnFaultStationQuantSCEC. friclaw<=2 (solveSWTW) never touches those
    slots or vini_*, so no creep term is added there -- matching the
    Fortran's solveSWTW, which never references FRIC_SLOT_VINI_*.
  - fric[TRACT_NORM]/[TRACT_STRIKE]/[TRACT_DIP] are the POST-friction-law
    nsdTractionVector(1:3) in BOTH branches (solveSWTW's clipped n,s,d at
    faulting.f90:183; solveRSF's final n,s,d at faulting.f90:273) -- read
    directly, no further transform needed.
  - fric[STATE]/[TP_TEMP]/[TP_NORM_TP]/[TP_PINI] are read directly (psi,
    temperature, pore pressure -- friclaw>=3 columns only).
  - Column 7 (n-slip) is computed by the Fortran but never appears in
    output_onfault_st's write list (library_output.f90:107-118/130-138) --
    confirmed by reading the write statement's argument order against the
    column-name line two lines above it -- so it is not recorded here.
  - The n-stress sign is `nStressOutSign * tnrm / 1e6`, the case's SCEC
    spec convention read from bGlobal.txt (board row 22a, PR #20); see
    driver.py's build_invariants (finv['st_sign']).

`write_offfault_stations` (output_offfault_st): idhist's dispOrVel is always
1 or 2 (eqdyna3d.f90:287-298's allocInitAfterMeshGen loop never writes 3), so
storeOffFaultStData's quantType==3 (nodalForceArr) branch is dead code and is
not reproduced -- only dispArr/velArr are read, exactly the Fortran's live
path.

Filename/header formatting reproduces Fortran's Iw.m / Fw.d edit descriptors
(_format_fortran_i / _format_fortran_f below), including their OVERFLOW
behavior (a value too wide for its field prints as literal asterisks) and,
for the on-fault strike field specifically, the item-22 fix that writes the
sign as a separate a1 literal so a 4-digit-with-sign value cannot overflow
the 3-digit magnitude sub-field the way a naive negative-i3.3 would.
`projectname`/`author` are globalvar.f90 compile-time literals ('San-Ti',
'Sophon'), not case parameters -- reproduced verbatim, not read from any
input file.
"""
import datetime

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


PROJECTNAME = 'San-Ti'   # globalvar.f90:88 -- compile-time literal, not input
AUTHOR = 'Sophon'         # globalvar.f90:89 -- compile-time literal, not input


def _fortran_nint(x):
    """Fortran's NINT: round half away from zero (not numpy's round-half-
    to-even). `x` is a plain Python/numpy scalar."""
    x = float(x)
    if x == 0.0:
        return 0
    return int(np.sign(x) * np.floor(np.abs(x) + 0.5))


def _format_fortran_i(value, width, mindigits):
    """Fortran's Iw.m integer edit descriptor: at least `mindigits` digits
    (zero-padded), a leading '-' for negative values (consuming one column
    of `width`), and `'*' * width` if the result does not fit -- Fortran's
    own overflow behavior, reproduced rather than guarded around."""
    value = int(value)
    digits = str(abs(value)).zfill(mindigits)
    s = ('-' + digits) if value < 0 else digits
    return '*' * width if len(s) > width else s


def _format_fortran_f(x, width, decimals):
    """Fortran's Fw.d fixed-point edit descriptor: `'*' * width` on overflow,
    else the plain fixed-point string (the callers below always adjustl+trim
    it afterward, exactly as the Fortran does, so left-padding is not
    reproduced -- only the digits and the overflow behavior matter)."""
    s = '%.*f' % (decimals, x)
    return '*' * width if len(s) > width else s


def onfault_filename(strike_m, depth_m, dip_rad):
    """library_output.f90:41-52's faultst<sss>dp<ddd>.txt name.

    `strike_m`/`depth_m` are xonfs's along-strike x and vertical-depth z, in
    meters (this port's bStations.txt units, already *1000 from km -- see
    readInputFiles.read_bstations). `dip_rad` is fltxyz(2,4,1): C_degen*pi/180
    for C_degen>3, else 90*pi/180 (meshgen.py build_fault_geometry's `fdip`,
    recomputed here verbatim rather than plumbed through, since it is a pure
    function of params['C_degen'])."""
    st_val = _fortran_nint(strike_m / 100.0)
    if st_val < 0:
        sttmp = '-' + _format_fortran_i(abs(st_val), 3, 3)
    else:
        sttmp = _format_fortran_i(st_val, 3, 3)
    dp_val = _fortran_nint(abs(depth_m) / np.sin(dip_rad) / 100.0)
    dptmp = _format_fortran_i(dp_val, 3, 3)
    return 'faultst%sdp%s.txt' % (sttmp, dptmp)


def onfault_location_stamp(strike_m, depth_m, dip_rad):
    """library_output.f90:56-69's '# location = on fault, ...' header line."""
    st_km = _format_fortran_f(strike_m / 1000.0, 5, 1).strip()
    dd_km = _format_fortran_f(abs(depth_m) / np.sin(dip_rad) / 1000.0, 5, 1).strip()
    return '# location = on fault, %s km along strike, %s km down-dip' % (st_km, dd_km)


def offfault_filename(x_m, y_m, z_m):
    """library_output.f90:155-165's body<yyy>st<xxx>dp<zzz>.txt name.
    x_m/y_m/z_m: x4nds's along-strike/off-fault/depth coordinates, meters."""
    b = _format_fortran_i(_fortran_nint(y_m / 100.0), 4, 3)
    s = _format_fortran_i(_fortran_nint(x_m / 100.0), 4, 3)
    d = _format_fortran_i(_fortran_nint(abs(z_m) / 100.0), 4, 3)
    return 'body%sst%sdp%s.txt' % (b, s, d)


def offfault_location_stamp(x_m, y_m, z_m):
    """library_output.f90:167-185's '# location = ... km off fault' line."""
    b_km = _format_fortran_f(y_m / 1000.0, 5, 1).strip()
    s_km = _format_fortran_f(x_m / 1000.0, 5, 1).strip()
    d_km = _format_fortran_f(abs(z_m) / 1000.0, 5, 1).strip()
    return '# location = %s km off fault, %s km along strike, %s km depth' % (b_km, s_km, d_km)


def _date_line():
    """library_output.f90:77-80/193-196's ' # date = mm/dd/yyyy hh:mm:ss'
    line. Per this mission's rule 23 scope, only the date/time VALUE is
    exempt from byte parity (it is wall-clock, not reproducible); the line's
    SHAPE matches the Fortran's write statement."""
    now = datetime.datetime.now()
    return ' # date = %d/%d/%d %d:%d:%d' % (
        now.month, now.day, now.year, now.hour, now.minute, now.second)


def _scec_row(t, rest, first_width=21, first_dec=13, rest_width=16, rest_dec=7):
    """One data row of output_onfault_st/output_offfault_st: column 1 in
    Fortran's 'E21.13' (no explicit exponent width -> default 2-digit
    exponent, UNLIKE write_frt's 'E18.7E4'), the remaining columns in
    'E16.7', same default 2-digit exponent. No leading blank (the Fortran
    format strings here carry no '1x', unlike output_frt's)."""
    fields = [format_fortran_e(t, first_width, first_dec, 2)]
    fields += [format_fortran_e(v, rest_width, rest_dec, 2) for v in rest]
    return ''.join(fields)


def write_onfault_stations(case_dir, S, on_st_hist):
    """Port of output_onfault_st. `on_st_hist`: (n_on, ncols, nstep) float
    array from driver.run's/run_mpi's carry (columns already in the exact
    output order -- see this module's docstring). Writes nothing (matches
    Fortran's `if(numOfOnFaultStCount>0)` guard) when there are no on-fault
    stations."""
    import os
    n_on = int(S['st_on_idx'].shape[0])
    if n_on == 0:
        return []
    friclaw = S['friclaw']; dip_rad = S['fault_dip_rad']
    ncols = 11 if friclaw >= 3 else 8
    if on_st_hist.shape[0] != n_on or on_st_hist.shape[1] != ncols:
        raise ValueError('write_onfault_stations: on_st_hist shape %r does not '
                          'match (n_on=%d, ncols=%d)' % (on_st_hist.shape, n_on, ncols))
    # nstep is the array's own step count (a shortened test run may pass
    # nsteps < S['nstep']; the header's num_time_steps line follows the data,
    # matching how nstep flows consistently through the Fortran).
    nstep = on_st_hist.shape[2]
    names8 = 't h-slip h-slip-rate h-shear-stress v-slip v-slip-rate v-shear-stress n-stress'
    names11 = names8 + ' psi temperature pressure'
    paths = []
    for i in range(n_on):
        strike_m = float(S['st_on_strike_m'][i]); depth_m = float(S['st_on_depth_m'][i])
        fname = onfault_filename(strike_m, depth_m, dip_rad)
        path = os.path.join(case_dir, fname)
        with open(path, 'w') as f:
            f.write(onfault_location_stamp(strike_m, depth_m, dip_rad) + '\n')
            f.write('# Project=%s\n' % PROJECTNAME)
            f.write('# Author=%s\n' % AUTHOR)
            f.write(_date_line() + '\n')
            f.write('# code = EQdyna\n')
            f.write('# element_size = %s\n' % float(S['dx']))
            f.write('# time_step =%8.4f s\n' % S['dt'])
            f.write('# num_time_steps =%6d\n' % nstep)
            f.write('# Column #1 = Time (s)\n')
            f.write('# Column #2 = horizontal slip (m)\n')
            f.write('# Column #3 = horizontal slip rate (m/s)\n')
            f.write('# Column #4 = horizontal shear stress (MPa)\n')
            f.write('# Column #5 = down-dip slip (m)\n')
            f.write('# Column #6 = down-dip slip rate (m/s)\n')
            f.write('# Column #7 = down-dip shear stress (MPa)\n')
            f.write('# Column #8 = normal stress (MPa)\n')
            if ncols == 11:
                f.write('# Column #9 = state variable psi (dimensionless)\n')
                f.write('# Column #10 = Temperature (degrees Kelvin)\n')
                f.write('# Column #11 = Pore pressure (MPa)\n')
                f.write('# Time series in 11 columns; column 1 in format E21.13, '
                        'columns 2-11 in format E16.7\n')
                f.write('# The line below lists the names of the data fields:\n')
                f.write(' ' + names11 + '\n')
            else:
                f.write('# Time series in 8 columns; column 1 in format E21.13, '
                        'columns 2-8 in format E16.7\n')
                f.write('# The line below lists the names of the data fields:\n')
                f.write(' ' + names8 + '\n')
            for j in range(nstep):
                row = on_st_hist[i, :, j]
                f.write(_scec_row(row[0], row[1:]) + '\n')
        paths.append(path)
    return paths


def write_offfault_stations(case_dir, S, off_st_hist):
    """Port of output_offfault_st. `off_st_hist`: (n_off, 7, nstep). Writes
    nothing when there are no off-fault stations (Fortran's
    `if(numOfOffFaultStCount>0)` guard)."""
    import os
    n_off = int(S['st_off_idx'].shape[0])
    if n_off == 0:
        return []
    if off_st_hist.ndim != 3 or off_st_hist.shape[:2] != (n_off, 7):
        raise ValueError('write_offfault_stations: off_st_hist shape %r does not '
                         'match (%d stations, 7 columns, nstep)'
                         % (off_st_hist.shape, n_off))
    nstep = off_st_hist.shape[2]
    paths = []
    for i in range(n_off):
        # Row 94 (owner ruling 2026-09-24): the FILE NAME is derived from the
        # REQUESTED coordinate (x_m/y_m/z_m, from bStations.txt) -- the
        # stable, resolution-independent identifier testsys/matrix.py's
        # GATE_STATIONS names files by -- while the header LOCATION STAMP
        # records the ACTUAL matched node (x_actual_m/y_actual_m/z_actual_m),
        # which can legitimately differ now that depth snaps to the nearest
        # node instead of requiring an exact grid-plane match. Mirrors
        # library_output.f90's output_offfault_st exactly (bodytmp/sttmp/
        # dptmp computed twice there too, once from each source).
        x_m = float(S['st_off_x_m'][i]); y_m = float(S['st_off_y_m'][i])
        z_m = float(S['st_off_z_m'][i])
        x_actual_m = float(S['st_off_x_actual_m'][i]); y_actual_m = float(S['st_off_y_actual_m'][i])
        z_actual_m = float(S['st_off_z_actual_m'][i])
        fname = offfault_filename(x_m, y_m, z_m)
        path = os.path.join(case_dir, fname)
        with open(path, 'w') as f:
            f.write(offfault_location_stamp(x_actual_m, y_actual_m, z_actual_m) + '\n')
            f.write('# Project=%s\n' % PROJECTNAME)
            f.write('# Author=%s\n' % AUTHOR)
            f.write(_date_line() + '\n')
            f.write('# code = EQdyna\n')
            f.write('# element_size = %s\n' % float(S['dx']))
            f.write('# time_step=%8.4f s\n' % S['dt'])
            f.write('# num_time_steps=%6d\n' % nstep)
            f.write('# Time series in 7 columns; column 1 in format E21.13, '
                    'columns 2-7 in format E16.7\n')
            f.write('# Column #1 = Time (s)\n')
            f.write('# Column #2 = horizontal displacement (m)\n')
            f.write('# Column #3 = horizontal velocity (m/s)\n')
            f.write('# Column #4 = vertical displacement (m)\n')
            f.write('# Column #5 = vertical velocity (m/s)\n')
            f.write('# Column #6 = normal displacement (m)\n')
            f.write('# Column #7 = normal velocity (m/s)\n')
            f.write('#\n')
            f.write('# The line below lists the names of the data fields:\n')
            f.write(' t h-disp h-vel v-disp v-vel n-disp n-vel\n')
            for j in range(nstep):
                row = off_st_hist[i, :, j]
                f.write(_scec_row(row[0], row[1:]) + '\n')
        paths.append(path)
    return paths


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
