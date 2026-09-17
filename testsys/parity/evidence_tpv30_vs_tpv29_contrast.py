#! /usr/bin/env python3
"""
TPV29 (elastic) vs TPV30 (Drucker-Prager viscoplastic) contrast, at matching
dx, from two COMPLETED fortran e2e runs of the two gated cases.

REPORT-ONLY, same pattern as evidence_tpv29_scec_comparison.py and
evidence_drv_a6_chaos.py in this directory: never a gate, never wired into
testsys/run.py, invoked by hand. Prints the metrics and an
acausal-rupture-time check, and always exits 0 -- a DRIFTED/CONTAMINATED
line is a prompt to tell a human, not a test failure.

Promoted from scratch/tpv30/compareTpv29Tpv30.py (rule 17 step 6 -- TPV29's
own cross-code validation was prose-only with a gitignored baseline, which is
exactly what this script and evidence_tpv29_scec_comparison.py exist to stop
happening again). Rewritten against the CURRENT compsets, not the old
scratch/tpv30/case_input_draft draft: no more G6 half-traction language (that
blocker does not reproduce -- see case_input/test.tpv30/README.md and
pathway_forward item 24(a)), and it reads the real par.dx / par.faultGeometry*
declarations from each case's OWN user_defined_params.py rather than
hardcoding the fault frame.

WHAT IS COMPARED
-----------------
Two completed EQdyna runs of test.tpv29 and test.tpv30, at the SAME dx and
rank count (default: the e2e sweep's own `test/test.tpv29` and
`test/test.tpv30` fortran-cell output directories -- this script does NOT
launch a run itself, rule 9; produce them with
`python3 testsys/e2e/run_e2e.py --cases test.tpv29,test.tpv30 --backends
fortran` if they are not already there). Reads each case's own
`frt.txt<rank>` files (one per MPI rank) directly -- the same raw per-node
on-fault output the sweep's own frt_canonical.py canonicalizes, but gridded
onto the (nz, nx) fault mesh here for the spatial contrast plots.

METRICS (identical formulas to the old scratch script, so this is a
promotion, not a rewrite of the physics):
  * ruptured area (rupture-time field AND independently, final-slip > 1 cm)
  * max / mean final slip, max |along-strike slip|, max peak slip rate
  * seismic moment and Mw, mu = rho*Vs^2 with the case's OWN rho (2670) --
    NOT scripts/lib.py's old hardcoded 2800 (fixed separately by item 24(e),
    lib.shearModulusFromPar; this script computes its own mu directly so it
    does not depend on that fix either way)
  * acausal rupture times: nodes whose recorded rupture time is EARLIER than
    a P wave leaving the hypocenter could possibly have reached them. No
    physical rupture process produces these; any non-zero count means the
    rupture-time (cplot) field is contaminated by a spontaneous failure that
    did not originate at the hypocenter -- this is what the old scratch
    draft's G6 blocker showed up as, and it is checked here on the CURRENT
    code so a future regression of the same shape is caught, not re-derived.

Usage:
    python3 testsys/parity/evidence_tpv30_vs_tpv29_contrast.py
    python3 testsys/parity/evidence_tpv30_vs_tpv29_contrast.py \\
        --tpv29-dir /path/to/tpv29/run --tpv30-dir /path/to/tpv30/run
    python3 testsys/parity/evidence_tpv30_vs_tpv29_contrast.py --no-plots

Writes two PNGs (rupture-time/slip contour, slip profiles) and a timestamped
JSON snapshot to testsys/parity/evidence_output/ (gitignored -- these numbers
are tied to whatever completed runs currently sit in test/, not golden
reference data; rule 4/rule 7), unless --no-plots / --no-json is given.
"""
import argparse
import datetime
import glob
import json
import os
import socket
import subprocess
import sys

import numpy as np

TESTSYS = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(TESTSYS))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)
OUT_DIR = os.path.join(TESTSYS, 'evidence_output')

from testsys import matrix  # noqa: E402

RUN_DIR_DEFAULT = {
    'test.tpv29': os.path.join(REPO_ROOT, 'test', 'test.tpv29'),
    'test.tpv30': os.path.join(REPO_ROOT, 'test', 'test.tpv30'),
}

RHO, VS, VP = 2670.0, 3464.0, 6000.0   # both cases' own material (spec p.11)
MU = RHO * VS ** 2
NORUPT = 9.0e4                          # EQdyna's "never ruptured" sentinel is 99999


def _dx_and_nproc(case):
    """par.dx and the fortran rank count this case is gated at, read from the
    case's own compset -- not hardcoded, so a future re-gate at a different dx
    changes this script's expectations automatically."""
    case_input = os.path.join(REPO_ROOT, 'case_input', case)
    sys.path.insert(0, case_input)
    try:
        saved = {k: sys.modules.pop(k) for k in
                  ('user_defined_params', 'defaultParameters', 'lib')
                  if k in sys.modules}
        from user_defined_params import par  # noqa
        dx = float(par.dx)
        fxmin, fxmax = float(par.fxmin), float(par.fxmax)
        fzmin, fzmax = float(par.fzmin), float(par.fzmax)
    finally:
        sys.path.remove(case_input)
        for k, v in saved.items():
            sys.modules[k] = v
        for k in ('user_defined_params',):
            sys.modules.pop(k, None)
    return dx, fxmin, fxmax, fzmin, fzmax, matrix.FORTRAN_RANKS[case]


def load_trial(run_dir, nproc, dx, fxmin, fxmax, fzmin, fzmax):
    """Grid frt.txt<rank> onto the (nz, nx) fault mesh. Returns dict of fields."""
    files = [os.path.join(run_dir, 'frt.txt%d' % i) for i in range(nproc)]
    missing = [f for f in files if not os.path.isfile(f)]
    if missing:
        raise SystemExit(
            '%s: missing %s -- this run is not complete or was not built at '
            '%d ranks.' % (run_dir, missing[0], nproc))
    a = np.vstack([np.loadtxt(f) for f in files])
    nx = round((fxmax - fxmin) / dx) + 1
    nz = round((fzmax - fzmin) / dx) + 1
    out = {k: np.full((nz, nx), np.nan) for k in
           ('rupt', 'slip', 'slipS', 'slipD', 'peakSr', 'tnrm', 'tstk')}
    ii = np.rint((a[:, 0] - fxmin) / dx).astype(int)
    jj = np.rint((a[:, 2] - fzmin) / dx).astype(int)
    out['rupt'][jj, ii] = a[:, 3]
    out['slip'][jj, ii] = np.hypot(a[:, 4], a[:, 5])
    out['slipS'][jj, ii] = a[:, 4]
    out['slipD'][jj, ii] = a[:, 5]
    out['peakSr'][jj, ii] = a[:, 9]
    out['tnrm'][jj, ii] = a[:, 11] / 1.0e6
    out['tstk'][jj, ii] = a[:, 12] / 1.0e6
    out['x'] = fxmin + np.arange(nx) * dx
    out['z'] = fzmin + np.arange(nz) * dx
    out['dx'] = dx
    if np.isnan(out['rupt']).any():
        raise SystemExit('%s: %d fault nodes missing from frt.txt<rank> -- '
                          'gridding is incomplete (wrong dx/nproc?)'
                          % (run_dir, int(np.isnan(out['rupt']).sum())))
    return out


def spurious_mask(f, xsource, zsource):
    """Acausal failures: nodes whose recorded rupture time is earlier than a P
    wave leaving the hypocenter could possibly have reached them."""
    X, Z = np.meshgrid(f['x'], f['z'])
    r = np.hypot(X - xsource, Z - zsource)
    return (f['rupt'] < NORUPT) & (f['rupt'] < r / VP)


def metrics(f, xsource, zsource, slip_thresh=0.01):
    dx = f['dx']
    cell = dx * dx
    rup = f['rupt'] < NORUPT
    sl = f['slip'] > slip_thresh
    mom = float(np.nansum(f['slip'][sl])) * cell * MU
    return dict(
        area_rupt_km2=float(rup.sum()) * cell / 1.0e6,
        area_slip_km2=float(sl.sum()) * cell / 1.0e6,
        max_slip_m=float(np.nanmax(f['slip'])),
        mean_slip_m=float(np.nanmean(f['slip'][sl])) if sl.any() else 0.0,
        max_slip_strike_m=float(np.nanmax(np.abs(f['slipS']))),
        max_peak_slip_rate_ms=float(np.nanmax(f['peakSr'])),
        moment_Nm=mom,
        mw=float(2.0 / 3.0 * np.log10(mom * 1.0e7) - 10.7) if mom > 0 else None,
        t_max_s=float(np.nanmax(f['rupt'][rup])) if rup.any() else None,
        t_min_s=float(np.nanmin(f['rupt'])),
        n_acausal=int(spurious_mask(f, xsource, zsource).sum()),
        n_ruptured=int(rup.sum()),
        n_nodes=int(f['rupt'].size),
    )


def hypocenter(case):
    """(xsource, zsource) from the case's own compset, m."""
    case_input = os.path.join(REPO_ROOT, 'case_input', case)
    sys.path.insert(0, case_input)
    try:
        saved = {k: sys.modules.pop(k) for k in
                  ('user_defined_params', 'defaultParameters', 'lib')
                  if k in sys.modules}
        from user_defined_params import par  # noqa
        xsource, zsource = float(par.xsource), float(par.zsource)
    finally:
        sys.path.remove(case_input)
        for k, v in saved.items():
            sys.modules[k] = v
        sys.modules.pop('user_defined_params', None)
    return xsource, zsource


def provenance():
    sha = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], cwd=REPO_ROOT,
                          capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(['git', 'status', '--porcelain'], cwd=REPO_ROOT,
                                 capture_output=True, text=True).stdout.strip())
    return dict(date_utc=datetime.datetime.utcnow().isoformat() + 'Z',
                host=socket.gethostname(), sha=sha, dirty=dirty)


def fig_compare(a, b, dx, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    X, Z = a['x'] / 1e3, a['z'] / 1e3
    smax = max(np.nanmax(a['slip']), np.nanmax(b['slip']))
    lev = np.arange(0.0, 20.1, 0.5)
    fig, axes = plt.subplots(2, 1, figsize=(11, 8.6), sharex=True, sharey=True)
    for ax, f, name in zip(axes, (a, b),
                            ('TPV29 (elastic)', 'TPV30 (Drucker-Prager viscoplastic)')):
        im = ax.pcolormesh(X, Z, f['slip'], cmap='magma', vmin=0, vmax=smax,
                           shading='nearest')
        rt = np.ma.masked_greater_equal(f['rupt'], NORUPT)
        cs = ax.contour(X, Z, rt, levels=lev, colors='w', linewidths=0.6)
        ax.clabel(cs, cs.levels[::4], fmt='%.0f', fontsize=7, colors='w')
        ax.set_ylabel('depth (km)')
        ax.set_title('%s   dx = %.0f m' % (name, dx), loc='left', fontsize=11)
        fig.colorbar(im, ax=ax, label='final slip (m)', pad=0.015)
    axes[1].set_xlabel('along-strike distance (km)')
    fig.suptitle('TPV29 vs TPV30: final slip with 0.5 s rupture-time contours',
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out, dpi=170)
    plt.close(fig)


def print_report(ma, mb, dx, nproc):
    rows = [
        ('ruptured area (rupture-time field, km2)', 'area_rupt_km2', '%9.1f'),
        ('ruptured area (final slip > 1 cm, km2)', 'area_slip_km2', '%9.1f'),
        ('max final slip (m)', 'max_slip_m', '%9.3f'),
        ('mean slip over slipping area (m)', 'mean_slip_m', '%9.3f'),
        ('max |along-strike slip| (m)', 'max_slip_strike_m', '%9.3f'),
        ('max peak slip rate (m/s)', 'max_peak_slip_rate_ms', '%9.3f'),
        ('seismic moment (N m, mu=rho*Vs^2)', 'moment_Nm', '%9.4g'),
        ('last rupture time on fault (s)', 't_max_s', '%9.3f'),
        ('first rupture time on fault (s)', 't_min_s', '%9.3f'),
        ('ruptured nodes (of total)', 'n_ruptured', '%9d'),
        ('acausal rupture times (nodes)', 'n_acausal', '%9d'),
    ]
    print('\n=== TPV29 vs TPV30 contrast at dx = %.0f m, %d ranks ===' % (dx, nproc))
    print('%-46s%10s%10s%13s' % ('quantity', 'TPV29', 'TPV30', 'TPV30/TPV29'))
    for label, key, fmt in rows:
        va, vb = ma[key], mb[key]
        if va is None or vb is None:
            print('%-46s%10s%10s%13s' % (label, 'n/a', 'n/a', 'n/a'))
            continue
        ratio = (vb / va) if (isinstance(va, (int, float)) and va) else float('nan')
        print(('%-46s' + fmt + fmt + '%13.3f') % (label, va, vb, ratio))
    if ma['mw'] is not None and mb['mw'] is not None:
        print('%-46s%10.3f%10.3f%13.3f' % ('Mw', ma['mw'], mb['mw'], mb['mw'] - ma['mw']))

    print('\nAcausal-rupture-time check (defect signature, per module docstring):')
    for name, m in (('TPV29', ma), ('TPV30', mb)):
        verdict = 'CLEAN' if m['n_acausal'] == 0 else 'CONTAMINATED'
        print('  %-6s: %d acausal node(s) of %d ruptured -- %s'
              % (name, m['n_acausal'], m['n_ruptured'], verdict))

    print('\nExpected physical direction (spec: plasticity dissipates energy,'
          ' so TPV30 should show REDUCED slip/moment relative to TPV29):')
    if mb['moment_Nm'] is not None and ma['moment_Nm']:
        direction = ('as expected (TPV30 < TPV29)'
                     if mb['moment_Nm'] < ma['moment_Nm'] else
                     'UNEXPECTED (TPV30 >= TPV29) -- worth a human look')
        print('  moment ratio TPV30/TPV29 = %.3f -- %s'
              % (mb['moment_Nm'] / ma['moment_Nm'], direction))

    print('\nThis script is REPORT-ONLY (module docstring): it never asserts '
          'and always exits 0.')


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--tpv29-dir', default=None,
                     help='completed test.tpv29 fortran-cell run directory '
                          '(default: %s)' % RUN_DIR_DEFAULT['test.tpv29'])
    ap.add_argument('--tpv30-dir', default=None,
                     help='completed test.tpv30 fortran-cell run directory '
                          '(default: %s)' % RUN_DIR_DEFAULT['test.tpv30'])
    ap.add_argument('--no-plots', action='store_true')
    ap.add_argument('--no-json', action='store_true')
    args = ap.parse_args()

    tpv29_dir = args.tpv29_dir or RUN_DIR_DEFAULT['test.tpv29']
    tpv30_dir = args.tpv30_dir or RUN_DIR_DEFAULT['test.tpv30']
    for case, d in (('test.tpv29', tpv29_dir), ('test.tpv30', tpv30_dir)):
        if not os.path.isdir(d):
            raise SystemExit(
                '%s: no such directory. This script does not launch a run '
                'itself (rule 9) -- produce it with\n'
                '  python3 testsys/e2e/run_e2e.py --cases %s --backends fortran\n'
                'or pass --%s-dir /path/to/completed/run.'
                % (d, case, case.split('.')[-1]))

    dx29, fxmin29, fxmax29, fzmin29, fzmax29, nproc29 = _dx_and_nproc('test.tpv29')
    dx30, fxmin30, fxmax30, fzmin30, fzmax30, nproc30 = _dx_and_nproc('test.tpv30')
    if (dx29, fxmin29, fxmax29, fzmin29, fzmax29, nproc29) != \
       (dx30, fxmin30, fxmax30, fzmin30, fzmax30, nproc30):
        raise SystemExit(
            'test.tpv29 and test.tpv30 compsets do not agree on dx/fault '
            'extent/rank count -- got (%r) vs (%r). This script only compares '
            'matched-resolution runs (spec: material properties are the only '
            'difference).' % ((dx29, fxmin29, fxmax29, fzmin29, fzmax29, nproc29),
                               (dx30, fxmin30, fxmax30, fzmin30, fzmax30, nproc30)))

    prov = provenance()
    print('==== provenance ====')
    for k, v in prov.items():
        print('  %-9s: %s' % (k, v))

    a = load_trial(tpv29_dir, nproc29, dx29, fxmin29, fxmax29, fzmin29, fzmax29)
    b = load_trial(tpv30_dir, nproc30, dx30, fxmin30, fxmax30, fzmin30, fzmax30)
    xsrc29, zsrc29 = hypocenter('test.tpv29')
    xsrc30, zsrc30 = hypocenter('test.tpv30')
    ma = metrics(a, xsrc29, zsrc29)
    mb = metrics(b, xsrc30, zsrc30)
    print_report(ma, mb, dx29, nproc29)

    if not args.no_plots:
        os.makedirs(OUT_DIR, exist_ok=True)
        png = os.path.join(OUT_DIR, 'tpv30_vs_tpv29_dx%d.png' % round(dx29))
        fig_compare(a, b, dx29, png)
        print('\nwrote %s' % png)

    if not args.no_json:
        os.makedirs(OUT_DIR, exist_ok=True)
        ts = datetime.datetime.utcnow().strftime('%Y%m%dT%H%M%SZ')
        out_path = os.path.join(OUT_DIR, 'evidence_tpv30_contrast_%s.json' % ts)
        with open(out_path, 'w') as f:
            json.dump(dict(provenance=prov, dx=dx29, nproc=nproc29,
                            tpv29_dir=tpv29_dir, tpv30_dir=tpv30_dir,
                            tpv29=ma, tpv30=mb), f, indent=2, default=str)
        print('wrote %s' % out_path)
    return 0


if __name__ == '__main__':
    sys.exit(main())
