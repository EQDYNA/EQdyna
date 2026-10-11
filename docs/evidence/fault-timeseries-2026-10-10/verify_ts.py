"""Verify src_ts* against every independent output in a finished case dir.
usage: python3 verify_ts.py CASE_DIR"""
import glob, os, re, sys
import numpy as np
import netCDF4 as nc
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'scripts'))
from readFaultTimeSeries import read_case, read_fault_ts

d = sys.argv[1]
os.chdir(d)
tags = sorted({re.match(r'src_ts_hdr(ft\d+_)?\d+$', f).group(1) or '' for f in os.listdir('.') if f.startswith('src_ts_hdr')})
TS = {tag: read_case('.', tag) for tag in tags}
res = []

def rep(name, a, b, digits=None):
    """digits: significant digits of the E-format text b was printed with;
    ok means every |a-b| is within half a unit in b's last printed digit
    (digits=None: just report; the caller judges)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    diff = np.abs(a - b)
    mx = diff.max() if diff.size else 0.0
    scale = np.abs(b).max() if b.size else 0.0
    if digits is None:
        ok = '-'
    else:
        with np.errstate(divide='ignore'):
            e = np.floor(np.log10(np.where(b == 0, 1e-300, np.abs(b)))) + 1
        ok = bool(np.all(diff <= 0.5 * 10.0 ** (e - digits) * (1 + 1e-9)))
    res.append((name, a.size, mx, scale, ok))

# 1. src_evol (only when outputGroundMotion=1): must equal field 0 bit for bit
for f in sorted(glob.glob('src_evol*')):
    m = re.match(r'src_evol(ft\d+_)?(\d+)$', f)
    tag, rank = m.group(1) or '', m.group(2)
    ts = read_fault_ts('src_ts_hdr%s%s' % (tag, rank))
    ev = np.fromfile(f, '<f8').reshape(-1, ts['meta']['nnode'])
    rep('src_evol%s%s == sliprate (bitwise)' % (tag, rank), ev, ts['sliprate'])
    res[-1] = res[-1][:4] + (np.array_equal(ev, ts['sliprate']),)

# 2. last record == frt.txt* (E18.7: 8 significant digits)
nfield = {tag: [fd['name'] for fd in TS[tag]['fields']] for tag in tags}
frt_cols = [(4, 'slip_strike'), (5, 'slip_dip'), (6, 'slip_normal'), (7, 'sliprate_strike'),
            (8, 'sliprate_dip'), (9, 'sliprate_peak'), (10, 'sliprate'), (11, 'traction_normal'),
            (12, 'traction_strike'), (13, 'traction_dip'),
            (14, 'vel_master_x'), (15, 'vel_master_y'), (16, 'vel_master_z'),
            (17, 'vel_slave_x'), (18, 'vel_slave_y'), (19, 'vel_slave_z'),
            (20, 'rsf_state'), (21, 'theta_pc')]
# frt.txt<r> holds rank r's nodes, fault by fault, in the same node order
acc = {}
for f in sorted(glob.glob('frt.txt*')):
    r = f[len('frt.txt'):]
    rows = np.loadtxt(f, ndmin=2)
    parts = [read_fault_ts('src_ts_hdr%s%s' % (t, r)) for t in tags if os.path.exists('src_ts_hdr%s%s' % (t, r))]
    xyz = np.vstack([p['xyz'] for p in parts])
    assert xyz.shape[0] == rows.shape[0] and np.abs(xyz - rows[:, :3]).max() <= 1e-7 * np.abs(xyz).max(), f
    for col, name in frt_cols:
        if all(name in p for p in parts):
            a = np.concatenate([p[name][-1] for p in parts])
            acc.setdefault(name, ([], []))
            acc[name][0].append(a); acc[name][1].append(rows[:, col])
for col, name in frt_cols:
    if name in acc:
        rep('frt last step %s' % name, np.concatenate(acc[name][0]), np.concatenate(acc[name][1]), 7)

# 3. last record == fault.dyna.r.nc (gridded from frt.txt*, so same rounding)
ds = nc.Dataset('fault.dyna.r.nc')
ncmap = [('shear_strike', 'traction_strike'), ('shear_dip', 'traction_dip'),
         ('effective_normal', 'traction_normal'), ('slip_rate', 'sliprate'),
         ('state_variable', 'rsf_state'), ('state_normal', 'theta_pc'),
         ('vxm', 'vel_master_x'), ('vym', 'vel_master_y'), ('vzm', 'vel_master_z'),
         ('vxs', 'vel_slave_x'), ('vys', 'vel_slave_y'), ('vzs', 'vel_slave_z')]
for vname, name in ncmap:
    a, b = [], []
    for k, t in enumerate(tags):
        if name not in TS[t]:
            continue
        xyz = TS[t]['xyz']
        fx, fz = np.unique(np.round(xyz[:, 0], 6)), np.unique(np.round(xyz[:, 2], 6))
        ii = np.searchsorted(fx, np.round(xyz[:, 0], 6)); jj = np.searchsorted(fz, np.round(xyz[:, 2], 6))
        v = ds[vname][:]
        v = np.asarray(v[k] if v.ndim == 3 else v)
        a.append(TS[t][name][-1]); b.append(v[jj, ii])
    if a:
        rep('nc last step %s' % vname, np.concatenate(a), np.concatenate(b), 7)

# 4. station histories faultst<tag><strike>dp<depth>.txt (E16.7 data, E21.13 time)
sign = int(open('bGlobal.txt').read().split('\n')[-3])
stat = {}
for f in sorted(glob.glob('faultst*.txt')):
    m = re.match(r'faultst(ft\d+_)?(-?\d+)dp(-?\d+)\.txt$', f)
    tag, xs, zs = m.group(1) or '', int(m.group(2)) * 100.0, -int(m.group(3)) * 100.0
    if tag not in TS:
        continue
    ts = TS[tag]
    k = np.argmin((ts['xyz'][:, 0] - xs) ** 2 + (ts['xyz'][:, 2] - zs) ** 2)
    lines = [l for l in open(f) if not l.lstrip().startswith('#') and l.strip() and not l.strip()[0].isalpha()]
    st = np.array([[float(x) for x in l.split()] for l in lines])
    n = min(len(st), len(ts['t']))
    cols = [(0, ts['t'], 1.0, 'time'), (1, ts['slip_strike'][:, k], 1.0, 'slip_strike'),
            (2, ts['sliprate_strike'][:, k], 1.0, 'sliprate_strike'),
            (3, ts['traction_strike'][:, k], 1e-6, 'traction_strike'),
            (4, ts['slip_dip'][:, k], -1.0, 'slip_dip'), (5, ts['sliprate_dip'][:, k], -1.0, 'sliprate_dip'),
            (6, ts['traction_dip'][:, k], -1e-6, 'traction_dip'),
            (7, ts['traction_normal'][:, k], sign * 1e-6, 'traction_normal')]
    if 'rsf_state' in ts:
        cols.append((8, ts['rsf_state'][:, k], 1.0, 'rsf_state'))
    if 'tp_temperature' in ts:
        cols.append((9, ts['tp_temperature'][:, k], 1.0, 'tp_temperature'))
    if ts['meta']['friclaw'] >= 3:
        # the station files add the background creep rate VINI_X along
        # strike (slot 26 = par.creep_slip_rate) to slip and slip rate
        sys.path.insert(0, '.')
        from user_defined_params import par
        vini = par.creep_slip_rate
        cols[1] = (1, ts['slip_strike'][:, k] + vini * ts['t'], 1.0, 'slip_strike (+VINI*t)')
        cols[2] = (2, ts['sliprate_strike'][:, k] + vini, 1.0, 'sliprate_strike (+VINI)')
    if 'tp_pore_change' in ts:
        # station column 11 = (pore-pressure change + initial pore pressure)/1e6;
        # the initial value is a constant, so take it from the first record
        p0 = st[0, 10] * 1e6 - ts['tp_pore_change'][0, k]
        cols.append((10, ts['tp_pore_change'][:, k] + p0, 1e-6, 'tp_pore (+ initial, MPa)'))
    for c, v, fac, name in cols:
        stat.setdefault(name, ([], []))
        stat[name][0].append(v[:n] * fac); stat[name][1].append(st[:n, c])
for name, (a, b) in stat.items():
    rep('station %s (%d files)' % (name, len(a)), np.concatenate(a), np.concatenate(b), 7 if name != 'time' else 13)

# 5. integrals over time
for t in tags:
    ts = TS[t]; dt = ts['dt']
    if ts['meta']['friclaw'] < 3:
        # friclaw 1/2: slip_path is the running sum of sliprate * dt
        rep('fault%s slip_path == cumsum(sliprate)*dt' % (t or '1'), ts['slip_path'], np.cumsum(ts['sliprate'] * dt, axis=0))
    for comp in ('strike', 'dip'):
        integ = np.cumsum(ts['sliprate_' + comp] * dt, axis=0)
        rep('fault%s slip_%s vs cumsum(sliprate_%s)*dt' % (t or '1', comp, comp), ts['slip_' + comp], integ)
    rep('fault%s traction_shear vs hypot' % (t or '1'), ts['traction_shear'], np.hypot(ts['traction_strike'], ts['traction_dip']))

print('%-55s %9s %12s %12s %s' % ('check', 'n', 'max|diff|', 'max|ref|', 'ok'))
for name, n, mx, sc, ok in res:
    print('%-55s %9d %12.3e %12.3e %s' % (name, n, mx, sc, ok))
for t in tags:
    m = TS[t]['meta']
    print('fault%s: %d nodes (all ranks), %d fields, %d steps; %d bytes per field per node per step'
          % (t or '1', TS[t]['xyz'].shape[0], m['nfield'], TS[t]['t'].size,
             sum(os.path.getsize(f) for f in glob.glob('src_ts%s[0-9]*' % t)) // (TS[t]['xyz'].shape[0] * m['nfield'] * TS[t]['t'].size)))
