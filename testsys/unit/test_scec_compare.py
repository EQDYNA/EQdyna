"""Both-ways tests (rule 14a) of the ONE SCEC overlay tool,
scripts/figures/scec_compare.py and its reader layer scec_readers.py.

Every on-disk layout the tool reads gets a synthetic fixture built here (no
archive, no shared dataset, no run directory needed), and every fixture
asserts a property that would come out WRONG if that layout's reader broke:
a masked-node count, a coordinate, a resolved file name, a station key. The
`test_mutant_*` tests then break a reader on purpose and assert the matching
fixture test's property no longer holds, so none of these can be green
while testing nothing.
"""
import math
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, '..', '..', 'scripts', 'figures')
sys.path.insert(0, os.path.abspath(FIG))

import scec_readers as R    # noqa: E402
import scec_compare         # noqa: E402

# a 9 x 4 node fault, 500 m spacing, hypocentre at (0, 750)
STRIKE = np.arange(-2000.0, 2001.0, 500.0)
DIP = np.arange(0.0, 1501.0, 500.0)
N_UNRUPT = 3                     # the last three nodes never rupture


def grid(sentinel):
    s, d = np.meshgrid(STRIKE, DIP, indexing='ij')
    s, d = s.ravel(), d.ravel()
    t = np.hypot(s, d - 750.0) / 3000.0
    t[-N_UNRUPT:] = sentinel
    return s, d, t


def write_cplot(path, sentinel, second_header=True):
    s, d, t = grid(sentinel)
    with open(path, 'w') as fh:
        fh.write('# problem = TPVxx\n# element_size = 100 m\n')
        fh.write('j k t\n')                    # the unparseable field line
        if second_header:
            fh.write('# a second comment block after the field line\n')
        for row in zip(s, d, t):
            fh.write('%14.6e %14.6e %14.6e\n' % row)


def write_frt(path, x, y, z, t):
    with open(path, 'w') as fh:
        for row in zip(x, y, z, t):
            fh.write('%14.6e %14.6e %14.6e %14.6e\n' % row)


def write_station(path, ncol=8, declared=None, lost_e=False, nt=20):
    tt = np.linspace(0.0, 1.9, nt)
    a = np.column_stack([tt] + [np.sin(tt + k) * (k + 1) for k in range(1, ncol)])
    with open(path, 'w') as fh:
        fh.write('# problem = TPVxx\n')
        fh.write('# Time series in %d columns of format e15.7\n'
                 % (declared if declared is not None else ncol))
        fh.write('t h-slip h-slip-rate h-shear-stress v-slip v-slip-rate '
                 'v-shear-stress n-stress\n')
        for i, row in enumerate(a):
            line = ' '.join('%15.7e' % v for v in row)
            if lost_e and i == 3:
                toks = line.split()
                toks[1] = '0.1341312-114'      # Fortran e15.7, 3-digit exponent
                line = ' '.join(toks)
            fh.write(line + '\n')
    return a


# ------------------------------------------------------------ sentinels

@pytest.mark.parametrize('sentinel', [1000.0, 99999.0, 1e9, 1e10])
def test_cplot_sentinel_masked(tmp_path, sentinel):
    write_cplot(tmp_path / 'cplot', sentinel)
    m = R.Model(str(tmp_path), 0)
    assert m.reader == 'cvws-cplot' and m.kind == 'archive'
    f = m.load_rupture()
    assert f['sentinel'] == sentinel
    assert int(np.isnan(f['t']).sum()) == N_UNRUPT
    assert f['t'].size == STRIKE.size * DIP.size
    assert f['dx_m'] == 500.0


def test_sentinel_not_invented_when_every_node_ruptured():
    s, d, t = grid(0.0)
    t[-N_UNRUPT:] = 1.0                # genuine late arrivals, repeated
    assert R.detect_sentinel(t)[0] is None


def test_mutant_greater_than_1e4_mask_misses_1000(tmp_path):
    """The bug the per-file detection exists for: a fixed `> 1e4` mask lets
    the 2015-2020 archives' 1000 s through as a real rupture time."""
    _, _, t = grid(1000.0)
    assert int((t > 1e4).sum()) == 0
    assert R.detect_sentinel(t)[0] == 1000.0


# ------------------------------------------------------------ our layouts

def _frt_vertical(y=0.0, sentinel=99999.0):
    s, d, t = grid(sentinel)
    return s, np.full_like(s, y), -d, t


def test_frt_canonical(tmp_path):
    x, y, z, t = _frt_vertical()
    write_frt(tmp_path / 'frt.canonical.txt', x, y, z, t)
    m = R.Model(str(tmp_path), 0)
    assert (m.reader, m.kind) == ('eqdyna-frt-canonical', 'eqdyna')
    f = m.load_rupture()
    assert np.allclose(np.sort(f['downdip']), np.sort(-z))   # down-dip = -z
    assert f['downdip'].min() >= 0.0
    assert int(np.isnan(f['t']).sum()) == N_UNRUPT
    assert 'vertical' in f['dipwise']


def test_frt_rank_files_concatenated(tmp_path):
    x, y, z, t = _frt_vertical()
    h = x.size // 2
    write_frt(tmp_path / 'frt.txt0', x[:h], y[:h], z[:h], t[:h])
    write_frt(tmp_path / 'frt.txt1', x[h:], y[h:], z[h:], t[h:])
    m = R.Model(str(tmp_path), 0)
    assert m.reader == 'eqdyna-frt-ranks'
    f = m.load_rupture()
    assert f['t'].size == x.size and f['source'] == '2 x frt.txt<rank>'


def test_probe_order_canonical_wins_over_ranks(tmp_path):
    x, y, z, t = _frt_vertical()
    write_frt(tmp_path / 'frt.canonical.txt', x, y, z, t)
    write_frt(tmp_path / 'frt.txt0', x[:3], y[:3], z[:3], t[:3])
    assert R.Model(str(tmp_path), 0).reader == 'eqdyna-frt-canonical'


def test_scec_rupture_time(tmp_path):
    write_cplot(tmp_path / 'SCECRuptureTime.txt', 99999.0)
    m = R.Model(str(tmp_path), 0)
    assert (m.reader, m.kind) == ('eqdyna-scec-rt', 'eqdyna')
    assert int(np.isnan(m.load_rupture()['t']).sum()) == N_UNRUPT


def test_dipping_frt_down_dip_is_hypot(tmp_path):
    s, d, t = grid(99999.0)
    dip = math.radians(60.0)
    write_frt(tmp_path / 'frt.canonical.txt', s, d * math.cos(dip),
              -d * math.sin(dip), t)
    f = R.Model(str(tmp_path), 0).load_rupture()
    assert 'dipping' in f['dipwise']
    assert np.allclose(np.unique(np.round(f['downdip'], 3)), DIP)
    assert abs(f['sin_dip'] - math.sin(dip)) < 1e-6


# ------------------------------------------------------------ station headers

def test_station_header_wrong_column_count_trusts_data(tmp_path):
    p = tmp_path / 'faultst000dp000'
    a = write_station(p, ncol=8, declared=7)
    got, note = R.load_station(str(p), 'fault')
    assert got.shape == a.shape
    assert note is not None and 'declares 7' in note and 'has 8' in note
    p2 = tmp_path / 'faultst000dp005'
    write_station(p2, ncol=8)
    assert R.load_station(str(p2), 'fault')[1] is None


def test_station_too_few_columns_refused(tmp_path):
    p = tmp_path / 'faultst000dp000'
    write_station(p, ncol=5)
    with pytest.raises(SystemExit):
        R.load_station(str(p), 'fault')


def test_lost_e_exponent_repaired(tmp_path):
    p = tmp_path / 'faultst000dp000.txt'
    a = write_station(p, lost_e=True)
    got, _ = R.load_station(str(p), 'fault')
    assert got.shape == a.shape
    assert got[3, 1] == pytest.approx(0.1341312e-114, rel=1e-12)


def test_mutant_lost_e_unrepaired_fails(tmp_path, monkeypatch):
    p = tmp_path / 'faultst000dp000.txt'
    write_station(p, lost_e=True)
    monkeypatch.setattr(R, '_floats', lambda s: [float(v) for v in s.split()])
    with pytest.raises(ValueError):
        R.load_station(str(p), 'fault')


def test_station_discovery_archive_and_ours(tmp_path):
    arc, ours = tmp_path / 'arc', tmp_path / 'ours'
    arc.mkdir(), ours.mkdir()
    for n in ('faultst-020dp000', 'faultst000dp015', 'body-030st000dp000'):
        write_station(arc / n)
    for n in ('faultst-020dp000.txt', 'faultst000dp015.txt',
              'faultstft2_010dp000.txt', 'body-030st000dp000.txt'):
        write_station(ours / n)
    fa, _ = R.discover_stations(str(arc), 'fault')
    fo, _ = R.discover_stations(str(ours), 'fault')
    assert set(fa) == set(fo) == {(-2000.0, 0.0), (0.0, 1500.0)}   # ft2 excluded
    bo, _ = R.discover_stations(str(ours), 'body')
    assert set(bo) == {(-3000.0, 0.0, 0.0)}


# ------------------------------------------------------------ multi-fault CVWS

def test_multifault_cvws_layouts(tmp_path):
    barall, kaneko, branch = (tmp_path / n for n in ('barall', 'kaneko', 'branch'))
    for d in (barall, kaneko, branch):
        d.mkdir()
    write_cplot(barall / 'cplot_1', 1e9)
    write_cplot(barall / 'cplot_2', 1e9)
    write_cplot(kaneko / 'cplot_1.txt', 1e10)
    write_cplot(kaneko / 'cplot_2.txt', 1e10)
    write_cplot(branch / 'cplot_main', 1e9)
    write_cplot(branch / 'cplot_branch', 1e9)
    for n in ('fault1st-100dp050', 'fault2st050dp050'):
        write_station(barall / n)
    for n in ('faultst-020dp000', 'branchst020dp000'):
        write_station(branch / n)
    want = {(barall, 1): 'cplot_1', (barall, 2): 'cplot_2',
            (kaneko, 1): 'cplot_1.txt', (kaneko, 2): 'cplot_2.txt',
            (branch, 1): 'cplot_main', (branch, 2): 'cplot_branch'}
    for (d, f), name in want.items():
        m = R.Model(str(d), 0, fault=f)
        assert os.path.basename(m.rupture_file) == name, (d, f)
        assert int(np.isnan(m.load_rupture()['t']).sum()) == N_UNRUPT
    assert set(R.Model(str(barall), 0, fault=1).station_files('fault')) == {(-10000.0, 5000.0)}
    assert set(R.Model(str(barall), 0, fault=2).station_files('fault')) == {(5000.0, 5000.0)}
    assert set(R.Model(str(branch), 0, fault=2).station_files('fault')) == {(2000.0, 0.0)}
    assert set(R.Model(str(branch), 0, fault=1).station_files('fault')) == {(-2000.0, 0.0)}


def test_multifault_without_fault_flag_hints(tmp_path):
    write_cplot(tmp_path / 'cplot_1', 1e9)
    with pytest.raises(SystemExit, match='--fault'):
        R.Model(str(tmp_path), 0)


def test_ours_ft2_station_tag(tmp_path):
    x, y, z, t = _frt_vertical()
    write_frt(tmp_path / 'frt.canonical.txt', x, y, z, t)
    write_station(tmp_path / 'faultst000dp000.txt')
    write_station(tmp_path / 'faultstft2_050dp050.txt')
    assert set(R.Model(str(tmp_path), 0, fault=2).station_files('fault')) == {(5000.0, 5000.0)}
    assert set(R.Model(str(tmp_path), 0, fault=1).station_files('fault')) == {(0.0, 0.0)}


# ------------------------------------------------------------ frt fault split

def test_frt_split_parallel_fault():
    x1, y1, z1, t1 = _frt_vertical(0.0)
    x2, y2, z2, t2 = _frt_vertical(5000.0)
    a = np.column_stack([np.r_[x1, x2], np.r_[y1, y2], np.r_[z1, z2], np.r_[t1, t2]])
    rows1, s1, _ = R.select_fault(a, 1)
    rows2, s2, rule = R.select_fault(a, 2)
    assert rows1.shape[0] == rows2.shape[0] == x1.size
    assert np.array_equal(s2, x2) and 'parallel' in rule


def test_frt_split_branch_strike_from_junction():
    x1, y1, z1, t1 = _frt_vertical(0.0)
    ang, xj = math.radians(30.0), 1000.0
    r = np.arange(500.0, 4001.0, 500.0)
    xb, yb = xj + r * math.cos(ang), r * math.sin(ang)
    a = np.column_stack([np.r_[x1, xb], np.r_[y1, yb],
                         np.r_[z1, np.zeros_like(r)], np.r_[t1, r / 3000.0]])
    rows2, s2, rule = R.select_fault(a, 2)
    assert rows2.shape[0] == r.size and 'BRANCH' in rule
    assert np.allclose(s2, r)


# ------------------------------------------------------------ TPV12 depth-named dp

def _dipping_case(tmp_path, names):
    s, d, t = grid(99999.0)
    dip = math.radians(60.0)
    write_frt(tmp_path / 'frt.canonical.txt', s, d * math.cos(dip),
              -d * math.sin(dip), t)
    for n in names:
        write_station(tmp_path / n)
    return R.Model(str(tmp_path), 0)


def test_depth_named_dp_converted(tmp_path):
    # down-dip 1500 m on a 60 deg fault is 1.3 km depth: written as dp013
    m = _dipping_case(tmp_path, ['faultst000dp013.txt', 'faultst-020dp000.txt'])
    assert set(m.station_files('fault')) == {(0.0, 1500.0), (-2000.0, 0.0)}
    assert any('DEPTH' in n for n in m.notes)


def test_correct_down_dip_names_untouched(tmp_path):
    m = _dipping_case(tmp_path, ['faultst000dp015.txt', 'faultst-020dp000.txt'])
    assert set(m.station_files('fault')) == {(0.0, 1500.0), (-2000.0, 0.0)}
    assert not any('DEPTH' in n for n in m.notes)


# ------------------------------------------------------------ mutants of the registry

@pytest.mark.parametrize('drop,make', [
    ('cvws-cplot', lambda d: write_cplot(d / 'cplot', 1000.0)),
    ('eqdyna-frt-canonical', lambda d: write_frt(d / 'frt.canonical.txt', *_frt_vertical())),
    ('eqdyna-scec-rt', lambda d: write_cplot(d / 'SCECRuptureTime.txt', 99999.0)),
    ('eqdyna-frt-ranks', lambda d: write_frt(d / 'frt.txt0', *_frt_vertical())),
])
def test_mutant_registry_entry_removed(tmp_path, monkeypatch, drop, make):
    """Each fixture is read by exactly its own registered reader: removing
    that entry makes the same directory unreadable."""
    make(tmp_path)
    assert R.Model(str(tmp_path), 0).reader == drop
    monkeypatch.setattr(R, 'RUPTURE_READERS',
                        [r for r in R.RUPTURE_READERS if r[0] != drop])
    with pytest.raises(SystemExit):
        R.Model(str(tmp_path), 0)


# ------------------------------------------------------------ every plot mode

def _two_models(tmp_path):
    ours, arc = tmp_path / 'ours', tmp_path / 'barall-faultmod-100m-2013'
    ours.mkdir(), arc.mkdir()
    write_frt(ours / 'frt.canonical.txt', *_frt_vertical())
    write_cplot(arc / 'cplot', 1e9)
    for n in ('faultst-020dp000', 'faultst000dp010', 'body-030st000dp000'):
        write_station(ours / (n + '.txt'))
        write_station(arc / n)
    return str(ours), str(arc)


@pytest.mark.parametrize('mode', ['cplot', 'res-series', 'ts-fault', 'ts-body'])
def test_plot_mode_smoke(tmp_path, mode, capsys):
    ours, arc = _two_models(tmp_path)
    out = tmp_path / ('%s.png' % mode)
    rc = scec_compare.main(['--plot', mode, '--models', ours, arc,
                            '--dpi', '30', '--note', '5 s term',
                            '--out', str(out)])
    assert rc == 0 and out.stat().st_size > 1000
    text = capsys.readouterr().out
    assert 'reader=eqdyna-frt-canonical' in text and 'reader=cvws-cplot' in text
    assert 'CVWS Barall (2013), FaultMod 100 m' in text
