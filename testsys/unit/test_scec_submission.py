"""Round-trip tests (rule 14a: each can go red) of scripts/scec_submission.py,
the SCEC CVWS submission-zip packager.

The committed test.tpv26 reference (frt.canonical.txt + stations/) is
packaged, the zip is extracted, and both are read back through the SAME
reader layer the overlay tool uses (scripts/figures/scec_readers.py). Every
station array and the rupture-time field must be EQUAL to the source --
not close: the packager copies data rows verbatim and writes cplot at the
source's own precision. The mutant tests prove the equality check can fail.
"""
import glob
import os
import re
import shutil
import subprocess
import sys
import zipfile

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
sys.path.insert(0, os.path.join(ROOT, 'scripts', 'figures'))

import scec_readers as R        # noqa: E402
import scec_submission as P     # noqa: E402

REF = os.path.join(ROOT, 'test.reference.results', 'test.tpv26')
REF22 = os.path.join(ROOT, 'test.reference.results', 'test.tpv22')
SHA = '0123456789abcdef0123456789abcdef01234567'
DATE = '2026/10/12'


def _package(src, tmp, name='sub.zip', tpv='26'):
    # a committed reference set has no run log: the override is explicit
    out = os.path.join(str(tmp), name)
    r = P.package(src, tpv, out, commit=SHA, date=DATE, allow_no_log=True)
    ext = os.path.join(str(tmp), name + '.d')
    with zipfile.ZipFile(out) as z:
        z.extractall(ext)
    return r, out, os.path.join(ext, 'tpv' + tpv, r['label'])


def _sorted_field(model_dir):
    f = R.Model(model_dir, 0).load_rupture()
    o = np.lexsort((f['downdip'], f['strike']))
    return f['strike'][o], f['downdip'][o], f['t'][o]


def _assert_roundtrip(src_stations, src_field_dir, sub):
    found_src, _ = R.discover_stations(src_stations, 'fault')
    found_sub, _ = R.discover_stations(sub, 'fault')
    body_src, _ = R.discover_stations(src_stations, 'body')
    body_sub, _ = R.discover_stations(sub, 'body')
    assert set(found_src) == set(found_sub) and len(found_src) == 12
    assert set(body_src) == set(body_sub) and len(body_src) == 6
    nrows = 0
    for which, a, b in (('fault', found_src, found_sub), ('body', body_src, body_sub)):
        for k in a:
            x, _ = R.load_station(a[k], which)
            y, _ = R.load_station(b[k], which)
            assert x.shape == y.shape and x.shape[0] > 1
            np.testing.assert_array_equal(x, y, err_msg=f'{which} {k}')
            nrows += x.shape[0]
    s0, d0, t0 = _sorted_field(src_field_dir)
    s1, d1, t1 = _sorted_field(sub)
    np.testing.assert_array_equal(s0, s1)
    np.testing.assert_array_equal(d0, d1)
    np.testing.assert_array_equal(np.isnan(t0), np.isnan(t1))
    np.testing.assert_array_equal(t0[~np.isnan(t0)], t1[~np.isnan(t1)])
    return nrows, int(np.isnan(t0).sum()), int((~np.isnan(t0)).sum())


def test_roundtrip_reference_equals_source(tmp_path):
    r, out, sub = _package(REF, tmp_path)
    nrows, n_unrupt, n_rupt = _assert_roundtrip(
        os.path.join(REF, 'stations'), REF, sub)
    # a property of the compared content, so this cannot pass on nothing
    print(f'round-trip: 18 stations, {nrows} rows, cplot {n_rupt} ruptured + '
          f'{n_unrupt} never-ruptured nodes, all equal')
    assert nrows > 1000 and n_rupt > 100 and n_unrupt > 100
    with zipfile.ZipFile(out) as z:
        names = sorted(os.path.basename(n) for n in z.namelist())
    assert len(names) == 20 and 'cplot' in names and 'PROVENANCE.md' in names
    assert not any(n.endswith('.txt') for n in names)


def test_headers_carry_code_version_commit_dx_dt(tmp_path):
    r, _, sub = _package(REF, tmp_path)
    version = open(os.path.join(ROOT, 'VERSION')).read().strip()
    for f in sorted(os.listdir(sub)):
        txt = open(os.path.join(sub, f)).read()
        assert 'Sophon' not in txt and 'San-Ti' not in txt, f
        if f == 'PROVENANCE.md':
            assert SHA in txt and version in txt
            assert 'UNVERIFIED' in txt and '--no-run-log' in txt
            assert '/' + 'home' not in txt          # no paths in the zip
            continue
        hdr = [l for l in txt.splitlines() if l.startswith('#')]
        assert '# code=EQdyna' in hdr
        assert f'# code_version={version}' in hdr
        assert f'# commit={SHA}' in hdr
        assert '# problem=TPV26' in hdr
        assert f'# author={P.AUTHOR_PLACEHOLDER}' in hdr
        assert '# element_size=500 m' in hdr
        if f == 'cplot':
            assert txt.splitlines()[len(hdr)] == 'j k t'
            continue
        dt = float(next(l for l in hdr if l.startswith('# time_step=')).split('=')[1])
        a = R.read_numeric(os.path.join(sub, f))
        assert dt == pytest.approx(np.median(np.diff(a[:, 0])), rel=1e-9)
        nt = int(next(l for l in hdr if l.startswith('# num_time_steps=')).split('=')[1])
        assert nt == a.shape[0]
        assert any(l.startswith('# location=') for l in hdr)


def test_cplot_uses_spec_never_ruptured_value(tmp_path):
    _, _, sub = _package(REF, tmp_path)
    a = R.read_numeric(os.path.join(sub, 'cplot'))
    src = R.read_numeric(os.path.join(REF, 'frt.canonical.txt'))
    assert (a[:, 2] == 1.0e9).sum() == (src[:, 3] == 99999.0).sum() > 0
    assert a[:, 1].min() >= 0.0                      # down-dip positive down


def test_sentinel_maps_only_the_exact_frt_value(tmp_path):
    d = _copy_ref(tmp_path)
    f = os.path.join(d, 'frt.canonical.txt')
    a = R.read_numeric(f)
    k = np.flatnonzero(a[:, 3] == 99999.0)[:2]
    a[k[0], 3] = 99998.0                             # NOT the sentinel
    np.savetxt(f, a, fmt='%19.7E')
    _, _, sub = _package(d, tmp_path)
    c = R.read_numeric(os.path.join(sub, 'cplot'))
    assert (c[:, 2] == 99998.0).sum() == 1
    assert (c[:, 2] == 1.0e9).sum() == (a[:, 3] == 99999.0).sum()


def test_cplot_hypocentre_is_earliest_at_spec_location(tmp_path):
    # independent of the reader layer: TPV26_27 v13 Part 1 puts the
    # hypocentre 5 km along strike on the negative side (station
    # faultst-050dp100, "hypocenter") at 10 km depth, so j = -5000, k = +10000
    # must hold the earliest rupture time. Tests the j sign and k = -z.
    _, _, sub = _package(REF, tmp_path)
    c = R.read_numeric(os.path.join(sub, 'cplot'))
    tmin = c[:, 2].min()
    at = c[(c[:, 0] == -5000.0) & (c[:, 1] == 10000.0)]
    assert at.shape[0] == 1 and at[0, 2] == tmin
    early = c[c[:, 2] == tmin]
    print(f'hypocentre t = {tmin:g} s, {early.shape[0]} node(s) at that time')
    assert np.allclose(early[:, :2].mean(axis=0), [-5000.0, 10000.0])
    for bad in ((5000.0, 10000.0), (-5000.0, -10000.0)):   # j or k flipped
        m = (c[:, 0] == bad[0]) & (c[:, 1] == bad[1])
        assert not m.any() or c[m, 2].min() > tmin


def test_multifault_tpv22_names_and_cplots(tmp_path):
    r, out, sub = _package(REF22, tmp_path, tpv='22')
    names = sorted(os.listdir(sub))
    src = sorted(os.listdir(os.path.join(REF22, 'stations')))
    print(f'tpv22: {len(src)} source stations -> {names}')
    assert r['nft'] == 2
    fault = [n for n in names if n.startswith('fault')]
    assert fault and all(re.match(r'^fault[12]st-?\d+dp-?\d+$', n) for n in fault)
    assert not any(n.startswith('faultst') for n in names)
    n1 = sum(1 for f in src if re.match(r'^faultst-?\d', f))
    n2 = sum(1 for f in src if f.startswith('faultstft2_'))
    assert n1 > 0 and n2 > 0
    assert sum(n.startswith('fault1st') for n in names) == n1
    assert sum(n.startswith('fault2st') for n in names) == n2
    assert 'fault2st050dp050' in names and 'fault1st000dp000' in names
    assert 'cplot_1' in names and 'cplot_2' in names and 'cplot' not in names
    frt = R.read_numeric(os.path.join(REF22, 'frt.canonical.txt'))
    on1 = np.abs(frt[:, 1]) < 0.5
    for ft, nn in ((1, int(on1.sum())), (2, int((~on1).sum()))):
        txt = open(os.path.join(sub, f'cplot_{ft}')).read()
        assert f'# fault={ft}' in txt.splitlines()
        assert R.read_numeric(os.path.join(sub, f'cplot_{ft}')).shape[0] == nn > 0


def test_refuses_fault_count_mismatch(tmp_path):
    d = os.path.join(str(tmp_path), 'src')
    shutil.copytree(REF22, d)
    for f in glob.glob(os.path.join(d, 'stations', 'faultstft2_*')):
        os.remove(f)
    with pytest.raises(SystemExit, match='fault count'):
        P.package(d, '22', os.path.join(str(tmp_path), 'x.zip'), SHA, DATE,
                  allow_no_log=True)


def test_refuses_no_fault_stations(tmp_path):
    d = _copy_ref(tmp_path)
    for f in glob.glob(os.path.join(d, 'stations', 'faultst*')):
        os.remove(f)
    with pytest.raises(SystemExit, match='no on-fault station'):
        P.package(d, '26', os.path.join(str(tmp_path), 'x.zip'), SHA, DATE,
                  allow_no_log=True)


def test_refuses_missing_run_log_without_override(tmp_path):
    with pytest.raises(SystemExit, match='no-run-log'):
        P.package(REF, '26', os.path.join(str(tmp_path), 'x.zip'), SHA, DATE)


def test_refuses_no_stamp_without_commit(tmp_path):
    with pytest.raises(SystemExit, match='--commit'):
        P.package(REF, '26', os.path.join(str(tmp_path), 'x.zip'), None, DATE,
                  allow_no_log=True)


def _banner(d, stamp):
    version = open(os.path.join(ROOT, 'VERSION')).read().strip()
    open(os.path.join(d, 'run.log'), 'w').write(
        f' ==================   Welcome to EQdyna {version} =======\n'
        f' ==================   (src {stamp})   ========\n')


def test_commit_verified_against_banner_stamp(tmp_path):
    head = subprocess.check_output(['git', '-C', ROOT, 'rev-parse', 'HEAD'],
                                   text=True).strip()
    stamp = P.stamp_at_commit(head)
    assert re.fullmatch(r'[0-9a-f]{64}', stamp)
    d = _copy_ref(tmp_path)
    _banner(d, stamp[:12])
    r = P.package(d, '26', os.path.join(str(tmp_path), 'ok.zip'), head, DATE)
    assert r['commit'] == head
    _banner(d, '0' * 12)
    with pytest.raises(SystemExit, match='source stamp'):
        P.package(d, '26', os.path.join(str(tmp_path), 'x.zip'), head, DATE)


def _edit_station(d, fn):
    f = os.path.join(d, 'stations', 'faultst000dp100.txt')
    txt = open(f).read()
    open(f, 'w').write(fn(txt))


def test_refuses_missing_header_time_step(tmp_path):
    d = _copy_ref(tmp_path)
    _edit_station(d, lambda t: re.sub(r'.*time_step.*\n', '', t))
    with pytest.raises(SystemExit, match='no header time_step'):
        P.package(d, '26', os.path.join(str(tmp_path), 'x.zip'), SHA, DATE,
                  allow_no_log=True)


def test_refuses_single_row_station(tmp_path):
    d = _copy_ref(tmp_path)
    def one_row(t):
        ls = t.splitlines(True)
        data = [i for i, l in enumerate(ls) if re.match(r'\s*[-0-9.]', l)]
        return ''.join(l for i, l in enumerate(ls) if i not in data[1:])
    _edit_station(d, one_row)
    with pytest.raises(SystemExit, match='data row'):
        P.package(d, '26', os.path.join(str(tmp_path), 'x.zip'), SHA, DATE,
                  allow_no_log=True)


def test_dt_tolerance_is_relative():
    assert P._dt_rel_tol('0.0167') == pytest.approx(0.5e-4 / 0.0167, rel=1e-3)
    assert P._dt_rel_tol('8.3333e-03') == pytest.approx(0.5e-7 / 8.3333e-3, rel=1e-3)
    # 1/60 s passes as '0.0167'; 1/50 s does not
    assert abs(0.0167 - 1 / 60) <= P._dt_rel_tol('0.0167') * 0.0167
    assert abs(0.0167 - 1 / 50) > P._dt_rel_tol('0.0167') * 0.0167


def test_column_labels_use_spec_wording(tmp_path):
    _, _, sub = _package(REF, tmp_path)
    f = open(os.path.join(sub, 'faultst-050dp100')).read().splitlines()
    assert '# Column #5 = vertical slip (m)' in f
    assert '# Column #8 = normal stress (MPa)' in f
    assert not any('down-dip slip' in l for l in f)
    b = open(os.path.join(sub, 'body030st050dp000')).read().splitlines()
    assert '# Column #4 = vertical displacement (m)' in b
    assert sum(l.startswith('# Column #') for l in b) == 7


def _copy_ref(tmp_path):
    d = os.path.join(str(tmp_path), 'src')
    shutil.copytree(REF, d)
    return d


def test_mutant_changed_value_is_caught(tmp_path):
    d = _copy_ref(tmp_path)
    f = os.path.join(d, 'stations', 'faultst000dp100.txt')
    lines = open(f).read().splitlines(True)
    i = max(i for i, l in enumerate(lines) if not l.lstrip().startswith('#'))
    cols = lines[i].split()
    cols[1] = '%.7E' % (float(cols[1]) + 1.0)
    lines[i] = '  ' + '   '.join(cols) + '\n'
    open(f, 'w').writelines(lines)
    _, _, sub = _package(d, tmp_path)
    with pytest.raises(AssertionError):
        _assert_roundtrip(os.path.join(REF, 'stations'), REF, sub)


def test_mutant_changed_rupture_time_is_caught(tmp_path):
    d = _copy_ref(tmp_path)
    f = os.path.join(d, 'frt.canonical.txt')
    a = R.read_numeric(f)
    k = int(np.argmax(a[:, 3] < 99999.0))
    a[k, 3] += 0.5
    np.savetxt(f, a, fmt='%19.7E')
    _, _, sub = _package(d, tmp_path)
    with pytest.raises(AssertionError):
        _assert_roundtrip(os.path.join(REF, 'stations'), REF, sub)


def test_lost_e_exponent_is_repaired(tmp_path):
    d = _copy_ref(tmp_path)
    f = os.path.join(d, 'stations', 'faultst000dp100.txt')
    txt = open(f).read()
    lines = txt.splitlines(True)
    i = max(i for i, l in enumerate(lines) if not l.lstrip().startswith('#'))
    cols = lines[i].split()
    cols[2] = '0.1341312-114'
    lines[i] = '  ' + '   '.join(cols) + '\n'
    open(f, 'w').writelines(lines)
    _, _, sub = _package(d, tmp_path)
    out = open(os.path.join(sub, 'faultst000dp100')).read()
    assert '0.1341312E-114' in out and not re.search(r'\d-114', out)


def test_refuses_version_mismatch(tmp_path):
    d = _copy_ref(tmp_path)
    open(os.path.join(d, 'run.log'), 'w').write(
        ' ==================   Welcome to EQdyna 0.0.1 =======\n'
        ' ==================   (src abcdef123456)   ========\n')
    with pytest.raises(SystemExit, match='VERSION'):
        P.package(d, '26', os.path.join(str(tmp_path), 'x.zip'), SHA, DATE)


def test_refuses_dx_disagreeing_with_frt_spacing(tmp_path):
    d = _copy_ref(tmp_path)
    for f in os.listdir(os.path.join(d, 'stations')):
        p = os.path.join(d, 'stations', f)
        txt = re.sub(r'element_size =\s*500\.0+', 'element_size = 100.0', open(p).read())
        open(p, 'w').write(txt)
    with pytest.raises(SystemExit, match='spacing'):
        P.package(d, '26', os.path.join(str(tmp_path), 'x.zip'), SHA, DATE,
                  allow_no_log=True)


# ---- TPV103/104: emit exactly the upload spec's 9 on-fault columns ----------
REF104 = os.path.join(ROOT, 'test.reference.results', 'test.tpv104')
SPEC104 = 't h-slip h-slip-rate h-shear-stress v-slip v-slip-rate v-shear-stress n-stress psi'


def _data_rows(path):
    _, _, data = P.split_station(path)
    return data


def test_tpv104_on_fault_files_carry_exactly_the_spec_columns(tmp_path):
    r, out, sub = _package(REF104, tmp_path, tpv='104')
    src = sorted(glob.glob(os.path.join(REF104, 'stations', 'faultst*.txt')))
    assert len(src) == 3
    nrows = 0
    for f in src:
        name = os.path.basename(f)[:-4]
        txt = open(os.path.join(sub, name)).read().splitlines()
        assert SPEC104 in txt
        assert '# Time series in 9 columns' in txt
        assert '# Column #9 = state variable psi (dimensionless)' in txt
        assert sum(l.startswith('# Column #') for l in txt) == 9
        assert not any('emperature' in l or 'ressure' in l for l in txt)
        got = [l for l in txt if not l.startswith('#') and l != SPEC104]
        want = _data_rows(f)
        assert len(got) == len(want) > 1
        for g, w in zip(got, want):
            cells = re.findall(r'\s*\S+', w)
            assert len(cells) == 11
            # the first 9 source cells, byte for byte; psi is source column 9
            assert g == ''.join(cells[:9])
        nrows += len(got)
    print(f'tpv104: 3 on-fault files, {nrows} rows, 9 spec columns each')
    with zipfile.ZipFile(out) as z:
        prov = [n for n in z.namelist() if n.endswith('PROVENANCE.md')][0]
        p = z.read(prov).decode()
    assert 'columns dropped' in p and 'pressure, temperature' in p


def test_tpv104_body_files_are_unchanged(tmp_path):
    _, _, sub = _package(REF104, tmp_path, tpv='104')
    for f in glob.glob(os.path.join(REF104, 'stations', 'body*.txt')):
        got = open(os.path.join(sub, os.path.basename(f)[:-4])).read().splitlines()
        want = _data_rows(f)
        assert got[-len(want):] == want


def test_extra_columns_still_refused_outside_tpv103_104(tmp_path):
    with pytest.raises(SystemExit, match='no spec wording'):
        P.package(REF104, '26', os.path.join(str(tmp_path), 'x.zip'), SHA,
                  DATE, allow_no_log=True)


def test_tpv104_refuses_a_missing_spec_column():
    with pytest.raises(SystemExit, match=r"missing from the run: \['psi'\]"):
        P.select_spec_columns('104', 'fault', SPEC104.replace(' psi', ' temperature'),
                              ['  1.0' * 9])


def test_no_columns_dropped_for_other_tpvs(tmp_path):
    _, out, _ = _package(REF, tmp_path)
    with zipfile.ZipFile(out) as z:
        prov = [n for n in z.namelist() if n.endswith('PROVENANCE.md')][0]
        assert 'columns dropped' not in z.read(prov).decode()


# ---- TPV105-3D: the spec's 11 on-fault columns, thermal pressurization kept --
REF1053D = os.path.join(ROOT, 'test.reference.results', 'test.tpv1053d')
SPEC1053D = SPEC104 + ' temperature pressure'


def test_tpv1053d_on_fault_files_keep_all_11_spec_columns(tmp_path):
    r, out, sub = _package(REF1053D, tmp_path, tpv='105-3D')
    src = sorted(glob.glob(os.path.join(REF1053D, 'stations', 'faultst*.txt')))
    assert len(src) == 3
    nrows = 0
    for f in src:
        txt = open(os.path.join(sub, os.path.basename(f)[:-4])).read().splitlines()
        assert '# problem=TPV105-3D' in txt
        assert SPEC1053D in txt
        assert sum(l.startswith('# Column #') for l in txt) == 11
        assert '# Column #9 = state variable psi (dimensionless)' in txt
        assert '# Column #10 = temperature (K)' in txt
        assert '# Column #11 = pore pressure (MPa)' in txt
        got = [l for l in txt if not l.startswith('#') and l != SPEC1053D]
        want = _data_rows(f)
        assert len(got) == len(want) > 1
        assert got == want          # all 11 columns, byte for byte
        nrows += len(got)
    print(f'tpv105-3D: 3 on-fault files, {nrows} rows, 11 spec columns each')
    with zipfile.ZipFile(out) as z:
        prov = [n for n in z.namelist() if n.endswith('PROVENANCE.md')][0]
        assert 'columns dropped' not in z.read(prov).decode()


def test_tpv1053d_key_is_case_insensitive():
    a = P.select_spec_columns('105-3d', 'fault', SPEC1053D, ['  1.0' * 11])
    b = P.select_spec_columns('105-3D', 'fault', SPEC1053D, ['  1.0' * 11])
    assert a == b and a[0] == SPEC1053D and a[2] == []


def test_tp_columns_refused_for_a_tpv_whose_spec_lacks_them(tmp_path):
    # the TP labels exist for TPV105-3D only; a TPV26 package of the same
    # 11-column run must still refuse them
    with pytest.raises(SystemExit, match=r"no spec wording: \['temperature', 'pressure'\]"):
        P.package(REF1053D, '26', os.path.join(str(tmp_path), 'x.zip'), SHA,
                  DATE, allow_no_log=True)
