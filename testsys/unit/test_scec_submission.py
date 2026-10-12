"""Round-trip tests (rule 14a: each can go red) of scripts/scec_submission.py,
the SCEC CVWS submission-zip packager.

The committed test.tpv26 reference (frt.canonical.txt + stations/) is
packaged, the zip is extracted, and both are read back through the SAME
reader layer the overlay tool uses (scripts/figures/scec_readers.py). Every
station array and the rupture-time field must be EQUAL to the source --
not close: the packager copies data rows verbatim and writes cplot at the
source's own precision. The mutant tests prove the equality check can fail.
"""
import os
import re
import shutil
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
SHA = '0123456789abcdef0123456789abcdef01234567'
DATE = '2026/10/12'


def _package(src, tmp, name='sub.zip'):
    out = os.path.join(str(tmp), name)
    r = P.package(src, '26', out, commit=SHA, date=DATE)
    ext = os.path.join(str(tmp), name + '.d')
    with zipfile.ZipFile(out) as z:
        z.extractall(ext)
    return r, out, os.path.join(ext, 'tpv26', r['label'])


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
        P.package(d, '26', os.path.join(str(tmp_path), 'x.zip'), SHA, DATE)
