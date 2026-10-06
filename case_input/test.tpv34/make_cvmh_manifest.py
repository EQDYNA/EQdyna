"""Write <dataset_root>/MANIFEST.json (size + md5 per file under raw/cvmh)
for the CVM-H data fetch_cvmh.sh downloaded, then make raw/ read-only.

Usage: python3 make_cvmh_manifest.py <dataset_root> [--vx-version 11.9.0]
Run only after fetch_cvmh.sh printed FETCH_DONE.
"""
import argparse, datetime, hashlib, json, os, stat

ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
ap.add_argument('root')
ap.add_argument('--vx-version', default='11.9.0',
                help='version string the built vx_lite reports')
a = ap.parse_args()
raw = os.path.join(a.root, 'raw', 'cvmh')
files, total = [], 0
for dp, dn, fn in os.walk(raw):
    dn.sort()
    for f in sorted(fn):
        p = os.path.join(dp, f)
        h = hashlib.md5()
        with open(p, 'rb') as fh:
            for chunk in iter(lambda: fh.read(1 << 24), b''):
                h.update(chunk)
        n = os.path.getsize(p)
        total += n
        files.append(dict(name=os.path.relpath(p, raw), size=n, md5=h.hexdigest()))
if not files:
    raise SystemExit(f'no files under {raw}')
man = dict(
    source_url='https://g-3a9041.a78b8.36fe.data.globus.org/ucvm/models/CVMH/cvmh/',
    source_description=('SCEC CVM-H 15.1.1 model data files as distributed by UCVM '
                        '(file list and host from SCECcode/cvmh main branch, model/config). '
                        'GOCAD voxets (CVM_HR/LR/CM), SCEC 1D background (CVMSM), Wills '
                        'Vs30 GTL, topo/base/moho surfaces and tsurf/ triangulated surfaces.'),
    code_url='https://github.com/SCECcode/cvmh',
    code_note=(f'vx_lite from SCECcode/cvmh main (reports Version {a.vx_version}), run as '
               '`vx_lite -s -z dep -m <raw/cvmh>`. SCEC TPV34 (TPV34_Description_v10, '
               'Part 2) prescribes CVM-H 15.1.0; this is 15.1.1, the version UCVM '
               'distributes -- label any TPV34 extraction with that difference.'),
    fetch_date=datetime.date.today().isoformat(),
    fetch_tool='EQdyna case_input/test.tpv34/fetch_cvmh.sh (curl, resumable, per file)',
    n_files=len(files), total_bytes=total, files=files)
with open(os.path.join(a.root, 'MANIFEST.json'), 'w') as f:
    json.dump(man, f, indent=1)
for dp, dn, fn in os.walk(os.path.join(a.root, 'raw')):
    for f in fn:
        os.chmod(os.path.join(dp, f), stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)
print('MANIFEST_DONE', len(files), total)
