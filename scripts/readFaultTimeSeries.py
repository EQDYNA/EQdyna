#! /usr/bin/env python3
"""
Reads the whole-fault time series EQdyna writes when
par.outputFaultTimeSeries = 1: per MPI rank and fault, a text header
src_ts_hdr<tag><rank> and a binary src_ts<tag><rank> (float64, native byte
order, no record markers, layout [step][node][field]; the header's
byteorder line names the order). See
docs/user/outputs.md for the field list.

    from readFaultTimeSeries import read_fault_ts, read_case
    ts = read_fault_ts('src_ts_hdr0')    # one rank / fault
    ts['slip_strike']                    # (nstep, nnode) array
    ts['xyz'], ts['t'], ts['dt']         # (nnode, 3) coords, (nstep,) times
    allts = read_case('.', tag='')       # every rank of fault 1, concatenated

Command line: python3 readFaultTimeSeries.py [case_dir] prints a summary of
every src_ts header in the directory.
"""
import glob
import os
import re
import sys
import warnings

import numpy as np


def read_header(hdr_path):
    """Parse a src_ts_hdr file: returns (meta dict, field list, xyz array).
    Each field is a dict with index, name, unit, desc."""
    meta, fields = {}, []
    with open(hdr_path) as f:
        lines = f.readlines()
    for line in lines:
        if not line.startswith('#'):
            continue
        body = line[1:].strip()
        if body.startswith('field '):
            parts = body.split(None, 4)
            fields.append(dict(index=int(parts[1]), name=parts[2], unit=parts[3],
                               desc=parts[4] if len(parts) > 4 else ''))
        elif ' = ' in body:
            key, val = body.split(' = ', 1)
            meta[key.strip()] = val.strip()
    for key in ('fault', 'rank', 'friclaw', 'nnode', 'nfield', 'nstep'):
        if key in meta:
            meta[key] = int(meta[key])
    meta['dt'] = float(meta['dt'])
    if len(fields) != meta['nfield'] or [fd['index'] for fd in fields] != list(range(meta['nfield'])):
        raise ValueError('%s: header lists %d field lines, nfield = %d'
                         % (hdr_path, len(fields), meta['nfield']))
    xyz = np.loadtxt(hdr_path, comments='#', ndmin=2)
    if xyz.shape != (meta['nnode'], 3):
        raise ValueError('%s: %r coordinate rows, expected (%d, 3)'
                         % (hdr_path, xyz.shape, meta['nnode']))
    return meta, fields, xyz


def read_fault_ts(hdr_path, fields=None, mmap=False):
    """Read one src_ts header + its binary (named in the header, looked up
    next to it). Returns a dict: one (nstep_written, nnode) array per field
    name, plus 'xyz', 't', 'dt', 'meta', 'fields'. `fields` limits which
    field arrays are returned; mmap=True maps the binary instead of reading
    it (the arrays are then strided views into it). A trailing partial
    record (a run killed mid-write) is dropped with a warning; the complete
    records before it are returned."""
    meta, flist, xyz = read_header(hdr_path)
    bin_path = os.path.join(os.path.dirname(hdr_path), meta['binary'])
    nnode, nfield = meta['nnode'], meta['nfield']
    dtype = {'little': '<f8', 'big': '>f8'}[meta.get('byteorder', 'little')]
    rec = nnode * nfield
    nbytes = os.path.getsize(bin_path)
    nrec = nbytes // (8 * rec)
    if nbytes % (8 * rec):
        warnings.warn('%s: %d bytes is not a whole number of %d-value records; '
                      'reading the %d complete records and dropping the partial last one'
                      % (bin_path, nbytes, rec, nrec))
    if mmap:
        raw = np.memmap(bin_path, dtype=dtype, mode='r', shape=(nrec * rec,))
    else:
        raw = np.fromfile(bin_path, dtype=dtype, count=nrec * rec)
    data = raw.reshape(nrec, nnode, nfield)
    out = dict(xyz=xyz, dt=meta['dt'], t=(np.arange(nrec) + 1) * meta['dt'],
               meta=meta, fields=flist)
    want = None if fields is None else set(fields)
    for fd in flist:
        if want is None or fd['name'] in want:
            out[fd['name']] = data[:, :, fd['index']]
    return out


def read_case(case_dir='.', tag='', fields=None):
    """Read every rank's file for one fault (tag '' for fault 1, 'ft2_' for
    fault 2, ...) and concatenate along the node axis, in rank order."""
    pat = re.compile(r'^src_ts_hdr' + re.escape(tag) + r'(\d+)$')
    hdrs = sorted((int(m.group(1)), os.path.join(case_dir, n))
                  for n in os.listdir(case_dir) for m in [pat.match(n)] if m)
    if not hdrs:
        raise FileNotFoundError('no src_ts_hdr%s<rank> in %s' % (tag, case_dir))
    parts = [read_fault_ts(h, fields=fields) for _, h in hdrs]
    nrec = {p['t'].size for p in parts}
    if len(nrec) != 1:
        raise ValueError('ranks wrote different step counts: %r' % sorted(nrec))
    out = dict(parts[0])
    out['xyz'] = np.concatenate([p['xyz'] for p in parts])
    out['ranks'] = np.concatenate([np.full(p['meta']['nnode'], r) for (r, _), p in zip(hdrs, parts)])
    for fd in parts[0]['fields']:
        if fd['name'] in parts[0]:
            out[fd['name']] = np.concatenate([p[fd['name']] for p in parts], axis=1)
    return out


def main(argv):
    case_dir = argv[1] if len(argv) > 1 else '.'
    hdrs = sorted(glob.glob(os.path.join(case_dir, 'src_ts_hdr*')))
    if not hdrs:
        print('no src_ts_hdr* files in %s' % case_dir)
        return 1
    for h in hdrs:
        ts = read_fault_ts(h, mmap=True)
        m = ts['meta']
        print('%s: fault %d rank %d, %d nodes x %d fields, %d of %d steps, dt = %g s'
              % (os.path.basename(h), m['fault'], m['rank'], m['nnode'], m['nfield'],
                 ts['t'].size, m['nstep'], m['dt']))
        for fd in ts['fields']:
            a = ts[fd['name']]
            print('  %2d %-16s %-5s min %+.4e max %+.4e'
                  % (fd['index'], fd['name'], fd['unit'], a.min(), a.max()))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
