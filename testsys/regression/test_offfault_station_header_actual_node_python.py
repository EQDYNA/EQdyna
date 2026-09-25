#! /usr/bin/env python3
"""
Regression guard (row 94 audit finding 4) for `write_offfault_stations`
(`src/python/eqdyna/library_output.py`): the header LOCATION STAMP must
record the ACTUAL matched node (S['st_off_*_actual_m']), not the REQUESTED
coordinate (S['st_off_*_m']) -- mirrors the Fortran-side
test_offfault_station_header_actual_node.py, and the same gap: the existing
row-114 test's fixtures use a matched node equal to the request, so a
revert of the stamp source could not be caught there.

Cheap (rule 9): calls write_offfault_stations directly with a synthetic S
dict; no case, no mesh, no subprocess.
"""
import os
import re
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PYSRC = os.path.join(ROOT, 'src', 'python')
sys.path.insert(0, PYSRC)


def main():
    import numpy as np
    from eqdyna import library_output

    S = dict(
        st_off_idx=np.array([6], dtype=np.int64),
        st_off_x_m=np.array([500.0]), st_off_y_m=np.array([-2000.0]), st_off_z_m=np.array([-300.0]),
        st_off_x_actual_m=np.array([600.0]), st_off_y_actual_m=np.array([-2100.0]),
        st_off_z_actual_m=np.array([-500.0]),
        dx=100.0, dt=0.01,
    )
    off_st_hist = np.ones((1, 7, 1))

    fails = []
    with tempfile.TemporaryDirectory() as tmp:
        paths = library_output.write_offfault_stations(tmp, S, off_st_hist)
        if [os.path.basename(p) for p in paths] != ['body-020st005dp003.txt']:
            fails.append('filename %r, expected [\'body-020st005dp003.txt\'] (request-derived)'
                         % [os.path.basename(p) for p in paths])
        else:
            with open(paths[0]) as f:
                text = f.read()
            m = re.search(r'# location = (\S+) km off fault, (\S+) km along strike, (\S+) km depth', text)
            if not m:
                fails.append('no location stamp line; file was:\n%s' % text)
            else:
                got = [float(v) for v in m.groups()]
                want = [-2.1, 0.6, 0.5]
                if got != want:
                    fails.append('location stamp = %r, expected %r (the ACTUAL matched node) -- '
                                 'got the REQUESTED coordinate instead if this shows '
                                 '-20.0/5.0/3.0' % (got, want))

    if fails:
        print('FAIL test_offfault_station_header_actual_node_python')
        for f in fails:
            print('  -', f)
        return 1
    print('SUCCESS test_offfault_station_header_actual_node_python '
          '(filename request-derived, header stamp actual-node-derived, and they differ)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
