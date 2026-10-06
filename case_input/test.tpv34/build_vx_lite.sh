#!/bin/bash
# Build vx_lite (the CVM-H query tool extract_cvmh_grid.py calls) from
# SCECcode/cvmh. Needs gcc, make, autoconf, automake, libtool.
# Usage: bash build_vx_lite.sh <work_dir>
#   clones https://github.com/SCECcode/cvmh into <work_dir>/cvmh (if absent)
#   and leaves <work_dir>/cvmh/src/vx_lite; pass that to extract_cvmh_grid.py
#   --vx-lite, with LD_LIBRARY_PATH reaching its libvxapi.so if needed.
set -u
[ $# -eq 1 ] || { echo "usage: $0 <work_dir>" >&2; exit 2; }
mkdir -p "$1" && cd "$1" || exit 1
[ -d cvmh ] || git clone --depth 1 https://github.com/SCECcode/cvmh.git cvmh || exit 1
cd cvmh || exit 1
libtoolize >/dev/null 2>&1
aclocal -I m4 && autoconf && automake --add-missing -f || exit 1
./configure --prefix="$PWD/install" || exit 1
make -C gctpc/source || exit 1
make -C src vx_lite || exit 1
ls -l src/vx_lite
