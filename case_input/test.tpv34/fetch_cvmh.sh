#!/bin/bash
# Fetch the SCEC CVM-H 15.1.1 model data files (the UCVM distribution; file
# list and host from SCECcode/cvmh main, model/config) into <dataset_root>/raw/cvmh.
# Resumable: re-run to finish an interrupted download.
# Usage: bash fetch_cvmh.sh <dataset_root>      e.g. ~/shared_dataset/scec_cvmh.15.1.1
# Then:  python3 make_cvmh_manifest.py <dataset_root>
set -u
[ $# -eq 1 ] || { echo "usage: $0 <dataset_root>" >&2; exit 2; }
DEST="$1/raw/cvmh"
BASE=https://g-3a9041.a78b8.36fe.data.globus.org/ucvm/models/CVMH/cvmh
mkdir -p "$DEST/tsurf"
cd "$DEST" || exit 1
fail=0
for f in 'base@@' 'BASE.gts' 'BATO.gts' 'CVM_CM_TAG@@' 'CVM_CM.vo' 'CVM_CM_VP@@' 'CVM_CM_VS@@' \
         'CVM_HR_TAG@@' 'CVM_HR.vo' 'CVM_HR_VP@@' 'CVM_HR_VS@@' 'CVM_LR.vo' 'CVMSM_flags@@' \
         'CVMSM_tag66@@' 'CVMSM_vp66@@' 'CVMSM_vs66@@' 'cvm_vs30_wills.hdr' 'cvm_vs30_wills.mdl' \
         'interfaces.vo' 'model_top@@' 'moho@@' 'MOHO.gts' 'topo_dem@@' \
         'tsurf/CMxVM_Model3D_CalMex_BATO.ts' 'tsurf/CMxVM_Model3D_CM_BASE_Folded.dxf' \
         'tsurf/CMxVM_Model3D_CM_BASE_Folded.ts' 'tsurf/CVMH_Basement64.ts' \
         'tsurf/CVMH_CalMex_BATO.ts' 'tsurf/CVMH_Moho64.ts' 'tsurf/CVMH_Moho.ts'; do
  curl -sSL -C - --retry 20 --retry-delay 5 -o "$f" "$BASE/$f" || { echo "FAILED $f"; fail=1; }
done
echo "$(find . -type f | wc -l) files, $(du -sh . | cut -f1)"   # expect 30 files, ~1.5G
[ $fail -eq 0 ] && echo FETCH_DONE || { echo FETCH_INCOMPLETE; exit 1; }
