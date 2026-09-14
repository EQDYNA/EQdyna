#! /usr/bin/env python3
"""
Exact per-rank element count for EQdyna's internal structured mesher.

Line-by-line replica of the only two routines that decide how many elements a
rank owns:

  src/meshgen.f90:getLocalOneDimCoorArrAndSize  -> global node-line size per
      dimension, then the 1-D block split across npx/npy/npz with a one-node
      halo overlap between neighbours
  src/countMeshEntities.f90                     -> elements counted where
      ix>=2 .and. iy>=2 .and. iz>=2, i.e. (nxl-1)*(nyl-1)*(nzl-1)

Why this file exists -- what the code prints, and what is wrong with each:

  * `library.f90:memory_estimate` at HEAD (dabacbf) prints only
        1.54 * totalNumOfElements/1e6 * npx*npy*npz  GB
    `totalNumOfElements` is rank 0's LOCAL count, so the printed value is a
    GLOBAL memory estimate (rank-0 count x rank count). It is near-invariant
    with rank count at fixed dx BY CONSTRUCTION -- that is the line working as
    designed, not a bug -- and it cannot be inverted to a per-rank count
    without separately knowing rank 0's share.
  * The working tree's rewritten `memory_estimate` does print
    "Cells per rank (rank 0)", which is a true per-rank count -- but only
    RANK 0's. Rank 0 is (mex,mey,mez)=(0,0,0), which takes the SMALLEST slice
    whenever a dimension does not divide evenly (meshgen.f90:476-480 gives the
    extra node line to the high-id ranks). Per-step cost is gated by the
    SLOWEST rank, i.e. the maximum, and "Cells total" computed as rank0 x
    nranks understates the global count by the imbalance. Measured here:
    TPV29 dx=500, 4x1x4 -> rank 0 has 53,176 but the max rank has 57,960
    (+9.0%) and the true global is 876,024, not 53,176x16 = 850,816 (-2.9%).
  * `scripts/case.setup:estimate_HPC_resource` is wrong twice over: it counts
    only the +/-nuni_y cells in y (dropping the geometric coarsening zone and
    the PML) and then multiplies a global cell count by the rank count.

The authoritative runtime number is `totalNumOfElements` written per rank into
`compTime<me>` by `library_output.f90:output_timeanalysis` (needs
globalvar.f90's `writeCompTime = 1`). This module reproduces it exactly and
gives max/min/mean without needing a run.

Validation: reproduces, with no run, the numbers a real 1024-rank set-up was
sized against (scratch/tpv29/hpc50m/README.md) -- node lines 985 x 247 x 493,
119,095,488 global elements, 119,164 max elem/rank, +2.5% imbalance -- and
matches the per-rank `totalNumOfElements` written by every rank at 1, 4, 16 and
48 ranks in testsys/perf/elem_scaling_last.json (`pred_per_rank_match`).

Usage:
    python3 elem_per_rank.py --case-dir <dir with bGlobal.txt/bModelGeometry.txt>
    python3 elem_per_rank.py --case-dir <dir> --np 8 1 6
or import: counts(params, npx, npy, npz) -> dict
"""
import argparse
import os
import sys

NPML = 6            # src/globalvar.f90:164
NP_LOOP = 1000000   # src/globalvar.f90:232 (loop bound only)


def _line_size(grid_size, front_edge, back_edge, min_coor, max_coor,
               n_uniform, rat, dim_id):
    """globalOneDimCoorArrSize for one dimension (meshgen.f90:453-470)."""
    coor, gs = front_edge, grid_size
    for i in range(1, NP_LOOP + 1):
        gs *= rat
        coor -= gs
        if coor <= min_coor:
            break
    front_edge_node_id = i + NPML

    coor, gs = back_edge, grid_size
    for j in range(1, NP_LOOP + 1):
        gs *= rat
        coor += gs
        if coor >= max_coor:
            break
    if dim_id == 3:
        j = -NPML
    return n_uniform + front_edge_node_id + j + NPML


def _local_size(n_global, n_mpi, mpi_id):
    """Local node-line size on one rank (meshgen.f90:473-480)."""
    per = (n_global + n_mpi - 1) // n_mpi
    residual = (n_global + n_mpi - 1) - per * n_mpi
    return per if mpi_id < (n_mpi - residual) else per + 1


def global_line_sizes(p):
    """(Nx, Ny, Nz) global node-line sizes. `p` is a dict of case parameters."""
    nx_uni = round((p['fxmax'] - p['fxmin']) / p['dx']) + 1
    ny_uni = p['dis4uniF'] + p['dis4uniB'] + 1
    nz_uni = round((p['fzmax'] - p['fzmin']) / p['dz']) + 1
    rat = p['rat']
    return (
        _line_size(p['dx'], p['fxmin'], p['fxmax'], p['xmin'], p['xmax'],
                   nx_uni, rat, 1),
        _line_size(p['dy'], -p['dis4uniF'] * p['dy'], p['dis4uniB'] * p['dy'],
                   p['ymin'], p['ymax'], ny_uni, rat, 2),
        _line_size(p['dz'], p['fzmin'], p['fzmax'], p['zmin'], p['zmax'],
                   nz_uni, rat, 3),
    )


def counts(p, npx, npy, npz):
    """Per-rank element counts for a decomposition. Returns a summary dict."""
    ngx, ngy, ngz = global_line_sizes(p)
    per_rank = []
    # meshgen.f90:calcXyzMPIId -- me = mex*npy*npz + mey*npz + mez
    for mex in range(npx):
        for mey in range(npy):
            for mez in range(npz):
                nxl = _local_size(ngx, npx, mex)
                nyl = _local_size(ngy, npy, mey)
                nzl = _local_size(ngz, npz, mez)
                per_rank.append((nxl - 1) * (nyl - 1) * (nzl - 1))
    n = npx * npy * npz
    return dict(
        node_lines=(ngx, ngy, ngz),
        global_elements=(ngx - 1) * (ngy - 1) * (ngz - 1),
        ranks=n,
        per_rank=per_rank,
        elem_max=max(per_rank),
        elem_min=min(per_rank),
        elem_mean=sum(per_rank) / n,
        imbalance=max(per_rank) / (sum(per_rank) / n) - 1.0,
        halo_inflation=sum(per_rank) / ((ngx - 1) * (ngy - 1) * (ngz - 1)),
    )


def read_case(case_dir):
    """Parse bGlobal.txt + bModelGeometry.txt + bFaultGeometry.txt of a set-up case."""
    def toks(name):
        with open(os.path.join(case_dir, name)) as f:
            return [ln.split() for ln in f]
    g = toks('bGlobal.txt')
    m = toks('bModelGeometry.txt')
    ft = toks('bFaultGeometry.txt')
    p = {}
    p['npx'], p['npy'], p['npz'] = (int(v) for v in g[13])
    p['term'], p['dt'] = float(g[15][0]), float(g[16][0])
    p['xmin'], p['xmax'] = (float(v) for v in m[0])
    p['ymin'], p['ymax'] = (float(v) for v in m[1])
    p['zmin'], p['zmax'] = (float(v) for v in m[2])
    p['dis4uniF'], p['dis4uniB'] = (int(v) for v in m[4])
    p['rat'] = float(m[5][0])
    p['dx'], p['dy'], p['dz'] = (float(v) for v in m[6])
    p['fxmin'], p['fxmax'] = (float(v) for v in ft[1])
    p['fzmin'], p['fzmax'] = (float(v) for v in ft[3])
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--case-dir', default='.')
    ap.add_argument('--np', type=int, nargs=3, metavar=('NPX', 'NPY', 'NPZ'),
                    help='override the decomposition in bGlobal.txt')
    a = ap.parse_args()
    p = read_case(a.case_dir)
    npx, npy, npz = a.np if a.np else (p['npx'], p['npy'], p['npz'])
    c = counts(p, npx, npy, npz)
    nstep = round(p['term'] / p['dt'])
    print(f"case         {os.path.abspath(a.case_dir)}")
    print(f"dx,dy,dz     {p['dx']:.0f} {p['dy']:.0f} {p['dz']:.0f}")
    print(f"decomp       {npx} x {npy} x {npz} = {c['ranks']} ranks")
    print(f"node lines   {c['node_lines']}")
    print(f"global elem  {c['global_elements']:,}")
    print(f"elem/rank    max {c['elem_max']:,}  mean {c['elem_mean']:,.0f}  "
          f"min {c['elem_min']:,}  imbalance {c['imbalance']:+.1%}")
    print(f"halo infl.   {c['halo_inflation']:.4f}x (sum of local counts / global)")
    print(f"nstep        {nstep}  (term {p['term']} / dt {p['dt']:.6g})")


if __name__ == '__main__':
    sys.exit(main())
