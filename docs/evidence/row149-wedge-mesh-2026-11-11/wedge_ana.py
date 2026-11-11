#!/usr/bin/env python3
"""Row 149 wedge-mesh scratch test: trace-station normal stress vs grid size, wedge mesh against insertion mesh.

usage: wedge_ana.py <insertion_root> <insertion_500_root> <wedge_root> <out_dir> [--allow-missing]

  insertion_root      holds t12_dx200, t12_dx150   (row149-force-budget runs, insertion mesh)
  insertion_500_root  holds t12_dx500_diag         (row149-trace-traction run, insertion mesh)
  wedge_root          holds w12_dx500, w12_dx500_p212, w12_dx200, w12_dx150 (this evidence)
Writes wedge_vs_insertion.csv, same_problem.csv, trace_profile.csv, rupture_time_field.csv and fig_dp000_vs_dx.png into out_dir.
A missing run directory is an error unless --allow-missing is given (then it is printed as NOT RUN).

Methods are those of the row149-force-budget-2026-11-10 evidence (budget_extra.py): the station value is column 8 of
faultst000dp000.txt at the sample nearest 4 s; the fixed-length mean is the nodal normal traction at x = 0, piecewise
linear in down-dip distance, averaged over the top L metres, plus the initial traction (linear down dip, C_elastic = 1).
"""
import sys, os, csv, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

T_EVAL = 4.0
DIP = 60.0
L_MEAN = 500.0
T0_SLOPE = -0.00739001      # MPa per metre down dip, TPV12 initial normal traction (T0 column of the force-budget at4s.csv)
FAULTMOD_100M = -3.0050     # TPV12 dp000 n-stress at 4 s, from the row149-trace-traction evidence
RATE_ONSET = 1.0e-3         # m/s, slip-rate threshold for the station rupture time
CHECK_STATIONS = ("faultst000dp015", "faultst000dp030", "faultst000dp075", "faultst000dp120", "faultst045dp075", "faultst120dp000")
NO_NODE = "no_node_at_this_dx"   # the solver writes a station file only where a fault node sits at the station
RUNS = [  # root index, dir, mesh path, nominal dx, partition
    (1, "t12_dx500_diag", "insertion", 500, "2,2,1"),
    (0, "t12_dx200", "insertion", 200, "2,2,1"),
    (0, "t12_dx150", "insertion", 150, "2,2,1"),
    (2, "w12_dx500", "wedge", 500, "2,2,1"),
    (2, "w12_dx500_p212", "wedge", 500, "2,1,2"),
    (2, "w12_dx200", "wedge", 200, "2,2,1"),
    (2, "w12_dx150", "wedge", 150, "2,2,1"),
]


def station(path):
    rows = []
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) < 8:
                continue
            try:
                rows.append([float(x) for x in p[:8]])
            except ValueError:
                continue
    if not rows:
        raise SystemExit("no station rows in " + path)
    return np.array(rows)   # t, strike slip, strike rate, strike shear, dip slip, dip rate, dip shear, n-stress


def spacing(rd):
    with open(os.path.join(rd, "bModelGeometry.txt")) as f:
        for line in f:
            p = line.split()
            if len(p) == 3:
                return [float(x) for x in p]
    raise SystemExit("no dx dy dz line in bModelGeometry.txt of " + rd)


def trace(rd):
    """x = 0 column of fault node pairs at the dump nearest 4 s: down-dip distance, traction change, its parts."""
    files = [f for f in sorted(glob.glob(os.path.join(rd, "diagTrace.*.txt"))) if os.path.getsize(f)]
    if not files:
        raise SystemExit("no diagTrace rows in " + rd)
    a = np.vstack([np.loadtxt(f, ndmin=2) for f in files])
    t = np.round(a[:, 1], 9)
    z = np.round(a[:, 2], 3) + 0.0
    tu = np.unique(t)
    tt = tu[np.argmin(np.abs(tu - T_EVAL))]
    arn, mS, mM, un = a[:, 3], a[:, 4], a[:, 5], a[:, 6:9]
    wM, wS = mS/(mS + mM), mM/(mS + mM)          # weight of the master-side and slave-side force in the traction
    per = 1.0/(arn*1e6)                          # nodal force (N) to traction (MPa)
    part = dict(ku_M=(a[:, 15:18]*un).sum(1)*per, hg_M=(a[:, 21:24]*un).sum(1)*per,
                ku_S=-(a[:, 12:15]*un).sum(1)*per, hg_S=-(a[:, 18:21]*un).sum(1)*per)
    c = wM*(part["ku_M"] + part["hg_M"]) + wS*(part["ku_S"] + part["hg_S"])
    zs = np.array(sorted(np.unique(z), reverse=True))
    idx = [np.flatnonzero(np.logical_and(z == zz, t == tt))[0] for zz in zs]
    top = dict((k, v[idx[0]]) for k, v in part.items())
    return tt, -zs/np.sin(np.radians(DIP)), c[idx], top, (wM[idx[0]], wS[idx[0]])


def frt(rd):
    """rupture-time field keyed by rounded (x, z); 99999 marks a node that never ruptured."""
    files = sorted(glob.glob(os.path.join(rd, "frt.txt*")))
    if not files:
        raise SystemExit("no frt files in " + rd)
    a = np.vstack([np.loadtxt(f, ndmin=2)[:, :4] for f in files])
    return dict(((round(r[0], 1), round(r[2], 1)), r[3]) for r in a)


def write_csv(path, head, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(head)
        w.writerows(rows)
    with open(path, newline="") as f:
        back = list(csv.reader(f))
    if len(back) != len(rows) + 1 or any(len(b) != len(head) for b in back):
        raise SystemExit("csv re-read mismatch: " + path)
    print("wrote %s: %d rows x %d cols (re-read ok)" % (path, len(rows), len(head)))


def main():
    args = [a for a in sys.argv[1:] if a != "--allow-missing"]
    allow = "--allow-missing" in sys.argv[1:]
    roots, out = args[0:3], args[3]
    head = ["run", "mesh_path", "partition", "dx_m", "dy_m", "dz_m", "along_dip_spacing_m", "dt_s", "station_rows", "t_station_s",
            "dp000_nstress_MPa", "trace_node_from_parts_MPa", "mean_total_top%dm_MPa" % L_MEAN,
            "stressdiv_hangingwall_MPa", "stressdiv_footwall_MPa", "hourglass_hangingwall_MPa", "hourglass_footwall_MPa",
            "weight_hangingwall", "weight_footwall"]
    head += ["dp000_nstress_min_3p5_to_4p5s_MPa", "dp000_nstress_max_3p5_to_4p5s_MPa", "dp000_dip_slip_4s_m"]
    rows, chk, pro = [], [], []
    for ri, name, path, dxn, part in RUNS:
        rd = os.path.join(roots[ri], name)
        if not os.path.isfile(os.path.join(rd, "DONE")) and ri == 2 or not os.path.isfile(os.path.join(rd, "faultst000dp000.txt")):
            if allow:
                print("NOT RUN / not finished:", name)
                continue
            raise SystemExit("missing or unfinished run " + rd)
        dx, dy, dz = spacing(rd)
        if abs(dx - dxn) > 1e-6:
            raise SystemExit("dx of %s is %g, expected %g" % (name, dx, dxn))
        st = station(os.path.join(rd, "faultst000dp000.txt"))
        i = int(np.argmin(np.abs(st[:, 0] - T_EVAL)))
        tt, sd, prof, parts, w = trace(rd)
        s = np.linspace(0.0, L_MEAN, 2001)
        mean = np.trapezoid(np.interp(s, sd, prof), s)/L_MEAN + T0_SLOPE*L_MEAN/2.0
        rows.append([name, path, part, "%.3f" % dx, "%.3f" % dy, "%.3f" % dz, "%.3f" % (sd[1] - sd[0]), "%.7f" % (st[1, 0] - st[0, 0]), len(st),
                     "%.6f" % st[i, 0], "%.4f" % st[i, 7], "%.4f" % prof[0], "%.4f" % mean,
                     "%.4f" % parts["ku_M"], "%.4f" % parts["ku_S"], "%.4f" % parts["hg_M"], "%.4f" % parts["hg_S"], "%.4f" % w[0], "%.4f" % w[1]])
        late = np.logical_and(st[:, 0] >= 3.5, st[:, 0] <= 4.5)
        rows[-1] += ["%.4f" % st[late, 7].min(), "%.4f" % st[late, 7].max(), "%.4f" % st[i, 4]]
        for k in range(len(sd)):
            pro.append([name, path, "%.0f" % dx, k, "%.3f" % sd[k], "%.4f" % prof[k], "%.4f" % (prof[k] + T0_SLOPE*sd[k])])
        for stn in CHECK_STATIONS:
            if not os.path.isfile(os.path.join(rd, stn + ".txt")):
                chk.append([name, path, part, "%.0f" % dx, stn] + [NO_NODE]*5)
                continue
            q = station(os.path.join(rd, stn + ".txt"))
            j = int(np.argmin(np.abs(q[:, 0] - T_EVAL)))
            rate = np.hypot(q[:, 2], q[:, 5])
            on = np.flatnonzero(rate > RATE_ONSET)
            chk.append([name, path, part, "%.0f" % dx, stn, "%.4f" % q[on[0], 0] if len(on) else "none", "%.4f" % q[j, 4], "%.4f" % q[j, 1],
                        "%.4f" % q[j, 7], "%.4f" % q[j, 6]])
    os.makedirs(out, exist_ok=True)
    write_csv(os.path.join(out, "wedge_vs_insertion.csv"), head, rows)
    chead = ["run", "mesh_path", "partition", "dx_m", "station", "rupture_time_s", "dip_slip_4s_m", "strike_slip_4s_m", "nstress_4s_MPa", "dip_shear_4s_MPa"]
    write_csv(os.path.join(out, "same_problem.csv"), chead, chk)
    write_csv(os.path.join(out, "trace_profile.csv"), ["run", "mesh_path", "dx_m", "node_row", "down_dip_m", "nstress_change_4s_MPa", "nstress_total_4s_MPa"], pro)
    d = [dict(zip(head, r)) for r in rows]
    fmt = "%-15s %-9s %-6s %7s %7s %8s %8s %9s %9s %9s : %8s %8s %8s %8s"
    print(fmt % ("run", "mesh", "part", "dx", "dy", "dz", "dip_sp", "dp000", "parts", "top500", "kuHW", "kuFW", "hgHW", "hgFW"))
    for r in d:
        print(fmt % (r["run"], r["mesh_path"], r["partition"], r["dx_m"], r["dy_m"], r["dz_m"], r["along_dip_spacing_m"],
              r["dp000_nstress_MPa"], r["trace_node_from_parts_MPa"], r[head[12]], r[head[13]], r[head[14]], r[head[15]], r[head[16]]))
    print("dp000 n-stress over 3.5-4.5 s (min, max) and dip slip at 4 s; then total n-stress at the top node rows, x = 0, 4 s:")
    for r in d:
        top = [q[6] for q in pro if q[0] == r["run"]][:6]
        print("  %-15s min %9s max %9s dslip %8s : %s" % (r["run"], r[head[19]], r[head[20]], r[head[21]], "  ".join(top)))
    print("same-problem check (rupture time s, dip slip at 4 s m, n-stress and dip shear at 4 s MPa):")
    for c in chk:
        print("  %-15s %-9s dx %4s %-16s trup %7s  dslip %8s  nstress %9s  dshear %9s" % (c[0], c[1], c[3], c[4], c[5], c[6], c[8], c[9]))
    print("rupture-time field, wedge against insertion at the same dx (fault nodes matched on x, z):")
    fr = []
    for wname, iname, iroot in (("w12_dx500", "t12_dx500_diag", 1), ("w12_dx500_p212", "t12_dx500_diag", 1), ("w12_dx200", "t12_dx200", 0), ("w12_dx150", "t12_dx150", 0)):
        if not any(r["run"] == wname for r in d):
            continue
        fw, fi = frt(os.path.join(roots[2], wname)), frt(os.path.join(roots[iroot], iname))
        keys = sorted(set(fw) & set(fi))
        tw, ti = np.array([fw[k] for k in keys]), np.array([fi[k] for k in keys])
        both = np.logical_and(tw < 1000.0, ti < 1000.0)
        dt = np.abs(tw - ti)[both]
        fr.append([wname, iname, len(fw), len(fi), len(keys), int((tw < 1000.0).sum()), int((ti < 1000.0).sum()), int(both.sum()),
                   "%.4f" % np.median(dt), "%.4f" % np.percentile(dt, 95), "%.4f" % dt.max()])
        print("  %-15s vs %-15s nodes %d / %d matched %d  ruptured %d / %d both %d  |dt| median %s p95 %s max %s s" % tuple(fr[-1]))
    write_csv(os.path.join(out, "rupture_time_field.csv"), ["wedge_run", "insertion_run", "wedge_nodes", "insertion_nodes", "matched_nodes", "wedge_ruptured",
              "insertion_ruptured", "both_ruptured", "abs_dt_median_s", "abs_dt_p95_s", "abs_dt_max_s"], fr)
    for path in ("insertion", "wedge"):
        v = dict((int(float(r["dx_m"])), float(r["dp000_nstress_MPa"])) for r in d if r["mesh_path"] == path and r["partition"] == "2,2,1")
        if 200 in v and 150 in v:
            print("VERDICT INPUT %-9s dp000 at 4 s: 200 m %.4f MPa, 150 m %.4f MPa, rise 200 -> 150 m = %+.4f MPa" % (path, v[200], v[150], v[150] - v[200]))
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    for path, col, mk in (("insertion", "tab:blue", "o"), ("wedge", "tab:red", "s")):
        pts = sorted((float(r["dx_m"]), float(r["dp000_nstress_MPa"])) for r in d if r["mesh_path"] == path and r["partition"] == "2,2,1")
        if pts:
            ax.plot([p[0] for p in pts], [p[1] for p in pts], "-", marker=mk, color=col, label="EQdyna, %s mesh" % path)
    ax.plot([100], [FAULTMOD_100M], "*", color="k", ms=12, label="FaultMod 100 m (reference marker)")
    ax.axhline(0, color="0.7", lw=0.5)
    ax.set_xlabel("dx = along-strike and along-dip node spacing (m)")
    ax.set_ylabel("normal stress at faultst000dp000, 4 s (MPa)")
    ax.set_title("TPV12 trace-station normal stress vs grid size:\ninsertion mesh and wedge-degeneration mesh", fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fp = os.path.join(out, "fig_dp000_vs_dx.png")
    fig.savefig(fp, dpi=120)
    plt.close(fig)
    print("wrote", fp)


if __name__ == "__main__":
    main()
