#!/usr/bin/env python3
"""Row 149 force budget, second pass: fixed-length means of the trace normal traction,
the viscous-hourglass spike against the baseline, and a text summary of the 4 s budget.
#
usage: budget_extra.py <runs_root> <csv_dir> <fig_dir> [<morning_root>]
#
Fixed-length mean (method of the row149-trace-traction evidence, ana.py): nodal normal
traction at x = 0, piecewise linear in down-dip distance, averaged over the top L metres.
  change = stress-divergence + hourglass parts (what the diagTrace dump holds)
  total  = change + initial normal traction. C_elastic = 1: the initial traction is linear in
           down-dip distance, slope read from the T0 column of at4s.csv. C_elastic = 0: the
           initial stress is already inside the element stress, so total = change.
"""
import sys, os, csv, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
#
T_EVAL = 4.0
LENGTHS = (250.0, 500.0, 1000.0)
RUNS = [  # dir, case key, label, dx, dip (deg), C_elastic, source
    ("t12_dx200", "t12", "TPV12", 200, 60.0, 1, "budget build"),
    ("t12_dx150", "t12", "TPV12", 150, 60.0, 1, "budget build"),
    ("t12hg2_dx200", "t12hg2", "TPV12 viscous-hourglass spike", 200, 60.0, 1, "budget build + C_hg = 2"),
    ("t12hg2_dx150", "t12hg2", "TPV12 viscous-hourglass spike", 150, 60.0, 1, "budget build + C_hg = 2"),
    ("t13_dx200", "t13", "TPV13", 200, 60.0, 0, "budget build"),
    ("t13_dx150", "t13", "TPV13", 150, 60.0, 0, "budget build"),
]
MORNING = [
    ("t12_dx500_diag", "t12", "TPV12", 500, 60.0, 1, "morning build"),
    ("t12_dx200_diag", "t12", "TPV12", 200, 60.0, 1, "morning build"),
    ("t12_dx150_diag", "t12", "TPV12", 150, 60.0, 1, "morning build"),
    ("t13_dx500_diag", "t13", "TPV13", 500, 60.0, 0, "morning build"),
    ("t13_dx250_diag", "t13", "TPV13", 250, 60.0, 0, "morning build"),
    ("t13_dx200_diag", "t13", "TPV13", 200, 60.0, 0, "morning build"),
    ("t13_dx150_diag", "t13", "TPV13", 150, 60.0, 0, "morning build"),
]
FAULTMOD_100M = {"t12": -3.0050, "t13": 1.5193}   # dp000 n-stress at 4 s, from the row149-trace-traction evidence
#
#
def station(path):
    t, s = [], []
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) < 8:
                continue
            try:
                v = [float(x) for x in p[:8]]
            except ValueError:
                continue
            t.append(v[0])
            s.append(v[7])
    if not t:
        raise SystemExit("no station rows in " + path)
    return np.array(t), np.array(s)
#
#
def trace_profile(rd, dip):
    rows = [np.loadtxt(f, ndmin=2) for f in sorted(glob.glob(os.path.join(rd, "diagTrace.*.txt"))) if os.path.getsize(f)]
    if not rows:
        raise SystemExit("no diagTrace rows in " + rd)
    a = np.vstack(rows)
    t = np.round(a[:, 1], 9)
    z = np.round(a[:, 2], 3) + 0.0
    tu = np.unique(t)
    tt = tu[np.argmin(np.abs(tu - T_EVAL))]
    arn, mS, mM, un = a[:, 3], a[:, 4], a[:, 5], a[:, 6:9]
    fS = a[:, 12:15] + a[:, 18:21]
    fM = a[:, 15:18] + a[:, 21:24]
    c = (mS * (fM * un).sum(1) - mM * (fS * un).sum(1)) / ((mS + mM) * arn) / 1e6
    zs = np.array(sorted(np.unique(z), reverse=True))
    prof = np.array([c[np.where((z == zz) & (t == tt))[0][0]] for zz in zs])
    return tt, -zs / np.sin(np.radians(dip)), prof
#
#
def read(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))
#
#
def t0_slopes(at):
    """initial normal traction per metre of down-dip distance (MPa/m), per case key, from node rows 2 and 3."""
    out = {}
    for name, key, label, dx, dip, cel, src in RUNS:
        if not cel:
            out.setdefault(key, 0.0)
            continue
        sl = []
        for r in at:
            if r["run"] == name and r["row"] in ("2", "3"):
                sl.append(float(r["T0_n"]) / (-float(r["z_m"]) / np.sin(np.radians(dip))))
        if not sl:
            continue
        if max(sl) - min(sl) > 1e-4 * abs(sl[0]):
            raise SystemExit("initial normal traction is not linear down dip for " + name)
        if key in out and abs(out[key] - sl[0]) > 1e-4 * abs(sl[0]):
            raise SystemExit("initial slope differs between grids for " + key)
        out[key] = sl[0]
    if "t12" in out:
        out.setdefault("t12hg2", out["t12"])
    return out
#
#
def means(runs_root, morning_root, at):
    slopes = t0_slopes(at)
    rows = []
    todo = [(runs_root, r) for r in RUNS]
    if morning_root:
        todo += [(morning_root, r) for r in MORNING]
    for root, (name, key, label, dx, dip, cel, src) in todo:
        rd = os.path.join(root, name)
        if not glob.glob(os.path.join(rd, "diagTrace.*.txt")):
            print("NOT RUN / no dump:", name)
            continue
        st_t, st_s = station(os.path.join(rd, "faultst000dp000.txt"))
        i = int(np.argmin(np.abs(st_t - T_EVAL)))
        tt, sd, prof = trace_profile(rd, dip)
        if key not in slopes:
            raise SystemExit("no initial-traction slope for " + key)
        row = [name, key, label, dx, src, "%.6f" % st_t[i], "%.4f" % st_s[i], len(st_t), "%.4f" % prof[0]]
        for L in LENGTHS:
            s = np.linspace(0.0, L, 2001)
            ch = np.trapezoid(np.interp(s, sd, prof), s) / L
            row += ["%.4f" % ch, "%.4f" % (ch + slopes[key] * L / 2.0)]
        rows.append(row)
    head = ["run", "case", "label", "dx_m", "source", "t_station_s", "dp000_nstress_MPa", "station_rows", "trace_node_change_MPa"]
    for L in LENGTHS:
        head += ["mean_change_top%dm_MPa" % L, "mean_total_top%dm_MPa" % L]
    return head, rows
#
#
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
#
#
def summary(at, es, cl):
    print("== closure (max over node rows): run, parts, inertia, |Ttrial-Tcode| MPa")
    for run in sorted(set(r["run"] for r in cl)):
        q = [r for r in cl if r["run"] == run]
        print("  %-13s %.1e %.1e %.1e" % (run, max(float(r["clos_parts_rel"]) for r in q), max(float(r["clos_inertia_rel"]) for r in q),
                                        max(float(r["max_abs_Ttrial_minus_Tcode_MPa"]) for r in q)))
    print("== 4 s budget, MPa: per side  BtS  hg  damp  grav  net  applied | kin  T0  Tcode")
    for r in at:
        if r["row"] not in ("1", "2"):
            continue
        for d in ("n", "s", "d"):
            txt = []
            for side in ("M", "S"):
                damp = sum(float(r["t%s_%s_%s" % (side, p, d)]) for p in ("mdamp", "kdamp", "hgdamp"))
                txt.append("%s %8.3f %8.3f %7.3f %7.3f %7.3f %8.3f" % (side, float(r["t%s_BtS_%s" % (side, d)]), float(r["t%s_hg_%s" % (side, d)]), damp,
                           float(r["t%s_grav_%s" % (side, d)]), float(r["t%snet_%s" % (side, d)]), float(r["t%sapp_%s" % (side, d)])))
            print("  %-13s row %s %s | %s | %s | %9.3f %8.3f %8.3f" % (r["run"], r["row"], d, txt[0], txt[1], float(r["kin_" + d]), float(r["T0_" + d]), float(r["Tcode_" + d])))
    print("== element stress beside the trace, MPa (mean over the x entries): run row side when n | nn nd zn zd zz | rjust pstrain")
    keys = sorted(set((r["run"], r["row"], r["side"], r["step"]) for r in es if r["row"] in ("1", "2")), key=lambda k: (k[0], k[1], k[2], float(k[3])))
    for run, row, side, step in keys:
        q = [r for r in es if (r["run"], r["row"], r["side"], r["step"]) == (run, row, side, step)]
        m = lambda k: float(np.mean([float(r[k]) for r in q]))
        print("  %-13s row %s side %2s t %7.4f n %d | %8.3f %8.3f %8.3f %8.3f %8.3f | %.5f %.3e" % (run, row, side, float(q[0]["t_s"]), len(q), m("sig_nn"), m("sig_nd"),
              m("sig_zn"), m("sig_zd"), m("sig_zz"), m("rjust"), m("pstrain")))
#
#
def figure(head, rows, fig_dir):
    d = [dict(zip(head, r)) for r in rows]
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.8))
    ax = axs[0]
    for key, src, col, mk, lab in (("t12", "budget build", "tab:blue", "o", "TPV12, stiffness hourglass (as shipped)"),
                                   ("t12hg2", "budget build + C_hg = 2", "tab:red", "s", "TPV12, viscous hourglass (scratch spike)"),
                                   ("t12", "morning build", "tab:blue", "x", "TPV12, uninstrumented build"),
                                   ("t13", "budget build", "tab:green", "o", "TPV13, as shipped"),
                                   ("t13", "morning build", "tab:green", "x", "TPV13, uninstrumented build")):
        pts = sorted((float(r["dx_m"]), float(r["dp000_nstress_MPa"])) for r in d if r["case"] == key and r["source"] == src)
        if pts:
            ax.plot([p[0] for p in pts], [p[1] for p in pts], "-" if mk != "x" else ":", marker=mk, color=col, label=lab)
    for key, col in (("t12", "tab:blue"), ("t13", "tab:green")):
        ax.plot([100], [FAULTMOD_100M[key]], "*", color=col, ms=11, mec="k", label="FaultMod 100 m, %s" % ("TPV12" if key == "t12" else "TPV13"))
    ax.axhline(0, color="0.7", lw=0.5)
    ax.set_xlabel("dx (m)")
    ax.set_ylabel("normal stress at the trace station, 4 s (MPa)")
    ax.set_title("Point value at the trace vs grid size", fontsize=9)
    ax.legend(fontsize=6.5)
    lo, hi = 0.0, 0.0
    for ax, key, title in ((axs[1], "t12", "TPV12: mean normal traction over a fixed down-dip length"),
                           (axs[2], "t13", "TPV13: mean normal traction over a fixed down-dip length")):
        for L, col in zip(LENGTHS, ("tab:orange", "tab:purple", "k")):
            for kk, src, ls, mk, tag in ((key, "budget build", "-", "o", "as shipped"), (key, "morning build", ":", "x", "uninstrumented build"),
                                         (key + "hg2", "budget build + C_hg = 2", "--", "s", "viscous hourglass spike")):
                pts = sorted((float(r["dx_m"]), float(r["mean_total_top%dm_MPa" % L])) for r in d if r["case"] == kk and r["source"] == src and float(r["dx_m"]) <= 2 * L)
                if pts:
                    ax.plot([p[0] for p in pts], [p[1] for p in pts], ls, marker=mk, color=col, mfc="none", label="top %d m, %s" % (L, tag))
                    lo = min(lo, min(p[1] for p in pts))
        pts = sorted((float(r["dx_m"]), float(r["dp000_nstress_MPa"])) for r in d if r["case"] == key and r["source"] == "morning build")
        if pts:
            ax.plot([p[0] for p in pts], [p[1] for p in pts], "-", marker="d", color="0.5", label="point value at the trace")
            lo = min(lo, min(p[1] for p in pts))
            hi = max(hi, max(p[1] for p in pts))
        ax.axhline(0, color="0.7", lw=0.5)
        ax.set_xlabel("dx (m)")
        ax.set_title(title, fontsize=9)
        ax.legend(fontsize=6, ncol=2)
    for ax in axs[1:]:
        ax.set_ylim(lo - 2.0, hi + 9.0)
    axs[1].set_ylabel("normal traction at 4 s (MPa)")
    fig.suptitle("Trace normal stress at 4 s: point value, viscous-hourglass spike, and fixed-length means (panels 2 and 3 share one scale)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = os.path.join(fig_dir, "fig_spike_and_means.png")
    fig.savefig(out, dpi=110)
    plt.close(fig)
    print("wrote", out)
#
#
if __name__ == "__main__":
    runs_root, csv_dir, fig_dir = sys.argv[1:4]
    morning_root = sys.argv[4] if len(sys.argv) > 4 else None
    at = read(os.path.join(csv_dir, "at4s.csv"))
    head, rows = means(runs_root, morning_root, at)
    write_csv(os.path.join(csv_dir, "means.csv"), head, rows)
    for r in rows:
        print("  " + "  ".join(str(x) for x in r))
    summary(at, read(os.path.join(csv_dir, "elem_stress.csv")), read(os.path.join(csv_dir, "closure.csv")))
    figure(head, rows, fig_dir)
