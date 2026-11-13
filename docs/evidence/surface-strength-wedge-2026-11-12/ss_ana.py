#!/usr/bin/env python3
"""Surface-strength wedge tests: station metrics, force parts at the trace node, kill-criterion inputs, figures.
usage: ss_ana.py SS_ROOT W12_BASE_ROOT ROW149_EVIDENCE_DIR WEDGE_TRACE_EVIDENCE_DIR OUT_DIR [--allow-missing]
  SS_ROOT                  this evidence's runs (w12c4_*, i12c4_*, i12c02_*, t36c4_*, t36c02_*)
  W12_BASE_ROOT            row-149 wedge baseline runs w12_dx500/200/150 (TPV12 wedge, 0.2 MPa), read-only
  ROW149_EVIDENCE_DIR      committed row-149 wedge evidence (insertion 0.2 MPa numbers whose raw runs are deleted)
  WEDGE_TRACE_EVIDENCE_DIR committed tpv36 wedge-trace evidence (4 MPa baseline at 500/250/200 m and own 50 m archive)
Station columns: t, h-slip, h-rate, h-shear, v-slip, v-rate, v-shear, n-stress (MPa, positive = tension).
Strength = c - mu(slip) * min(n, 0), linear slip weakening (faulting.f90 NewtonRaphson, friclaw 1).
A missing run is an error unless --allow-missing (then printed NOT RUN)."""
import csv, glob, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

args = [a for a in sys.argv[1:] if a != "--allow-missing"]
ALLOW = "--allow-missing" in sys.argv[1:]
SS, WB, R149, WTR, OUT = args[:5]
os.makedirs(OUT, exist_ok=True)
CASE = {  # dip, t_eval, mu_s, mu_d, D0
    "tpv12": (60.0, 4.0, 0.70, 0.10, 0.50),
    "tpv36": (15.0, 13.0, 0.575, 0.45, 0.18),
}
RUNS = [  # label, root, dir, case, mesh, trace cohesion MPa, dx
    ("wedge c=0.2 (baseline)", WB, "w12_dx500", "tpv12", "wedge", 0.2, 500),
    ("wedge c=0.2 (baseline)", WB, "w12_dx200", "tpv12", "wedge", 0.2, 200),
    ("wedge c=0.2 (baseline)", WB, "w12_dx150", "tpv12", "wedge", 0.2, 150),
    ("wedge c=0.2 (16-rank decomposition control)", SS, "w12c02_dx500_p422", "tpv12", "wedge", 0.2, 500),
    ("wedge c=4 uniform (degenerate)", SS, "w12u4_dx500", "tpv12", "wedge", 4.0, 500),
    ("wedge c=4 taper (test 1)", SS, "w12c4_dx500", "tpv12", "wedge", 4.0, 500),
    ("wedge c=4 taper (test 1)", SS, "w12c4_dx200", "tpv12", "wedge", 4.0, 200),
    ("wedge c=4 taper (test 1)", SS, "w12c4_dx150", "tpv12", "wedge", 4.0, 150),
    ("insertion c=4 taper", SS, "i12c4_dx500", "tpv12", "insertion", 4.0, 500),
    ("insertion c=4 taper", SS, "i12c4_dx200", "tpv12", "insertion", 4.0, 200),
    ("insertion c=4 taper", SS, "i12c4_dx150", "tpv12", "insertion", 4.0, 150),
    ("insertion c=0.2 (rerun)", SS, "i12c02_dx150", "tpv12", "insertion", 0.2, 150),
    ("wedge c=4 (baseline rerun)", SS, "t36c4_dx500", "tpv36", "wedge", 4.0, 500),
    ("wedge c=4 (baseline rerun)", SS, "t36c4_dx250", "tpv36", "wedge", 4.0, 250),
    ("wedge c=0.2 (test 2)", SS, "t36c02_dx500", "tpv36", "wedge", 0.2, 500),
    ("wedge c=0.2 (test 2)", SS, "t36c02_dx250", "tpv36", "wedge", 0.2, 250),
    ("wedge c=0.2 (test 2)", SS, "t36c02_dx200", "tpv36", "wedge", 0.2, 200),
]


def fl(x):
    if "E" not in x and re.search(r"\d[-+]\d{3}$", x):
        x = x[:-4] + "E" + x[-4:]
    return float(x)


def station(p):
    rows = []
    for line in open(p):
        s = line.split()
        if len(s) < 8 or s[0].startswith("#") or s[0] == "t":
            continue
        rows.append([fl(x) for x in s[:8]])
    if not rows:
        raise SystemExit("no station rows in " + p)
    return np.array(rows)


UNIFORM = {"w12u4_dx500"}   # uniform 4 MPa (degenerate, no nucleation); every other tpv12 c=4 run is the surface taper


def cohesion(case, c_trace, dd, uniform=False):
    if case == "tpv12":   # test 1: c_trace at the trace, linear to the unchanged 0.2 MPa at 8 km down dip (set_case.py), or uniform
        if uniform or c_trace <= 0.2:
            return c_trace
        return 0.2 + (c_trace - 0.2) * (0.0 if dd >= 8000.0 else 1.0 - dd / 8000.0)
    return c_trace / 8000.0 * (8000.0 - dd) if dd <= 8000.0 else 0.0   # tpv36: slope*(8e3 - dd), slope = c_trace/8e3


def metrics(a, case, c_trace, dd, uniform=False):
    dip, te, mus, mud, d0 = CASE[case]
    rate = np.hypot(a[:, 2], a[:, 5])
    on = np.flatnonzero(rate > 1e-3)
    k0 = on[0] if len(on) else len(a)
    i = int(np.argmin(np.abs(a[:, 0] - te))) if case == "tpv12" else int(np.searchsorted(a[:, 0], te, side="right") - 1)
    slip = np.hypot(a[:, 1], a[:, 4])
    mu = mus - (mus - mud) * np.minimum(slip, d0) / d0
    c = cohesion(case, c_trace, dd, uniform)
    strength = c - mu * np.minimum(a[:, 7], 0.0)
    shear = np.hypot(a[:, 3], a[:, 6])
    ratio = shear / np.maximum(strength, 1e-9)
    post = slice(k0, i + 1)
    ten = a[post, 7] > 0.0
    pinned = np.logical_and(ten, ratio[post] >= 0.99)
    return dict(arrival_s=a[k0, 0] if len(on) else float("nan"), t_eval_s=a[i, 0], dip_slip_m=a[i, 4], slip_m=slip[i],
                n_MPa=a[i, 7], n_max_post_arrival_MPa=a[post, 7].max() if k0 <= i else float("nan"), n_min_MPa=a[:i + 1, 7].min(),
                n_max_MPa=a[:i + 1, 7].max(), shear_MPa=shear[i], strength_MPa=strength[i], shear_over_strength=ratio[i],
                cohesion_MPa=c, frac_post_arrival_tensile=ten.mean() if k0 <= i else 0.0,
                frac_post_arrival_tensile_and_pinned=pinned.mean() if k0 <= i else 0.0, rows=len(a))


def parts(rd, dip, te):
    """trace-node traction change and its per-side parts from the instrumented diagTrace dump (method of row-149 wedge_ana.py)."""
    files = [f for f in sorted(glob.glob(os.path.join(rd, "diagTrace.*.txt"))) if os.path.getsize(f)]
    if not files:
        return None
    a = np.vstack([np.array([[fl(x) for x in l.split()] for l in open(f) if l.strip()]) for f in files])   # fl: Fortran 3-digit exponent without E
    t = np.round(a[:, 1], 9)
    z = np.round(a[:, 2], 3) + 0.0
    tu = np.unique(t)
    tt = tu[np.argmin(np.abs(tu - te))]
    arn, mS, mM, un = a[:, 3], a[:, 4], a[:, 5], a[:, 6:9]
    wM, wS = mS / (mS + mM), mM / (mS + mM)
    per = 1.0 / (arn * 1e6)
    p = dict(ku_M=(a[:, 15:18] * un).sum(1) * per, hg_M=(a[:, 21:24] * un).sum(1) * per,
             ku_S=-(a[:, 12:15] * un).sum(1) * per, hg_S=-(a[:, 18:21] * un).sum(1) * per)
    tot = wM * (p["ku_M"] + p["hg_M"]) + wS * (p["ku_S"] + p["hg_S"])
    j = np.flatnonzero(np.logical_and(z == z.max(), t == tt))[0]
    return dict(parts_t_s=tt, parts_sum_MPa=tot[j], stressdiv_hw_MPa=p["ku_M"][j], stressdiv_fw_MPa=p["ku_S"][j],
                hourglass_hw_MPa=p["hg_M"][j], hourglass_fw_MPa=p["hg_S"][j])


KEYS = ["arrival_s", "t_eval_s", "dip_slip_m", "slip_m", "n_MPa", "n_max_post_arrival_MPa", "n_min_MPa", "n_max_MPa", "shear_MPa",
        "strength_MPa", "shear_over_strength", "cohesion_MPa", "frac_post_arrival_tensile", "frac_post_arrival_tensile_and_pinned", "rows"]
PK = ["parts_t_s", "parts_sum_MPa", "stressdiv_hw_MPa", "stressdiv_fw_MPa", "hourglass_hw_MPa", "hourglass_fw_MPa"]
table, series = {"tpv12": [], "tpv36": []}, {}
for lab, root, nm, case, mesh, c0, dx in RUNS:
    rd = os.path.join(root, nm)
    if not os.path.isfile(os.path.join(rd, "DONE")):
        if ALLOW:
            print("NOT RUN / not finished:", nm)
            continue
        raise SystemExit("missing or unfinished run " + rd)
    pr = parts(rd, CASE[case][0], CASE[case][1])
    stns = sorted(os.path.basename(p)[7:-4] for p in glob.glob(os.path.join(rd, "faultst000dp*.txt")))
    for s in stns:
        a = station(os.path.join(rd, "faultst%s.txt" % s))
        dd = float(s[5:]) * 100.0
        m = metrics(a, case, c0, dd, nm in UNIFORM)
        series[(nm, s)] = a
        row = [nm, lab, mesh, "%g" % c0, dx, s, "%g" % dd] + ["%.6g" % m[k] for k in KEYS]
        row += ["%.6g" % pr[k] if (pr and s == "000dp000") else "" for k in PK]
        table[case].append(row)
head = ["run", "label", "mesh", "trace_cohesion_MPa", "dx_m", "station", "down_dip_m"] + KEYS + PK
for case, fn in (("tpv12", "test1_tpv12.csv"), ("tpv36", "test2_tpv36.csv")):
    p = os.path.join(OUT, fn)
    with open(p, "w", newline="") as f:
        csv.writer(f).writerows([head] + table[case])
    back = list(csv.reader(open(p)))
    assert len(back) == len(table[case]) + 1 and all(len(b) == len(head) for b in back), p
    print("wrote %s: %d rows x %d cols (re-read ok)" % (p, len(table[case]), len(head)))


def get(case, nm, s, k):
    for r in table[case]:
        if r[0] == nm and r[5] == s:
            return float(r[head.index(k)])
    return float("nan")


print("\nBASELINE REPRODUCTION (committed numbers in brackets)")
for nm, ref in (("w12_dx500", 30.0603), ("w12_dx200", 84.9702), ("w12_dx150", 110.5110), ("i12c02_dx150", -2.2977)):
    print("  %-13s dp000 n at 4 s %+.4f [%+.4f]  dip slip %.4f" % (nm, get("tpv12", nm, "000dp000", "n_MPa"), ref, get("tpv12", nm, "000dp000", "dip_slip_m")))
for nm, s, ref in (("w12_dx150", "000dp015", 5.94), ("w12_dx150", "000dp030", 7.95), ("i12c02_dx150", "000dp015", 8.38), ("i12c02_dx150", "000dp030", 9.19)):
    print("  %-13s %s dip slip at 4 s %.4f [%.2f]" % (nm, s, get("tpv12", nm, s, "dip_slip_m"), ref))
wt = list(csv.DictReader(open(os.path.join(WTR, "trace_summary.csv"))))
for r in wt:
    if r["case"] == "tpv36" and r["station"] in ("000dp000", "000dp010") and r["dx_m"] in ("500", "250"):
        nm = "t36c4_dx" + r["dx_m"]
        print("  %-13s %s n_max %s [%s] slip13 %s [%s] arrival %s [%s]" % (nm, r["station"], "%.5g" % get("tpv36", nm, r["station"], "n_max_MPa"), r["n_max_MPa"],
              "%.5g" % get("tpv36", nm, r["station"], "slip_m"), r["slip_at_t_eval_m"], "%.5g" % get("tpv36", nm, r["station"], "arrival_s"), r["arrival_s"]))
nm0 = "w12c02_dx500_p422"
print("\nDECOMPOSITION CONTROL (16 ranks 4x2x2 vs baseline 4 ranks 2x2x1, wedge c=0.2, 500 m)")
for s_ in ("000dp000", "000dp015", "000dp075"):
    print("  %s n %+.5f vs %+.5f  dip slip %.5f vs %.5f" % (s_, get("tpv12", nm0, s_, "n_MPa"), get("tpv12", "w12_dx500", s_, "n_MPa"),
          get("tpv12", nm0, s_, "dip_slip_m"), get("tpv12", "w12_dx500", s_, "dip_slip_m")))
print("\nKILL-CRITERION INPUTS")
for dx in (200, 150):
    w, i = get("tpv12", "w12c4_dx%d" % dx, "000dp000", "n_MPa"), get("tpv12", "i12c4_dx%d" % dx, "000dp000", "n_MPa")
    sw, si = get("tpv12", "w12c4_dx%d" % dx, "000dp030", "dip_slip_m"), get("tpv12", "i12c4_dx%d" % dx, "000dp030", "dip_slip_m")
    live = sw > 0.1 and si > 0.1   # a run with no rupture would pass C1 vacuously
    print("  C1 dx %d: wedge c=4 dp000 n %+.4f, insertion c=4 %+.4f, excess %+.4f MPa (dp030 slip %.3f / %.3f m) -> %s" % (
        dx, w, i, w - i, sw, si, ("PASS" if (w <= 0 or w - i <= 5.0) else "FAIL") if live else "DEGENERATE (no rupture)"))
for s in ("000dp015", "000dp030", "000dp000"):
    w, i, i2 = (get("tpv12", n, s, "dip_slip_m") for n in ("w12c4_dx150", "i12c4_dx150", "i12c02_dx150"))
    print("  C2 dx 150 %s: wedge c=4 %.4f m, insertion c=4 %.4f m (loss %.1f %%), insertion c=0.2 %.4f m (loss %.1f %%)%s" % (
        s, w, i, 100 * (1 - w / i), i2, 100 * (1 - w / i2), ("  -> %s" % ("PASS" if 1 - w / i < 0.05 else "FAIL")) if s == "000dp015" else ""))
for dx in (500, 250, 200):
    nm = "t36c02_dx%d" % dx
    nmax, pin = get("tpv36", nm, "000dp000", "n_max_post_arrival_MPa"), get("tpv36", nm, "000dp000", "frac_post_arrival_tensile_and_pinned")
    print("  C3 dx %d: tpv36 c=0.2 dp000 n max after arrival %+.4f MPa, tensile-and-pinned fraction %.3f -> %s" % (
        dx, nmax, pin, "TENSION/PINNING" if (nmax > 1.0 or pin > 0) else "none"))

# figures: one per test, against its baseline
for case, fn, runs, panels in (
        ("tpv12", "fig_test1_tpv12.png", ["w12_dx200", "w12_dx150", "w12c4_dx200", "w12c4_dx150", "i12c4_dx150", "i12c02_dx150"],
         [("000dp000", 7, "trace dp000 n-stress (MPa, +tension)"), ("000dp000", 4, "trace dp000 dip slip (m)"), ("000dp015", 4, "dp015 (1.5 km) dip slip (m)")]),
        ("tpv36", "fig_test2_tpv36.png", ["t36c4_dx500", "t36c4_dx250", "t36c02_dx500", "t36c02_dx250", "t36c02_dx200"],
         [("000dp000", 7, "trace dp000 n-stress (MPa, +tension)"), ("000dp000", -1, "trace dp000 slip (m)"), ("000dp010", -1, "dp010 (1 km) slip (m)")])):
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    for nm in runs:
        lab = [r for r in RUNS if r[2] == nm][0]
        for x, (s, col, title) in zip(ax, panels):
            if (nm, s) not in series:
                continue
            a = series[(nm, s)]
            y = np.hypot(a[:, 1], a[:, 4]) if col < 0 else a[:, col]
            x.plot(a[:, 0], y, "--" if "c=0.2" in lab[0] and case == "tpv12" or "c=4" in lab[0] and case == "tpv36" else "-",
                   label="%s %s %d m" % (lab[0], lab[4] if case == "tpv12" else "", lab[6]))
            x.set_title(title, fontsize=9)
    for x in ax:
        x.set_xlabel("t (s)"); x.grid(alpha=.3)
    ax[0].axhline(0, color="k", lw=.5); ax[0].legend(fontsize=6)
    fig.suptitle("%s wedge mesh: surface strength %s (dashed = baseline)" % (case, "0.2 -> 4 MPa at trace, tapered to 0.2 MPa at 8 km down dip" if case == "tpv12" else "4 -> 0.2 MPa"), fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, fn), dpi=120); plt.close(fig)
    print("wrote", os.path.join(OUT, fn))
