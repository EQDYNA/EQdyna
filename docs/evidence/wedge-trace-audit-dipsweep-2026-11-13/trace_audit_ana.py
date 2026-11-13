"""Trace-node audit analysis. usage: trace_audit_ana.py <run_dir> <dip_deg> <label> [npairs=4]
Reads audit_pairs.0000.txt / audit_elems.0000.txt written by the scratch diagAudit dump (serial run).
Prints: per-pair table, per-id attachment table, side verdict. Writes <label>_pairs.csv and <label>_attach.csv to cwd."""
import sys, math, csv
import numpy as np
d, dip, label = sys.argv[1], float(sys.argv[2]), sys.argv[3]
npairs = int(sys.argv[4]) if len(sys.argv) > 4 else 4
t = math.tan(math.radians(dip)); c = math.cos(math.radians(dip))
hdr = open(d + "/audit_pairs.0000.txt").readline().split()[1:]
P = np.loadtxt(d + "/audit_pairs.0000.txt")
E = np.loadtxt(d + "/audit_elems.0000.txt")
print(label, "header:", dict(zip(["dx","dy","dz","tol","faultDegenAngle","fymin","nftnd"], hdr[8:] if False else open(d+"/audit_pairs.0000.txt").readline().split()[8:])))
# columns: i ns nm x y z arn fnms_s fnms_m un(3) us(3) ud(3)
order = np.argsort(-P[:, 5])
sel = P[order]
signed = lambda y, z: (z + y * t) * c          # >0: hanging wall (above plane, y+ wall)
pr = []
for r in sel[:npairs + 6]:
    pr.append(dict(i=int(r[0]), slave=int(r[1]), master=int(r[2]), y=r[4], z=r[5], arn=r[6], mass_s=r[7], mass_m=r[8],
                   un=r[9:12].tolist(), us=r[12:15].tolist(), ud=r[15:18].tolist()))
ref = pr[npairs + 3]
with open(label + "_pairs.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["pair","slave","master","y","z","arn","arn_ratio_to_ref","mass_slave","mass_master","mass_s_ratio","mass_m_ratio","un_1","un_2","un_3","us_1","us_2","us_3","ud_1","ud_2","ud_3"])
    for k, p in enumerate(pr):
        w.writerow([k, p["slave"], p["master"], p["y"], p["z"], p["arn"], p["arn"]/ref["arn"], p["mass_s"], p["mass_m"], p["mass_s"]/ref["mass_s"], p["mass_m"]/ref["mass_m"]] + p["un"] + p["us"] + p["ud"])
        if k < npairs + 4: print("pair%d y=%.1f z=%.1f arn=%.4g (x%.3f ref) mass s/m=%.4g/%.4g (x%.3f / x%.3f ref) un=%s us=%s ud=%s" % (k, p["y"], p["z"], p["arn"], p["arn"]/ref["arn"], p["mass_s"], p["mass_m"], p["mass_s"]/ref["mass_s"], p["mass_m"]/ref["mass_m"], np.round(p["un"],1), np.round(p["us"],1), np.round(p["ud"],1)))
print("reference pair for ratios: pair%d (z=%.1f)" % (npairs + 3, ref["z"]))
# attachments
ids = {}
for p in pr[:npairs]:
    ids[p["slave"]] = "slave"; ids[p["master"]] = "master"
rows = []
for e in E:
    eid, ty = int(e[0]), int(e[1]); nodes = e[2:10].astype(int); yz = e[10:26].reshape(8, 2)
    dist = np.array([signed(y, z) for y, z in yz]); nf = np.isin(nodes, list(ids))
    # off-plane nodes (not on the fault plane within 1 m): decide side
    off = np.abs(dist) > 1.0
    s = dist[off]
    side = "plane" if s.size == 0 else ("HW" if (s > 0).all() else ("FW" if (s < 0).all() else "MIXED"))
    cy, cz = yz[:, 0].mean(), yz[:, 1].mean()
    for n in sorted(set(nodes[nf].tolist())):   # one row per (node id, element), even if the element repeats the id (collapsed wedges)
        rows.append(dict(node=int(n), kind=ids[int(n)], pair=[q["slave"] if ids[int(n)]=="slave" else q["master"] for q in pr[:npairs]].index(int(n)), elem=eid, type=ty, cy=cy, cz=cz, side=side, n_off=int(off.sum())))
with open(label + "_attach.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print("--- attachments (pair, id kind, #elements by type/side)")
from collections import Counter
for k in range(npairs):
    for kind in ("slave", "master"):
        rr = [r for r in rows if r["pair"] == k and r["kind"] == kind]
        cnt = Counter((r["type"], r["side"]) for r in rr)
        print("pair%d %-6s n_elem=%d  %s" % (k, kind, len(rr), dict(sorted(cnt.items()))))
bad_m = [r for r in rows if r["kind"] == "master" and r["side"] in ("FW", "MIXED")]
bad_s = [r for r in rows if r["kind"] == "slave" and r["side"] in ("HW", "MIXED")]
tot_m = sum(1 for r in rows if r["kind"] == "master"); tot_s = sum(1 for r in rows if r["kind"] == "slave")
print("VERDICT %s: master-attached elems on FW/MIXED side: %d of %d; slave-attached on HW/MIXED side: %d of %d; distinct (node id, element) pairs %d; elements read %d" % (label, len(bad_m), tot_m, len(bad_s), tot_s, len(rows), len(E)))
for r in bad_m + bad_s: print("  MISWIRE?", r)
