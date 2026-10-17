"""连通性判定：韧带附着节点与骨面节点的距离分布 —— 决定用合并节点还是 tied contact。"""
import numpy as np
from pathlib import Path

T = Path(r"D:\Project\climbing_fall_analysis\temp\thums")
z = np.load(T / "thums_lowerlimb.npz", allow_pickle=True)
nid, xyz = z["node_ids"], z["node_xyz"]
i_of = {int(v): i for i, v in enumerate(nid)}

BONES = {"TIBIA": [81000700, 81000600], "FIBULA": [81000900, 81000800],
         "TALUS": [81001100, 81001000], "CALCANEUS": [81001300, 81001200],
         "NAVICULAR": [81001400, 81001401]}
bone_pts, bone_nid = {}, {}
for nm, pids in BONES.items():
    pts, ids = [], []
    for pid in pids:
        for kind in ("solid", "shell"):
            k = f"{pid:08d}_{kind}"
            if k in z:
                u = sorted({int(v) for v in z[k].ravel()})
                pts.append(xyz[[i_of[n] for n in u]])
                ids += u
    bone_pts[nm] = np.vstack(pts)
    bone_nid[nm] = set(ids)

LIG = [(81100701, "PTFL 距腓后"), (81100801, "CFL 跟腓"),
       (81100901, "ATFL 距腓前"), (81101001, "PTCL 距跟后"),
       (81101101, "下胫腓"), (81101201, "内侧1 胫距后深"),
       (81101202, "内侧2 胫跟"), (81101203, "内侧3 胫舟"),
       (81101204, "内侧4 胫距前深"), (81101401, "跖长")]

print(f"{'pid':>9} {'名称':16s} {'节点':>4} {'重合<1e-6':>9} {'<0.1mm':>7} "
      f"{'<0.5mm':>7} {'<2mm':>6}  最近骨")
print("-" * 100)
for pid, nm in LIG:
    e = z[f"{pid:08d}_shell"]
    nodes = sorted({int(n) for n in e.ravel()})
    P = xyz[[i_of[n] for n in nodes]]
    best = min(bone_pts, key=lambda b: np.min(np.linalg.norm(
        bone_pts[b][None, :, :] - P[:, None, :], axis=2)))
    D = np.min(np.linalg.norm(
        bone_pts[best][None, :, :] - P[:, None, :], axis=2), axis=1)
    n0 = int((D < 1e-6).sum())
    n1 = int((D < 0.1).sum())
    n5 = int((D < 0.5).sum())
    n2 = int((D < 2.0).sum())
    print(f"{pid:>9} {nm:16s} {len(nodes):>4} {n0:>9} {n1:>7} {n5:>7} {n2:>6}  {best}")

print("\n\n=== 吴恺四条：几何先验 vs 吴恺排序 ===")
EA = {81101202: 510, 81101203: 510, 81100801: 510, 81101101: 75.4}
NAME = {81101202: "胫跟", 81101203: "胫舟", 81100801: "跟腓",
        81101101: "下胫腓前"}
sc = []
for pid in EA:
    e = z[f"{pid:08d}_shell"]
    nodes = sorted({int(n) for n in e.ravel()})
    P = xyz[[i_of[n] for n in nodes]]
    area = 0.0
    for q in e:
        p = xyz[[i_of[int(n)] for n in q]]
        area += 0.5 * np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0]))
        area += 0.5 * np.linalg.norm(np.cross(p[2] - p[0], p[3] - p[0]))
    c = P - P.mean(0)
    _, _, vt = np.linalg.svd(c, full_matrices=False)
    pr = c @ vt[0]
    L = pr.max() - pr.min()
    sc.append((EA[pid] * area / L, NAME[pid], EA[pid], area, L))
print("  吴恺:   胫跟 > 胫舟 > 跟腓 > 下胫腓前")
print("  先验:   " + " > ".join(nm for _, nm, *_ in sorted(sc, reverse=True)))
for k, nm, ea, a, L in sorted(sc, reverse=True):
    print(f"    {nm:10s} k={k:8.1f} N/mm  (EA={ea}, A={a:6.2f}, L={L:5.1f})")
