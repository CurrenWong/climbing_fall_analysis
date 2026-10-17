"""诊断：关节面节点的「法向 vs 指向对侧骨」夹角分布 —— 判断低对齐率是几何真相还是法向 bug。
同时对比两种节点法向：Newell 多边形法向 vs 三点法向。
"""
from collections import Counter
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

NPZ = Path(r"D:\Project\climbing_fall_analysis\temp\thums\thums_lowerlimb.npz")
TET_FACES = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
HEX_FACES = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
             (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
WEDGE_FACES = [(0, 2, 1), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)]
GROUPS = {"tibia": ["81000700_solid", "81000600_solid", "81000601_solid"],
          "talus": ["81001100_solid", "81001000_solid"]}

d = np.load(NPZ)
nid, xyz = d["node_ids"], d["node_xyz"]
idx = {int(v): i for i, v in enumerate(nid)}


def boundary(keys):
    acc = {}
    for k in keys:
        for row in d[k]:
            uniq = sorted(set(int(x) for x in row))
            tmpl = {4: TET_FACES, 8: HEX_FACES, 6: WEDGE_FACES}.get(len(uniq))
            if not tmpl:
                continue
            P = xyz[[idx[v] for v in uniq]]
            cen = P.mean(0)
            for f in tmpl:
                if max(f) >= len(uniq):
                    continue
                verts = [uniq[z] for z in f]
                pts = xyz[[idx[v] for v in verts]]
                n = np.zeros(3)
                m = len(pts)
                for i in range(m):
                    n += np.cross(pts[i], pts[(i + 1) % m])
                if np.dot(n, pts.mean(0) - cen) < 0:
                    n = -n
                    verts = verts[::-1]
                key = tuple(sorted(verts))
                if key in acc:
                    acc[key]["cnt"] += 1
                else:
                    acc[key] = {"verts": tuple(verts), "n": n, "cnt": 1, "pts": pts}
    return [v for v in acc.values() if v["cnt"] == 1]


fA = boundary(GROUPS["tibia"])
fB = boundary(GROUPS["talus"])
nA = set().union(*[set(f["verts"]) for f in fA])
nB = set().union(*[set(f["verts"]) for f in fB])
idsA = np.array(sorted(nA))
PA = xyz[[idx[int(v)] for v in idsA]]
PB = xyz[[idx[int(v)] for v in sorted(nB)]]
extra = np.array([f["pts"].mean(0) for f in fB])
QB = np.vstack([PB, extra])

# 两种法向
def normals_newell():
    acc = {v: np.zeros(3) for v in nA}
    for f in fA:
        n = f["n"]
        ln = np.linalg.norm(n)
        if ln < 1e-12:
            continue
        unit = n / ln
        for v in f["verts"]:
            if v in acc:
                acc[v] += unit * ln            # ln ∝ 2×面积
    return {v: (s / np.linalg.norm(s) if np.linalg.norm(s) > 1e-12 else np.array([0., 0., 1.]))
            for v, s in acc.items()}

def normals_tri():
    acc = {v: np.zeros(3) for v in nA}
    for f in fA:
        p = f["pts"]
        c = np.cross(p[1] - p[0], p[2] - p[0]) * 0.5
        for v in f["verts"]:
            if v in acc:
                acc[v] += c
    return {v: (s / np.linalg.norm(s) if np.linalg.norm(s) > 1e-12 else np.array([0., 0., 1.]))
            for v, s in acc.items()}

tq = cKDTree(QB)
gap, nq = tq.query(PA, k=1)
dirv = QB[nq] - PA
dirv /= np.maximum(np.linalg.norm(dirv, axis=1, keepdims=True), 1e-12)

for name, fn in (("Newell", normals_newell()), ("三点", normals_tri())):
    N = np.array([fn[int(v)] for v in idsA])
    al = np.einsum("ij,ij->i", N, dirv)
    sel = gap < 6.0
    print(f"\n=== {name} 法向 ===")
    print(f"  gap<6mm 的 {int(sel.sum())} 个节点，align 分布：")
    a = al[sel]
    print(f"    min {a.min():.3f}  p25 {np.percentile(a,25):.3f}  中位 "
          f"{np.median(a):.3f}  p75 {np.percentile(a,75):.3f}  max {a.max():.3f}")
    for thr in (0.0, 0.3, 0.5, 0.7):
        print(f"    align > {thr}: {int((a>thr).sum()):>4d}")
    print("  最接近的 12 个节点的 gap / align / 坐标:")
    order = np.argsort(gap)[:12]
    for i in order:
        v = int(idsA[i])
        p = PA[i]
        print(f"    gap {gap[i]:7.3f}  align {al[i]:7.3f}  "
              f"xyz ({p[0]/10:7.2f},{p[1]/10:7.2f},{p[2]/10:7.2f}) cm")
