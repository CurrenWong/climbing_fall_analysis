"""CARTFIX 验收 ②：真·零穿透检查（对侧骨的真实外表面 + 射线奇偶内测）。

对侧骨表面 = 该骨部件(CORT∪SPON)单元面中**只出现 1 次**的面（与 cartilage.py 同法）。
对每个软骨顶面节点：
   d = min 到对侧骨表面的精确点-三角距离（无符号，必定 ≥0）
   inside = 射线奇偶（3 个方向多数表决）
判据：d ≥ 0（距离恒真）+ **inside 必须全 0**（这才是真穿透）。
"""
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
NPZ = ROOT / "temp/thums/thums_lowerlimb.npz"
db = np.load(NPZ, allow_pickle=True)
xyz = db["node_xyz"]
idx = {int(v): i for i, v in enumerate(db["node_ids"])}

TALUS = ["81001100_solid", "81001000_solid"]
TIBIA = ["81000700_solid", "81000600_solid", "81000601_solid"]


def surface_tris(keys):
    cnt = {}
    for k in keys:
        for row in db[k]:
            v = [int(z) for z in row if int(z) not in (0, -1)]
            v = list(dict.fromkeys(v))
            if len(v) < 3:
                continue
            if len(v) == 4:            # 四面体: 4 个面
                faces = [(v[0], v[1], v[2]), (v[0], v[1], v[3]),
                         (v[0], v[2], v[3]), (v[1], v[2], v[3])]
            elif len(v) == 8:          # 六面体 → 6 个四边形面 → 各拆 2 三角
                H = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4),
                     (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
                faces = []
                for q in H:
                    a, b, c, d = (v[q[0]], v[q[1]], v[q[2]], v[q[3]])
                    faces += [(a, b, c), (a, c, d)]
            else:
                faces = [(v[0], v[1], v[2])]
            for f in faces:
                cnt[tuple(sorted(f))] = cnt.get(tuple(sorted(f)), 0) + 1
    return [f for f, c in cnt.items() if c == 1]


def tris_xyz(faces):
    out = []
    for f in faces:
        a, b, c = f
        if a in idx and b in idx and c in idx:
            out.append([xyz[idx[a]], xyz[idx[b]], xyz[idx[c]]])
    return np.array(out)


def pt_tri_dist(p, T):
    """p:(3,) T:(n,3,3) -> (n,) 精确点-三角距离（Ericson 区域法，向量化）。"""
    a, b, c = T[:, 0], T[:, 1], T[:, 2]
    ab, ac = b - a, c - a
    ap = p - a
    d1 = np.einsum("ij,ij->i", ab, ap)
    d2 = np.einsum("ij,ij->i", ac, ap)
    bp = p - b
    d3 = np.einsum("ij,ij->i", ab, bp)
    d4 = np.einsum("ij,ij->i", ac, bp)
    cp = p - c
    d5 = np.einsum("ij,ij->i", ab, cp)
    d6 = np.einsum("ij,ij->i", ac, cp)
    q = np.empty_like(a)
    # 默认：面内
    vc = d1 * d4 - d3 * d2
    vb = d5 * d2 - d1 * d6
    va = d3 * d6 - d5 * d4
    den = 1.0 / (va + vb + vc)
    q = a + ab * (vb * den)[:, None] + ac * (vc * den)[:, None]
    m = (d1 <= 0) & (d2 <= 0)
    q[m] = a[m]
    m = (d3 >= 0) & (d4 <= d3)
    q[m] = b[m]
    m = (vc <= 0) & (d1 >= 0) & (d3 <= 0)
    q[m] = a[m] + (d1 / (d1 - d3))[m, None] * ab[m]
    m = (d6 >= 0) & (d5 <= d6)
    q[m] = c[m]
    m = (vb <= 0) & (d2 >= 0) & (d6 <= 0)
    q[m] = a[m] + (d2 / (d2 - d6))[m, None] * ac[m]
    m = (va <= 0) & ((d4 - d3) >= 0) & ((d5 - d6) >= 0)
    q[m] = b[m] + ((d4 - d3) / ((d4 - d3) + (d5 - d6)))[m, None] * (c[m] - b[m])
    return np.linalg.norm(q - p, axis=1)


def ray_inside(p, T, Tn, d=None):
    """射线奇偶：返回 True=inside。用 3 个方向多数表决。"""
    cnt = 0
    for k, direc in enumerate(((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))):
        dd = np.array(direc) + 0.0003 * np.array((k * 0.11 + 0.05, k * 0.07 + 0.03, 0.017))
        dd = dd / np.linalg.norm(dd)
        e1 = T[:, 1] - T[:, 0]
        e2 = T[:, 2] - T[:, 0]
        h = np.cross(dd, e2)
        a = np.einsum("ij,ij->i", e1, h)
        ok = np.abs(a) > 1e-12
        f = 1.0 / np.where(ok, a, 1.0)
        s = p - T[:, 0]
        u = f * np.einsum("ij,ij->i", s, h)
        q = np.cross(s, e1)
        v = f * (q @ dd)
        t = f * np.einsum("ij,ij->i", e2, q)
        hit = ok & (u >= 0) & (u <= 1) & (v >= 0) & (u + v <= 1) & (t > 1e-9)
        if int(hit.sum()) % 2 == 1:
            cnt += 1
    return cnt >= 2


TB = tris_xyz(surface_tris(TALUS))
TT = tris_xyz(surface_tris(TIBIA))
print(f"对侧骨表面三角: talus {len(TB)}  tibia {len(TT)}")

for tag, path in (("BASELINE cart_m1", ROOT / "temp/thums/cart_m1/cartilage_mesh.npz"),
                  ("FIXED   cart_m1_fix", ROOT / "temp/thums/cart_m1_fix/cartilage_mesh.npz")):
    d = np.load(path, allow_pickle=True)
    print("=" * 88)
    print(tag)
    for j, Bt, nm in (("tibiotalar_tibia", TB, "talus"), ("tibiotalar_talus", TT, "tibia")):
        N = np.asarray(d[f"{j}_nodes"], float)
        off = np.asarray(d[f"{j}_is_offset"], bool)
        pts = N[off]
        dmin = np.array([pt_tri_dist(p, Bt).min() for p in pts])
        nins = sum(1 for p in pts if ray_inside(p, Bt, None))
        print(f"  {j:20s} 顶面 {len(pts):3d} 节点 → {nm} 骨面: "
              f"距离 min {dmin.min():+.4f} 中位 {np.median(dmin):+.4f} mm  "
              f"(<0.05mm 的 {int((dmin<0.05).sum())})  射线内测命中(真穿透) {nins}")
