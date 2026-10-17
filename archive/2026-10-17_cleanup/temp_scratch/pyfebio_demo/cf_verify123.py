"""CARTFIX 验收 ①②③：完整采样 det(J) / 零穿透（有符号） / 几何厚度。

穿透判据（有符号）：把 A 侧软骨的**顶面节点**对 B 侧软骨的**底面（= 骨面）三角网**
求最近点 + 有向距离 d = n̂·(p − q)，n̂ 取「B 侧底面的法向指向关节腔」那一侧。
d < 0 ⇒ A 侧顶面越过了 B 侧骨面 ⇒ 真穿透。
"""
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(r"D:\Project\climbing_fall_analysis")
RS = (0.0, 1 / 6, 1 / 3, 1 / 2, 2 / 3, 1.0)
TS = (0.0, 0.25, 0.5, 0.75, 1.0)
GS = [(r, s, t) for r in RS for s in RS if r + s <= 1 + 1e-9 for t in TS]
GG = np.array([np.array([
    [-(1 - t), -(1 - t), -(1 - r - s)],
    [1 - t, 0.0, -r],
    [0.0, 1 - t, -s],
    [-t, -t, (1 - r - s)],
    [t, 0.0, r],
    [0.0, t, s],
]) for (r, s, t) in GS])


def dets(P6):
    return np.linalg.det(np.einsum("nic,mia->nmca", P6, GG))


def ix(x):
    x = int(x)
    return (-x - 1) if x < 0 else (x - 1)


def load(p):
    d = np.load(p, allow_pickle=True)
    out = {}
    for j in ("tibiotalar_tibia", "tibiotalar_talus"):
        N = np.asarray(d[f"{j}_nodes"], float)
        E = np.asarray(d[f"{j}_elems"], np.int64)
        off = np.asarray(d[f"{j}_is_offset"], bool)
        C = np.array([[ix(x) for x in row] for row in E])
        out[j] = dict(N=N, E=E, off=off, C=C, B=C[:, :3], T=C[:, 3:])
    return out


def closest_on_tri(p, a, b, c):
    """点 p 到三角形 (a,b,c) 的最近点（Ericson 算法）。"""
    ab, ac, ap = b - a, c - a, p - a
    d1, d2 = ab @ ap, ac @ ap
    if d1 <= 0 and d2 <= 0:
        return a
    bp = p - b
    d3, d4 = ab @ bp, ac @ bp
    if d3 >= 0 and d4 <= d3:
        return b
    vc = d1 * d4 - d3 * d2
    if vc <= 0 and d1 >= 0 and d3 <= 0:
        return a + (d1 / (d1 - d3)) * ab
    cp = p - c
    d5, d6 = ab @ cp, ac @ cp
    if d6 >= 0 and d5 <= d6:
        return c
    vb = d5 * d2 - d1 * d6
    if vb <= 0 and d2 >= 0 and d6 <= 0:
        return a + (d2 / (d2 - d6)) * ac
    va = d3 * d6 - d5 * d4
    if va <= 0 and (d4 - d3) >= 0 and (d5 - d6) >= 0:
        return b + ((d4 - d3) / ((d4 - d3) + (d5 - d6))) * (c - b)
    den = 1.0 / (va + vb + vc)
    return a + ab * (vb * den) + ac * (vc * den)


def pen_check(src, dst, name):
    """src 的顶面节点 vs dst 的底面三角（= 骨面）。返回 (min_d, n_pen, n_lt_0p05)。"""
    Ns, Nd = src["N"], dst["N"]
    tb = np.unique(dst["B"])
    tri = dst["B"]
    P = Nd[tb]
    Nrm = np.cross(Nd[tri[:, 1]] - Nd[tri[:, 0]], Nd[tri[:, 2]] - Nd[tri[:, 0]])
    ctr = Nd[tri].mean(1)
    # 法向定向：指向 dst 的顶面（关节腔）
    up = Nd[dst["T"]].mean(1) - Nd[tri].mean(1)
    sgn = np.sign(np.einsum("ij,ij->i", Nrm, up))
    Nrm = Nrm * sgn[:, None]
    Nrm /= (np.linalg.norm(Nrm, axis=1, keepdims=True) + 1e-15)
    # 采样点云加速（顶点 + 边中点 + 质心）
    samp = np.vstack([Nd[tri].reshape(-1, 3),
                      ((Nd[tri] + Nd[tri[:, [1, 2, 0]]]) / 2).reshape(-1, 3),
                      ctr])
    tree = cKDTree(samp)
    tops = np.where(src["off"])[0]
    ds = []
    for p in Ns[tops]:
        _, i = tree.query(p, k=1)
        e = int(i) % len(tri)
        q = closest_on_tri(p, Nd[tri[e, 0]], Nd[tri[e, 1]], Nd[tri[e, 2]])
        ds.append(float(Nrm[e] @ (p - q)))
    ds = np.array(ds)
    print(f"  {name}: 顶面节点 {len(ds)} 个  有符号距 {name.split('->')[1]} 骨面: "
          f"min {ds.min():+.4f}mm  中位 {np.median(ds):+.4f}   "
          f"穿透(<0) {int((ds < 0).sum())} 个   贴近(<0.05) {int((ds < 0.05).sum())} 个")
    return ds


def thickness(D):
    th = np.linalg.norm(D["N"][D["T"]].mean(1) - D["N"][D["B"]].mean(1), axis=1)
    return th


for tag, path in (("BASELINE cart_m1", ROOT / "temp/thums/cart_m1/cartilage_mesh.npz"),
                  ("FIXED   cart_m1_fix", ROOT / "temp/thums/cart_m1_fix/cartilage_mesh.npz")):
    print("=" * 92)
    print(tag, " ", path.name)
    D = load(path)
    print("① worst-det(J) 完整采样（105 点）")
    for j in ("tibiotalar_tibia", "tibiotalar_talus"):
        dt = dets(D[j]["N"][D[j]["C"]]).min(axis=1)
        print(f"   {j:20s} n={len(dt):4d}  min {dt.min():+9.4f}   ≤0 的 {int((dt<=0).sum()):3d}  "
              f"≤1e-2 的 {int((dt<=1e-2).sum()):3d}")
    print("③ 几何厚度（底面质心→顶面质心）")
    for j in ("tibiotalar_tibia", "tibiotalar_talus"):
        th = thickness(D[j])
        bad = int(((th < 0.3) | (th > 1.6)).sum())
        print(f"   {j:20s} min {th.min():.4f}  中位 {np.median(th):.4f}  max {th.max():.4f}  "
              f"越界(<0.3 或 >1.6) {bad} 个")
    print("② 有符号穿透检查（A 顶面 vs B 底面骨面）")
    pen_check(D["tibiotalar_tibia"], D["tibiotalar_talus"],
              "tibia_top->talus_bone")
    pen_check(D["tibiotalar_talus"], D["tibiotalar_tibia"],
              "talus_top->tibia_bone")
