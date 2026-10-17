"""CARTFIX 诊断 7：决定性判据 —— 逐顶点「抬升方向 ⇄ 底面 facet 法向」的夹角。

楔形要不开裂/不翻转，必要条件：3 个顶点的 disp 都投影到底面法向的**同一侧**。
若某顶点 u_i·n̂_facet ≤ 0 ⇒ 该角上"抬升方向切着面走" ⇒ 该角 detJ 翻负。
"""
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
CM1 = np.load(ROOT / "temp/thums/cart_m1/cartilage_mesh.npz", allow_pickle=True)

RS = (0.0, 1 / 6, 1 / 3, 1 / 2, 2 / 3, 1.0)
TS = (0.0, 0.25, 0.5, 0.75, 1.0)
GS = [(r, s, t) for r in RS for s in RS if r + s <= 1 + 1e-9 for t in TS]


def dets(P6):
    out = np.full(len(P6), np.inf)
    for (r, s, t) in GS:
        G = np.array([
            [-(1 - t), -(1 - t), -(1 - r - s)],
            [1 - t, 0.0, -r],
            [0.0, 1 - t, -s],
            [-t, -t, (1 - r - s)],
            [t, 0.0, r],
            [0.0, t, s],
        ])
        out = np.minimum(out, np.linalg.det(np.einsum("nia,ic->nac", P6, G)))
    return out


def ix(x):
    x = int(x)
    return (-x - 1) if x < 0 else (x - 1)


for j in ["tibiotalar_tibia", "tibiotalar_talus"]:
    N = np.asarray(CM1[f"{j}_nodes"], float)
    E = np.asarray(CM1[f"{j}_elems"], np.int64)
    C = np.array([[ix(x) for x in row] for row in E])
    B, T = C[:, :3], C[:, 3:]
    P6 = N[C]
    W = dets(P6)
    nb = np.cross(N[B[:, 1]] - N[B[:, 0]], N[B[:, 2]] - N[B[:, 0]])
    nbu = nb / (np.linalg.norm(nb, axis=1, keepdims=True) + 1e-15)
    d = N[T] - N[B]                                        # (n,3,3)
    dm = d.mean(axis=1)
    # 依 mean disp 定向法向
    sgn = np.sign(np.einsum("ij,ij->i", nbu, dm))
    # 逐顶点 disp 的投影（正=与顶面同侧）与夹角
    proj = np.einsum("nij,nj->ni", d, sgn[:, None] * nbu)
    ang = np.degrees(np.arccos(np.clip(np.einsum("nij,nj->ni", d / np.linalg.norm(d, axis=2, keepdims=True),
                                                 sgn[:, None] * nbu), -1, 1)))
    pmin = proj.min(axis=1)
    amin = ang.max(axis=1)
    print("=" * 84)
    print(f"{j}: elems {len(E)}  坏(≤1e-2) {int((W<=1e-2).sum())}")
    print(f"  逐顶点投影 min_i(d_i·n̂):  坏 中位 {np.median(pmin[W<=1e-2]):+.4f} "
          f"(min {pmin[W<=1e-2].min():+.4f})  |  好 中位 {np.median(pmin[W>1e-2]):+.4f} "
          f"(min {pmin[W>1e-2].min():+.4f})")
    print(f"  逐顶点夹角 max_i(u_i,n̂):  坏 中位 {np.median(amin[W<=1e-2]):.1f}° "
          f"(max {amin[W<=1e-2].max():.1f}°)  |  好 中位 {np.median(amin[W>1e-2]):.1f}° "
          f"(max {amin[W>1e-2].max():.1f}°)")
    for thr in (0.0, 0.2, 0.4, 0.6, 0.8):
        m = pmin < thr
        tot = int(m.sum())
        hit = int((W[m] <= 1e-2).sum())
        miss = int(((W <= 1e-2) & ~m).sum())
        fa = tot - hit
        print(f"    判据 min_i(d_i·n̂) < {thr:.1f}mm: 命中 {tot:>4} 个，真坏 {hit:>3}，"
              f"漏报 {miss:>3}，误报 {fa:>4}")
    print(f"  整体 min_i(d_i·n̂) 分位: " +
          " ".join(f"p{p}={np.percentile(pmin,p):+.3f}" for p in (1, 5, 10, 25, 50, 75, 100)))
    np.savez(ROOT / f"temp/pyfebio_demo/cf_p7_{j}.npz", W=W, pmin=pmin, amin=amin, dm=dm)
