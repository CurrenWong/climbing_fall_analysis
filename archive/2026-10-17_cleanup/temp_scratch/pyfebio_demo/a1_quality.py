"""A1：subtalar 两侧软骨的网格质量对比（只测）——定位为何 calcaneus 侧一接就翻。
指标：厚度（底→顶）、完整采样 det(J)、面内尺寸/厚度比、剔除单元的空间聚集。"""
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
cg = np.load(ROOT / "temp/thums/cart_m1/cartilage_mesh.npz", allow_pickle=True)
PENTA_FACES = [(0, 1, 2), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)]


def cnidx(x):
    x = int(x)
    return (-x - 1) if x < 0 else (x - 1)


def penta6_dN(r, s, t):
    L = (1 - r - s, r, s)
    dL = ((-1.0, -1.0), (1.0, 0.0), (0.0, 1.0))
    G = np.zeros((6, 3))
    for i in range(3):
        G[i] = ((1 - t) / 2 * dL[i][0], (1 - t) / 2 * dL[i][1], -L[i] / 2)
        G[i + 3] = ((1 + t) / 2 * dL[i][0], (1 + t) / 2 * dL[i][1], L[i] / 2)
    return G


def detJ(X, idx, full):
    P = X[list(idx)]
    RS = (0.0, 1 / 6, 1 / 3, 1 / 2, 2 / 3) + ((1.0,) if full else ())
    TS = (0.0, 0.25, 0.5, 0.75, 1.0)
    w = None
    for r in RS:
        for s in RS:
            if r + s > 1.0 + 1e-9:
                continue
            for t in TS:
                d = float(np.linalg.det(P.T @ penta6_dN(r, s, t)))
                w = d if w is None else min(w, d)
    return w


def sv(P):
    dt = (P[3] + P[4] + P[5] - P[0] - P[1] - P[2]) / 3.0
    return float(np.linalg.det(np.stack([P[1] - P[0], P[2] - P[0], dt])))


for j in ("subtalar_talus", "subtalar_calcaneus"):
    X = np.asarray(cg[f"{j}_nodes"], float)
    isoff = np.asarray(cg[f"{j}_is_offset"], bool)
    alive, dead = [], []
    for row in cg[f"{j}_elems"]:
        raw = [int(x) for x in row if int(x) != 0]
        if len(raw) != 6:
            continue
        idx = [cnidx(x) for x in raw]
        ok = sv(X[idx]) > 0 and detJ(X, idx, False) > 0
        if not ok:
            alt = idx[3:] + idx[:3]
            ok = sv(X[alt]) > 0 and detJ(X, alt, False) > 0
            idx = alt if ok else idx
        (alive if ok else dead).append(idx)

    def stats(idxs, tag):
        th, fdet, ar, cen = [], [], [], []
        for idx in idxs:
            P = X[idx]
            # 厚度 = 顶面质心 - 底面质心 的距离（不管朝向）
            th.append(np.linalg.norm(P[3:].mean(0) - P[:3].mean(0)))
            fdet.append(detJ(X, idx, True))
            # 面内特征边长 / 厚度
            base = P[:3]
            e = [np.linalg.norm(base[(i + 1) % 3] - base[i]) for i in range(3)]
            ar.append(max(e) / max(th[-1], 1e-12))
            cen.append(P.mean(0))
        th, fdet, ar = map(np.array, (th, fdet, ar))
        print(f"  [{tag}] n={len(idxs)}  厚度 min/med/max = "
              f"{th.min():.4f}/{np.median(th):.4f}/{th.max():.4f} mm")
        print(f"        完整采样 det(J): min {fdet.min():.3e} | ≤0 {int((fdet<=0).sum())}"
              f" | 0<x≤1e-2 {int(((fdet>0)&(fdet<=1e-2)).sum())}"
              f" | 1e-2<x≤1e-1 {int(((fdet>1e-2)&(fdet<=1e-1)).sum())}")
        print(f"        长细比 max/厚 min/med = {ar.max():.1f} / {np.median(ar):.1f}")
        return np.array(cen), fdet

    ca, da = stats(alive, f"{j} 存活")
    if dead:
        cd, dd = stats(dead, f"{j} 剔除")
        # 剔除单元是否聚集：最近邻（与存活单元质心）
        from scipy.spatial import cKDTree
        t = cKDTree(ca)
        d, _ = t.query(cd)
        print(f"        剔除单元到最近存活单元的距离: med {np.median(d):.2f} "
              f"max {d.max():.2f} mm（大 ⇒ 成簇/成洞）")
    # 最差 5 个存活单元的位置
    order = np.argsort([detJ(X, i, True) for i in alive])[:5]
    print("        最差存活单元质心 (x,y,z) / det_full:")
    for k in order:
        print(f"          {np.round(ca[k],1)}  {detJ(X, alive[k], True):.3e}")
