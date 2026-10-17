"""CARTFIX 实验 4：每个底面节点，是否存在使「所有相邻 facet」都正投影的方向？

对每个节点做球面暴力搜索（Fibonacci 球 ~4000 点）：
   f(u) = min_{e∋i}  (u · n̂_e^U)
若 max_u f(u) ≥ 0 则可行。
同时报告「可行节点里，把方向取成 argmax f」后的 detJ（含 6 个"正投影也坏"的单元）。
"""
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
CM1 = np.load(ROOT / "temp/thums/cart_m1/cartilage_mesh.npz", allow_pickle=True)
RS = (0.0, 1 / 6, 1 / 3, 1 / 2, 2 / 3, 1.0)
TS = (0.0, 0.25, 0.5, 0.75, 1.0)
GS = [(r, s, t) for r in RS for s in RS if r + s <= 1 + 1e-9 for t in TS]
NSPH = 6000


def fib_sphere(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    th = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], axis=1)


SPH = fib_sphere(NSPH)


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
    bt = sorted(set(B.ravel().tolist()))
    d0 = {}
    for r in range(len(E)):
        for a in range(3):
            d0[int(B[r, a])] = N[int(T[r, a])] - N[int(B[r, a])]
    t_i = {b: np.linalg.norm(d0[b]) for b in bt}
    U = np.sum([d0[b] for b in bt], axis=0)
    U = U / np.linalg.norm(U)
    fac = {}
    for r in range(len(E)):
        P = N[B[r]]
        cr = np.cross(P[1] - P[0], P[2] - P[0])
        nv = cr / (np.linalg.norm(cr) + 1e-15)
        if np.dot(nv, U) < 0:
            nv = -nv
        for q in B[r]:
            fac.setdefault(int(q), []).append(nv)
    print("=" * 90)
    print(f"{j}:  底节点 {len(bt)}  U={np.round(U,3)}")
    best_u, fmax = {}, {}
    for b in bt:
        M = np.array(fac[b])                       # (k,3)
        vals = SPH @ M.T                           # (NSPH, k)
        f = vals.min(axis=1)
        k = int(np.argmax(f))
        best_u[b] = SPH[k]
        fmax[b] = float(f[k])
    fm = np.array([fmax[b] for b in bt])
    print(f"   max_u min_e(u·n̂) 分布: min {fm.min():+.4f}  中位 {np.median(fm):+.4f}  "
          f"<0 的节点 {int((fm < 0).sum())}  <0.02 的 {int((fm < 0.02).sum())}")
    worst = sorted(zip(bt, fm), key=lambda z: z[1])[:10]
    print("   最差 10 个节点 (node, max-min投影, 相邻facet数):")
    for b, v in worst:
        print(f"      {b:>4}  {v:+.4f}   k={len(fac[b])}")
    # 用最优方向重建（保持模长）
    N2 = N.copy()
    for r in range(len(E)):
        for a in range(3):
            b = int(B[r, a])
            N2[int(T[r, a])] = N[b] + t_i[b] * best_u[b]
    W = dets(N2[C])
    th = np.linalg.norm(N2[T].mean(1) - N2[B].mean(1), axis=1)
    print(f"   用 argmax 方向: bad {int((W<=1e-2).sum())} neg {int((W<=0).sum())} "
          f"min detJ {W.min():+.4f}  厚度 [{th.min():.4f},{th.max():.4f}] "
          f"位移 {max(np.linalg.norm(N2[T]-N[T],axis=2).ravel()):.3f}mm")
    # 再叠一层：对仍坏的单元做「顶节点梯度爬山」（仅移动其 3 个顶节点）
    for it in range(300):
        W = dets(N2[C])
        badi = np.where(W <= 1e-2)[0]
        if not len(badi):
            break
        grad = np.zeros_like(N2)
        for e in badi:
            idxs = C[e]
            P6 = N2[idxs]
            r, s, t = min(GS, key=lambda g: _dk(P6, g))
            G = np.array([
                [-(1 - t), -(1 - t), -(1 - r - s)],
                [1 - t, 0.0, -r],
                [0.0, 1 - t, -s],
                [-t, -t, (1 - r - s)],
                [t, 0.0, r],
                [0.0, t, s],
            ])
            J = P6.T @ G
            det = np.linalg.det(J)
            Jinv = np.linalg.inv(J)
            for i in range(6):
                g = det * (G[i] @ Jinv)          # d(det)/dP_i
                grad[idxs[i]] += g
        step = 0.02
        nrm = np.linalg.norm(grad, axis=1, keepdims=True)
        N2 = N2 + step * grad / (nrm + 1e-12) * (nrm > 0)
        # 约束：只动顶节点；模长投影到 [0.3,1.6]
        for r in range(len(E)):
            for a in range(3):
                b = int(B[r, a])
                v = N2[int(T[r, a])] - N[b]
                L = np.linalg.norm(v)
                if L > 1.6:
                    N2[int(T[r, a])] = N[b] + v * (1.6 / L)
                elif L < 0.3:
                    N2[int(T[r, a])] = N[b] + v * (0.3 / max(L, 1e-9))
    W = dets(N2[C])
    th = np.linalg.norm(N2[T].mean(1) - N2[B].mean(1), axis=1)
    print(f"   + 顶节点梯度爬山(≤300步): bad {int((W<=1e-2).sum())} neg {int((W<=0).sum())} "
          f"min detJ {W.min():+.4f}  厚度 [{th.min():.4f},{th.max():.4f}] "
          f"位移 {max(np.linalg.norm(N2[T]-N[T],axis=2).ravel()):.3f}mm")
