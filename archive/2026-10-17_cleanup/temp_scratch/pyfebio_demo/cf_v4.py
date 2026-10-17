"""CARTFIX 快速验证 v4：球面约束 Adam —— **只转方向，厚度模长逐个节点严格保持**。

d_i = t_i · u_i， t_i ≡ |d_raw_i| 固定 ⇒ 厚度场（逐节点）完全不变。
优化 u_i ∈ S²，F = Σ_e Σ_{det<τ}(τ-det)² + λ_sm Σ|d_i-d_j|²。
"""
import os
import time
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
CM1 = np.load(ROOT / "temp/thums/cart_m1/cartilage_mesh.npz", allow_pickle=True)
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


def dets_of(P6):
    return np.linalg.det(np.einsum("nic,mia->nmca", P6, GG))


def ix(x):
    x = int(x)
    return (-x - 1) if x < 0 else (x - 1)


def solve(j, tau=0.6, lam_sm=0.001, iters=2000, lr=0.05, verbose=True):
    t0 = time.time()
    N = np.asarray(CM1[f"{j}_nodes"], float)
    E = np.asarray(CM1[f"{j}_elems"], np.int64)
    C = np.array([[ix(x) for x in row] for row in E])
    B, T = C[:, :3], C[:, 3:]
    bt = sorted(set(B.ravel().tolist()))
    bi = {b: i for i, b in enumerate(bt)}
    bidx = np.vectorize(bi.__getitem__)(B)
    d_raw = np.zeros((len(bt), 3))
    for a in range(3):
        d_raw[bidx[:, a]] = N[T[:, a]] - N[B[:, a]]
    tt = np.linalg.norm(d_raw, axis=1)                     # 逐节点厚度（固定）
    U = d_raw / tt[:, None]
    edges = np.array(sorted({(min(bidx[e, a], bidx[e, (a + 1) % 3]),
                             max(bidx[e, a], bidx[e, (a + 1) % 3]))
                             for e in range(len(E)) for a in range(3)}))

    def rebuild(Uu):
        D = tt[:, None] * Uu
        N2 = N.copy()
        for a in range(3):
            N2[T[:, a]] = N[B[:, a]] + D[bidx[:, a]]
        return N2, D

    def gradf(N2, dt, sel):
        g = np.zeros((len(bt), 3))
        for k, e in enumerate(sel):
            P6 = N2[C[e]]
            act = np.where(dt[k] < tau)[0]
            for mm in act:
                gm = GG[mm]
                J = P6.T @ gm
                det = float(np.linalg.det(J))
                Jinv = np.linalg.inv(J + (1e-9 if det <= 1e-12 else 0.0) * np.eye(3))
                w = -2.0 * (tau - det) * det
                for i in range(3, 6):
                    g[bidx[e, i - 3]] += w * (gm[i] @ Jinv)
        return g

    m = np.zeros_like(U); v = np.zeros_like(U)
    best = (None, None)
    for it in range(iters):
        N2, D = rebuild(U)
        dt = dets_of(N2[C])
        F = float(np.sum(np.maximum(0.0, tau - dt) ** 2))
        if lam_sm:
            dd = D[edges[:, 0]] - D[edges[:, 1]]
            F += lam_sm * float(np.sum(dd ** 2))
        nneg = int((dt.min(axis=1) <= 0).sum())
        cur = (nneg, F)
        if best[0] is None or cur < best[0]:
            best = (cur, U.copy())
        if verbose and it % 400 == 0:
            print(f"     it={it:>5} F={F:9.4f} min {dt.min():+9.4f} "
                  f"bad {int((dt.min(axis=1)<=1e-2).sum()):>3} neg {nneg:>3} "
                  f"({time.time()-t0:.1f}s)", flush=True)
        sel = np.where(dt.min(axis=1) < tau)[0]
        if not len(sel):
            if verbose:
                print(f"     it={it}: 全部达标，提前结束", flush=True)
            break
        gd = gradf(N2, dt, sel)                            # ∂F/∂d
        if lam_sm:
            np.add.at(gd, edges[:, 0], 2 * lam_sm * (D[edges[:, 0]] - D[edges[:, 1]]))
            np.add.at(gd, edges[:, 1], 2 * lam_sm * (D[edges[:, 1]] - D[edges[:, 0]]))
        gu = tt[:, None] * gd                              # ∂F/∂u = t ∂F/∂d
        gu = gu - np.sum(gu * U, axis=1, keepdims=True) * U   # 切空间投影
        m = 0.9 * m + 0.1 * gu
        v = 0.999 * v + 0.001 * gu ** 2
        mh = m / (1 - 0.9 ** (it + 1))
        vh = v / (1 - 0.999 ** (it + 1))
        U = U - lr * mh / (np.sqrt(vh) + 1e-12)
        U = U / (np.linalg.norm(U, axis=1, keepdims=True) + 1e-15)
    U = best[1] if best[1] is not None else U
    N2, D = rebuild(U)
    dt = dets_of(N2[C])
    mins = dt.min(axis=1)
    th = np.linalg.norm(N2[T].mean(1) - N2[B].mean(1), axis=1)
    mv = float(np.max(np.linalg.norm(N2[T] - N[T], axis=2)))
    print(f"  ► {j}: bad(≤1e-2) {int((mins<=1e-2).sum())} neg {int((mins<=0).sum())} "
          f"min detJ {mins.min():+.4f}", flush=True)
    print(f"    单元厚度 [min {th.min():.4f} 中位 {np.median(th):.4f} max {th.max():.4f}]  "
          f"节点厚度不变: {np.allclose(np.linalg.norm(D,axis=1), tt)}  "
          f"顶面位移 {mv:.3f}mm  用时 {time.time()-t0:.1f}s", flush=True)
    np.savez(ROOT / f"temp/pyfebio_demo/cf_v4_{j}.npz", U=U, tt=tt, bt=np.array(bt))
    return U, mins, th, mv, N2, dict(N=N, T=T, B=B, bidx=bidx, tt=tt)


for j in ["tibiotalar_tibia", "tibiotalar_talus"]:
    for tau, lr in ((0.6, 0.05), (1.2, 0.08)):
        print("=" * 70, flush=True)
        print(f"{j}  tau={tau} lr={lr}", flush=True)
        solve(j, tau=tau, lr=lr, iters=2500)
