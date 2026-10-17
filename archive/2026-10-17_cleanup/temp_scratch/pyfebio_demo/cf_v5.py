"""CARTFIX 快速验证 v5：v3 的稳定 Adam + 「保模长」正则。

F = Σ_e Σ_{det<τ}(τ-det)² + λ_t Σ_i (|d_i|-t_i^raw)² + λ_sm Σ_{i~j}|d_i-d_j|²
λ_t 越大 ⇒ 逐节点厚度越接近原值。
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


def solve(j, tau=0.6, lam_t=0.05, lam_sm=0.001, iters=1400, lr=0.010,
          tau2=0.30, iters2=400, lr2=0.006, verbose=True):
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
    tt = np.linalg.norm(d_raw, axis=1)
    edges = np.array(sorted({(min(bidx[e, a], bidx[e, (a + 1) % 3]),
                             max(bidx[e, a], bidx[e, (a + 1) % 3]))
                             for e in range(len(E)) for a in range(3)}))

    def rebuild(D):
        N2 = N.copy()
        for a in range(3):
            N2[T[:, a]] = N[B[:, a]] + D[bidx[:, a]]
        return N2

    def gradf(N2, dt, tau_, sel):
        g = np.zeros_like(D)
        for k, e in enumerate(sel):
            P6 = N2[C[e]]
            act = np.where(dt[k] < tau_)[0]
            for mm in act:
                gm = GG[mm]
                J = P6.T @ gm
                det = float(np.linalg.det(J))
                Jinv = np.linalg.inv(J + (1e-9 if det <= 1e-12 else 0.0) * np.eye(3))
                w = -2.0 * (tau_ - det) * det
                for i in range(3, 6):
                    g[bidx[e, i - 3]] += w * (gm[i] @ Jinv)
        return g

    D = d_raw.copy()
    m = np.zeros_like(D); v = np.zeros_like(D)
    best = (None, None, None)
    for st, (tau_, it_n, lr_) in enumerate(((tau, iters, lr), (tau2, iters2, lr2))):
        for it in range(it_n):
            N2 = rebuild(D)
            dt = dets_of(N2[C])
            F = float(np.sum(np.maximum(0.0, tau_ - dt) ** 2))
            F += lam_t * float(np.sum((np.linalg.norm(D, axis=1) - tt) ** 2))
            if lam_sm:
                dd = D[edges[:, 0]] - D[edges[:, 1]]
                F += lam_sm * float(np.sum(dd ** 2))
            mins = dt.min(axis=1)
            cur = (int((mins <= 0).sum()), F)
            if best[0] is None or cur < best[0]:
                best = (cur, D.copy(), float(mins.min()))
            if verbose and it % 400 == 0:
                print(f"     st{st} it={it:>5} F={F:9.4f} min {mins.min():+9.4f} "
                      f"bad {int((mins<=1e-2).sum()):>3} neg {cur[0]:>3} "
                      f"({time.time()-t0:.1f}s)", flush=True)
            sel = np.where(mins < tau_)[0]
            if not len(sel):
                if verbose:
                    print(f"     st{st} it={it}: 达标，提前结束", flush=True)
                break
            g = gradf(N2, dt, tau_, sel)
            L = np.linalg.norm(D, axis=1)
            u = D / np.maximum(L, 1e-12)[:, None]
            Fr = np.sum((L - tt) ** 2)
            g += 2 * lam_t * (L - tt)[:, None] * u
            if lam_sm:
                np.add.at(g, edges[:, 0], 2 * lam_sm * (D[edges[:, 0]] - D[edges[:, 1]]))
                np.add.at(g, edges[:, 1], 2 * lam_sm * (D[edges[:, 1]] - D[edges[:, 0]]))
            m = 0.9 * m + 0.1 * g
            v = 0.999 * v + 0.001 * g ** 2
            mh = m / (1 - 0.9 ** (it + 1))
            vh = v / (1 - 0.999 ** (it + 1))
            D = D - lr_ * mh / (np.sqrt(vh) + 1e-12)
            L = np.linalg.norm(D, axis=1)
            D = D * (np.clip(L, 0.30, 1.60) / np.maximum(L, 1e-12))[:, None]
        if best[1] is not None:
            D = best[1]
    D = best[1] if best[1] is not None else D
    N2 = rebuild(D)
    dt = dets_of(N2[C])
    mins = dt.min(axis=1)
    th = np.linalg.norm(N2[T].mean(1) - N2[B].mean(1), axis=1)
    nodeL = np.linalg.norm(D, axis=1)
    mv = float(np.max(np.linalg.norm(N2[T] - N[T], axis=2)))
    print(f"  ► {j} (λ_t={lam_t}): bad(≤1e-2) {int((mins<=1e-2).sum())} "
          f"neg {int((mins<=0).sum())} min detJ {mins.min():+.4f}", flush=True)
    print(f"    单元厚度 [{th.min():.4f}, {np.median(th):.4f}, {th.max():.4f}]  "
          f"节点厚度 偏差 max|Δ|={np.abs(nodeL-tt).max():.4f}mm  顶面位移 {mv:.3f}mm  "
          f"用时 {time.time()-t0:.1f}s", flush=True)
    return (D, mins, th, mv, N2)


if __name__ == "__main__":
    for lam_t in (0.05, 0.5):
        for j in ["tibiotalar_tibia", "tibiotalar_talus"]:
            print("=" * 70, flush=True)
            print(f"{j}  lam_t={lam_t}", flush=True)
            solve(j, lam_t=lam_t)
