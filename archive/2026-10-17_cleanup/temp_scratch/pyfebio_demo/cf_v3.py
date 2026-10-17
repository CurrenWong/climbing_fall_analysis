"""CARTFIX 快速验证 v3：Adam 解缠（无需线搜索）+ 分阶段（先大 tau 后收敛）。

阶段 A: tau=0.6, lam_reg=0.002, lam_sm=0.001, 2000 it
阶段 B: tau=0.25（只求脱离负值）, lam_reg=0.0005, 800 it
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


def build(j):
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
    edges = np.array(sorted({(min(bidx[e, a], bidx[e, (a + 1) % 3]),
                             max(bidx[e, a], bidx[e, (a + 1) % 3]))
                             for e in range(len(E)) for a in range(3)}))
    return dict(N=N, E=E, C=C, B=B, T=T, bt=bt, bidx=bidx, d_raw=d_raw, edges=edges)


def solve(j, stages, verbose=True):
    P = build(j)
    N, E, C, B, T, bidx, d_raw, edges = (P["N"], P["E"], P["C"], P["B"], P["T"],
                                         P["bidx"], P["d_raw"], P["edges"])
    t0 = time.time()

    def rebuild(D):
        N2 = N.copy()
        for a in range(3):
            N2[T[:, a]] = N[B[:, a]] + D[bidx[:, a]]
        return N2

    # 顶点 -> 其所有相邻单元（底面顶点行号）
    nb_of_e = [set(bidx[e].tolist()) for e in range(len(E))]
    node2e = {}
    for e in range(len(E)):
        for b in nb_of_e[e]:
            node2e.setdefault(b, []).append(e)

    def gradf(D, N2, dt, tau, sel):
        g = np.zeros_like(D)
        for k, e in enumerate(sel):
            P6 = N2[C[e]]
            act = np.where(dt[k] < tau)[0]
            if not len(act):
                continue
            for mm in act:
                gm = GG[mm]
                J = P6.T @ gm
                det = float(np.linalg.det(J))
                Jinv = np.linalg.inv(J + (1e-9 if det <= 1e-12 else 0.0) * np.eye(3))
                w = -2.0 * (tau - det) * det
                for i in range(3, 6):
                    g[bidx[e, i - 3]] += w * (gm[i] @ Jinv)
        return g

    D = d_raw.copy()
    for st in stages:
        tau, lam_reg, lam_sm, iters, lr = (st["tau"], st["lam_reg"], st["lam_sm"],
                                           st["iters"], st["lr"])
        m = np.zeros_like(D)
        v = np.zeros_like(D)
        best = (1e30, None, None)
        allsel = np.arange(len(E))
        for it in range(iters):
            N2 = rebuild(D)
            dt = dets_of(N2[C])
            F = float(np.sum(np.maximum(0.0, tau - dt) ** 2))
            F += lam_reg * float(np.sum((D - d_raw) ** 2))
            if lam_sm:
                dd = D[edges[:, 0]] - D[edges[:, 1]]
                F += lam_sm * float(np.sum(dd ** 2))
            nneg = int((dt.min(axis=1) <= 0).sum())
            bad = int((dt.min(axis=1) <= 1e-2).sum())
            cur = (nneg, F)
            if best[1] is None or cur < best[0]:
                best = (cur, D.copy(), dt.min(axis=1).min())
            if it % 400 == 0 and verbose:
                print(f"     it={it:>5} F={F:9.4f} min {dt.min():+9.4f} "
                      f"bad {bad:>3} neg {nneg:>3}  ({time.time()-t0:.1f}s)", flush=True)
            g = gradf(D, N2, dt, tau, allsel)
            g += 2 * lam_reg * (D - d_raw)
            if lam_sm:
                np.add.at(g, edges[:, 0], 2 * lam_sm * (D[edges[:, 0]] - D[edges[:, 1]]))
                np.add.at(g, edges[:, 1], 2 * lam_sm * (D[edges[:, 1]] - D[edges[:, 0]]))
            m = 0.9 * m + 0.1 * g
            v = 0.999 * v + 0.001 * g ** 2
            mh = m / (1 - 0.9 ** (it + 1))
            vh = v / (1 - 0.999 ** (it + 1))
            D = D - lr * mh / (np.sqrt(vh) + 1e-12)
            L = np.linalg.norm(D, axis=1)
            D = D * (np.clip(L, 0.30, 1.60) / np.maximum(L, 1e-12))[:, None]
        if best[1] is not None:
            D = best[1]
        if verbose:
            print(f"   阶段(tau={tau}) 结束: 最佳 neg {best[0][0]} min {best[2]:+.4f}", flush=True)
    # 最终只求「全正」的抛光阶段：仅当仍有 ≤1e-2 的单元时，用逐单元块优化
    N2 = rebuild(D)
    dt = dets_of(N2[C])
    mins = dt.min(axis=1)
    for _ in range(40):
        bad = np.where(mins <= 1e-2)[0]
        if not len(bad):
            break
        for e in bad:
            sel = np.array(sorted({x for b in bidx[e] for x in node2e[int(b)]}))
            # 仅对这些节点的相邻单元做 20 步 Adam，目标 tau=0.15，无 reg
            m2 = np.zeros_like(D); v2 = np.zeros_like(D)
            for k in range(40):
                N2 = rebuild(D)
                d2 = dets_of(N2[C[sel]])
                g = gradf(D, N2, d2, 0.15, sel)
                m2 = 0.9 * m2 + 0.1 * g
                v2 = 0.999 * v2 + 0.001 * g ** 2
                mh = m2 / (1 - 0.9 ** (k + 1)); vh = v2 / (1 - 0.999 ** (k + 1))
                D = D - 0.004 * mh / (np.sqrt(vh) + 1e-12)
                L = np.linalg.norm(D, axis=1)
                D = D * (np.clip(L, 0.30, 1.60) / np.maximum(L, 1e-12))[:, None]
        N2 = rebuild(D)
        dt = dets_of(N2[C])
        mins = dt.min(axis=1)
    th = np.linalg.norm(N2[T].mean(1) - N2[B].mean(1), axis=1)
    mv = float(np.max(np.linalg.norm(N2[T] - N[T], axis=2)))
    print(f"  ► {j}: bad(≤1e-2) {int((mins<=1e-2).sum())} neg {int((mins<=0).sum())} "
          f"min detJ {mins.min():+.4f}", flush=True)
    print(f"    厚度 [{th.min():.4f}, {np.median(th):.4f}, {th.max():.4f}] "
          f"顶面位移 {mv:.3f}mm  用时 {time.time()-t0:.1f}s", flush=True)
    np.savez(ROOT / f"temp/pyfebio_demo/cf_v3_{j}.npz", D=D, bt=np.array(P["bt"]))
    return D, mins, th, mv, N2, P


STAGES = [
    dict(tau=0.6, lam_reg=0.0020, lam_sm=0.0010, iters=1500, lr=0.010),
    dict(tau=0.30, lam_reg=0.0008, lam_sm=0.0005, iters=800, lr=0.006),
]
if __name__ == "__main__":
    for j in ["tibiotalar_tibia", "tibiotalar_talus"]:
        print("=" * 70, flush=True)
        print(j, flush=True)
        solve(j, STAGES)
