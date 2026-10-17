"""CARTFIX 快速验证 v6：球面（严格保模长）+ 回溯线搜索 + 分阶段 tau。

d_i = t_i·u_i，t_i ≡ |d_raw_i| 固定（逐节点厚度严格不变）。
阶段 A: tau=0.5, 1500 it;  阶段 B: tau=0.18, 1200 it（只求脱离负值）
"""
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


def solve(j, stages, lam_sm=0.001, verbose=True):
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

    def F_of(N2, D, dt, tau_):
        F = float(np.sum(np.maximum(0.0, tau_ - dt) ** 2))
        if lam_sm:
            dd = D[edges[:, 0]] - D[edges[:, 1]]
            F += lam_sm * float(np.sum(dd ** 2))
        return F

    def gradf(N2, dt, tau_, sel):
        g = np.zeros((len(bt), 3))
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

    best = (None, None, None)
    for st, (tau_, it_n, step0) in enumerate(stages):
        step = step0
        stall = 0
        for it in range(it_n):
            N2, D = rebuild(U)
            dt = dets_of(N2[C])
            mins = dt.min(axis=1)
            F = F_of(N2, D, dt, tau_)
            cur = (int((mins <= 0).sum()), F)
            if best[0] is None or cur < best[0]:
                best = (cur, U.copy(), float(mins.min()))
            if verbose and it % 300 == 0:
                print(f"     st{st} it={it:>5} τ={tau_} F={F:10.4f} min {mins.min():+9.4f} "
                      f"bad {int((mins<=1e-2).sum()):>3} neg {cur[0]:>3} step={step:.2e} "
                      f"({time.time()-t0:.1f}s)", flush=True)
            sel = np.where(mins < tau_)[0]
            if not len(sel):
                if verbose:
                    print(f"     st{st} it={it}: 达标", flush=True)
                break
            gd = gradf(N2, dt, tau_, sel)
            if lam_sm:
                np.add.at(gd, edges[:, 0], 2 * lam_sm * (D[edges[:, 0]] - D[edges[:, 1]]))
                np.add.at(gd, edges[:, 1], 2 * lam_sm * (D[edges[:, 1]] - D[edges[:, 0]]))
            gu = tt[:, None] * gd
            gu = gu - np.sum(gu * U, axis=1, keepdims=True) * U     # 切空间
            gmax = np.abs(gu).max()
            if gmax < 1e-14:
                break
            ok = False
            for al in (1.0, 0.5, 0.25, 0.1, 0.03, 0.01, 0.003):
                Un = U - al * step * gu / gmax
                Un = Un / (np.linalg.norm(Un, axis=1, keepdims=True) + 1e-15)
                N2n, Dn = rebuild(Un)
                dtn = dets_of(N2n[C])
                Fn = F_of(N2n, Dn, dtn, tau_)
                if Fn < F - 1e-15:
                    U, step = Un, min(step, 0.05)
                    ok = True
                    break
            if not ok:
                step *= 0.5
                stall += 1
                if stall > 12 or step < 1e-7:
                    if verbose:
                        print(f"     st{st} it={it}: 线搜索停滞（step→{step:.2e}），"
                              f"进入下一阶段", flush=True)
                    break
            else:
                stall = 0
        if best[1] is not None:
            U = best[1]
    U = best[1] if best[1] is not None else U
    N2, D = rebuild(U)
    dt = dets_of(N2[C])
    mins = dt.min(axis=1)
    th = np.linalg.norm(N2[T].mean(1) - N2[B].mean(1), axis=1)
    mv = float(np.max(np.linalg.norm(N2[T] - N[T], axis=2)))
    print(f"  ► {j}: bad(≤1e-2) {int((mins<=1e-2).sum())} neg {int((mins<=0).sum())} "
          f"min detJ {mins.min():+.4f}", flush=True)
    print(f"    单元厚度 [{th.min():.4f}, {np.median(th):.4f}, {th.max():.4f}]  "
          f"逐节点厚度严格不变: {np.allclose(np.linalg.norm(D,axis=1), tt)}  "
          f"顶面位移 {mv:.3f}mm  用时 {time.time()-t0:.1f}s", flush=True)
    np.savez(ROOT / f"temp/pyfebio_demo/cf_v6_{j}.npz", U=U, tt=tt, bt=np.array(bt))
    return U, mins, th, mv, N2


STAGES = [(0.5, 1500, 0.06), (0.18, 1200, 0.03)]
if __name__ == "__main__":
    for j in ["tibiotalar_tibia", "tibiotalar_talus"]:
        print("=" * 70, flush=True)
        print(j, flush=True)
        solve(j, STAGES)
