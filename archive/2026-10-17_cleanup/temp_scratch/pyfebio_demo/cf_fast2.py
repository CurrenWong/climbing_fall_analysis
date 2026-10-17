"""CARTFIX 快速验证 v2：子集目标（只算选中单元的 detJ）+ 时间预算。

全局解缠 ⇄ 残留单元局部块修复（交替）。
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
TAU = float(os.environ.get("CF_TAU", "1.0"))
LAM_REG = float(os.environ.get("CF_LAMREG", "0.002"))
LAM_SM = float(os.environ.get("CF_LAMSM", "0.001"))
BUDGET = float(os.environ.get("CF_BUDGET", "420"))       # 每个 joint 的秒数预算


def dets_of(P6):
    return np.linalg.det(np.einsum("nic,mia->nmca", P6, GG))


def ix(x):
    x = int(x)
    return (-x - 1) if x < 0 else (x - 1)


def run(j):
    t0 = time.time()
    N = np.asarray(CM1[f"{j}_nodes"], float)
    E = np.asarray(CM1[f"{j}_elems"], np.int64)
    C = np.array([[ix(x) for x in row] for row in E])
    B, T = C[:, :3], C[:, 3:]
    bt = sorted(set(B.ravel().tolist()))
    bi = {b: i for i, b in enumerate(bt)}
    bidx = np.vectorize(bi.__getitem__)(B)                       # (nE,3) 顶点行号
    occ = {b: [] for b in bt}
    for e in range(len(E)):
        for a in range(3):
            occ[int(B[e, a])].append((e, a))
    d_raw = np.zeros((len(bt), 3))
    for b in bt:
        e, a = occ[b][0]
        d_raw[bi[b]] = N[int(T[e, a])] - N[b]
    edges = np.array(sorted({(min(bidx[e, a], bidx[e, (a + 1) % 3]),
                             max(bidx[e, a], bidx[e, (a + 1) % 3]))
                             for e in range(len(E)) for a in range(3)}))
    e2nodes = bidx                                          # (nE,3)
    nb_of_e = [set(bidx[e].tolist()) for e in range(len(E))]
    node2e = {b: [e for e in range(len(E)) if b in nb_of_e[e]] for b in range(len(bt))}

    def rebuild(D):
        N2 = N.copy()
        for a in range(3):
            N2[T[:, a]] = N[B[:, a]] + D[bidx[:, a]]
        return N2

    def obj(D, sel):
        N2 = rebuild(D)
        dt = dets_of(N2[C[sel]])
        F = float(np.sum(np.maximum(0.0, TAU - dt) ** 2))
        F += LAM_REG * float(np.sum((D - d_raw) ** 2))
        if LAM_SM:
            dd = D[edges[:, 0]] - D[edges[:, 1]]
            F += LAM_SM * float(np.sum(dd ** 2))
        return F, dt, N2

    def gvec(D, N2, dt, sel):
        g = np.zeros_like(D)
        for k, e in enumerate(sel):
            P6 = N2[C[e]]
            m = np.where(dt[k] < TAU)[0]
            if not len(m):
                continue
            for mm in m:
                gm = GG[mm]
                J = P6.T @ gm
                det = float(np.linalg.det(J))
                Jinv = np.linalg.inv(J + (1e-9 if det <= 1e-12 else 0.0) * np.eye(3))
                w = -2.0 * (TAU - det) * det
                for i in range(3, 6):
                    g[bidx[e, i - 3]] += w * (gm[i] @ Jinv)
        g += 2 * LAM_REG * (D - d_raw)
        if LAM_SM:
            np.add.at(g, edges[:, 0], 2 * LAM_SM * (D[edges[:, 0]] - D[edges[:, 1]]))
            np.add.at(g, edges[:, 1], 2 * LAM_SM * (D[edges[:, 1]] - D[edges[:, 0]]))
        return g

    def descend(D, sel, iters, scale, tcap):
        sel = np.asarray(sel)
        F, dt, N2 = obj(D, sel)
        for _ in range(iters):
            if time.time() > tcap:
                break
            act = np.where(dt.min(axis=1) < TAU)[0]
            if not len(act):
                break
            g = gvec(D, N2, dt, sel[act])
            gmax = np.abs(g).max()
            if gmax < 1e-13:
                break
            for al in (1.0, 0.5, 0.25, 0.1, 0.03, 0.01, 0.003):
                Dn = D - al * scale * g / gmax
                L = np.linalg.norm(Dn, axis=1)
                Dn = Dn * (np.clip(L, 0.30, 1.60) / np.maximum(L, 1e-12))[:, None]
                Fn, dtn, N2n = obj(Dn, sel)
                if Fn < F - 1e-15:
                    D, F, dt, N2 = Dn, Fn, dtn, N2n
                    break
            else:
                break
        return D, dt, N2

    tcap = t0 + BUDGET
    allsel = np.arange(len(E))
    D = d_raw.copy()
    for rnd in range(8):
        D, dt, N2 = descend(D, allsel, 1500, 0.03, tcap)
        mins = dt.min(axis=1)
        bad = np.where(mins < TAU)[0]
        print(f"   round {rnd}: bad(≤1e-2) {int((mins<=1e-2).sum())} "
              f"min {mins.min():+.4f}  bad(≤0) {int((mins<=0).sum())} "
              f"({time.time()-t0:.1f}s)", flush=True)
        if not len(bad):
            break
        for e in bad:
            fwd = sorted({x for b in bidx[e] for x in node2e[int(b)]})
            D, dt, N2 = descend(D, np.array(fwd), 120, 0.06, tcap)
        if time.time() > tcap:
            print("   [预算耗尽]", flush=True)
            break
    dt_all = dets_of(N2[C])
    mins = dt_all.min(axis=1)
    th = np.linalg.norm(N2[T].mean(1) - N2[B].mean(1), axis=1)
    mv = float(np.max(np.linalg.norm(N2[T] - N[T], axis=2)))
    print(f"  ► {j}: bad(≤1e-2) {int((mins<=1e-2).sum())}  neg {int((mins<=0).sum())}  "
          f"min detJ {mins.min():+.4f}", flush=True)
    print(f"    厚度 [{th.min():.4f}, {np.median(th):.4f}, {th.max():.4f}]  "
          f"顶面位移 {mv:.3f}mm  用时 {time.time()-t0:.1f}s", flush=True)
    np.savez(ROOT / f"temp/pyfebio_demo/cf_fast2_{j}.npz", D=D, bt=np.array(bt))


for j in ["tibiotalar_tibia", "tibiotalar_talus"]:
    print("=" * 70, flush=True)
    print(j, flush=True)
    run(j)
