"""CARTFIX 快速验证：全局解缠 ⇄ 残留单元局部修复（交替），目标 bad=0。

全局：F = Σ_e Σ_{det<τ}(τ-det)² + λ_reg|d-d_raw|² + λ_sm Σ|d_i-d_j|²
局部：只对「残留坏单元」的 3 个底节点做块优化（目标 = 这些节点的所有相邻单元都不翻）
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


def all_det(P6):
    return np.linalg.det(np.einsum("nic,mia->nmca", P6, GG))


def ix(x):
    x = int(x)
    return (-x - 1) if x < 0 else (x - 1)


def prep(j):
    N = np.asarray(CM1[f"{j}_nodes"], float)
    E = np.asarray(CM1[f"{j}_elems"], np.int64)
    C = np.array([[ix(x) for x in row] for row in E])
    B, T = C[:, :3], C[:, 3:]
    bt = sorted(set(B.ravel().tolist()))
    bi = {b: i for i, b in enumerate(bt)}
    occ = {b: [] for b in bt}
    for e in range(len(E)):
        for a in range(3):
            occ[int(B[e, a])].append((e, a))
    d_raw = np.zeros((len(bt), 3))
    for b in bt:
        e, a = occ[b][0]
        d_raw[bi[b]] = N[int(T[e, a])] - N[b]
    edges = set()
    for e in range(len(E)):
        for a in range(3):
            x, y = bi[int(B[e, a])], bi[int(B[e, (a + 1) % 3])]
            edges.add((min(x, y), max(x, y)))
    return dict(N=N, E=E, C=C, B=B, T=T, bt=bt, bi=bi, d_raw=d_raw,
                edges=np.array(sorted(edges)))


def run(j, nit=3000):
    P = prep(j)
    N, E, C, B, T, bi, d_raw, edges = (P["N"], P["E"], P["C"], P["B"], P["T"],
                                       P["bi"], P["d_raw"], P["edges"])

    def rebuild(D):
        N2 = N.copy()
        for e in range(len(E)):
            for a in range(3):
                N2[int(T[e, a])] = N[int(B[e, a])] + D[bi[int(B[e, a])]]
        return N2

    def obj(D, sel=None):
        N2 = rebuild(D)
        dt = all_det(N2[C])
        sel = np.arange(len(E)) if sel is None else sel
        F = float(np.sum(np.maximum(0.0, TAU - dt[sel]) ** 2))
        F += LAM_REG * float(np.sum((D - d_raw) ** 2))
        if LAM_SM:
            dd = D[edges[:, 0]] - D[edges[:, 1]]
            F += LAM_SM * float(np.sum(dd ** 2))
        return F, dt, N2

    def grad_of(D, N2, dt, sel):
        g = np.zeros_like(D)
        for e in sel:
            P6 = N2[C[e]]
            for m in np.where(dt[e] < TAU)[0]:
                gm = GG[m]
                J = P6.T @ gm
                det = float(np.linalg.det(J))
                Jinv = np.linalg.inv(J + (1e-9 if det <= 1e-12 else 0.0) * np.eye(3))
                w = -2.0 * (TAU - det) * det
                for i in range(3, 6):
                    g[bi[int(B[e, i - 3])]] += w * (gm[i] @ Jinv)
        g += 2 * LAM_REG * (D - d_raw)
        if LAM_SM:
            np.add.at(g, edges[:, 0], 2 * LAM_SM * (D[edges[:, 0]] - D[edges[:, 1]]))
            np.add.at(g, edges[:, 1], 2 * LAM_SM * (D[edges[:, 1]] - D[edges[:, 0]]))
        return g

    def descend(D, sel, iters, scale):
        F, dt, N2 = obj(D, sel)
        for _ in range(iters):
            bad = np.array([e for e in sel if dt[e].min() < TAU])
            if not len(bad):
                break
            g = grad_of(D, N2, dt, bad)
            gmax = np.abs(g).max()
            if gmax < 1e-13:
                break
            for al in (1.0, 0.5, 0.25, 0.1, 0.03, 0.01, 0.003):
                Dn = D - al * scale * g / gmax
                L = np.linalg.norm(Dn, axis=1)
                Dn = Dn * (np.clip(L, 0.30, 1.60) / np.maximum(L, 1e-12))[:, None]
                Fn, dtn, N2n = obj(Dn, sel)
                if Fn < F - 1e-14:
                    D, F, dt, N2 = Dn, Fn, dtn, N2n
                    break
            else:
                break
        return D, dt, N2

    t0 = time.time()
    D = d_raw.copy()
    for rnd in range(6):
        D, dt, N2 = descend(D, np.arange(len(E)), nit, 0.02)
        mins = dt.min(axis=1)
        bad = np.where(mins < TAU)[0]
        if not len(bad):
            break
        # 局部块修复：逐个残留单元的 3 个底节点（含其全部相邻单元）
        nadj = {b: [e for e in range(len(E))
                    if any(bi[int(B[e, a])] == b for a in range(3))]
                for b in bi.values()}
        for e in bad:
            nodes = [bi[int(B[e, a])] for a in range(3)]
            fwd = sorted(set(sum((nadj[b] for b in nodes), [])))
            D, dt, N2 = descend(D, np.array(fwd), 400, 0.05)
        mins = dt.min(axis=1)
        print(f"   round {rnd}: bad {int((mins<=1e-2).sum())} min {mins.min():+.4f} "
              f"({time.time()-t0:.1f}s)")
    th = np.linalg.norm(N2[T].mean(1) - N2[B].mean(1), axis=1)
    mv = float(np.max(np.linalg.norm(N2[T] - N[T], axis=2)))
    return D, dt, th, mv, N2, P


print(f"TAU={TAU} LAM_REG={LAM_REG} LAM_SM={LAM_SM}")
for j in ["tibiotalar_tibia", "tibiotalar_talus"]:
    print("=" * 70)
    print(j)
    D, dt, th, mv, N2, P = run(j)
    print(f"  ► bad(≤1e-2) {int((dt.min(1)<=1e-2).sum())}  neg {int((dt.min(1)<=0).sum())}  "
          f"min detJ {dt.min():+.4f}")
    print(f"    厚度 [{th.min():.4f}, {np.median(th):.4f}, {th.max():.4f}]  "
          f"顶面位移 {mv:.3f}mm")
    np.savez(ROOT / f"temp/pyfebio_demo/cf_fast_{j}.npz", D=D, bt=np.array(P["bt"]))
