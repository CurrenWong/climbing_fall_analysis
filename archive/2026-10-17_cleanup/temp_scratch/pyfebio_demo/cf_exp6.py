"""CARTFIX 实验 6：加强版解缠（全活跃角点二次罚 + 场光滑 + 回溯线搜索）。

F = Σ_e Σ_{角点 det<τ} (τ-det)²  +  λ_reg Σ_i |d_i-d_i^raw|²  +  λ_sm Σ_{i~j} |d_i-d_j|²
约束 0.3 ≤ |d_i| ≤ 1.6。
"""
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
CM1 = np.load(ROOT / "temp/thums/cart_m1/cartilage_mesh.npz", allow_pickle=True)
RS = (0.0, 1 / 6, 1 / 3, 1 / 2, 2 / 3, 1.0)
TS = (0.0, 0.25, 0.5, 0.75, 1.0)
GS = [(r, s, t) for r in RS for s in RS if r + s <= 1 + 1e-9 for t in TS]
GG = []
for (r, s, t) in GS:
    GG.append(np.array([
        [-(1 - t), -(1 - t), -(1 - r - s)],
        [1 - t, 0.0, -r],
        [0.0, 1 - t, -s],
        [-t, -t, (1 - r - s)],
        [t, 0.0, r],
        [0.0, t, s],
    ]))
GG = np.array(GG)                      # (nS, 6, 3)


def all_det(P6):
    return np.linalg.det(np.einsum("nic,mia->nmca", P6, GG))


def ix(x):
    x = int(x)
    return (-x - 1) if x < 0 else (x - 1)


def optimize(j, tau=0.5, lam_reg=0.002, lam_sm=0.001, nit=1500, verbose=True):
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
    # 边（底面节点图）
    edges = set()
    for e in range(len(E)):
        for a in range(3):
            x = bi[int(B[e, a])]
            y = bi[int(B[e, (a + 1) % 3])]
            edges.add((min(x, y), max(x, y)))
    edges = np.array(sorted(edges))
    ii, jj, kk = C[:, :3], C[:, 3:], np.array([bi[int(q)] for q in C[:, :3].ravel()]).reshape(-1, 3)

    def rebuild(D):
        N2 = N.copy()
        for e in range(len(E)):
            for a in range(3):
                N2[int(T[e, a])] = N[int(B[e, a])] + D[bi[int(B[e, a])]]
        return N2

    def obj(D):
        N2 = rebuild(D)
        dt = all_det(N2[C])
        act = np.maximum(0.0, tau - dt)
        F = float(np.sum(act ** 2)) + lam_reg * float(np.sum((D - d_raw) ** 2))
        if lam_sm:
            dd = D[edges[:, 0]] - D[edges[:, 1]]
            F += lam_sm * float(np.sum(dd ** 2))
        return F, dt, N2

    D = d_raw.copy()
    F, dt, N2 = obj(D)
    if verbose:
        print(f"   init bad {int((dt.min(1)<=1e-2).sum())} min {dt.min():+.4f} F={F:.4f}")
    for it in range(nit):
        mins = dt.min(axis=1)
        bad = np.where(mins < tau)[0]
        if not len(bad):
            break
        grad = np.zeros_like(D)
        for e in bad:
            P6 = N2[C[e]]
            dts = dt[e]
            idx = np.where(dts < tau)[0]
            for m in idx:
                g = GG[m]
                J = P6.T @ g
                det = float(np.linalg.det(J))
                if det > 1e-12:
                    Jinv = np.linalg.inv(J)
                else:
                    Jinv = np.linalg.inv(J + 1e-9 * np.eye(3))
                w = -2.0 * (tau - det) * det
                for i in range(3, 6):
                    b = bi[int(B[e, i - 3])]
                    grad[b] += w * (g[i] @ Jinv)
        grad += 2 * lam_reg * (D - d_raw)
        if lam_sm:
            np.add.at(grad, edges[:, 0], 2 * lam_sm * (D[edges[:, 0]] - D[edges[:, 1]]))
            np.add.at(grad, edges[:, 1], 2 * lam_sm * (D[edges[:, 1]] - D[edges[:, 0]]))
        gn = np.linalg.norm(grad)
        if gn < 1e-13:
            break
        step = 0.6
        for _ in range(40):
            Dn = D - step * 0.25 * grad / gn
            L = np.linalg.norm(Dn, axis=1)
            Dn = Dn * (np.clip(L, 0.30, 1.60) / np.maximum(L, 1e-12))[:, None]
            Fn, dtn, N2n = obj(Dn)
            if Fn < F - 1e-13:
                D, F, dt, N2 = Dn, Fn, dtn, N2n
                break
            step *= 0.6
        else:
            break
    th = np.linalg.norm(N2[T].mean(1) - N2[B].mean(1), axis=1)
    mv = float(np.max(np.linalg.norm(N2[T] - N[T], axis=2)))
    return dict(bad=int((dt.min(1) <= 1e-2).sum()), neg=int((dt.min(1) <= 0).sum()),
                wmin=float(dt.min()), th=(float(th.min()), float(np.median(th)), float(th.max())),
                mv=mv, D=D, dt=dt, N2=N2, C=C, B=B, T=T, N=N, bt=bt)


if __name__ == "__main__":
    for j in ["tibiotalar_tibia", "tibiotalar_talus"]:
        print("=" * 84)
        print(j)
        r = optimize(j)
        print(f"   ► bad {r['bad']} neg {r['neg']} min detJ {r['wmin']:+.4f}  "
              f"厚度 {r['th']}  顶面位移 {r['mv']:.3f}mm")
        np.savez(ROOT / f"temp/pyfebio_demo/cf_opt_{j}.npz", D=r["D"], bt=np.array(r["bt"]))
