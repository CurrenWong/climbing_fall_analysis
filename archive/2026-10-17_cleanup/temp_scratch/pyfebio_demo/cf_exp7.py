"""CARTFIX 实验 7：残留坏单元诊断 + 优化器参数扫描（目标 bad=0 两边都达到）。"""
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
    edges = np.array(sorted(edges))
    return N, E, C, B, T, bt, bi, d_raw, edges


def opt(j, tau=0.5, lam_reg=0.002, lam_sm=0.001, nit=4000, scale=0.02, verbose=False):
    N, E, C, B, T, bt, bi, d_raw, edges = prep(j)

    def rebuild(D):
        N2 = N.copy()
        for e in range(len(E)):
            for a in range(3):
                N2[int(T[e, a])] = N[int(B[e, a])] + D[bi[int(B[e, a])]]
        return N2

    def obj(D):
        N2 = rebuild(D)
        dt = all_det(N2[C])
        F = float(np.sum(np.maximum(0.0, tau - dt) ** 2))
        F += lam_reg * float(np.sum((D - d_raw) ** 2))
        if lam_sm:
            dd = D[edges[:, 0]] - D[edges[:, 1]]
            F += lam_sm * float(np.sum(dd ** 2))
        return F, dt, N2

    D = d_raw.copy()
    F, dt, N2 = obj(D)
    for it in range(nit):
        bad = np.where(dt.min(axis=1) < tau)[0]
        if not len(bad):
            break
        grad = np.zeros_like(D)
        for e in bad:
            P6 = N2[C[e]]
            for m in np.where(dt[e] < tau)[0]:
                g = GG[m]
                J = P6.T @ g
                det = float(np.linalg.det(J))
                Jinv = np.linalg.inv(J + (1e-9 if det <= 1e-12 else 0.0) * np.eye(3))
                w = -2.0 * (tau - det) * det
                for i in range(3, 6):
                    grad[bi[int(B[e, i - 3])]] += w * (g[i] @ Jinv)
        grad += 2 * lam_reg * (D - d_raw)
        if lam_sm:
            np.add.at(grad, edges[:, 0], 2 * lam_sm * (D[edges[:, 0]] - D[edges[:, 1]]))
            np.add.at(grad, edges[:, 1], 2 * lam_sm * (D[edges[:, 1]] - D[edges[:, 0]]))
        gmax = np.abs(grad).max()
        if gmax < 1e-13:
            break
        for al in (1.0, 0.5, 0.25, 0.1, 0.03, 0.01, 0.003):
            Dn = D - al * scale * grad / gmax
            L = np.linalg.norm(Dn, axis=1)
            Dn = Dn * (np.clip(L, 0.30, 1.60) / np.maximum(L, 1e-12))[:, None]
            Fn, dtn, N2n = obj(Dn)
            if Fn < F - 1e-14:
                D, F, dt, N2 = Dn, Fn, dtn, N2n
                break
        else:
            break
    return D, dt, N2


print("=== 残留坏单元诊断（原网格）===")
for j in ["tibiotalar_tibia", "tibiotalar_talus"]:
    N, E, C, B, T, bt, bi, d_raw, edges = prep(j)
    dt0 = all_det(N[C]).min(axis=1)
    k = int(np.argmin(dt0))
    print(f"{j}: worst #{k} detJ={dt0[k]:.4f}; 该单元 3 顶点 disp:")
    for a in range(3):
        b = int(B[k, a])
        print(f"    col{a}: bottom节点 {b}  d={np.round(d_raw[bi[b]],4).tolist()} "
              f"|d|={np.linalg.norm(d_raw[bi[b]]):.4f}  所属单元数={sum(1 for e in range(len(E)) for q in range(3) if int(B[e,q])==b)}")

print()
print("=== 参数扫描 ===")
for tau in (0.5, 1.0, 2.0):
    for lam_sm in (0.0, 0.001, 0.01):
        line = f"tau={tau} lam_sm={lam_sm}: "
        for j in ["tibiotalar_tibia", "tibiotalar_talus"]:
            D, dt, N2 = opt(j, tau=tau, lam_sm=lam_sm, nit=4000)
            N, E, C, B, T, bt, bi, d_raw, edges = prep(j)
            th = np.linalg.norm(N2[T].mean(1) - N2[B].mean(1), axis=1)
            mv = float(np.max(np.linalg.norm(N2[T] - N[T], axis=2)))
            line += (f"[{j[11:]}: bad {int((dt.min(1)<=1e-2).sum())} min {dt.min():+.3f} "
                     f"th [{th.min():.3f},{th.max():.3f}] mv {mv:.2f}] ")
        print(line)
