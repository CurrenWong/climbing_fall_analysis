"""CARTFIX 实验 3：若干「抬升方向场」候选对比（底面/模长都不动）。

候选（都保持 t_i=|d_i| 不变，只改方向 u_i）：
  raw        : 原样
  navg       : 面积加权节点法向（按全局参考方向 U 定向）
  navg*S     : 再 Laplacian 平滑 S 次
  X+cone     : 再做法向锥可行性迭代（u ← normalize(u+Σ违反 n̂)）
"""
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
CM1 = np.load(ROOT / "temp/thums/cart_m1/cartilage_mesh.npz", allow_pickle=True)
RS = (0.0, 1 / 6, 1 / 3, 1 / 2, 2 / 3, 1.0)
TS = (0.0, 0.25, 0.5, 0.75, 1.0)
GS = [(r, s, t) for r in RS for s in RS if r + s <= 1 + 1e-9 for t in TS]
EPS = 0.05


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


def cone_fix(u, ns, nit=80):
    u = u / (np.linalg.norm(u) + 1e-15)
    for _ in range(nit):
        viol = [nv for nv in ns if np.dot(u, nv) <= EPS]
        if not viol:
            return u, True
        u = u + np.sum(viol, axis=0)
        u = u / (np.linalg.norm(u) + 1e-15)
    return u, False


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
    u_raw = {b: d0[b] / (t_i[b] + 1e-15) for b in bt}
    U = np.sum([d0[b] for b in bt], axis=0)
    U = U / np.linalg.norm(U)
    # facet normals（解析三角顺序 + 全局 U 定向）
    fac = {}
    for r in range(len(E)):
        P = N[B[r]]
        cr = np.cross(P[1] - P[0], P[2] - P[0])
        ar2 = np.linalg.norm(cr)
        nv = cr / (ar2 + 1e-15)
        if np.dot(nv, U) < 0:
            nv = -nv
        for q in B[r]:
            fac.setdefault(int(q), []).append((nv, ar2 / 2))
    adj = {b: set() for b in bt}
    for r in range(len(E)):
        for a in range(3):
            for c in range(3):
                if a != c:
                    adj[int(B[r, a])].add(int(B[r, c]))

    def build(u_new):
        N2 = N.copy()
        for r in range(len(E)):
            for a in range(3):
                b = int(B[r, a])
                N2[int(T[r, a])] = N[b] + t_i[b] * u_new[b]
        W = dets(N2[C])
        th = np.linalg.norm(N2[T].mean(1) - N2[B].mean(1), axis=1)
        moved = max(np.linalg.norm(N2[int(T[r, a])] - N[int(T[r, a])])
                    for r in range(len(E)) for a in range(3))
        # 逐顶点投影（用存储顺序法向，即 build 后由代码 newell 定向的那一个）
        nb = np.cross(N2[B[:, 1]] - N2[B[:, 0]], N2[B[:, 2]] - N2[B[:, 0]])
        nbu = nb / (np.linalg.norm(nb, axis=1, keepdims=True) + 1e-15)
        dd = N2[T] - N2[B]
        sg = np.sign(np.einsum("ij,ij->i", nbu, dd.mean(1)))
        pr = np.einsum("nij,nj->ni", dd, sg[:, None] * nbu)
        return W, th, moved, pr.min()

    print("=" * 92)
    print(f"{j}  U={np.round(U,3)}")
    cands = {}
    cands["raw"] = {b: u_raw[b].copy() for b in bt}
    navg = {}
    for b in bt:
        v = np.sum([nv * ar for nv, ar in fac[b]], axis=0)
        navg[b] = v / (np.linalg.norm(v) + 1e-15)
    cands["navg"] = navg
    for S in (1, 2, 4, 8):
        cur = {b: navg[b].copy() for b in bt}
        for _ in range(S):
            nxt = {}
            for b in bt:
                v = cur[b] + np.mean([cur[x] for x in adj[b]], axis=0)
                nxt[b] = v / (np.linalg.norm(v) + 1e-15)
            cur = nxt
        cands[f"navg+sm{S}"] = cur
    for nm, base in (("raw", {b: u_raw[b].copy() for b in bt}), ("navg", navg)):
        cur, nf = {}, 0
        for b in bt:
            u, ok = cone_fix(base[b].copy(), [nv for nv, _ in fac[b]])
            cur[b] = u
            nf += 0 if ok else 1
        cands[f"{nm}+cone"] = cur
        print(f"    [{nm}+cone] 法向锥不可行节点数 = {nf}")
    for nm, uu in cands.items():
        W, th, moved, pmin = build(uu)
        print(f"   {nm:12s} bad {int((W<=1e-2).sum()):>3}  neg {int((W<=0).sum()):>3}  "
              f"min detJ {W.min():+9.4f} | 厚度 [{th.min():.4f},{th.max():.4f}] 中位 {np.median(th):.4f} | "
              f"顶面位移 {moved:.3f}mm | min proj {pmin:+.3f}")
