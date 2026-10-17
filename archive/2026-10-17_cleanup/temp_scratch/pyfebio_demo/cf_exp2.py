"""CARTFIX 实验 2：把「抬升方向」限制在节点的法向锥内（保模长）。

判据（诊断 7）：单元翻转 ⟺ min_i d_i·n̂_facet < 0。
⇒ 修复目标：对每个底面节点，把偏移方向 u_i 调到「其所有相邻必 facet 法向的正半空间
   交集」（法向锥）之内；模长 t_i=|d_i| 保持不变（⇒ 厚度场不动）。
算法：从原始 u_i 出发，迭代 u ← normalize(u + Σ_{被违反的 n̂}) 直到全部正投影。
"""
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
CM1 = np.load(ROOT / "temp/thums/cart_m1/cartilage_mesh.npz", allow_pickle=True)

RS = (0.0, 1 / 6, 1 / 3, 1 / 2, 2 / 3, 1.0)
TS = (0.0, 0.25, 0.5, 0.75, 1.0)
GS = [(r, s, t) for r in RS for s in RS if r + s <= 1 + 1e-9 for t in TS]


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


EPS = float(__import__("os").environ.get("CF_EPS", "0.05"))     # 要求的最小正投影(mm)


def run(j, mode, nit=60):
    N = np.asarray(CM1[f"{j}_nodes"], float)
    E = np.asarray(CM1[f"{j}_elems"], np.int64)
    off = np.asarray(CM1[f"{j}_is_offset"], bool)
    C = np.array([[ix(x) for x in row] for row in E])
    B, T = C[:, :3], C[:, 3:]
    bt = sorted(set(B.ravel().tolist()))
    bi = {b: i for i, b in enumerate(bt)}
    # 逐底节点 原始 disp
    d0 = {}
    for r in range(len(E)):
        for a in range(3):
            d0[int(B[r, a])] = N[int(T[r, a])] - N[int(B[r, a])]
    t_i = {b: np.linalg.norm(d0[b]) for b in bt}
    u_i = {b: d0[b] / (t_i[b] + 1e-15) for b in bt}
    # 底面 facet 法向（依 mean 原始 disp 定向）
    fac = {}
    for r in range(len(E)):
        P = N[B[r]]
        nv = np.cross(P[1] - P[0], P[2] - P[0])
        nv = nv / (np.linalg.norm(nv) + 1e-15)
        dm = np.mean([d0[int(q)] for q in B[r]], axis=0)
        if np.dot(nv, dm) < 0:
            nv = -nv
        ar = 0.5 * np.linalg.norm(np.cross(P[1] - P[0], P[2] - P[0]))
        for q in B[r]:
            fac.setdefault(int(q), []).append((nv, ar))
    # 迭代：把方向推进法向锥
    u_new = dict(u_i)
    nchg = 0
    for b in bt:
        u = u_i[b].copy()
        for _ in range(nit):
            viol = [nv for nv, _ in fac[b] if np.dot(u, nv) <= EPS]
            if not viol:
                break
            u = u + np.sum(viol, axis=0)
            u = u / (np.linalg.norm(u) + 1e-15)
        u_new[b] = u
        if np.degrees(np.arccos(np.clip(np.dot(u, u_i[b]), -1, 1))) > 0.01:
            nchg += 1
    if mode == "raw":
        u_new = u_i
    # 重建顶面（内/边节点用重心插值 ⇒ 这里直接按单元顶点法线性插值，M=1 时等价）
    N2 = N.copy()
    for b in bt:
        for r in range(len(E)):
            for a in range(3):
                if int(B[r, a]) == b:
                    N2[int(T[r, a])] = N[b] + t_i[b] * u_new[b]
    P6 = N2[C]
    W = dets(P6)
    th = np.linalg.norm(N2[T].mean(1) - N2[B].mean(1), axis=1)
    moved = max(np.linalg.norm(N2[t] - N[t]) for t in T.ravel())
    # 逐顶点投影
    nb = np.cross(N2[B[:, 1]] - N2[B[:, 0]], N2[B[:, 2]] - N2[B[:, 0]])
    nbu = nb / (np.linalg.norm(nb, axis=1, keepdims=True) + 1e-15)
    d = N2[T] - N2[B]
    sgn = np.sign(np.einsum("ij,ij->i", nbu, d.mean(1)))
    proj = np.einsum("nij,nj->ni", d, sgn[:, None] * nbu)
    return dict(bad=int((W <= 1e-2).sum()), neg=int((W <= 0).sum()), wmin=float(W.min()),
                thmin=float(th.min()), thmax=float(th.max()), thmed=float(np.median(th)),
                moved=float(moved), nchg=nchg, pmin=float(proj.min()))


print(f"EPS={EPS}")
for j in ["tibiotalar_tibia", "tibiotalar_talus"]:
    for mode in ("raw", "cone"):
        r = run(j, mode)
        print(f"{j:20s} {mode:5s}: bad {r['bad']:>3} neg {r['neg']:>3} min detJ {r['wmin']:+9.4f} | "
              f"厚度 [{r['thmin']:.4f},{r['thmax']:.4f}] 中位 {r['thmed']:.4f} | "
              f"顶面最大位移 {r['moved']:.4f}mm | 改向节点 {r['nchg']} | min proj {r['pmin']:+.4f}")
