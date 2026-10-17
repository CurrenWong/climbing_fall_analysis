"""CARTFIX 实验 5：顶节点「解缠（untangling）」优化。

变量：每个底面节点的顶偏移 d_i（3 dof）。底面完全不动。
目标：F = Σ_e max(0, τ - worst detJ_e)² + λ Σ_i |d_i - d_i^raw|²
约束：0.3 ≤ |d_i| ≤ 1.6 mm
先用有限差分校验解析梯度（本项目有"符号翻车"前科，必须先验）。
"""
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
CM1 = np.load(ROOT / "temp/thums/cart_m1/cartilage_mesh.npz", allow_pickle=True)
RS = (0.0, 1 / 6, 1 / 3, 1 / 2, 2 / 3, 1.0)
TS = (0.0, 0.25, 0.5, 0.75, 1.0)
GS = [(r, s, t) for r in RS for s in RS if r + s <= 1 + 1e-9 for t in TS]
TAU = 0.5
LAM = 0.002
DLO, DHI = 0.30, 1.60


def G_of(r, s, t):
    return np.array([
        [-(1 - t), -(1 - t), -(1 - r - s)],
        [1 - t, 0.0, -r],
        [0.0, 1 - t, -s],
        [-t, -t, (1 - r - s)],
        [t, 0.0, r],
        [0.0, t, s],
    ])


def all_det(P6):
    """(n,6,3) -> (n, |GS|) 所有采样点的 detJ"""
    out = np.empty((len(P6), len(GS)))
    for m, (r, s, t) in enumerate(GS):
        out[:, m] = np.linalg.det(np.einsum("nia,ic->nac", P6, G_of(r, s, t)))
    return out


def ix(x):
    x = int(x)
    return (-x - 1) if x < 0 else (x - 1)


# ══ ① 梯度校验 ══
j = "tibiotalar_tibia"
N = np.asarray(CM1[f"{j}_nodes"], float)
E = np.asarray(CM1[f"{j}_elems"], np.int64)
C = np.array([[ix(x) for x in row] for row in E])
P6 = N[C[0]]
g, (r, s, t) = min(((G_of(*q), q) for q in GS), key=lambda z: np.linalg.det(P6.T @ z[0]))
J = P6.T @ g
det = float(np.linalg.det(J))
Jinv = np.linalg.inv(J)
print("解析梯度按 (i,c) 逐分量校验：")
num = np.zeros(6)
for i in range(6):
    for c in range(3):
        Pp = P6.copy(); Pp[i, c] += 1e-6
        Pm = P6.copy(); Pm[i, c] -= 1e-6
        # 注意：det 是标量，逐分量：
        num[i] += 0  # placeholder
# 逐分量真正校验
print("逐分量有限差分校验 (i,c): 解析 vs 数值")
ok = True
for i in range(6):
    for c in range(3):
        Pp = P6.copy(); Pp[i, c] += 1e-6
        Pm = P6.copy(); Pm[i, c] -= 1e-6
        dnum = (np.linalg.det(Pp.T @ g) - np.linalg.det(Pm.T @ g)) / 2e-6
        # 解析：d(det)/dP_{i,c} = det * Σ_b G[i,b] * Jinv[b,c]
        dana = det * float(g[i] @ Jinv[:, c])
        if abs(dnum - dana) > 1e-4 * max(1.0, abs(dnum)):
            ok = False
            print(f"   ✗ i={i} c={c}: 解析 {dana:.6f} vs 数值 {dnum:.6f}")
print("  梯度校验:", "✅ 通过" if ok else "❌ 不通过")


# ══ ② 解缠 ══
def untangle(j, init, nit=400, verbose=True):
    N = np.asarray(CM1[f"{j}_nodes"], float)
    E = np.asarray(CM1[f"{j}_elems"], np.int64)
    C = np.array([[ix(x) for x in row] for row in E])
    B, T = C[:, :3], C[:, 3:]
    bt = sorted(set(B.ravel().tolist()))
    bi = {b: i for i, b in enumerate(bt)}
    # 顶点 -> (elem, col) 映射
    occ = {b: [] for b in bt}
    for e in range(len(E)):
        for a in range(3):
            occ[int(B[e, a])].append((e, a))
    d_raw = np.zeros((len(bt), 3))
    for b in bt:
        e, a = occ[b][0]
        d_raw[bi[b]] = N[int(T[e, a])] - N[b]
    D = np.array([init[b] for b in bt], float)
    pos = N.copy()

    def rebuild(Dv):
        N2 = N.copy()
        for b in bt:
            for e, a in occ[b]:
                N2[int(T[e, a])] = N[b] + Dv[bi[b]]
        return N2

    def obj(Dv):
        N2 = rebuild(Dv)
        dt = np.min(all_det(N2[C]), axis=1)
        defc = np.maximum(0.0, TAU - dt)
        reg = np.sum((Dv - d_raw) ** 2)
        return float(np.sum(defc ** 2) + LAM * reg), dt, N2

    F, dt, N2 = obj(D)
    if verbose:
        print(f"   init: bad(≤1e-2) {int((dt<=1e-2).sum())} min detJ {dt.min():+.4f} F={F:.3f}")
    hist = []
    for it in range(nit):
        # 取出「违规单元」的采样点 argmin，累加梯度
        bad = np.where(dt < TAU)[0]
        if not len(bad):
            break
        grad = np.zeros_like(D)
        for e in bad:
            idxs = C[e]
            P6 = N2[idxs]
            dets = np.array([np.linalg.det(P6.T @ G_of(*q)) for q in GS])
            m = int(np.argmin(dets))
            dmin = dets[m]
            if dmin >= TAU:
                continue
            g = G_of(*GS[m])
            J = P6.T @ g
            det = dmin
            Jinv = np.linalg.inv(J)
            w = 2.0 * (TAU - det)               # ∂F/∂det * (-1)
            for i in range(6):
                pass
            for i in range(3, 6):               # 只动顶节点（i=3,4,5 对应 T）
                ci = idxs[i]
                b = int(B[e, i - 3])
                # ∂det/∂T_ic = det * Σ_b G[i,b] * Jinv[b,c]
                gg = -det * np.array([g[i] @ Jinv[:, c] for c in range(3)])
                grad[bi[b]] += w * gg
        grad += 2 * LAM * (D - d_raw)
        gn = np.linalg.norm(grad)
        if gn < 1e-12:
            break
        # 回溯线搜索
        step = 0.5
        improved = False
        for _ in range(30):
            Dn = D - step * grad / gn * 0.3
            L = np.linalg.norm(Dn, axis=1)
            Lc = np.clip(L, DLO, DHI)
            Dn = Dn * (Lc / np.maximum(L, 1e-12))[:, None]
            Fn, dtn, N2n = obj(Dn)
            if Fn < F - 1e-12:
                D, F, dt, N2 = Dn, Fn, dtn, N2n
                improved = True
                break
            step *= 0.5
        hist.append((it, F, float(dt.min()), int((dt <= 1e-2).sum())))
        if not improved:
            break
    if verbose and hist:
        for k in (0, len(hist) // 2, len(hist) - 1):
            it, F, mn, nb = hist[k]
            print(f"   it={it:>4}  F={F:9.4f}  min detJ {mn:+9.4f}  bad {nb}")
    th = np.linalg.norm(N2[T].mean(1) - N2[B].mean(1), axis=1)
    mv = float(np.max(np.linalg.norm(N2[T] - N[T], axis=2)))
    return dt, th, mv, D, N2


print()
for j in ["tibiotalar_tibia", "tibiotalar_talus"]:
    N = np.asarray(CM1[f"{j}_nodes"], float)
    E = np.asarray(CM1[f"{j}_elems"], np.int64)
    C = np.array([[ix(x) for x in row] for row in E])
    B, T = C[:, :3], C[:, 3:]
    bt = sorted(set(B.ravel().tolist()))
    raw, nrm = {}, {}
    for e in range(len(E)):
        for a in range(3):
            b = int(B[e, a])
            v = N[int(T[e, a])] - N[b]
            raw[b] = v
            nrm[b] = v / np.linalg.norm(v)
    # 顶点法向初值
    fac = {}
    U = np.sum([raw[b] for b in bt], axis=0)
    U = U / np.linalg.norm(U)
    for e in range(len(E)):
        P = N[B[e]]
        nv = np.cross(P[1] - P[0], P[2] - P[0])
        nv = nv / np.linalg.norm(nv)
        if np.dot(nv, U) < 0:
            nv = -nv
        for q in B[e]:
            fac.setdefault(int(q), []).append(nv)
    navg = {b: (lambda v: v / np.linalg.norm(v))(
        np.sum(fac[b], axis=0)) for b in bt}
    print("=" * 84)
    print(f"{j}")
    for nm, init in (("raw", dict(raw)), ("navg", {b: navg[b] * np.linalg.norm(raw[b]) for b in bt})):
        print(f"  ── 初值 {nm} ──")
        dt, th, mv, D, N2 = untangle(j, init)
        print(f"     ► 最终 bad(≤1e-2) {int((dt<=1e-2).sum())} neg {int((dt<=0).sum())}  "
              f"min detJ {dt.min():+.4f}  厚度 [{th.min():.4f},{th.max():.4f}] "
              f"中位 {np.median(th):.4f}  顶面位移 {mv:.3f}mm")
