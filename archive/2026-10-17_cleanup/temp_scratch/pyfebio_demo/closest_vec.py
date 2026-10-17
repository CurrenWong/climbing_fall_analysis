"""向量化 point→triangle 最近点（Ericson），带自检。返回距离/三角号/垂足/法向。"""
import numpy as np

EPS = 1e-30


def closest_vec(P, tris):
    """P:(n,3)  tris:(m,3,3) -> dist(n,), tidx(n,), foot(n,3)"""
    P = np.atleast_2d(np.asarray(P, float))
    n, m = len(P), len(tris)
    A = tris[:, 0][None]        # (1,m,3)
    B = tris[:, 1][None]
    C = tris[:, 2][None]
    AB = B - A                  # (1,m,3)
    AC = C - A
    Pp = P[:, None, :]          # (n,1,3)
    AP = Pp - A                 # (n,m,3)
    d1 = np.einsum("nmk,nmk->nm", AP, AB)
    d2 = np.einsum("nmk,nmk->nm", AP, AC)
    BP = Pp - B
    d3 = np.einsum("nmk,nmk->nm", BP, AB)
    d4 = np.einsum("nmk,nmk->nm", BP, AC)
    CP = Pp - C
    d5 = np.einsum("nmk,nmk->nm", CP, AB)
    d6 = np.einsum("nmk,nmk->nm", CP, AC)
    va = d3 * d6 - d5 * d4
    vb = d5 * d2 - d1 * d6
    vc = d1 * d4 - d3 * d2

    mAB = np.einsum("mk,mk->m", AB.reshape(-1, 3), AB.reshape(-1, 3))
    mAC = np.einsum("mk,mk->m", AC.reshape(-1, 3), AC.reshape(-1, 3))
    mAB = mAB[None]
    mAC = mAC[None]

    def proj_v(vtx):
        return np.broadcast_to(vtx, (n, m, 3))

    cands, dists = [], []
    # A
    cands.append(proj_v(A))
    dists.append(np.linalg.norm(AP, axis=2))
    # B
    cands.append(proj_v(B))
    dists.append(np.linalg.norm(BP, axis=2))
    # C
    cands.append(proj_v(C))
    dists.append(np.linalg.norm(CP, axis=2))
    # edge AB: v = d1/(d1-d3)
    den_ab = d1 - d3
    v = np.clip(d1 / np.where(np.abs(den_ab) < EPS, 1.0, den_ab)
                if False else d1 / np.where(np.abs(den_ab) < EPS, EPS, den_ab), 0, 1)[:, :, None]
    cands.append(A + v * AB)
    dists.append(np.linalg.norm(Pp - (A + v * AB), axis=2))
    # edge AC: w = d2/(d2-d6)
    den_ac = d2 - d6
    w = np.clip(d2 / np.where(np.abs(den_ac) < EPS, EPS, den_ac), 0, 1)[:, :, None]
    cands.append(A + w * AC)
    dists.append(np.linalg.norm(Pp - (A + w * AC), axis=2))
    # edge BC: w = (d4-d3)/((d4-d3)+(d5-d6))
    den_bc = (d4 - d3) + (d5 - d6)
    tb = np.clip((d4 - d3) / np.where(np.abs(den_bc) < EPS, EPS, den_bc),
                 0, 1)[:, :, None]
    cands.append(B + tb * (C - B))
    dists.append(np.linalg.norm(Pp - (B + tb * (C - B)), axis=2))
    # face
    den = va + vb + vc
    den_s = np.where(np.abs(den) < EPS, 1.0, den)
    vv = vb / den_s
    ww = vc / den_s
    Fpt = A + vv[:, :, None] * AB + ww[:, :, None] * AC
    cands.append(Fpt)
    dists.append(np.linalg.norm(Pp - Fpt, axis=2))

    # 正确的区域选择（严格按 Ericson 顺序，用掩码互斥）
    inside = np.zeros((n, m), bool)
    reg = np.full((n, m), -1, int)

    def take(mask, idx):
        nonlocal inside
        sel = mask & ~inside
        reg[sel] = idx
        inside = inside | sel

    take((d1 <= 0) & (d2 <= 0), 0)
    take((d3 >= 0) & (d4 <= d3), 1)
    take((vc <= 0) & (d1 >= 0) & (d3 <= 0), 3)
    take((d6 >= 0) & (d5 <= d6), 2)
    take((vb <= 0) & (d2 >= 0) & (d6 <= 0), 4)
    take((va <= 0) & (d4 - d3 >= 0) & (d5 - d6 >= 0), 5)
    take(np.ones((n, m), bool), 6)
    # 兜底：面投影若重心为负则退回最近边
    bad = reg == 6
    if bad.any():
        reg[bad & ((vv < 0) | (ww < 0) | (vv + ww > 1))] = 6  # 保持，下面用 cands

    D = np.stack(dists, 0)             # (7,n,m)
    sel = np.take_along_axis(D, reg[None], 0)[0]      # (n,m)
    j = np.argmin(sel, axis=1)
    i = np.arange(n)
    dist = sel[i, j]
    tidx = j
    # 取对应候选点
    CS = np.stack(cands, 0)            # (7,n,m,3)
    foot = CS[reg[i, j], i, j]
    nrm = np.cross(AB[0, j], AC[0, j])
    ln = np.linalg.norm(nrm, axis=1, keepdims=True)
    nrm = nrm / np.where(ln < EPS, 1, ln)
    return dist, tidx, foot, nrm


if __name__ == "__main__":
    import pttri
    rng = np.random.default_rng(7)
    A, B, C = np.array([0., 0, 0]), np.array([1., 0, 0]), np.array([0., 1, 0])
    tris = np.array([[A, B, C]])
    tests = [([0.3, 0.3, 2.0], 2.0), ([5.0, 0.0, 0.0], 4.0), ([-1., 0., 0.], 1.0),
             ([0., 5., 0.], 4.0), ([0.5, 0.5, 0.], 0.0), ([0.5, 0.6, 3.0],
              np.sqrt(9 + 0.005)), ([1., 1., 0.], np.sqrt(0.5))]
    ok = True
    for p, exp in tests:
        got = closest_vec(np.array([p], float), tris)[0][0]
        g = abs(got - exp) < 1e-9
        ok &= g
        print(f"  [{'OK' if g else 'XX'}] {str(p):22s} 期望 {exp:.6f} 实得 {got:.6f}")
    # 随机对照标量版
    err = 0.0
    for _ in range(6):
        T = rng.normal(size=(9, 3, 3)) * 2
        Q = rng.normal(size=(11, 3)) * 2
        dv = closest_vec(Q, T)[0]
        ds = pttri.dist_pt_tris(Q, T)
        err = max(err, float(np.abs(dv - ds).max()))
    print(f"\n  与标量版最大偏差 = {err:.3e}  [{'OK' if err < 1e-9 else 'XX'}]")
    print(f"\n自检 {'ALL PASS' if ok and err < 1e-9 else 'FAILED'}")
