"""A1 几何预检（只测不改）：subtalar_talus ↔ subtalar_calcaneus 软骨外表面间隙。

- 正对照 = tibiotalar 对（已知能正常求解的既有接触对）⇒ 检查器本身有正对照。
- 间隙定义：
    d_min  : A 面片质心 → B 面片的最近点距离（point-triangle）
    d_n    : 沿 A 外法向的间隙 (p - c)·n_A（>0 = 分离；<0 = 穿透）
- 法向相对：最近点对处 n_A·n_B（<0 = 大体相对/面对面；>0 = 大体同向）
用法: python temp/pyfebio_demo/a1_subtalar_gap.py
"""
import json
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
CART = ROOT / "temp/thums/cart_m1/cartilage_mesh.npz"
OUTJ = ROOT / "temp/pyfebio_demo/a1_subtalar_gap.json"

PENTA_FACES = [(0, 1, 2), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)]
TET_FACES = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
MIN_FACET_AREA = 1e-8


def cnidx(x):
    x = int(x)
    return (-x - 1) if x < 0 else (x - 1)


def penta6_dN(r, s, t):
    L = (1 - r - s, r, s)
    dL = ((-1.0, -1.0), (1.0, 0.0), (0.0, 1.0))
    G = np.zeros((6, 3))
    for i in range(3):
        G[i] = ((1 - t) / 2 * dL[i][0], (1 - t) / 2 * dL[i][1], -L[i] / 2)
        G[i + 3] = ((1 + t) / 2 * dL[i][0], (1 + t) / 2 * dL[i][1], L[i] / 2)
    return G


_PPTS = [(1.0, 0.0, -1.0), (0.0, 1.0, -1.0), (0.0, 0.0, -1.0),
         (1.0, 0.0, 1.0), (0.0, 1.0, 1.0), (0.0, 0.0, 1.0),
         (1 / 3, 1 / 3, 0.0), (0.25, 0.25, -0.5), (0.25, 0.25, 0.5),
         (0.5, 0.25, -0.5), (0.25, 0.5, -0.5), (0.5, 0.25, 0.5),
         (0.25, 0.5, 0.5), (0.4, 0.4, 0.0)]


def min_det_penta6(X, idx):
    P = X[list(idx)]
    worst = np.inf
    for (r, s, t) in _PPTS:
        d = float(np.linalg.det(P.T @ penta6_dN(r, s, t)))
        worst = min(worst, d)
    return worst


def sv_penta6(P):
    dt = (P[3] + P[4] + P[5] - P[0] - P[1] - P[2]) / 3.0
    return float(np.linalg.det(np.stack([P[1] - P[0], P[2] - P[0], dt])))


def outer_facets(cg, j, det_min=0.0, apply_gate=True):
    """软骨外表面 facet：全部顶点为 is_offset；按 ankle_lig_feb 同款
    detJ 闸门（DET_MIN，acceptance 用 0.0）决定单元存活；法向朝外。"""
    X = np.asarray(cg[f"{j}_nodes"], float)
    isoff = np.asarray(cg[f"{j}_is_offset"], bool)
    alive, gate_dropped = 0, 0
    for row in cg[f"{j}_elems"]:
        raw = [int(x) for x in row if int(x) != 0]
        if len(raw) != 6:
            continue
        idx = [cnidx(x) for x in raw]
        if apply_gate:
            d = min_det_penta6(X, idx)
            if not (sv_penta6(X[idx]) > 0 and d > det_min):
                alt = idx[3:] + idx[:3]
                if sv_penta6(X[alt]) > 0 and min_det_penta6(X, alt) > det_min:
                    idx = alt
                else:
                    gate_dropped += 1
                    continue
        alive += 1
        rest_all = [k for k in range(6)]
        for f in PENTA_FACES:
            fv = [idx[z] for z in f]
            if not all(isoff[k] for k in fv):
                continue
            if len(set(fv)) != len(fv):
                continue
            P = X[fv]
            a = 0.5 * np.linalg.norm(np.cross(P[1] - P[0], P[2] - P[0]))
            if a < MIN_FACET_AREA:
                continue
            rest = [idx[z] for z in rest_all if z not in f]
            cen = X[rest].mean(0)
            n = np.zeros(3)
            for i in range(3):
                n += np.cross(P[i], P[(i + 1) % 3])
            if np.dot(n, P.mean(0) - cen) < 0:
                fv = fv[::-1]
                P = X[fv]
                n = -n
            yield P, n
    # (generator 方向) 计数不便从生成器返回，alive/gate 在外面单独打印


def gate_stats(cg, j, det_min=0.0):
    X = np.asarray(cg[f"{j}_nodes"], float)
    alive = dropped = 0
    for row in cg[f"{j}_elems"]:
        raw = [int(x) for x in row if int(x) != 0]
        if len(raw) != 6:
            continue
        idx = [cnidx(x) for x in raw]
        ok = sv_penta6(X[idx]) > 0 and min_det_penta6(X, idx) > det_min
        if not ok:
            alt = idx[3:] + idx[:3]
            ok = sv_penta6(X[alt]) > 0 and min_det_penta6(X, alt) > det_min
        alive += 1 if ok else 0
        dropped += 0 if ok else 1
    return alive, dropped


def pt_tri_dist(p, T):
    """p:(3,) 到三角形 T(3,3) 的距离（Ericson 线段/面最近点）。"""
    a, b, c = T
    ab, ac, ap = b - a, c - a, p - a
    d1, d2 = ab @ ap, ac @ ap
    if d1 <= 0 and d2 <= 0:
        return np.linalg.norm(ap)
    bp = p - b
    d3, d4 = ab @ bp, ac @ bp
    if d3 >= 0 and d4 <= d3:
        return np.linalg.norm(bp)
    vc = d1 * d4 - d3 * d2
    if vc <= 0 and d1 >= 0 and d3 <= 0:
        v = d1 / (d1 - d3)
        return np.linalg.norm(a + v * ab - p)
    cp = p - c
    d5, d6 = ab @ cp, ac @ cp
    if d6 >= 0 and d5 <= d6:
        return np.linalg.norm(cp)
    vb = d5 * d2 - d1 * d6
    if vb <= 0 and d2 >= 0 and d6 <= 0:
        w = d2 / (d2 - d6)
        return np.linalg.norm(a + w * ac - p)
    va = d3 * d6 - d5 * d4
    if va <= 0 and (d4 - d3) >= 0 and (d5 - d6) >= 0:
        w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
        return np.linalg.norm(b + w * (c - b) - p)
    denom = 1.0 / (va + vb + vc)
    v, w = vb * denom, vc * denom
    return np.linalg.norm(a + ab * v + ac * w - p)


def _batch_pt_tri(p, T):
    """向量化 point→triangles 距离：p(3,), T(M,3,3) → (M,) 距离 + (M,3) 最近点。
    Ericson《Real-Time Collision Detection》 closest point on triangle，逐分支向量化。"""
    a, b, c = T[:, 0], T[:, 1], T[:, 2]
    ab, ac = b - a, c - a
    ap = p - a
    d1 = np.einsum("ij,ij->i", ab, ap)
    d2 = np.einsum("ij,ij->i", ac, ap)
    bp = p - b
    d3 = np.einsum("ij,ij->i", ab, bp)
    d4 = np.einsum("ij,ij->i", ac, bp)
    cp = p - c
    d5 = np.einsum("ij,ij->i", ab, cp)
    d6 = np.einsum("ij,ij->i", ac, cp)
    vc = d1 * d4 - d3 * d2
    vb = d5 * d2 - d1 * d6
    va = d3 * d6 - d5 * d4

    Q = np.empty_like(a)
    # 区域判定（与标量版同序）
    m_a = (d1 <= 0) & (d2 <= 0)
    m_b = (d3 >= 0) & (d4 <= d3)
    m_ab = (~m_a) & (~m_b) & (vc <= 0) & (d1 >= 0) & (d3 <= 0)
    m_c = (d6 >= 0) & (d5 <= d6)
    m_ac = (~m_a) & (~m_b) & (~m_ab) & (~m_c) & (vb <= 0) & (d2 >= 0) & (d6 <= 0)
    m_bc = (~m_a) & (~m_b) & (~m_ab) & (~m_c) & (~m_ac) & (va <= 0) & \
        ((d4 - d3) >= 0) & ((d5 - d6) >= 0)
    m_face = ~(m_a | m_b | m_ab | m_c | m_ac | m_bc)

    Q[:] = a
    Q[m_b] = b[m_b]
    with np.errstate(divide="ignore", invalid="ignore"):
        v = np.where(d1 - d3 != 0, d1 / np.where(d1 - d3 != 0, d1 - d3, 1), 0)
        w = np.where(d2 - d6 != 0, d2 / np.where(d2 - d6 != 0, d2 - d6, 1), 0)
        denom = va + vb + vc
        den = np.where(denom != 0, denom, 1)
        vv = np.where(denom != 0, vb / den, 0)
        ww = np.where(denom != 0, vc / den, 0)
    Q[m_ab] = (a + v[:, None] * ab)[m_ab]
    Q[m_ac] = (a + w[:, None] * ac)[m_ac]
    Q[m_bc] = (b + w[:, None] * (c - b))[m_bc]
    Q[m_face] = (a + vv[:, None] * ab + ww[:, None] * ac)[m_face]
    # 面区域的投影必须仍在三角形内，否则退回最近边（用标量校正距离）
    dist = np.linalg.norm(Q - p, axis=1)
    d_true = np.array([pt_tri_dist(p, TT) for TT in T]) if m_face.any() else dist
    dist = np.where(m_face, d_true, dist)
    # 最近点：对 face 区域重新取与 p 距离 == dist 的候选
    if m_face.any():
        idx = np.where(m_face)[0]
        for i in idx:
            TT = T[i]
            aa, bb, cc = TT
            cand = [aa, bb, cc]
            for (u, vv2) in ((aa, bb), (bb, cc), (cc, aa)):
                e = vv2 - u
                t = float(np.clip(np.dot(p - u, e) / max(np.dot(e, e), 1e-30), 0, 1))
                cand.append(u + t * e)
            n = np.cross(bb - aa, cc - aa)
            nn = float(np.dot(n, n))
            if nn > 1e-30:
                tt = float(np.dot(p - aa, n)) / nn
                cand.append(p - tt * n)
            for q in cand:
                if abs(np.linalg.norm(q - p) - dist[i]) <= 1e-9:
                    Q[i] = q
                    break
    return dist, Q


def closest_on_tri(p, T):
    """(dist, 最近点) —— 保持原签名。"""
    d, q = _batch_pt_tri(p, T[None, :, :])
    return float(d[0]), q[0]


def pair_report(name, FA, FB, nfacets):
    """FA/FB: list of (P(3,3), n(3))"""
    TA = np.array([f[0] for f in FA])
    NA = np.array([f[1] for f in FA])
    TB = np.array([f[0] for f in FB])
    NB = np.array([f[1] for f in FB])
    CA, CB = TA.mean(1), TB.mean(1)

    def one_way(C, N, TT, NN):
        # ★ N 是**未归一化的面积向量**（cart_outer: Σ cross = 2·area·n̂）
        #   ⇒ 投影/点积前必须归一化，否则 gap 被 |n| 放大（实测 54.9 mm 假值）。
        Nh = N / np.maximum(np.linalg.norm(N, axis=1, keepdims=True), 1e-30)
        dn, dm, dots, coss = [], [], [], []
        for i in range(len(C)):
            d, Q = _batch_pt_tri(C[i], TT)
            k = int(np.argmin(d))
            bq = Q[k]
            dm.append(float(d[k]))
            dn.append(float(np.dot(bq - C[i], Nh[i])))
            ni, nk = Nh[i], NN[k] / max(np.linalg.norm(NN[k]), 1e-30)
            dots.append(float(np.dot(ni, nk)))
            v = bq - C[i]
            coss.append(float(np.dot(ni, v) / max(np.linalg.norm(v), 1e-30)))
        return np.array(dm), np.array(dn), np.array(dots), np.array(coss)

    dm_a, dn_a, dot_a, cos_a = one_way(CA, NA, TB, NB)
    dm_b, dn_b, dot_b, cos_b = one_way(CB, NB, TA, NA)

    def stat(x):
        return dict(min=float(np.min(x)), med=float(np.median(x)),
                    max=float(np.max(x)), p05=float(np.percentile(x, 5)),
                    p95=float(np.percentile(x, 95)))

    def side(dm, dn, dot, cos):
        return {"d_min": stat(dm), "gap_n": stat(dn),
                "n_pen_gt_0.1": int((dn < -0.1).sum()),
                "n_pen_any": int((dn < 0).sum()),
                "facing_neg_frac": float((dot < 0).mean()),
                "cos_med": float(np.median(cos)),
                "n_small_gap_le_1mm": int((dm <= 1.0).sum())}

    rep = {
        "name": name,
        "n_facets_A": int(len(FA)), "n_facets_B": int(len(FB)),
        "n_facets_raw": nfacets,
        "A_to_B": side(dm_a, dn_a, dot_a, cos_a),
        "B_to_A": side(dm_b, dn_b, dot_b, cos_b),
        "facing_dot_med": float(np.median(np.concatenate([dot_a, dot_b]))),
    }
    return rep


def main():
    cg = np.load(CART, allow_pickle=True)
    res = {"cart": str(CART), "pairs": {}}
    JOINTS = ["tibiotalar_tibia", "tibiotalar_talus",
              "subtalar_talus", "subtalar_calcaneus"]
    for j in JOINTS:
        alive, dropped = gate_stats(cg, j, 0.0)
        res.setdefault("gate", {})[j] = {"alive": alive, "dropped_DET_MIN=0.0": dropped}
        print(f"[gate] {j}: 单元存活 {alive} / 被 DET_MIN=0.0 闸门剔除 {dropped}")
    for pair, a, b in (("tibiotalar(正对照)", "tibiotalar_tibia", "tibiotalar_talus"),
                       ("subtalar(目标)", "subtalar_talus", "subtalar_calcaneus")):
        FA = list(outer_facets(cg, a))
        FB = list(outer_facets(cg, b))
        nraw = len(list(outer_facets(cg, a, apply_gate=False)))
        rep = pair_report(pair, FA, FB, nraw)
        res["pairs"][pair] = rep
        print(f"\n=== {pair}: {a}({len(FA)}) ↔ {b}({len(FB)}) ===")
        for side in ("A_to_B", "B_to_A"):
            s = rep[side]
            print(f"  {side}: d_min min/med/max = "
                  f"{s['d_min']['min']:.4f}/{s['d_min']['med']:.4f}/{s['d_min']['max']:.4f} mm"
                  f" | 沿归一化法向间隙 gap_n min/med/max = "
                  f"{s['gap_n']['min']:.4f}/{s['gap_n']['med']:.4f}/{s['gap_n']['max']:.4f} mm"
                  f" | 穿透>0.1mm {s['n_pen_gt_0.1']} / 任何穿透 {s['n_pen_any']}"
                  f" | d_min≤1mm 面片 {s['n_small_gap_le_1mm']}"
                  f" | 法向相对(dot<0) {s['facing_neg_frac']*100:.1f}%"
                  f" | cos_med {s['cos_med']:.3f}")
        print(f"  最近点对法向点积中位 = {rep['facing_dot_med']:.3f}"
              f"（<0 = 面对面）")
    OUTJ.write_text(json.dumps(res, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n-> {OUTJ}")


if __name__ == "__main__":
    main()
