"""负 vs 正 hex：底面 quad 的平面外翘曲有多大？（薄板单元翻负的经典机制是翘曲）

薄板理论：矩形盒单元 det(J) 恒正（a·b·c/8），所以翻负必然来自**翘曲/扭曲**。
量化：warp = 角点到最佳拟合平面的最大偏离 / 对角线长。
"""
import numpy as np
from pathlib import Path

THUMS = Path(r"D:\Project\climbing_fall_analysis\temp\thums")
src = np.load(THUMS / "cartilage_mesh.npz", allow_pickle=True)
new = np.load(THUMS / "cartilage_mesh_hex.npz", allow_pickle=True)
bpos = {int(v): src["bone_xyz"][i] for i, v in enumerate(src["bone_node_ids"])}

G = 1 / np.sqrt(3)
GP = [(-G, -G, -G), (G, -G, -G), (G, G, -G), (-G, G, -G),
      (-G, -G, G), (G, -G, G), (G, G, G), (-G, G, G)]
REF = [(-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1),
       (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1)]


def detJ_min(P):
    w = 1e30
    for s in GP:
        dN = np.zeros((8, 3))
        for i in range(8):
            a, b, c = REF[i]
            dN[i, 0] = a * (1 + b * s[1]) * (1 + c * s[2]) / 8
            dN[i, 1] = b * (1 + a * s[0]) * (1 + c * s[2]) / 8
            dN[i, 2] = c * (1 + a * s[0]) * (1 + b * s[1]) / 8
        w = min(w, float(np.linalg.det(P.T @ dN)))
    return w


def newell(p):
    n = np.zeros(3)
    for i in range(len(p)):
        n += np.cross(p[i], p[(i + 1) % len(p)])
    return n


def warp4(p):
    """四边形的平面外翘曲：角点到最佳拟合平面的最大距离 / 对角线。"""
    c = p.mean(axis=0)
    n = newell(p)
    n = n / (np.linalg.norm(n) + 1e-15)
    d = np.abs((p - c) @ n)
    diag = 0.5 * (np.linalg.norm(p[2] - p[0]) + np.linalg.norm(p[3] - p[1]))
    return float(d.max() / (diag + 1e-15))


def angs4(p):
    a = []
    for k in range(4):
        u = p[(k - 1) % 4] - p[k]
        w = p[(k + 1) % 4] - p[k]
        c = np.dot(u, w) / (np.linalg.norm(u) * np.linalg.norm(w) + 1e-15)
        a.append(np.degrees(np.arccos(np.clip(c, -1, 1))))
    return a


for j in ["tibiotalar_tibia", "subtalar_calcaneus"]:
    nd = new[f"{j}_nodes"]
    E = new[f"{j}_elems"]

    def P_of(r):
        return nd[-r - 1] if r < 0 else bpos[r]

    negr, posr = [], []
    for row in E:
        blk = [int(x) for x in row if int(x) != 0]
        m = len(blk) // 2
        if m != 4:
            continue
        P = np.array([P_of(r) for r in blk])
        pb, pt = P[:4], P[4:]
        # ★ 先按挤出方向把底面调成 CCW（否则绕向会污染所有量）
        ext = pt.mean(axis=0) - pb.mean(axis=0)
        if np.dot(newell(pb), ext) < 0:
            P = P[[3, 2, 1, 0, 7, 6, 5, 4]]
            pb, pt = P[:4], P[4:]
        dj = detJ_min(P)
        diag = 0.5 * (np.linalg.norm(pb[2] - pb[0]) + np.linalg.norm(pb[3] - pb[1]))
        h = float(np.mean([np.linalg.norm(pt[k] - pb[k]) for k in range(4)]))
        rec = dict(wb=warp4(pb), wt=warp4(pt), diag=diag, h=h,
                   amin=min(angs4(pb)), dj=dj,
                   wr=warp4(pb) / max(h / diag, 1e-9))
        (negr if dj <= 0 else posr).append(rec)

    print(f"\n{'='*80}\n[{j}]  绕向修正后  负 {len(negr)} / 正 {len(posr)}\n{'='*80}")
    for tag, arr in (("负", negr), ("正", posr)):
        if not arr:
            continue
        f = lambda k: np.median([r[k] for r in arr])
        print(f"  {tag}: 底面翘曲/对角={f('wb'):.4f}  顶面翘曲={f('wt'):.4f}  "
              f"翘曲/厚度比={f('wr'):.3f}")
        print(f"      对角={f('diag'):.3f}mm  厚={f('h'):.3f}mm  "
              f"厚/对角={f('h')/f('diag'):.4f}  底面最小内角={f('amin'):.1f}°")
    # 判决指标
    if negr and posr:
        wn = np.median([r["wb"] for r in negr])
        wp = np.median([r["wb"] for r in posr])
        print(f"  ⇒ 底面翘曲 负/正 = {wn/max(wp,1e-12):.2f}×")
        rn = np.median([r["wr"] for r in negr])
        rp = np.median([r["wr"] for r in posr])
        print(f"  ⇒ (翘曲/厚度) 负/正 = {rn/max(rp,1e-12):.2f}×")
