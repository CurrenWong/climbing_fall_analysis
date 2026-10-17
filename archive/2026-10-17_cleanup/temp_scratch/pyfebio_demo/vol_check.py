"""精确体积校验：三线性 Jacobian 的 2x2x2 Gauss 积分（对 trilinear hex 是精确的）。
对比 原 npz vs 加密后 npz，只看 *_solid。
"""
from pathlib import Path
import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
# ★ 形状函数定义在 u,v,w ∈ [0,1] ⇒ Gauss 点也必须在 [0,1]：
#   2 点 Gauss on [0,1] 的节点是 0.5 ± 0.5/√3，权重各 0.5 ⇒ 每维 Σw = 1
G = 0.5 / np.sqrt(3.0)
GAUSS = [(a, b, c) for a in (0.5 - G, 0.5 + G)
         for b in (0.5 - G, 0.5 + G) for c in (0.5 - G, 0.5 + G)]
HEX_UVW = np.array([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                    (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)], float)


def dN_hex(u, v, w):
    """dN/du, dN/dv, dN/dw（FEBio 局部顺序）。"""
    du = np.array([-(1-v)*(1-w), (1-v)*(1-w), v*(1-w), -v*(1-w),
                   -(1-v)*w, (1-v)*w, v*w, -v*w], float)
    dv = np.array([-(1-u)*(1-w), -u*(1-w), u*(1-w), (1-u)*(1-w),
                   -(1-u)*w, -u*w, u*w, (1-u)*w], float)
    dw = np.array([-(1-u)*(1-v), -u*(1-v), -u*v, -(1-u)*v,
                   (1-u)*(1-v), u*(1-v), u*v, (1-u)*v], float)
    return du, dv, dw


def vol_hex(P):
    V = 0.0
    for u, v, w in GAUSS:
        du, dv, dw = dN_hex(u, v, w)
        J = np.stack([du @ P, dv @ P, dw @ P], axis=0)
        V += abs(np.linalg.det(J))
    return V / 8.0          # Gauss 权重 = 1，共 8 点 -> /8（参数域体积 1）


def vol_tet(P):
    return abs(np.dot(np.cross(P[1] - P[0], P[2] - P[0]), P[3] - P[0]) / 6.0)


def total(npz):
    d = np.load(npz, allow_pickle=True)
    ids = {int(v): i for i, v in enumerate(d["node_ids"])}
    X = d["node_xyz"].astype(float)
    tot = 0.0
    cnt = {"hex8": 0, "tet4": 0, "wedge6": 0}
    for k in d.files:
        if not k.endswith("_solid"):
            continue
        for row in d[k]:
            uq = sorted({int(x) for x in row if x > 0})
            if len(uq) not in (4, 6, 8):
                continue
            P = np.array([X[ids[x]] for x in uq])
            if len(uq) == 8:
                tot += vol_hex(P); cnt["hex8"] += 1
            elif len(uq) == 4:
                tot += vol_tet(P); cnt["tet4"] += 1
            else:
                cnt["wedge6"] += 1
    return tot, cnt


a, ca = total(ROOT / "temp/thums/thums_lowerlimb.npz")
b, cb = total(ROOT / "temp/thums/thums_lowerlimb_r1.npz")
print("=== 精确体积（三线性 Jacobian Gauss 积分） ===")
print(f"  原始   {a:,.4f} mm^3   计数 {ca}")
print(f"  加密后 {b:,.4f} mm^3   计数 {cb}")
print(f"  偏差   {(b-a)/a*100:+.6f}%")
print(f"\n  体积应完全守恒（红色加密是精确细分）。")
