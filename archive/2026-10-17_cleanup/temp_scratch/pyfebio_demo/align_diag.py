"""刚体对齐诊断：找使关节间隙平方和最小的相对刚体位姿。
若对齐后残余中位间隙→0 ⇒ 只需重新摆位（Anderson §3.3）；若仍大 ⇒ 形状真不匹配。"""
import re
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from closest_vec import closest_vec

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio")
AX = np.array([0.8153, 0.135, 0.5631])
AX = AX / np.linalg.norm(AX)

txt = (W / "_B.feb").read_text(encoding="latin-1")
XYZ = {}
for m in re.finditer(
        r'<node id="(\d+)">\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)', txt):
    XYZ[int(m.group(1))] = np.array([float(m.group(2)), float(m.group(3)),
                                     float(m.group(4))])
sf = {}
for m in re.finditer(r'<Surface name="([^"]+)"[^>]*>(.*?)</Surface>', txt, re.S):
    sf[m.group(1)] = [tuple(int(y) for y in s.split(",")) for s in
                      re.findall(r"<[a-z0-9]+ id=[^>]*>\s*([\d,]+)\s*</", m.group(2))]
els = {}
for m in re.finditer(
        r'<Elements\b[^>]*type="(tet4|hex8|penta6)"[^>]*name="([^"]+)"[^>]*>(.*?)</Elements>',
        txt, re.S):
    els[m.group(2)] = {int(y) for s in
                       re.findall(r"<elem\b[^>]*>\s*([\d,]+)\s*</elem>", m.group(3))
                       for y in s.split(",")}

P = np.array([XYZ[i] for i in sorted({x for f in sf["surf_tibiotalar_tibia"] for x in f})])
T = np.array([np.array([XYZ[i] for i in f]) for f in sf["surf_tibiotalar_talus"]])
talC = np.array([XYZ[i] for i in els["talus_cort__hex8"]]).mean(0)

# 三角朝向：背离距骨质心
nr = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0])
cen = T.mean(1)
flip = np.einsum("ij,ij->i", nr, cen - talC) < 0
T[flip] = T[flip][:, ::-1]
nr = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0])
nr /= np.linalg.norm(nr, axis=1, keepdims=True)


def rod(r):
    th = np.linalg.norm(r)
    if th < 1e-12:
        return np.eye(3)
    k = r / th
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * (K @ K)


def signed_gap(q):
    R = rod(q[:3])
    tr = T @ R.T + q[3:6]
    d, ti, foot, nm = closest_vec(P, tr)
    s = np.einsum("ij,ij->i", P - foot, nm)
    return s, d


def cost(q):
    s, _ = signed_gap(q)
    return float(np.sum(s ** 2))


print("=" * 78)
print("初始（未对齐）")
print("=" * 78)
s0, d0 = signed_gap(np.zeros(6))
print(f"  |gap| 中位 {np.median(np.abs(s0)):.4f}  max {np.abs(s0).max():.4f}  "
      f"负值(侵入) {int((s0 < 0).sum())}/{len(s0)}")
print(f"  惯性距离 中位 {np.median(d0):.4f}  <0.05mm {int((d0<0.05).sum())}  "
      f"<0.5mm {int((d0<0.5).sum())}")

print("\n" + "=" * 78)
print("优化相对刚体位姿（6 自由度，最小化 Σgap²）")
print("=" * 78)
best, bcost = None, np.inf
for seed in (np.zeros(6),
             np.array([0.0, 0, 0, 0, 0, 1.0]),
             np.array([0.05, 0, 0, 0, 0, -1.0]),
             np.array([-0.05, 0.0, 0.0, 0.0, 0.0, 2.0])):
    r = minimize(cost, seed, method="Powell",
                 options=dict(maxiter=4000, xtol=1e-6, ftol=1e-6))
    if r.fun < bcost:
        bcost, best = r.fun, r.x.copy()
q = best
R = rod(q[:3])
ang = np.degrees(np.linalg.norm(q[:3]))
s1, d1 = signed_gap(q)
print(f"  最优：旋转 {ang:.3f}°  轴 {np.round(q[:3]/(np.linalg.norm(q[:3]) or 1),3)}"
      f"  平移 {np.round(q[3:6],4)} mm")
print(f"  Σgap² : {cost(np.zeros(6)):.4f} → {bcost:.4f}")
print(f"\n{'':10s} {'中位|gap|':>11} {'p95|gap|':>10} {'max|gap|':>10} "
      f"{'<0.05mm':>9} {'<0.1mm':>8} {'<0.5mm':>8} {'侵入数':>7}")
for nm, s in (("未对齐", s0), ("对齐后", s1)):
    print(f"  {nm:10s} {np.median(np.abs(s)):>11.4f} "
          f"{np.percentile(np.abs(s),95):>10.4f} {np.abs(s).max():>10.4f} "
          f"{int((np.abs(s)<0.05).sum()):>9d} {int((np.abs(s)<0.1).sum()):>8d} "
          f"{int((np.abs(s)<0.5).sum()):>8d} {int((s<0).sum()):>7d}")

print("\n" + "=" * 78)
print("判据")
print("=" * 78)
res = np.median(np.abs(s1))
print(f"  对齐后残余中位 |gap| = {res:.4f} mm")
print(f"  {'⇒ ✅ 姿态问题（对齐即可闭合）' if res < 0.2 else '⇒ ❌ 形状真不匹配（对齐救不了）'}")
print(f"  参照：软骨厚度 1.5mm；Anderson 用 1.5~1.7mm")
