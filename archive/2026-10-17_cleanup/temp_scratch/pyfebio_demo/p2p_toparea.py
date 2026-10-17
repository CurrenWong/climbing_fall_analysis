"""P2' 判定：tibia_top 面的真实面积 vs 投影面积（对照 47.8 mm²）。"""
import re
from pathlib import Path

import numpy as np

R = Path(r"D:\Project\climbing_fall_analysis")
FEB = R / "temp/thums/febio_lig/ankle_lig_WK.feb"
txt = FEB.read_text(encoding="latin-1", errors="replace")

XYZ = {}
for m in re.finditer(r'<node id="(\d+)">\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,'
                     r'\s*([-\d.eE+]+)', txt):
    XYZ[int(m.group(1))] = np.array([float(m.group(2)), float(m.group(3)),
                                     float(m.group(4))])

AXIS = np.array([0.0, 0.0, 1.0])   # ALIGN_AXIS=1 已把加载轴转到 +z


def surf(name):
    m = re.search(r'<Surface name="' + name + r'"[^>]*>(.*?)</Surface>', txt, re.S)
    if not m:
        return []
    return [[int(x) for x in f.split(",")] for f in
            re.findall(r'<[a-z0-9]+ id="\d+">\s*([\d,]+)\s*</', m.group(1))]


print("=" * 96)
print("各 Surface：三角数 / 真实面积 / 投影到加载轴(z)的面积 / 平均法向·z")
print("=" * 96)
print(f"{'Surface':<30}{'面':>6}{'真实面积':>12}{'投影面积':>12}"
      f"{'|n·z|均值':>12}{'真实/投影':>11}")
for nm in ("tibia_top", "fibula_top", "surf_tibiotalar_tibia",
           "surf_tibiotalar_talus", "bonesurf_tibiotalar_tibia",
           "bonesurf_tibiotalar_talus"):
    fs = surf(nm)
    if not fs:
        continue
    A = 0.0
    Ap = 0.0
    cosn = []
    for f in fs:
        try:
            P = np.array([XYZ[n] for n in f], float)
        except KeyError:
            continue
        if len(P) < 3:
            continue
        # 三角化扇形
        for i in range(1, len(P) - 1):
            a, b, c = P[0], P[i], P[i + 1]
            n = np.cross(b - a, c - a)
            ar = np.linalg.norm(n) / 2
            if ar <= 0:
                continue
            A += ar
            c_ = abs(float(n @ AXIS)) / (2 * ar) if ar else 0
            cosn.append(c_)
            Ap += ar * c_
    print(f"{nm:<30}{len(fs):>6}{A:>12.1f}{Ap:>12.1f}"
          f"{np.mean(cosn):>12.3f}{(A/Ap if Ap else float('nan')):>11.2f}")

print("\n" + "=" * 96)
print("对照")
print("=" * 96)
print("  [模型输出] 顶面投影面积：胫骨 47.8 mm²，腓骨 12.7 mm²")
print("  [截面积参考] 胫骨远端截面 ≈ 250~450 mm²")
print()
fa = surf("tibia_top")
P = []
for f in fa:
    P.extend([XYZ[n] for n in f if n in XYZ])
P = np.array(P)
print(f"  tibia_top 节点 z 范围 {P[:,2].min():.3f} .. {P[:,2].max():.3f} mm")
print(f"  tibia_top 节点 xy 范围 "
      f"x {P[:,0].min():.1f}..{P[:,0].max():.1f}  "
      f"y {P[:,1].min():.1f}..{P[:,1].max():.1f}")
ext = (P[:, 0].max()-P[:, 0].min()) * (P[:, 1].max()-P[:, 1].min())
print(f"  xy 包围盒面积 = {ext:.1f} mm²（凸包会略小）")
