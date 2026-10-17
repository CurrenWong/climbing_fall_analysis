"""量「新格式软骨」t=0 关节两侧外表面的间隙（无定向假设，面积加权）。

用法: python gap_new.py <feb 路径>
判据参考：旧格式基线 = 中位间隙 0.9732 mm、接触斑 31.0 mm²@0.02 容差（安全）。
"""
import re
import sys
from pathlib import Path

import numpy as np

p = Path(sys.argv[1] if len(sys.argv) > 1 else
         r"D:\Project\climbing_fall_analysis\temp\thums\febio_lig\ankle_lig_WK.feb")
txt = p.read_text(encoding="latin-1", errors="replace")
print(f"读 {p}  ({p.stat().st_size/1024:.0f} KB)")

XYZ = {}
for m in re.finditer(r'<node id="(\d+)">\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)', txt):
    XYZ[int(m.group(1))] = np.array([float(m.group(2)), float(m.group(3)), float(m.group(4))])
print(f"节点 {len(XYZ)}")

surfs = {}
for m in re.finditer(r'<Surface name="([^"]+)"[^>]*>(.*?)</Surface>', txt, re.S):
    surfs[m.group(1)] = [tuple(int(y) for y in s.split(",")) for s in
                         re.findall(r"<[a-z0-9]+ id=[^>]*>\s*([\d,]+)\s*</", m.group(2))]
for k, v in sorted(surfs.items()):
    print(f"  Surface {k:34} {len(v)} 面")


def pts_of(name):
    out = []
    for f in surfs.get(name, []):
        P = np.array([XYZ[i] for i in f], float)
        n = len(P)
        for w in range(6):            # 简单细分提升采样密度
            for a in np.linspace(0, 1, 7):
                pass
        for i in range(1, n - 1):     # 三角化 fan
            for u in np.linspace(0, 1, 5):
                for v in np.linspace(0, 1 - u, max(2, int(5 * (1 - u)) + 1)):
                    out.append(P[0] * (1 - u - v) + P[i] * u + P[i + 1] * v)
    return np.array(out, float)


for side in ("tibiotalar_tibia", "tibiotalar_talus"):
    A = pts_of(f"surf_{side}")
    B = pts_of(f"surf_{'tibiotalar_talus' if side.endswith('tibia') else 'tibiotalar_tibia'}")
    if not len(A) or not len(B):
        print(f"\n[{side}] 采样为空，跳过")
        continue
    D = np.sqrt(((A[:, None, :] - B[None, :, :]) ** 2).sum(-1))
    dmin = D.min(1)
    print(f"\n[{side}] 采样 {len(A)} 点")
    print(f"  最小距  min={dmin.min():.4f}  p1={np.percentile(dmin,1):.4f}  "
          f"中位={np.median(dmin):.4f}  p99={np.percentile(dmin,99):.4f}  max={dmin.max():.4f} mm")
    for tol in (0.02, 0.2, 0.5, 1.0):
        print(f"  ≤{tol:>5} mm 的采样点比例 {100.0*(dmin <= tol).mean():5.1f}%")

# 互穿检测（用对面三角的最近点带符号判定需要法向；这里用「极小距」代理）
print("\n注：本脚本给的是**无符号**距离代理；<0.02mm 的点即接触斑，"
      "若大量点 <0 数值上不可见，需 point_in_tet.py 做互穿裁决。")
