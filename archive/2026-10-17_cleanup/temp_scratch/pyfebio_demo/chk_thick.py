"""检查 cart_outer 取的是不是软骨的**外表面**。

判据：surf_* 节点到自身底面(tiesurf_*)的距离 ≈ 软骨厚度(1~2mm)。
若 ≈ 0 ⇒ cart_outer 取错了面（取到骨侧底面）⇒ 该批节点同时受 tie 与接触约束 ⇒ 必然冲突。
"""
import re
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

p = Path(sys.argv[1] if len(sys.argv) > 1 else
         r"D:\Project\climbing_fall_analysis\temp\thums\febio_lig\ankle_lig_WK.feb")
txt = p.read_text(encoding="latin-1", errors="replace")
XYZ = {}
for m in re.finditer(r'<node id="(\d+)">\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)', txt):
    XYZ[int(m.group(1))] = np.array([float(m.group(2)), float(m.group(3)), float(m.group(4))])
surfs = {}
for m in re.finditer(r'<Surface name="([^"]+)"[^>]*>(.*?)</Surface>', txt, re.S):
    surfs[m.group(1)] = [tuple(int(y) for y in s.split(",")) for s in
                         re.findall(r"<[a-z0-9]+ id=[^>]*>\s*([\d,]+)\s*</", m.group(2))]
print(f"文件 {p.name}  ({p.stat().st_size/1024:.0f} KB)")
for k, v in sorted(surfs.items()):
    ids = sorted({x for f in v for x in f})
    print(f"  {k:34} {len(v):>4} 面 / {len(ids):>4} 节点")

for j in ("tibiotalar_tibia", "tibiotalar_talus"):
    top = sorted({x for f in surfs.get(f"surf_{j}", []) for x in f})
    bot = sorted({x for f in surfs.get(f"tiesurf_{j}", []) for x in f})
    bone = sorted({x for f in surfs.get(f"bonesurf_{j}", []) for x in f})
    print(f"\n=== {j} ===")
    if not top or not bot:
        print("  (面为空，跳过)")
        continue
    Pt = np.array([XYZ[i] for i in top], float)
    tree = cKDTree(np.array([XYZ[i] for i in bot], float))
    d, _ = tree.query(Pt)
    print(f"  surf 节点 -> 自身底面(tiesurf) 距离:  min={d.min():.4f}  中位={np.median(d):.4f}  "
          f"max={d.max():.4f} mm")
    print(f"  ⇒ 若中位 ≈ 软骨厚度(1~2mm) 则 surf 是正确的**外表面**；"
          f"若 ≈ 0 则取错了面")
    if bone:
        Pb = np.array([XYZ[i] for i in bone], float)
        tb = cKDTree(Pb)
        d2, _ = tb.query(Pt)
        print(f"  surf 节点 -> 骨面(bonesurf) 距离:   min={d2.min():.4f}  中位={np.median(d2):.4f} mm")
