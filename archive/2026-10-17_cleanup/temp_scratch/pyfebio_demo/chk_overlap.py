"""检查接触面(surf_*) 与 tie 面(tiesurf_*) 是否重叠 / 各面是否朝外。

若 surf_* 与 tiesurf_* 共用面 ⇒ 同一节点同时被 tie 和接触约束 ⇒ 冲突。
"""
import re
import sys
from pathlib import Path

import numpy as np

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

for j in ("tibiotalar_tibia", "tibiotalar_talus"):
    A = {tuple(sorted(f)) for f in surfs.get(f"surf_{j}", [])}
    B = {tuple(sorted(f)) for f in surfs.get(f"tiesurf_{j}", [])}
    inter = A & B
    na = {x for f in surfs.get(f"surf_{j}", []) for x in f}
    nb = {x for f in surfs.get(f"tiesurf_{j}", []) for x in f}
    print(f"\n=== {j} ===")
    print(f"  surf 面 {len(A)} / tiesurf 面 {len(B)} / **共用面 {len(inter)}**")
    print(f"  surf 节点 {len(na)} / tiesurf 节点 {len(nb)} / **共用节点 {len(na & nb)}**")
    # 节点到骨的距离检查：tie 面节点应贴骨面，surf 面节点应在外
    for nm in (f"surf_{j}", f"tiesurf_{j}"):
        ids = sorted({x for f in surfs.get(nm, []) for x in f})
        P = np.array([XYZ[i] for i in ids], float)
        print(f"    {nm:34} 节点 {len(ids):>5}  bbox "
              f"{np.round(P.min(0),1)} .. {np.round(P.max(0),1)}")
