"""正确的面片面积（shoelace，非位置向量叉积）+ 面片质量体检：
   tiesurf/bonesurf/surf 的最小面积、退化面、重复面。"""
import re
from pathlib import Path

import numpy as np

feb = Path(r"D:\Project\climbing_fall_analysis\temp/thums/febio_lig/ankle_lig_WK.feb")
t = feb.read_text(encoding="latin-1")
xyz = {}
for m in re.finditer(r'<node id="(\d+)">([-0-9.eE+,]+)</node>', t):
    xyz[int(m.group(1))] = np.array([float(x) for x in m.group(2).split(",")])


def area(f):
    P = [xyz[i] for i in f]
    a = np.zeros(3)
    for i in range(len(P)):
        a += np.cross(P[i], P[(i + 1) % len(P)])
    return 0.5 * float(np.linalg.norm(a))


for m in re.finditer(r'<Surface name="([^"]+)">(.*?)</Surface>', t, re.S):
    name, body = m.groups()
    ids = [[int(x) for x in g.split(",")]
           for g in re.findall(r'<(?:tri3|quad4) id="\d+">([0-9,]+)</(?:tri3|quad4)>',
                               body)]
    a = np.array([area(f) for f in ids]) if ids else np.array([0.0])
    dup = len(ids) - len({tuple(sorted(f)) for f in ids})
    degen = int((a < 1e-6).sum())
    tiny = int((a < 1e-3).sum())
    flag = "  🚨" if degen or tiny else ""
    print(f"{name:34s} n={len(ids):4d} area min={a.min():.3e} "
          f"p1={np.percentile(a,1):.3e} med={np.median(a):.3e} "
          f"dup={dup} <1e-6:{degen} <1e-3:{tiny}{flag}")
