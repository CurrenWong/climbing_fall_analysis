"""看 FEBio 首次翻转的元件 3/4/5（tibia_cort tet4）到底是哪些节点、在哪。"""
import re
from pathlib import Path

import numpy as np

feb = Path(r"D:\Project\climbing_fall_analysis\temp/thums/febio_lig/ankle_lig_WK.feb")
t = feb.read_text(encoding="latin-1")
xyz = {}
for m in re.finditer(r'<node id="(\d+)">([-0-9.eE+,]+)</node>', t):
    xyz[int(m.group(1))] = np.array([float(x) for x in m.group(2).split(",")])
# 面
node2el = {}
for m in re.finditer(r'<Elements type="(\w+)" name="([^"]+)">(.*?)</Elements>',
                     t, re.S):
    et, name, body = m.groups()
    for e in re.finditer(r'<elem id="(\d+)">([0-9,]+)</elem>', body):
        eid, conn = int(e.group(1)), [int(x) for x in e.group(2).split(",")]
        node2el.setdefault(name, {})[eid] = conn

for eid in (1, 2, 3, 4, 5, 6):
    conn = node2el["tibia_cort__tet4"][eid]
    P = np.array([xyz[c] for c in conn])
    det = np.linalg.det(np.stack([P[1] - P[0], P[2] - P[0], P[3] - P[0]]))
    print(f"elem {eid}: nodes {conn}")
    print(f"   坐标:\n{np.round(P, 3)}")
    print(f"   det = {det:.4f}   面积/边长: "
          f"{[round(float(np.linalg.norm(P[(i+1)%4]-P[i])),3) for i in range(4)]}")
    # 这些节点还被哪些单元用
    others = set()
    for dn, d in node2el.items():
        for e2, c2 in d.items():
            if set(c2) & set(conn) and dn != "tibia_cort__tet4":
                others.add(dn)
    print(f"   与之共享节点的其它域: {sorted(others)}")
