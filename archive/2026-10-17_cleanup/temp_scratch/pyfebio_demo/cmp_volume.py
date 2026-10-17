"""直接比两条离散路径切出的**软骨体积**——若不等，就是几何不一致。

tet4 路径：三角底行 → WEDGE_T 拆 3 tet；四边底行 → HEX_T 拆 6 tet
prism 路径：三角底行 → 1 penta6；四边底行 → 拆 2 penta6
两者应覆盖同一体积。
"""
import numpy as np
from pathlib import Path

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums")
cg = np.load(W / "cartilage_mesh.npz", allow_pickle=True)
nid = cg["bone_node_ids"]
xyz = np.asarray(cg["bone_xyz"], float)
idx = {int(v): i for i, v in enumerate(nid)}

WEDGE_T = [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]
HEX_T = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6),
         (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]


def tet_vol(pts):
    return abs(np.linalg.det(np.stack([pts[1] - pts[0], pts[2] - pts[0],
                                       pts[3] - pts[0]]))) / 6.0


for j in ("tibiotalar_tibia", "tibiotalar_talus"):
    E = cg[f"{j}_elems"]
    O = np.asarray(cg[f"{j}_offsets"], float)
    v_tet = v_pri = 0.0
    n3 = n4 = 0
    for row in E:
        vals = [int(x) for x in row if int(x) != 0]
        if not vals or vals[0] <= 0:
            continue
        m = len(vals) // 2
        if any(b not in idx for b in vals[:m]):
            continue
        B = [xyz[idx[b]] for b in vals[:m]]
        T = [O[-x - 1] for x in vals[m:]]
        full = np.array(B + T)
        if m == 3:
            n3 += 1
            for t in WEDGE_T:
                v_tet += tet_vol(full[list(t)])
            v_pri += abs(np.linalg.det(np.stack(
                [full[1] - full[0], full[2] - full[0],
                 (full[3] + full[4] + full[5] - full[0] - full[1]
                  - full[2]) / 3.0]))) / 2.0
        else:
            n4 += 1
            for t in HEX_T:
                v_tet += tet_vol(full[list(t)])
            for tri in ((0, 1, 2), (0, 2, 3)):
                P = np.array([full[z] for z in tri]
                             + [full[z + m] for z in tri])
                v_pri += abs(np.linalg.det(np.stack(
                    [P[1] - P[0], P[2] - P[0],
                     (P[3] + P[4] + P[5] - P[0] - P[1] - P[2]) / 3.0]))) / 2.0
    print(f"\n{'='*70}\n{j}（三角底 {n3}，四边底 {n4}）")
    print(f"  tet4 路径体积   = {v_tet:10.3f} mm³")
    print(f"  prism 路径体积  = {v_pri:10.3f} mm³")
    print(f"  差              = {v_pri - v_tet:+10.3f} mm³  "
          f"({(v_pri/v_tet-1)*100:+.2f}%)")
