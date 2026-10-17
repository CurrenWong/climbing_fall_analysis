"""P2② 定位坏单元：在哪、是否在关节接触面/tie 面上（决定可否剔除）。"""
import re
from pathlib import Path

import numpy as np

R = Path(r"D:\Project\climbing_fall_analysis")
FEB = R / "temp/thums/febio_lig/ankle_lig_WK.feb"
txt = FEB.read_text(encoding="latin-1", errors="replace")
GP2 = 1 / np.sqrt(3)


def dNw(r, s, t):
    return np.array([-(1-t), -(1-t), -(1-r-s),
                     (1-t), 0, -r,
                     0, (1-t), -s,
                     -t, -t, (1-r-s),
                     t, 0, r,
                     0, t, s]).reshape(6, 3)


PTS = [(1/6, 1/6, .5), (2/3, 1/6, .5), (1/6, 2/3, .5),
       (1/3, 1/3, 0.), (1/3, 1/3, 1.), (1/3, 1/3, .5),
       (0., 0., 0.), (1., 0., 0.), (0., 1., 0.),
       (0., 0., 1.), (1., 0., 1.), (0., 1., 1.)]


def worst_det(P):
    return min(float(np.linalg.det(dNw(*p).T @ P)) for p in PTS)


XYZ, ELEMS = {}, {}
for m in re.finditer(r'<node id="(\d+)">\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,'
                     r'\s*([-\d.eE+]+)', txt):
    XYZ[int(m.group(1))] = np.array([float(m.group(2)), float(m.group(3)),
                                     float(m.group(4))])
for m in re.finditer(r'<Elements\s+type="(\w+)"[^>]*?name="([^"]*)"[^>]*>'
                     r'(.*?)</Elements>', txt, re.S):
    kind, nm, body = m.group(1), m.group(2), m.group(3)
    conn = [[int(x) for x in c.split(",")] for c in
            re.findall(r'<elem[^>]*>\s*([\d,\s]+?)\s*</elem>', body)]
    ELEMS[nm] = (kind, conn)

# 关节面 / tie 面的节点集合
SURF = {}
for m in re.finditer(r'<Surface name="([^"]+)"[^>]*>(.*?)</Surface>', txt, re.S):
    fs = [[int(x) for x in f.split(",")] for f in
          re.findall(r'<[a-z0-9]+ id="\d+">\s*([\d,]+)\s*</', m.group(2))]
    SURF[m.group(1)] = fs
print("Surfaces:", {k: len(v) for k, v in SURF.items()})

for nm in ("tibiotalar_talus__penta6", "tibiotalar_tibia__penta6"):
    if nm not in ELEMS:
        continue
    kind, conn = ELEMS[nm]
    bad = []
    for idx, c in enumerate(conn):
        P = np.array([XYZ[n] for n in c], float)
        w = worst_det(P)
        ctr = P.mean(0)
        if w <= 1e-2:
            bad.append((idx, w, ctr))
    print("\n" + "=" * 96)
    print(f"{nm}：{len(conn)} 单元，worst-det(J) ≤ 1e-2 的 {len(bad)} 个")
    print("=" * 96)
    for idx, w, ctr in sorted(bad, key=lambda x: x[1]):
        print(f"  单元#{idx:<4} worst-det(J)={w:>12.5e}  中心=({ctr[0]:8.2f},"
              f"{ctr[1]:8.2f},{ctr[2]:8.2f})")
    # 是否在接触面/tie 面上
    for sname, fs in SURF.items():
        if "tibiotalar" not in sname:
            continue
        sset = sorted({n for f in fs for n in f})
        hit = 0
        for idx, w, _ in bad:
            c = conn[idx]
            if any(n in sset for n in c):
                hit += 1
        print(f"    → 与 {sname} 共节点的坏单元数：{hit}/{len(bad)}")

# 体积亏损估算
print("\n" + "=" * 96)
print("整体质量：所有软骨单元的 worst-det(J) 分布")
print("=" * 96)
allw = []
for nm in ("tibiotalar_talus__penta6", "tibiotalar_tibia__penta6"):
    kind, conn = ELEMS[nm]
    for c in conn:
        P = np.array([XYZ[n] for n in c], float)
        allw.append(worst_det(P))
allw = np.array(allw)
for q in (0, 0.5, 1, 5, 25, 50, 100):
    print(f"  p{q:<4} = {np.percentile(allw, q):>12.5e}")
print(f"  ≤0 : {(allw<=0).sum()}   ≤1e-2 : {(allw<=1e-2).sum()}   "
      f"≤0.1 : {(allw<=0.1).sum()}   ≤0.5 : {(allw<=0.5).sum()}")
