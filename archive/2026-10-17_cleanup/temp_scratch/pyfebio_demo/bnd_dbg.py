"""定位软骨的 18 条 OCC 边界边: 坐标、长度、以及和 VTK 判定的差异。"""

import os
import sys
from collections import Counter
from pathlib import Path

import numpy as np

SDK = Path(r"D:\Program\gmsh-sdk\gmsh-4.15.2-Windows64-sdk")
LIB = SDK / "lib"
os.environ["PATH"] = f"{LIB};{SDK/'bin'};{os.environ['PATH']}"
sys.path.insert(0, str(LIB))

ROOT = Path(__file__).resolve().parents[2]
GEOM = ROOT / "temp" / "ankle_geom"
STL = GEOM / "tibial_talus_cartilage_r.stl"

import gmsh  # noqa: E402
import pyvista as pv  # noqa: E402

# --- VTK 侧精确统计 ---
raw = pv.read(str(STL))
print(f"VTK: {raw.n_points} 点 {raw.n_cells} 面")
faces = raw.faces.reshape(-1, 4)[:, 1:]
e = np.sort(np.vstack([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1)
u, c = np.unique(e, axis=0, return_counts=True)
print(f"VTK 边计数: {dict(sorted(Counter(c.tolist()).items()))}  非2次的: {(c != 2).sum()}")

# --- OCC 侧 ---
gmsh.initialize(str(LIB))
gmsh.option.setNumber("General.Terminal", 0)
gmsh.option.setNumber("General.Verbosity", 0)
gmsh.open(str(STL))
et, _, en = gmsh.model.mesh.getElements(2)
tri = next(np.array(n).reshape(-1, 3) for x, n in zip(et, en) if x == 2)
ntags, ncoord, _ = gmsh.model.mesh.getNodes()
ncoord = np.array(ncoord, dtype=np.float64).reshape(-1, 3)
tag2row = {int(x): i for i, x in enumerate(ntags)}
rows = np.vectorize(tag2row.get)(tri)

gmsh.model.add("chk")
occ = gmsh.model.occ
pts = [occ.addPoint(*ncoord[i], tag=int(i + 1)) for i in range(len(ntags))]

# 用 VTK 的行号建边, 直接统计
rows_edges = np.sort(
    np.vstack([rows[:, [0, 1]], rows[:, [1, 2]], rows[:, [2, 0]]]), axis=1
)
uu, cc = np.unique(rows_edges, axis=0, return_counts=True)
print(f"\n行号空间边计数: {dict(sorted(Counter(cc.tolist()).items()))}  非2次: {(cc!=2).sum()}")

e2c = {}
surfs = []
for a, b, c_ in rows:
    ta, tb, tc = pts[a], pts[b], pts[c_]
    cur = []
    for u_, v_ in ((ta, tb), (tb, tc), (tc, ta)):
        k = (min(u_, v_), max(u_, v_))
        t = e2c.get(k)
        if t is None:
            t = occ.addLine(u_, v_)
            e2c[k] = t
        cur.append(t)
    surfs.append(occ.addPlaneSurface([occ.addCurveLoop(cur)]))
occ.synchronize()
bnd = gmsh.model.getBoundary([(2, x) for x in surfs])
print(f"\nOCC 边界边: {len(bnd)} 条")

# OCC 的 getBoundary 返回 (dim, tag); 拿到对应 curve -> 两端点
print("\n边界边的几何位置:")
lens = []
for dim, tag in bnd:
    try:
        bpts = gmsh.model.getBoundary([(dim, tag)], combined=False, oriented=False, recursive=True)
    except Exception:
        bpts = []
    # 直接用 getValue 拿 curve 端点
    try:
        lo = gmsh.model.getValue(dim, tag, [])  # 可能不适用
    except Exception:
        lo = None
    lens.append(tag)

# 用 OCC curve 的参数范围拿端点坐标
print(f"  curve tags: {lens}")
for tag in lens[:20]:
    try:
        xmin = gmsh.model.getParametrization(dim, tag)
    except Exception:
        pass

# 直接: curve 1D 的两个端点 = 其 boundary points
print("\n逐条看边界 curve 的两端点:")
for dim, tag in bnd:
    try:
        ends = gmsh.model.getBoundary([(dim, tag)], combined=True, oriented=False, recursive=True)
        # ends 是 (1, p1tag) 和 (1, p2tag)
        pts_tag = [t for d, t in ends if d == 0]
        if len(pts_tag) == 2:
            a = occ.getMass(0, pts_tag[0])  # 面积=0, 不好取
    except Exception:
        pass
    # 简化: 用 model.getEntity 拿 boundbox
    try:
        bb = gmsh.model.getBoundingBox(0, int(tag))
        print(f"    curve {tag}: {np.round(bb, 5)}")
    except Exception as e:
        print(f"    curve {tag}: err {str(e)[:30]}")

gmsh.finalize()
