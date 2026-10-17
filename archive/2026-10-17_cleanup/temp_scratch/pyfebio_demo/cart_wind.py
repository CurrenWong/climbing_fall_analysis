"""软骨: 诊断+修复 三角形绕向一致性, 再走 bones 验证过的 plane-surface 路径成实体。

假设: OCC 判 18 条边界边, 不是"开口"也不是"非流形", 而是**相邻三角形绕向不一致**
      (OCC 要求共享边被两个面反向遍历)。VTK clean 只合并点, 不修绕向。
"""

import os
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

SDK = Path(r"D:\Program\gmsh-sdk\gmsh-4.15.2-Windows64-sdk")
LIB = SDK / "lib"
os.environ["PATH"] = f"{LIB};{SDK / 'bin'};{os.environ['PATH']}"
sys.path.insert(0, str(LIB))

ROOT = Path(__file__).resolve().parents[2]
GEOM = ROOT / "temp" / "ankle_geom"
STL = GEOM / "tibial_talus_cartilage_r.stl"

import gmsh  # noqa: E402
import pyvista as pv  # noqa: E402

T0 = time.time()


def log(m):
    print(f"[{time.time()-T0:5.1f}s] {m}", flush=True)


def winding_report(points, faces, tag):
    """统计每条共享边的两个三角形是否反向遍历 (一致绕向的必要条件)。"""
    d = {}
    for tri in faces:
        for a, b in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
            k = (min(a, b), max(a, b))
            d.setdefault(k, []).append((int(a), int(b)))
    opp = same = multi = 0
    for k, dirs in d.items():
        if len(dirs) > 2:
            multi += 1
        elif len(dirs) == 2:
            if dirs[0] == (dirs[1][1], dirs[1][0]):
                opp += 1
            else:
                same += 1
    log(f"{tag}: 边 {len(d)}  反向配对 {opp}  同向(绕向冲突) {same}  非流形边 {multi}")
    return same, multi


# ---------- 1. 原始 ----------
raw = pv.read(str(STL))
faces = raw.faces.reshape(-1, 4)[:, 1:]
winding_report(raw.points, faces, "原始")

# ---------- 2. clean (合点+删退化) ----------
cleaned = raw.clean(point_merging=True, tolerance=1e-6, absolute=True)
f2 = cleaned.faces.reshape(-1, 4)[:, 1:]
winding_report(cleaned.points, f2, "clean")

# ---------- 3. compute_normals 一致性修正 ----------
try:
    oriented = cleaned.compute_normals(
        consistent_normals=True, auto_orient_normals=True, split_vertices=False, inplace=False
    )
    f3 = oriented.faces.reshape(-1, 4)[:, 1:]
    winding_report(oriented.points, f3, "compute_normals(consistent)")
    SRC = oriented
    src_name = "_cart_oriented.stl"
except Exception as e:
    log(f"compute_normals FAIL: {str(e)[:90]}")
    SRC = cleaned
    src_name = "_cart_fixed.stl"

out = GEOM / src_name
SRC.save(str(out))
log(f"存 {src_name}  ({SRC.n_points} 点 / {SRC.n_cells} 面, vol={SRC.volume*1e6:.4f} cm³)")

# ---------- 4. plane-surface 路径 (bones 验证过, 0.00% 误差) ----------
gmsh.initialize(str(LIB))
gmsh.option.setNumber("General.Terminal", 0)
gmsh.option.setNumber("General.Verbosity", 0)
gmsh.open(str(out))
et, _, en = gmsh.model.mesh.getElements(2)
tri = next(np.array(n).reshape(-1, 3) for x, n in zip(et, en) if x == 2)
ntags, ncoord, _ = gmsh.model.mesh.getNodes()
ncoord = np.array(ncoord, dtype=np.float64).reshape(-1, 3)
tag2row = {int(x): i for i, x in enumerate(ntags)}
rows = np.vectorize(tag2row.get)(tri)
log(f"gmsh 读到 {len(ntags)} 点 {len(tri)} 三角")

gmsh.model.add("cart")
occ = gmsh.model.occ
pts = [occ.addPoint(*ncoord[i], tag=int(i + 1)) for i in range(len(ntags))]
e2c = {}
surfs = []
for a, b, c in rows:
    ta, tb, tc = pts[a], pts[b], pts[c]
    cur = []
    for u, v in ((ta, tb), (tb, tc), (tc, ta)):
        k = (min(u, v), max(u, v))
        t = e2c.get(k)
        if t is None:
            t = occ.addLine(u, v)
            e2c[k] = t
        cur.append(t)
    surfs.append(occ.addPlaneSurface([occ.addCurveLoop(cur)]))
occ.synchronize()
bnd = gmsh.model.getBoundary([(2, t) for d, t in gmsh.model.getEntities(2)])
log(f"边界边 {len(bnd)}  {'✅ 封闭!' if len(bnd) == 0 else '❌'}")

if len(bnd) == 0:
    sl = occ.addSurfaceLoop([t for d, t in gmsh.model.getEntities(2)])
    occ.addVolume([sl])
    occ.synchronize()
    v3 = gmsh.model.getEntities(3)
    if v3:
        log(f"实体体积 = {float(occ.getMass(3, v3[0][1]))*1e6:.4f} cm³")
    for cl in (0.002,):
        gmsh.option.setNumber("Mesh.MeshSizeMax", cl)
        gmsh.option.setNumber("Mesh.MeshSizeMin", cl * 0.4)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.optimize("Netgen")
        gmsh.write(str(ROOT / "temp" / "ankle_mesh" / f"cart_{int(cl*1000)}.msh"))
        log(f"char_len={cl*1000:.0f}mm -> 已写 cart_{int(cl*1000)}.msh")
gmsh.finalize()
log("完成")
