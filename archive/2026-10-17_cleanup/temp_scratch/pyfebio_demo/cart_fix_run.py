"""用 VTK clean 修复非流形后的软骨 -> gmsh BREP -> 四面体。验证体积。"""

import os
import sys
import time
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

T0 = time.time()


def log(m):
    print(f"  [{time.time()-T0:6.1f}s] {m}", flush=True)


# --- 1. VTK clean 修非流形, 存成新 STL ---
raw = pv.read(str(STL))
cleaned = raw.clean(point_merging=True, tolerance=1e-6, absolute=True)
log(f"clean -> {cleaned.n_points} 点 {cleaned.n_cells} 面 vol={cleaned.volume*1e6:.4f} cm³")
out = GEOM / "_cart_fixed.stl"
cleaned.save(str(out))
log(f"存 {out.name} ({out.stat().st_size:,} B)")

# --- 2. 送 gmsh 建 BREP ---
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
gmsh.option.setNumber("Geometry.Tolerance", 1e-8)
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
log(f"BREP 面 {len(surfs)}")

bnd = gmsh.model.getBoundary([(2, t) for d, t in gmsh.model.getEntities(2)])
log(f"边界边 {len(bnd)}  {'✅ 封闭!' if len(bnd)==0 else '❌'}")

if len(bnd) == 0:
    sl = occ.addSurfaceLoop([t for d, t in gmsh.model.getEntities(2)])
    occ.addVolume([sl])
    occ.synchronize()
    v3 = gmsh.model.getEntities(3)
    if v3:
        vol = float(occ.getMass(3, v3[0][1]))
        log(f"✅ 实体体积 = {vol*1e6:.4f} cm³  (VTK clean 体积 {cleaned.volume*1e6:.4f})")

    # 四面体网格化
    for cl in (0.002, 0.001):
        gmsh.option.setNumber("Mesh.MeshSizeMax", cl)
        gmsh.option.setNumber("Mesh.MeshSizeMin", cl * 0.4)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)
        gmsh.model.mesh.clear()
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.optimize("Netgen")
        gt, gc, _ = gmsh.model.mesh.getNodes()
        gc = np.array(gc, dtype=np.float64).reshape(-1, 3)
        t3, tg3, en3 = gmsh.model.mesh.getElements(3)
        nt = 0
        tets = None
        rmap = {int(x): i for i, x in enumerate(gt)}
        for e_, n_ in zip(t3, en3):
            if e_ == 4:
                tets = np.vectorize(rmap.get)(np.array(n_).reshape(-1, 4)).astype(np.int64)
                nt = len(tets)
                break
        if tets is not None:
            q = gc[tets]
            tv = float(
                np.abs(
                    np.einsum("ij,ij->i", q[:, 1] - q[:, 0],
                              np.cross(q[:, 2] - q[:, 0], q[:, 3] - q[:, 0])) / 6.0
                ).sum()
            )
            err = abs(tv - vol) / vol * 100
            log(f"  char_len={cl*1000:.1f}mm -> 节点 {len(gc):,} 四面体 {nt:,} "
                f"vol={tv*1e6:.4f} cm³ 误差={err:.3f}% {'✅' if err<1 else '❌'}")
        else:
            log(f"  char_len={cl*1000:.1f}mm -> 无四面体")

gmsh.finalize()
log("完成")
