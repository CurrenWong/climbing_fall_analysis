"""软骨缝合: 只测 2 个最可能成功的策略, 加进度输出避免超时无反馈。"""

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

T0 = time.time()


def log(msg):
    print(f"  [{time.time()-T0:6.1f}s] {msg}", flush=True)


def build(tol=1e-8):
    gmsh.open(str(STL))
    et, _, en = gmsh.model.mesh.getElements(2)
    tri = next(np.array(n).reshape(-1, 3) for x, n in zip(et, en) if x == 2)
    ntags, ncoord, _ = gmsh.model.mesh.getNodes()
    ncoord = np.array(ncoord, dtype=np.float64).reshape(-1, 3)
    tag2row = {int(x): i for i, x in enumerate(ntags)}
    rows = np.vectorize(tag2row.get)(tri)

    gmsh.model.add("chk")
    gmsh.option.setNumber("Geometry.Tolerance", tol)
    gmsh.option.setNumber("Geometry.ToleranceBoolean", tol)
    occ = gmsh.model.occ
    pts = [occ.addPoint(*ncoord[i], tag=int(i + 1)) for i in range(len(ntags))]
    log(f"点建完 {len(pts)}")

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
    log(f"面建完 {len(surfs)}")
    occ.synchronize()
    log("synchronize 完")
    return occ, surfs


def check(tag):
    ents2 = gmsh.model.getEntities(2)
    bnd = gmsh.model.getBoundary([(2, t) for d, t in ents2])
    log(f"{tag}: 面 {len(ents2)}  边界边 {len(bnd)}")
    return len(bnd)


# ---------- 策略 1: 极宽松容差 ----------
print("\n[1] Geometry.Tolerance=1e-3")
gmsh.initialize(str(LIB))
gmsh.option.setNumber("General.Terminal", 0)
gmsh.option.setNumber("General.Verbosity", 0)
occ, surfs = build(tol=1e-3)
nb = check("baseline tol=1e-3")
if nb == 0:
    sl = occ.addSurfaceLoop([t for d, t in gmsh.model.getEntities(2)])
    occ.addVolume([sl])
    occ.synchronize()
    v3 = gmsh.model.getEntities(3)
    if v3:
        log(f"✅ 封闭! 体积 = {float(occ.getMass(3, v3[0][1]))*1e6:.4f} cm³")
gmsh.finalize()
log("finalize 完")

# ---------- 策略 2: 用 VTK 先把几何"焊死"再送 gmsh ----------
print("\n[2] VTK merge_points(tol=1e-6) 后再送 gmsh")
import pyvista as pv  # noqa: E402

raw = pv.read(str(STL))
merged = raw.merge_points(tolerance=1e-6) if hasattr(raw, "merge_points") else raw
log(f"merge 后 {merged.n_points} 点 / {merged.n_cells} 面")

faces = merged.faces.reshape(-1, 4)[:, 1:]
e = np.sort(np.vstack([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1)
u, c = np.unique(e, axis=0, return_counts=True)
from collections import Counter

log(f"边计数 {dict(sorted(Counter(c.tolist()).items()))}")

mergedstl = ROOT / "temp" / "ankle_geom" / "_cart_merged.stl"
merged.save(str(mergedstl))

gmsh.initialize(str(LIB))
gmsh.option.setNumber("General.Terminal", 0)
gmsh.option.setNumber("General.Verbosity", 0)
gmsh.open(str(mergedstl))
et, _, en = gmsh.model.mesh.getElements(2)
tri = next(np.array(n).reshape(-1, 3) for x, n in zip(et, en) if x == 2)
ntags, ncoord, _ = gmsh.model.mesh.getNodes()
ncoord = np.array(ncoord, dtype=np.float64).reshape(-1, 3)
tag2row = {int(x): i for i, x in enumerate(ntags)}
rows = np.vectorize(tag2row.get)(tri)
gmsh.model.add("chk")
occ = gmsh.model.occ
pts = [occ.addPoint(*ncoord[i], tag=int(i + 1)) for i in range(len(ntags))]
log(f"点建完 {len(pts)}")
e2c = {}
surfs = []
for a, b, cc in rows:
    ta, tb, tc = pts[a], pts[b], pts[cc]
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
log(f"面建完 {len(surfs)}")
bnd = gmsh.model.getBoundary([(2, t) for d, t in gmsh.model.getEntities(2)])
log(f"边界边 {len(bnd)}")
if len(bnd) == 0:
    sl = occ.addSurfaceLoop([t for d, t in gmsh.model.getEntities(2)])
    occ.addVolume([sl])
    occ.synchronize()
    v3 = gmsh.model.getEntities(3)
    if v3:
        log(f"✅ 封闭! 体积 = {float(occ.getMass(3, v3[0][1]))*1e6:.4f} cm³")
gmsh.finalize()
log("总耗时 %.1f s" % (time.time() - T0))
