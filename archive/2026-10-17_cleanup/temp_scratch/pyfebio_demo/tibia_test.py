"""降采样 + 四面体化: 验证大部件(胫骨 OBJ, 10万顶点)能否走通全链路。"""

import os
import sys
import time
from pathlib import Path

import numpy as np

SDK = Path(r"D:\Program\gmsh-sdk\gmsh-4.15.2-Windows64-sdk")
LIB = SDK / "lib"
os.environ["PATH"] = f"{LIB};{SDK / 'bin'};{os.environ['PATH']}"
sys.path.insert(0, str(LIB))

ROOT = Path(__file__).resolve().parents[2]
GEOM = ROOT / "temp" / "ankle_geom"
MESH = ROOT / "temp" / "ankle_mesh"

import gmsh  # noqa: E402
import pyvista as pv  # noqa: E402

T0 = time.time()


def log(m):
    print(f"[{time.time()-T0:6.1f}s] {m}", flush=True)


def occ_build(pts, faces, char_len, out_msh=None):
    """plane-surface BREP -> 四面体。返回 (边界边数, nodes, tets, occ_vol)。"""
    gmsh.initialize(str(LIB))
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("General.Verbosity", 0)
        gmsh.model.add("p")
        occ = gmsh.model.occ
        ptag = [occ.addPoint(float(p[0]), float(p[1]), float(p[2])) for p in pts]
        e2c = {}
        surfs = []
        for a, b, c in faces:
            cur = []
            for u, v in ((a, b), (b, c), (c, a)):
                k = (min(int(u), int(v)), max(int(u), int(v)))
                t = e2c.get(k)
                if t is None:
                    t = occ.addLine(ptag[u], ptag[v])
                    e2c[k] = t
                cur.append(t)
            surfs.append(occ.addPlaneSurface([occ.addCurveLoop(cur)]))
        occ.synchronize()
        bnd = gmsh.model.getBoundary([(2, t) for d, t in gmsh.model.getEntities(2)])
        if len(bnd) != 0:
            return len(bnd), None, None, None
        sl = occ.addSurfaceLoop([t for d, t in gmsh.model.getEntities(2)])
        occ.addVolume([sl])
        occ.synchronize()
        v3 = gmsh.model.getEntities(3)
        if not v3:
            return 0, None, None, None
        v_occ = float(occ.getMass(3, v3[0][1]))
        gmsh.option.setNumber("Mesh.MeshSizeMax", char_len)
        gmsh.option.setNumber("Mesh.MeshSizeMin", char_len * 0.4)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.optimize("Netgen")
        gt, gc, _ = gmsh.model.mesh.getNodes()
        gc = np.array(gc, dtype=float).reshape(-1, 3)
        rm = {int(t): i for i, t in enumerate(gt)}
        t3, _, en3 = gmsh.model.mesh.getElements(3)
        tets = None
        for et, en in zip(t3, en3):
            if et == 4:
                tets = np.vectorize(rm.get)(np.array(en).reshape(-1, 4)).astype(np.int64)
                break
        if out_msh and tets is not None:
            gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
            gmsh.write(str(out_msh))
        return 0, gc, tets, v_occ
    finally:
        try:
            gmsh.finalize()
        except Exception:
            pass


def tet_vol(tets, nodes):
    q = nodes[tets]
    return float(
        np.abs(
            np.einsum(
                "ij,ij->i", q[:, 1] - q[:, 0], np.cross(q[:, 2] - q[:, 0], q[:, 3] - q[:, 0])
            )
            / 6.0
        ).sum()
    )


src = GEOM / "OK_tibia.obj"
log(f"读 {src.name}")
m = pv.read(str(src))
bb = np.asarray(m.bounds, dtype=float)
log(f"  {m.n_points:,} 点 / {m.n_cells:,} 面  vol={m.volume*1e6:.3f} cm³")
log(f"  bbox = {np.round((bb[3:] - bb[:3]) * 100, 2)} cm")

# 连通体拆分 (OBJ 可能含胫骨+腓骨)
try:
    conn = m.connectivity("all")
    nreg = conn["RegionId"].max() + 1
    log(f"  连通体 {nreg} 个  各体积(cm³):", )
    for r in range(nreg):
        sub = m.extract_cells(np.where(conn["RegionId"] == r)[0])
        if sub.n_cells > 50:
            log(f"    region {r}: {sub.n_points:>6} 点 {sub.n_cells:>7} 面 vol={abs(sub.volume)*1e6:.3f}")
except Exception as e:
    log(f"  connectivity FAIL: {str(e)[:60]}")

# 降采样到 ~6000 面
for target_tri in (6000, 3000):
    red = 1.0 - target_tri / m.n_cells
    dm = m.decimate(red) if red > 0 else m
    dm = dm.clean(point_merging=True, tolerance=1e-7, absolute=True)
    f = dm.faces.reshape(-1, 4)[:, 1:]
    keep = (f[:, 0] != f[:, 1]) & (f[:, 1] != f[:, 2]) & (f[:, 0] != f[:, 2])
    f = f[keep]
    used = np.unique(f)
    rm2 = np.full(dm.n_points, -1, dtype=np.int64)
    rm2[used] = np.arange(len(used))
    f = rm2[f]
    pts = np.asarray(dm.points[used], dtype=float)
    log(f"降采样 -> 目标 {target_tri} 面, 实得 {len(f)} 面 / {len(pts)} 点")
    nb, nodes, tets, v_occ = occ_build(pts, f, char_len=0.004, out_msh=MESH / f"tibia_{target_tri}.msh")
    if nb != 0:
        log(f"  边界边 {nb} ❌ 非封闭")
        continue
    if tets is None:
        log("  无四面体")
        continue
    tv = tet_vol(tets, nodes)
    err = abs(tv - v_occ) / v_occ * 100
    log(
        f"  ✅ 节点 {len(nodes):,}  四面体 {len(tets):,}  "
        f"vol={tv*1e6:.3f} cm³  (OCC {v_occ*1e6:.3f}, 误差 {err:.3f}%)"
    )
    break
