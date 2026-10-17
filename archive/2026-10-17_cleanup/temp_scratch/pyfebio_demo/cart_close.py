"""软骨封口: 加大合并容差扫, 找能闭合的最小容差。只测 2 个软骨。"""

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

import gmsh  # noqa: E402
import pyvista as pv  # noqa: E402

T0 = time.time()


def log(m):
    print(f"[{time.time()-T0:5.1f}s] {m}", flush=True)


def occ_test(pts, faces):
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
        vol = None
        if len(bnd) == 0:
            sl = occ.addSurfaceLoop([t for d, t in gmsh.model.getEntities(2)])
            occ.addVolume([sl])
            occ.synchronize()
            v3 = gmsh.model.getEntities(3)
            if v3:
                vol = float(occ.getMass(3, v3[0][1]))
        return len(bnd), vol
    finally:
        try:
            gmsh.finalize()
        except Exception:
            pass


for name in ("tibial_talus_cartilage_r.stl", "talus_cartilage_r.stl"):
    p = GEOM / name
    s = pv.read(str(p))
    print(f"\n{name}  表面 {s.n_points} 点 {s.n_cells} 面 vol={s.volume*1e6:.4f} cm³")
    for tol in (1e-4, 3e-4, 6e-4, 1e-3):
        c = s.clean(point_merging=True, tolerance=tol, absolute=True)
        f = c.faces.reshape(-1, 4)[:, 1:]
        keep = (f[:, 0] != f[:, 1]) & (f[:, 1] != f[:, 2]) & (f[:, 0] != f[:, 2])
        f = f[keep]
        used = np.unique(f)
        rm = np.full(c.n_points, -1, dtype=np.int64)
        rm[used] = np.arange(len(used))
        f = rm[f]
        pts = np.asarray(c.points[used], dtype=float)
        nb, vol = occ_test(pts, f)
        dv = (vol - s.volume) / s.volume * 100 if vol else float("nan")
        r = f"✅ 封闭 vol={vol*1e6:.4f} cm³ (偏差 {dv:+.2f}%)" if nb == 0 and vol else "❌"
        log(f"  tol={tol:.0e} 点 {len(pts):>5} 三角 {len(f):>5} 边界边 {nb:>4}  {r}")
        if nb == 0:
            break
