"""软骨缝合: 在 add() 之后设 Geometry.Tolerance, 并用当次网格建点。"""

import os
import sys
from pathlib import Path

import numpy as np

SDK = Path(r"D:\Program\gmsh-sdk\gmsh-4.15.2-Windows64-sdk")
LIB = SDK / "lib"
os.environ["PATH"] = f"{LIB};{SDK/'bin'};{os.environ['PATH']}"
sys.path.insert(0, str(LIB))

ROOT = Path(__file__).resolve().parents[2]
GEOM = ROOT / "temp" / "ankle_geom"

import gmsh  # noqa: E402

STL = GEOM / "tibial_talus_cartilage_r.stl"


def build(tol, fix_edges=0, fix_faces=0):
    gmsh.initialize(str(LIB))
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("General.Verbosity", 0)
        gmsh.open(str(STL))

        # --- 1. 先在只读模型里拿到 STL 的三角网 (不设容差, 避免 open 时就变) ---
        et, _, en = gmsh.model.mesh.getElements(2)
        tri = next(np.array(n).reshape(-1, 3) for e, n in zip(et, en) if e == 2)
        ntags, ncoord, _ = gmsh.model.mesh.getNodes()
        ncoord = np.array(ncoord, dtype=np.float64).reshape(-1, 3)
        tag2row = {int(x): i for i, x in enumerate(ntags)}
        rows = np.vectorize(tag2row.get)(tri)  # (n,3) 行号, 与 ncoord 同序

        # --- 2. 切到干净模型, 这时才设几何容差 ---
        gmsh.model.add("chk")
        gmsh.option.setNumber("Geometry.Tolerance", tol)
        gmsh.option.setNumber("Geometry.ToleranceBoolean", tol)
        gmsh.option.setNumber("Geometry.OCCFixSmallEdges", fix_edges)
        gmsh.option.setNumber("Geometry.OCCFixSmallFaces", fix_faces)

        occ = gmsh.model.occ
        pts = []
        for i in range(len(ntags)):
            q = ncoord[i]
            pts.append(occ.addPoint(float(q[0]), float(q[1]), float(q[2]), tag=int(i + 1)))

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

        # 注意: 此刻 nsurf 可能因合并而变少
        surfs = [t for d, t in gmsh.model.getEntities(2)]
        bnd = gmsh.model.getBoundary([(2, x) for x in surfs])
        vol = None
        if len(bnd) == 0:
            sl = occ.addSurfaceLoop(surfs)
            occ.addVolume([sl])
            occ.synchronize()
            v3 = gmsh.model.getEntities(3)
            if v3:
                vol = float(occ.getMass(3, v3[0][1]))
        return len(surfs), len(bnd), vol
    finally:
        gmsh.finalize()


print("=== 扫 Geometry.Tolerance (在 add() 之后设) ===")
for tol in (1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3):
    try:
        nf, nb, vol = build(tol)
        flag = "OK 封闭" if nb == 0 else f"边界边 {nb}"
        v = f"  vol={vol*1e6:7.3f} cm3" if vol is not None else ""
        print(f"  tol={tol:.0e}  面={nf:>6,}  {flag}{v}")
    except Exception as e:
        print(f"  tol={tol:.0e}  FAIL {type(e).__name__}: {str(e)[:60]}")

print("\n=== tol=1e-5 + OCCFixSmall* ===")
for fe, ff in ((1, 0), (1, 1), (0, 1)):
    try:
        nf, nb, vol = build(1e-5, fix_edges=fe, fix_faces=ff)
        flag = "OK 封闭" if nb == 0 else f"边界边 {nb}"
        v = f"  vol={vol*1e6:7.3f} cm3" if vol is not None else ""
        print(f"  fixE={fe} fixF={ff}  面={nf:>6,}  {flag}{v}")
    except Exception as e:
        print(f"  fixE={fe} fixF={ff}  FAIL {str(e)[:60]}")
