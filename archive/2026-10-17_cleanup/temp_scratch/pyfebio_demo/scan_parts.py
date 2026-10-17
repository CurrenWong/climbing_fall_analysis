"""扫所有右踝部件: 不同去退化/合并策略下, OCC 能否判封闭。

思路: gmsh 会丢弃退化三角形(合并点后两顶点重合的), 丢掉后就出现小洞。
      所以要主动删掉退化三角形, 并把因此产生的悬挂边一并处理。
"""

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


PARTS = [
    "talus_r.stl",
    "calcaneus_r.stl",
    "tibial_talus_cartilage_r.stl",
    "talus_cartilage_r.stl",
    "midfoot_r.stl",
]


def prepare(stl: Path, tol: float, drop_degen: bool) -> tuple[np.ndarray, np.ndarray] | None:
    """返回 (points, faces) 或 None。"""
    s = pv.read(str(stl))
    c = s.clean(point_merging=True, tolerance=tol, absolute=True)
    f = c.faces.reshape(-1, 4)[:, 1:]
    if drop_degen:
        good = (
            (f[:, 0] != f[:, 1]) & (f[:, 1] != f[:, 2]) & (f[:, 0] != f[:, 2])
        )
        f = f[good]
        # 重新压缩点索引
        used = np.unique(f)
        remap = np.full(c.n_points, -1, dtype=np.int64)
        remap[used] = np.arange(len(used))
        f = remap[f]
        pts = c.points[used]
    else:
        pts = c.points
    if len(f) < 4:
        return None
    return np.asarray(pts, dtype=float), np.asarray(f, dtype=np.int64)


def occ_boundary(pts: np.ndarray, faces: np.ndarray):
    """建 plane-surface BREP, 返回 (边界边数, 实体体积 or None)。"""
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
    except Exception as e:
        return -1, str(e)[:40]
    finally:
        try:
            gmsh.finalize()
        except Exception:
            pass


print(f"{'部件':<32}{'tol':>7} {'删退化':>6} {'三角':>7} {'边界边':>7}  结果")
for name in PARTS:
    p = GEOM / name
    if not p.exists():
        log(f"[!] 缺 {name}")
        continue
    s = pv.read(str(p))
    print(f"  {name}  (表面 {s.n_points} 点 {s.n_cells} 面 vol={s.volume*1e6:.3f} cm³)")
    for tol in (1e-6, 1e-5):
        for dd in (False, True):
            r = prepare(p, tol, dd)
            if r is None:
                continue
            pts, faces = r
            nb, vol = occ_boundary(pts, faces)
            if isinstance(vol, str):
                res = f"ERR {vol}"
                vol = None
            elif nb == 0 and vol:
                res = f"✅ 封闭 vol={vol*1e6:.4f}"
            else:
                res = "❌"
            print(f"  {'':<30}{tol:>7.0e} {str(dd):>6} {len(faces):>7} {nb:>7}  {res}")
    log(f"  ^ {name} 完成")
