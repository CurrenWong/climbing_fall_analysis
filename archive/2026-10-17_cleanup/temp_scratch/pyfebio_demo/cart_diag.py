"""诊断踝关节软骨 STL 为何被 OCC 判为非封闭。"""

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
import pyvista as pv  # noqa: E402

TARGETS = ["tibial_talus_cartilage_r.stl", "talus_cartilage_r.stl", "talus_r.stl"]

for t in TARGETS:
    p = GEOM / t
    if not p.exists():
        print(f"\n[!] 缺 {t}")
        continue

    print(f"\n{'='*70}\n{t}")

    # --- VTK 侧判断 ---
    surf = pv.read(str(p))
    raw = surf
    clean = surf.clean(point_merging=True, tolerance=1e-9)
    print(f"  VTK 原始 : {raw.n_points:>6} 点 {raw.n_cells:>6} 面 vol={raw.volume*1e6:>7.3f} cm³")
    print(f"  VTK clean: {clean.n_points:>6} 点 {clean.n_cells:>6} 面 vol={clean.volume*1e6:>7.3f} cm³")
    print(f"  VTK 开放边: {clean.n_open_edges}")

    # 三角形边统计: 每条边应被恰好 2 个三角形共享
    faces = clean.faces.reshape(-1, 4)[:, 1:]
    edges = np.vstack([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    edges_sorted = np.sort(edges, axis=1)
    uniq, counts = np.unique(edges_sorted, axis=0, return_counts=True)
    from collections import Counter

    c = Counter(counts.tolist())
    print(f"  边计数分布: {dict(sorted(c.items()))}")
    bad = uniq[counts != 2]
    print(f"  非2次边: {len(bad)} 条")
    if len(bad):
        pts = clean.points
        for e in bad[:5]:
            mid = pts[e].mean(0)
            print(f"      {tuple(e)}  中点={np.round(mid*1000, 2)} mm")
        # 这些边的长度分布
        if len(bad) > 1:
            lens = np.linalg.norm(pts[bad[:, 0]] - pts[bad[:, 1]], axis=1)
            print(f"      非2次边长度: min={lens.min()*1e6:.3f} max={lens.max()*1e6:.3f} um")

    # --- OCC 侧判断 ---
    gmsh.initialize(str(LIB))
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("General.Verbosity", 0)
        gmsh.open(str(p))
        et, etg, en = gmsh.model.mesh.getElements(2)
        tri = None
        for e, nodes in zip(et, en):
            if e == 2:
                tri = np.array(nodes).reshape(-1, 3)
                break
        ntags, ncoord, _ = gmsh.model.mesh.getNodes()
        ncoord = np.array(ncoord).reshape(-1, 3)
        tag2row = {int(x): i for i, x in enumerate(ntags)}
        tt = np.vectorize(tag2row.get)(tri)

        gmsh.model.add("chk")
        occ = gmsh.model.occ
        for i, tg in enumerate(ntags):
            q = ncoord[i]
            occ.addPoint(float(q[0]), float(q[1]), float(q[2]), tag=int(tg))
        e2c = {}
        surfs = []
        for a, b, c_ in tt:
            ta, tb, tc = int(ntags[a]), int(ntags[b]), int(ntags[c_])
            cur = []
            for u, v in ((ta, tb), (tb, tc), (tc, ta)):
                k = (min(u, v), max(u, v))
                tg2 = e2c.get(k)
                if tg2 is None:
                    tg2 = occ.addLine(u, v)
                    e2c[k] = tg2
                cur.append(tg2)
            surfs.append(occ.addPlaneSurface([occ.addCurveLoop(cur)]))
        occ.synchronize()
        bnd = gmsh.model.getBoundary([(2, x) for x in surfs])
        print(f"  OCC 三角={len(tt):,} 点={len(ntags):,} 面={len(surfs):,}")
        print(f"  OCC getBoundary 边界边数 = {len(bnd)}  {'✅封闭' if len(bnd)==0 else '❌非封闭'}")
    finally:
        gmsh.finalize()
