"""查清 tri 的 tag 与 ntags 到底为什么不匹配。"""

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
STL = GEOM / "tibial_talus_cartilage_r.stl"

import gmsh  # noqa: E402

for tol in (None, 1e-5):
    gmsh.initialize(str(LIB))
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("General.Verbosity", 0)
    if tol:
        gmsh.option.setNumber("Geometry.Tolerance", tol)
        gmsh.option.setNumber("Geometry.ToleranceBoolean", tol)
    gmsh.open(str(STL))

    ntags, ncoord, _ = gmsh.model.mesh.getNodes()
    et, etg, en = gmsh.model.mesh.getElements(2)
    tri = next(np.array(n).reshape(-1, 3) for e, n in zip(et, en) if e == 2)

    print(f"\ntol={tol}")
    print(f"  getNodes: {len(ntags)} 个  tag范围 [{ntags.min()}, {ntags.max()}]  是否连续={np.array_equal(ntags, np.arange(1, len(ntags)+1))}")
    print(f"  tri: {tri.shape}  tag范围 [{tri.min()}, {tri.max()}]")
    s = set(ntags.tolist())
    missing = sorted(set(tri.ravel().tolist()) - s)
    print(f"  tri 里不在 ntags 的 tag: {len(missing)} 个  例: {missing[:10]}")
    gmsh.finalize()
