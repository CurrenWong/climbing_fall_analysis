"""查 gmsh Geometry 容差选项 与 pyvista 点合并 API 的正确名字。"""

import os
import sys
from pathlib import Path

SDK = Path(r"D:\Program\gmsh-sdk\gmsh-4.15.2-Windows64-sdk")
LIB = SDK / "lib"
os.environ["PATH"] = f"{LIB};{SDK/'bin'};{os.environ['PATH']}"
sys.path.insert(0, str(LIB))

import gmsh  # noqa: E402
import pyvista as pv  # noqa: E402

gmsh.initialize(str(LIB))
gmsh.option.setNumber("General.Terminal", 0)
gmsh.option.setNumber("General.Verbosity", 0)
gmsh.model.add("x")

print("=== Geometry.* 容差/缝合选项 ===")
for o in [
    "Geometry.Tolerance",
    "Geometry.ToleranceBoolean",
    "Geometry.ToleranceBooleanPreserveNumbering",
    "Geometry.OCCFixSmallEdges",
    "Geometry.OCCFixSmallFaces",
    "Geometry.OCCSewFaces",
    "Geometry.ToleranceMergePoints",
    "Geometry.OCCUnionUnify",
    "Geometry.OCCImportLabels",
]:
    try:
        print(f"  OK  {o} = {gmsh.option.getNumber(o)}")
    except Exception as e:
        print(f"  --  {o}  ({str(e)[:40]})")

print("\n=== occ 缝合/修复相关函数 ===")
print("  ", [a for a in dir(gmsh.model.occ) if any(k in a.lower() for k in ("sew", "fix", "heal"))])

gmsh.finalize()

print("\n=== pyvista 点合并 API ===")
print("  pv 模块级:", [n for n in dir(pv) if any(k in n.lower() for k in ("merge", "clean", "unique"))])
print("  PolyData:", [n for n in dir(pv.PolyData) if any(k in n.lower() for k in ("merge", "clean", "unique"))])
print("  PolyDataFilters:", [n for n in dir(pv.PolyDataFilters) if any(k in n.lower() for k in ("merge", "clean", "unique"))])
