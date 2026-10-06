"""踝关节几何处理: STL -> 四面体网格 (官方 gmsh SDK + classifySurfaces 缝合)。

对应 ANKLE_FE_TODO.md 步骤 3-4。

为什么必须用官方 SDK:
    pip 装的 gmsh wheel 只有 Python 绑定, occ API 残缺
    (缺 addDiscreteSurface / reclassify / geo.addTriangle),
    且 merge(STL) 只得 "Discrete surface", 无法成实体。
    官方 Windows64 SDK (gmsh.info/bin/Windows/gmsh-*-Windows64-sdk.zip) 完整。

STL -> 实体的正确流程 (依据官方 examples/api/glue_and_remesh_stl.py):
    merge(stl)                          读离散三角网
    mesh.removeDuplicateNodes()         合并重合点(容差 Geometry.Tolerance)
    mesh.classifySurfaces(angle)        关键! 把共面三角形归类成曲面实体
    mesh.createGeometry()               为离散曲面建参数化几何
    geo.addSurfaceLoop(所有曲面)          成环
    geo.addVolume([loop])               封闭成实体
    mesh.generate(3)                    四面体网格化

踩过的坑:
    - 逐三角形 occ.addPlaneSurface 建的面互相独立, 即使共面 OCC 也不认,
      getBoundary 永远非空 -> 必须用 classifySurfaces
    - 医学 STL 有 float32 噪声造成的非流形边(VTK 边计数出现 4/6 而非 2),
      会让 OCC 缝合失败 -> 先用 pyvista clean(tol=1e-6, absolute=True) 修
    - DLL 在 SDK 的 lib/ 不在 bin/ -> gmsh.initialize(str(LIB)) + PATH
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

SDK = Path(r"D:\Program\gmsh-sdk\gmsh-4.15.2-Windows64-sdk")
LIB = SDK / "lib"
os.environ["PATH"] = f"{LIB};{SDK / 'bin'};{os.environ['PATH']}"
if str(LIB) not in sys.path:
    sys.path.insert(0, str(LIB))

ROOT = Path(__file__).resolve().parents[2]
GEOM = ROOT / "temp" / "ankle_geom"
MESH = ROOT / "temp" / "ankle_mesh"
MESH.mkdir(parents=True, exist_ok=True)


def fix_nonmanifold(stl: Path, tol: float = 1e-6) -> Path:
    """用 pyvista clean 修非流形边, 返回清理后的 STL 路径。"""
    import pyvista as pv

    surf = pv.read(str(stl))
    cl = surf.clean(point_merging=True, tolerance=tol, absolute=True)
    if (cl.n_points, cl.n_cells) == (surf.n_points, surf.n_cells):
        return stl  # 无需清理
    dst = stl.parent / f"_{stl.stem}_clean.stl"
    cl.save(str(dst))
    return dst


def stl_to_tet(
    stl: Path,
    char_len: float,
    angle: float = 1.5707963,  # pi/2, 官方默认
    out_msh: Path | None = None,
    verbose: bool = False,
):
    """STL -> 四面体网格。返回 (nodes (n,3) float64, tets (m,4) int64) 或 None。"""
    import gmsh

    src = fix_nonmanifold(stl)

    gmsh.initialize(str(LIB))
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("General.Verbosity", 0)

        # 1) 读离散 STL
        gmsh.merge(str(src))
        # 2) 合并重合点
        gmsh.option.setNumber("Geometry.Tolerance", 1e-6)
        try:
            gmsh.model.mesh.removeDuplicateNodes()
        except Exception:
            pass
        # 3) 关键: 把共面三角形归类成曲面实体
        gmsh.model.mesh.classifySurfaces(angle)
        # 4) 为离散曲面建参数化几何
        gmsh.model.mesh.createGeometry()
        gmsh.model.geo.synchronize()

        surfs = [t for d, t in gmsh.model.getEntities(2)]
        if verbose:
            print(f"    classifySurfaces -> {len(surfs)} 个曲面")
        if not surfs:
            return None

        # 5) 成环 -> 实体
        try:
            loop = gmsh.model.geo.addSurfaceLoop(surfs)
            gmsh.model.geo.addVolume([loop])
            gmsh.model.geo.synchronize()
        except Exception as e:
            if verbose:
                print(f"    addVolume FAIL: {str(e)[:80]}")
            return None

        if not gmsh.model.getEntities(3):
            if verbose:
                print("    未成实体(非封闭)")
            return None

        # 6) 四面体网格化
        gmsh.option.setNumber("Mesh.MeshSizeMax", char_len)
        gmsh.option.setNumber("Mesh.MeshSizeMin", char_len * 0.4)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)  # Delaunay
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.optimize("Netgen")

        gtags, gcoord, _ = gmsh.model.mesh.getNodes()
        gcoord = np.array(gcoord, dtype=np.float64).reshape(-1, 3)
        rmap = {int(t): i for i, t in enumerate(gtags)}
        t3, _, en3 = gmsh.model.mesh.getElements(3)
        tets = None
        for et, en in zip(t3, en3):
            if et == 4:
                tets = np.vectorize(rmap.get)(np.array(en).reshape(-1, 4)).astype(np.int64)
                break

        if out_msh and tets is not None:
            gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
            gmsh.write(str(out_msh))
        return gcoord, tets
    finally:
        try:
            gmsh.finalize()
        except Exception:
            pass


def tet_volume(tets: np.ndarray, nodes: np.ndarray) -> float:
    q = nodes[tets]
    return float(
        np.abs(
            np.einsum(
                "ij,ij->i", q[:, 1] - q[:, 0], np.cross(q[:, 2] - q[:, 0], q[:, 3] - q[:, 0])
            )
            / 6.0
        ).sum()
    )


if __name__ == "__main__":
    import time

    import pyvista as pv

    targets = sys.argv[1:] or [
        "talus_r.stl",
        "calcaneus_r.stl",
        "tibial_talus_cartilage_r.stl",
        "talus_cartilage_r.stl",
    ]
    for t in targets:
        p = GEOM / t
        if not p.exists():
            print(f"[!] 缺 {t}")
            continue
        surf = pv.read(str(p))
        print(f"\n{t}   表面体积 {surf.volume*1e6:.3f} cm³")
        t0 = time.time()
        for cl in (0.006, 0.003):
            try:
                r = stl_to_tet(p, cl, out_msh=MESH / f"{p.stem}_{int(cl*1000)}.msh", verbose=True)
                if r is None:
                    print(f"  char_len={cl*1000:.0f}mm -> 失败")
                    continue
                nodes, tets = r
                if tets is None:
                    print(f"  char_len={cl*1000:.0f}mm -> 无四面体")
                    continue
                vol = tet_volume(tets, nodes)
                err = abs(vol - surf.volume) / surf.volume * 100
                print(
                    f"  char_len={cl*1000:.0f}mm -> 节点 {len(nodes):>7,} 四面体 {len(tets):>8,} "
                    f"vol={vol*1e6:7.3f} cm³ 误差={err:5.2f}% {'✅' if err < 2 else '❌'}"
                    f"  [{time.time()-t0:.0f}s]"
                )
            except Exception as e:
                print(f"  char_len={cl*1000:.0f}mm FAIL {type(e).__name__}: {str(e)[:70]}")
