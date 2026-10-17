"""用 VTK Delaunay3D 把闭合 STL 四面体化 (pyvista 4-node TETRA = celltype 10)。"""

from pathlib import Path

import numpy as np
import pyvista as pv

ROOT = Path(__file__).resolve().parents[2]
GEOM = ROOT / "temp" / "ankle_geom"

TET = int(pv.CellType.TETRA)  # pyvista 编号 = 10


def tet_from_stl(stl: Path, tol_frac: float = 0.02, offset: float = 2.0):
    """闭合 STL -> (surface, tets ndarray (n,4) 0-based, nodes ndarray (m,3))。"""
    surf = pv.read(str(stl)).clean(point_merging=True, tolerance=1e-9)
    b = np.asarray(surf.bounds, dtype=float)
    diag = float(np.linalg.norm(b[3:] - b[:3]))
    grid = surf.delaunay_3d(alpha=0.0, tol=tol_frac * diag, offset=offset)

    # 混合单元时 cells_dict 可能缺项, 直接从 celltypes 过滤线性四面体
    conn = grid.cells.reshape(-1, 5)[:, 1:]  # VTK: [npts, n0..nk]
    types = np.asarray(grid.celltypes)
    tmask = types == TET
    if not tmask.any():
        return surf, None, None, grid
    nodes = np.asarray(grid.points, dtype=float)
    return surf, conn[tmask].astype(np.int64), nodes, grid


if __name__ == "__main__":
    import sys
    from collections import Counter

    targets = sys.argv[1:] or [
        "talus_r.stl",
        "calcaneus_r.stl",
        "tibial_talus_cartilage_r.stl",
    ]
    for t in targets:
        p = GEOM / t
        if not p.exists():
            print(f"[!] 缺 {t}")
            continue
        print(f"\n{t}")
        for tol_f in (0.08, 0.04, 0.02):
            try:
                surf, tets, nodes, grid = tet_from_stl(p, tol_frac=tol_f)
                if tets is None:
                    ct = Counter(np.asarray(grid.celltypes).tolist())
                    print(f"  tol={tol_f*100:>4.1f}%  无线性四面体, celltypes={dict(ct)}")
                    continue
                # 四面体体积
                q = nodes[tets]
                vol = np.abs(
                    np.einsum(
                        "ij,ij->i", q[:, 1] - q[:, 0], np.cross(q[:, 2] - q[:, 0], q[:, 3] - q[:, 0])
                    )
                    / 6.0
                ).sum()
                v_surf = surf.volume
                print(
                    f"  tol={tol_f*100:>4.1f}%  tets={len(tets):>8,}  nodes={len(nodes):>8,}  "
                    f"vol={vol*1e6:>7.2f} cm³ (表面 {v_surf*1e6:.2f}, 误差 {abs(vol-v_surf)/max(v_surf,1e-12)*100:>5.1f}%)"
                )
            except Exception as e:
                print(f"  tol={tol_f*100:>4.1f}%  FAIL {str(e)[:70]}")
