"""用第三方库（VTK via pyvista）独立验算体积 —— 第三种算法，排除自写代码的自证偏差。

同时给凸包体积作为**上界**（任何实体的体积都不能超过其节点凸包）。
并把混用单元类型正确映射：4 唯一节点=TETRA(10)、6=WEDGE(13)、8=HEXAHEDRON(12)。
"""
from __future__ import annotations

import numpy as np
import pyvista as pv

NPZ = r"D:\Project\climbing_fall_analysis\temp\thums\thums_lowerlimb.npz"

VTK_TETRA, VTK_HEXAHEDRON, VTK_WEDGE = 10, 12, 13
# wedge 的标准 6 节点顺序（用于从坍塌 hex 里还原，这里只做启发式映射）
HEX_TETS = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6),
            (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]


def tet_vol(p):
    return np.einsum("ij,ij->i",
                     np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]),
                     p[:, 3] - p[:, 0]) / 6.0


def my_volume(coords, uniq_hist_keys):
    """我的算法：tets 直接算；hex 按 6-四面体分解。"""
    tot = 0.0
    nd = np.array([len(set(r)) for r in coords.reshape(len(coords), -1, 3)[:, :, 0].tolist()])
    for i, c in enumerate(coords):
        u = len(set(map(tuple, c)))
        if u == 4:
            tot += abs(float(tet_vol(c[None, :4, :])[0]))
        elif u == 8:
            tot += abs(float(sum(tet_vol(c[None, list(t), :])[0] for t in HEX_TETS)))
        elif u == 6:
            # 6 唯一节点 -> 当成 wedge，用 3 个四面体分解 (0,1,2,3) (1,2,3,4)... 简化
            tot += abs(float(tet_vol(c[None, [0, 1, 2, 4], :])[0])) + \
                   abs(float(tet_vol(c[None, [1, 2, 3, 4], :])[0])) + \
                   abs(float(tet_vol(c[None, [2, 3, 4, 5], :])[0]))
    return tot


d = np.load(NPZ)
nid = d["node_ids"]
xyz = d["node_xyz"]
idx = {int(v): i for i, v in enumerate(nid)}

TARGETS = ["81000900_solid", "81000801_solid", "81000800_solid",
           "81000700_solid", "81000600_solid", "81000601_solid",
           "81001100_solid", "81001000_solid",
           "81001300_solid", "81001200_solid"]

print(f"{'part':18s} {'tets':>7s} {'hex':>6s} {'wedge':>6s} "
      f"{'my_vol':>12s} {'vtk_vol':>12s} {'hull_vol':>12s} {'vtk/my':>8s}")
print("-" * 96)

for key in TARGETS:
    arr = d[key]
    flat = np.vectorize(lambda x: idx.get(int(x), -1))(arr)
    coords = xyz[flat]                                    # (n, w, 3)
    nuniq = np.array([len(set(r)) for r in arr.tolist()])

    # ── 构建 VTK 网格（tets 为主，hex 也交给 VTK）──
    cells, celltypes = [], []
    for i, u in enumerate(nuniq):
        base = i * arr.shape[1]
        if u == 4:
            cells.append(4); cells.extend([base + j for j in (0, 1, 2, 3)])
            celltypes.append(VTK_TETRA)
        elif u == 8:
            cells.append(8); cells.extend([base + j for j in range(8)])
            celltypes.append(VTK_HEXAHEDRON)
        elif u == 6:
            cells.append(6); cells.extend([base + j for j in (0, 1, 2, 3, 4, 5)])
            celltypes.append(VTK_WEDGE)
    pts = coords.reshape(-1, 3)
    grid = pv.UnstructuredGrid(np.array(cells, np.int64),
                               np.array(celltypes, np.uint8), pts)
    vtk_vol = float(grid.volume)

    my_vol = my_volume(coords, None)

    # ── 凸包上界 ──
    try:
        from scipy.spatial import ConvexHull
        hull = ConvexHull(np.unique(pts, axis=0))
        hull_vol = float(hull.volume)
    except Exception as exc:
        hull_vol = float("nan")

    n_hex = int((nuniq == 8).sum())
    n_we = int((nuniq == 6).sum())
    n_tet = int((nuniq == 4).sum())
    print(f"{key:18s} {n_tet:>7,} {n_hex:>6,} {n_we:>6,} "
          f"{my_vol:>12.1f} {vtk_vol:>12.1f} {hull_vol:>12.1f} "
          f"{vtk_vol/my_vol:>8.4f}"
          + ("  ⚠️ vtк>hull!" if vtk_vol > hull_vol else ""))

print("\n对照解剖量级（cm³，整块骨）:")
print("  距骨 24–30 | 跟骨 45–60 | 腓骨 30–45 | 胫骨 300–400")
