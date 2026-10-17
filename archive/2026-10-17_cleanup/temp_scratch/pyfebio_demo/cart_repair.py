"""修复软骨的 11 条 4-共享边 + 1 条 6-共享边(非流形), 重建干净的三角网。"""

import os
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

SDK = Path(r"D:\Program\gmsh-sdk\gmsh-4.15.2-Windows64-sdk")
LIB = SDK / "lib"
os.environ["PATH"] = f"{LIB};{SDK/'bin'};{os.environ['PATH']}"
sys.path.insert(0, str(LIB))

ROOT = Path(__file__).resolve().parents[2]
GEOM = ROOT / "temp" / "ankle_geom"
STL = GEOM / "tibial_talus_cartilage_r.stl"

import pyvista as pv  # noqa: E402

T0 = time.time()


def log(m):
    print(f"  [{time.time()-T0:6.1f}s] {m}", flush=True)


def edge_report(pts, faces):
    e = np.sort(np.vstack([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1)
    u, c = np.unique(e, axis=0, return_counts=True)
    return u, c


raw = pv.read(str(STL))
log(f"原始 {raw.n_points} 点 {raw.n_cells} 面")
merged = raw.merge_points(tolerance=1e-6)
log(f"merge(1e-6) {merged.n_points} 点 {merged.n_cells} 面")

faces = merged.faces.reshape(-1, 4)[:, 1:]
u, c = edge_report(merged.points, faces)
log(f"边计数 {dict(sorted(Counter(c.tolist()).items()))}")
bad = u[c > 2]
log(f"非流形边 {len(bad)} 条")

for e in bad:
    lens = np.linalg.norm(merged.points[e[0]] - merged.points[e[1]])
    # 哪些三角形用这条边
    used = np.where(
        ((faces == e[0]).any(1) & (faces == e[1]).any(1))
    )[0]
    log(f"  边 {tuple(e)}  长 {lens*1e6:.2f} um  被三角 {used.tolist()} 使用")

# ---- 用 vtkFeatureEdges / extract_feature_edges 找并处理 ----
log("\n尝试 VTK 的非流形修复:")

# 方法1: extract_feature_edges(boundary_edges=True, feature_edges=True) 看它怎么判
try:
    fe = merged.extract_feature_edges(
        boundary_edges=True, feature_edges=True, manifold_edges=False,
        non_manifold_edges=True, feature_angle=30,
    )
    log(f"extract_feature_edges -> {fe.n_lines} 条")
except Exception as e:
    log(f"  extract_feature_edges FAIL: {str(e)[:80]}")

# 方法2: vtkCleanPolyData 的 PointMerging (VTK 内部会拆非流形)
try:
    cl = merged.clean(point_merging=True, tolerance=1e-6, absolute=True)
    f2 = cl.faces.reshape(-1, 4)[:, 1:]
    u2, c2 = edge_report(cl.points, f2)
    log(f"clean(1e-6,absolute) -> {cl.n_points} 点 {cl.n_cells} 面  边计数 {dict(sorted(Counter(c2.tolist()).items()))}")
except Exception as e:
    log(f"  clean FAIL: {str(e)[:80]}")

# 方法3: vtkFeatureEdges + removeCells, 手动剔除问题三角形再补洞
log("\n方法3: 手动拆非流形 (逐三角形判断法向一致性)")
pts = merged.points
faces = np.array(faces)
# 计算每个三角形的法向
n = np.cross(pts[faces[:, 1]] - pts[faces[:, 0]], pts[faces[:, 2]] - pts[faces[:, 0]])
keep = np.ones(len(faces), dtype=bool)
dropped = 0
for e in bad:
    used = np.where(((faces == e[0]).any(1)) & ((faces == e[1]).any(1)))[0]
    if len(used) < 2:
        continue
    # 只保留法向"最一致"的前两个 (即同一闭合面上的)
    for grp in _group_by_edge(faces, used, e):
        if len(grp) <= 2:
            continue
        # 保留法向平均方向最接近的 2 个
        nn = n[grp]
        mean = nn.mean(0)
        mean /= np.linalg.norm(mean)
        cosang = (nn @ mean)
        order = np.argsort(-cosang)
        for idx in grp[order[2:]]:
            keep[idx] = False
            dropped += 1
log(f"丢弃 {dropped} 个三角形")

f3 = faces[keep]
u3, c3 = edge_report(pts, f3)
log(f"清理后 边计数 {dict(sorted(Counter(c3.tolist()).items()))}")


def _group_by_edge(faces, used, e):
    """把 used 里的三角形按'经过该边的次序'分组 (每组是一条连续折线)。"""
    groups = []
    remaining = list(used)
    # 简化: 每两个相邻三角形一组
    while len(remaining) >= 2:
        g = remaining[:2]
        remaining = remaining[2:]
        groups.append(g)
    if remaining:
        groups.append(remaining)
    return groups
