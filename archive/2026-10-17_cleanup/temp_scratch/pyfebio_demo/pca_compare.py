"""PCA 主轴长度对比：THUMS 距骨/跟骨 vs IITD。
bbox 在 THUMS 全局站立坐标系里会被旋转放大，PCA 主轴是旋转不变的，可比。
"""
from pathlib import Path
import numpy as np
import pyvista as pv
import trimesh
import gmsh
import os, sys

SDK = Path(r"D:\Program\gmsh-sdk\gmsh-4.15.2-Windows64-sdk")
os.environ["PATH"] = f"{SDK/'lib'};{SDK/'bin'};{os.environ['PATH']}"
sys.path.insert(0, str(SDK / "lib"))

NPZ = r"D:\Project\climbing_fall_analysis\temp\thums\thums_lowerlimb.npz"
IITD = Path(r"D:\Project\climbing_fall_analysis\temp\ankle_geom")
MSH = Path(r"D:\Project\climbing_fall_analysis\temp\ankle_mesh")


def pca_extents(P: np.ndarray) -> np.ndarray:
    Q = P - P.mean(0)
    _, _, Vt = np.linalg.svd(Q, full_matrices=False)
    return np.ptp(Q @ Vt.T, axis=0)


def srt(e):
    return np.sort(e)[::-1]


print("=== PCA 主轴范围（长度，旋转不变）===")
print(f"{'source':28s} {'axes (mm, 降序)':>34s} {'diag':>8s} {'vol(raw)':>10s}")
print("-" * 88)

d = np.load(NPZ)
nid = d["node_ids"]; xyz = d["node_xyz"]
idx = {int(v): i for i, v in enumerate(nid)}

for label, key in (("THUMS R_TALUS (cort+spon)", ["81001100_solid", "81001000_solid"]),
                   ("THUMS R_CALCANEUS (cort+spon)", ["81001300_solid", "81001200_solid"])):
    pts = []
    for k in key:
        arr = d[k]
        flat = np.vectorize(lambda x: idx.get(int(x), -1))(arr)
        pts.append(xyz[flat].reshape(-1, 3))
    P = np.unique(np.vstack(pts), axis=0)
    e = srt(pca_extents(P))
    print(f"{label:28s} {e[0]:9.2f}{e[1]:9.2f}{e[2]:9.2f}   "
          f"{np.linalg.norm(e):8.2f}   (见上表)")

# IITD
gmsh.initialize(str(SDK / "lib"))
gmsh.option.setNumber("General.Terminal", 0)
for name, f in (("IITD talus_r", "talus_r_8.msh"), ("IITD calcaneus_r", "calcaneus_r_8.msh")):
    gmsh.clear()
    gmsh.open(str(MSH / f))
    nt, nc, _ = gmsh.model.mesh.getNodes()
    P = np.asarray(nc, float).reshape(-1, 3) * 1000.0     # m -> mm
    e = srt(pca_extents(P))
    print(f"{name:28s} {e[0]:9.2f}{e[1]:9.2f}{e[2]:9.2f}   {np.linalg.norm(e):8.2f}")
gmsh.finalize()

# IITD 原 STL（更准，未被网格化扰动）
for name, f in (("IITD talus_r (原STL)", "talus_r.stl"),
                ("IITD calcaneus_r (原STL)", "calcaneus_r.stl")):
    m = pv.read(str(IITD / f))
    tri = m.extract_surface(algorithm="dataset_surface").triangulate()
    P = np.asarray(tri.points, float) * 1000.0
    tm = trimesh.Trimesh(P, np.asarray(tri.faces, np.int64).reshape(-1, 4)[:, 1:])
    e = srt(pca_extents(P))
    print(f"{name:28s} {e[0]:9.2f}{e[1]:9.2f}{e[2]:9.2f}   "
          f"{np.linalg.norm(e):8.2f}   vol={tm.volume:10.1f} mm³")

print("\n说明：若 THUMS 主轴 ≈ IITD 主轴，则体积差应 ≈ 1.0；"
      "若主轴比 ~1.2 而体积比 ~1.8，说明形状不同而非整体缩放。")
