"""Maharaj 自带足部 STL 是否被降采样？（决定它能否当 FE 几何源）"""
from pathlib import Path
import numpy as np
import pyvista as pv
import trimesh

MAH = Path(r"D:\Project\climbing_fall_analysis\model\FootwithinLowerlimbModel")
IITD = Path(r"D:\Project\climbing_fall_analysis\temp\ankle_geom")

PAIRS = [
    ("talus_r",    MAH / "FootModel_talus_r.stl",    IITD / "talus_r.stl"),
    ("calcn_r",    MAH / "FootModel_calcn_r.stl",    IITD / "calcaneus_r.stl"),
    ("midfoot_r",  MAH / "FootModel_midfoot_r.stl",  IITD / "midfoot_r.stl"),
]

def info(p):
    m = pv.read(str(p))
    tri = m.extract_surface(algorithm="dataset_surface").triangulate()
    pts = np.asarray(tri.points, float)
    faces = np.asarray(tri.faces, dtype=np.int64).reshape(-1, 4)[:, 1:]
    tm = trimesh.Trimesh(vertices=pts, faces=faces, process=True)
    ext = (pts.max(0) - pts.min(0)) * 100
    return len(pts), len(faces), tm.is_watertight, tm.volume * 1e6, ext

print(f"{'bone':10s} {'source':12s} {'bytes':>9s} {'pts':>7s} {'tris':>7s} "
      f"{'close':>6s} {'vol(cm3)':>9s}  ext x/y/z (cm)")
print("-" * 96)
for name, mp, ip in PAIRS:
    for tag, p in (("Maharaj", mp), ("IITD原版", ip)):
        if not p.exists():
            print(f"{name:10s} {tag:12s}  MISSING {p}")
            continue
        pts, tris, wt, vol, ext = info(p)
        print(f"{name:10s} {tag:12s} {p.stat().st_size:9d} {pts:7d} {tris:7d} "
              f"{str(wt):>6s} {vol:9.3f}  {ext[0]:5.2f}/{ext[1]:5.2f}/{ext[2]:5.2f}")
    print()
