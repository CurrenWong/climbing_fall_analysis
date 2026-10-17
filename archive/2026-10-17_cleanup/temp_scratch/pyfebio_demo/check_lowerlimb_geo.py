"""检查 Rajagopal vtp 骨几何能否用于 FE 网格化（修正单位：vtp/STL 均为米）。
"""
from pathlib import Path
import numpy as np
import pyvista as pv
import trimesh

GEO = Path(r"D:\Project\climbing_fall_analysis\model\opensim\FullBodyModel-4.0\Geometry")
IITD_GEO = Path(r"D:\Project\climbing_fall_analysis\model\Model\Geometry")
IITD = Path(r"D:\Project\climbing_fall_analysis\temp\ankle_geom")

FILES = [
    ("Raj r_tibia", GE0 := GEO / "r_tibia.vtp"),
    ("Raj r_fibula", GEO / "r_fibula.vtp"),
    ("Raj l_tibia", GEO / "l_tibia.vtp"),
    ("Raj r_femur", GEO / "r_femur.vtp"),
    ("Raj r_patella", GEO / "r_patella.vtp"),
    ("Raj r_talus", GEO / "r_talus.vtp"),
    ("Raj r_foot", GEO / "r_foot.vtp"),
    ("IITD talus_r", IITD / "talus_r.stl"),
    ("IITD calcaneus_r", IITD / "calcaneus_r.stl"),
    ("IITD OK_tibia.obj", IITD / "OK_tibia.obj"),
]

print("单位: 米 -> 下面全部换算成 cm / cm^3")
print(f"{'name':20s} {'raw_pts':>8s} {'raw_cell':>8s} {'tris':>7s} {'close':>6s} "
      f"{'vol(cm3)':>10s} {'ext x/y/z (cm)':>22s} {'euler':>6s}")
print("-" * 100)

for label, p in FILES:
    if not p.exists():
        print(f"{label:20s}  MISSING  {p}")
        continue
    try:
        m = pv.read(str(p))
        raw_pts, raw_cells = m.n_points, m.n_cells
        tri = m.extract_surface(algorithm="dataset_surface").triangulate()
        pts = np.asarray(tri.points, float)
        faces = np.asarray(tri.faces, dtype=np.int64).reshape(-1, 4)[:, 1:]
        tm = trimesh.Trimesh(vertices=pts, faces=faces, process=True)
        ext = (pts.max(0) - pts.min(0)) * 100  # m -> cm
        vol_cm3 = tm.volume * 1e6
        print(f"{label:20s} {raw_pts:8d} {raw_cells:8d} {len(faces):7d} "
              f"{str(tm.is_watertight):>6s} {vol_cm3:10.3f} "
              f"{ext[0]:7.2f}/{ext[1]:7.2f}/{ext[2]:7.2f} {tm.euler_number:6d}")
    except Exception as exc:
        print(f"{label:20s}  ERROR {type(exc).__name__}: {exc}")

print()
print("=== IITD model/Model/Geometry 内容 ===")
if IITD_GEO.exists():
    for f in sorted(IITD_GEO.iterdir()):
        print(f"  {f.stat().st_size:10d}  {f.name}")
else:
    print("  (不存在)")

print()
print("=== IITD model/Model 顶层 ===")
top = IITD_GEO.parent
for f in sorted(top.iterdir()):
    tag = "<DIR>" if f.is_dir() else f"{f.stat().st_size:>9d}"
    print(f"  {tag}  {f.name}")
