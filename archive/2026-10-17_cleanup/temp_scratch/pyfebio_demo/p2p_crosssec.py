"""P2' 关键：加载面是否只覆盖截面的一小片（47.8 mm² vs 真实胫骨截面）。"""
import re
from pathlib import Path

import numpy as np
import h5py

R = Path(r"D:\Project\climbing_fall_analysis")
H5 = R / "temp/thums/febio_lig/ankle_lig_WK_stress.h5"
FEB = R / "temp/thums/febio_lig/ankle_lig_WK.feb"

with h5py.File(H5, "r") as f:
    nd = f["meshes/0/nodes"][:]
    nid = np.asarray(nd["id"], int)
    X = np.column_stack([nd["x"], nd["y"], nd["z"]]).astype(float)
    print("=== nodesets ===")
    ns = {}
    for k in f["meshes/0/nodesets"].keys():
        v = f["meshes/0/nodesets"][k]
        a = np.asarray(v[:], int).ravel()
        ns[k] = a
        print(f"  {k:<24} {len(a):>5} 项  id {a.min()}..{a.max()}")
    domkeys = list(f["meshes/0/domains"].keys())

pos = {int(v): i for i, v in enumerate(nid)}

# .feb 里的 Surface（加载面）
txt = FEB.read_text(encoding="latin-1", errors="replace")
print("\n=== .feb 的 <Surface> ===")
for m in re.finditer(r'<Surface name="([^"]+)"[^>]*>(.*?)</Surface>', txt, re.S):
    nm, body = m.group(1), m.group(2)
    faces = re.findall(r'<[a-z0-9]+ id="\d+">\s*([\d,]+)\s*</', body)
    print(f"  {nm:<28} {len(faces)} 面")

# 域 -> 节点
DOMS = {}
for m in re.finditer(r'<Elements\s+type="(\w+)"[^>]*?name="([^"]*)"[^>]*>'
                     r'(.*?)</Elements>', txt, re.S):
    kind, nm, body = m.group(1), m.group(2), m.group(3)
    conn = [[int(x) for x in c.split(",")] for c in
            re.findall(r'<elem[^>]*>\s*([\d,\s]+?)\s*</elem>', body)]
    DOMS[nm] = (kind, [c for c in conn if len(c) >= 3])


def cross_area_at_z(dom, zlo):
    """域内 z>=zlo 的节点的凸性近似：用 xy 投影的凸包面积。"""
    nodes = np.unique([n for c in dom[1] for n in c])
    idx = [pos[n] for n in nodes if n in pos]
    P = X[idx]
    sel = P[P[:, 2] >= zlo]
    if len(sel) < 3:
        return 0.0, len(sel)
    xy = sel[:, :2]
    try:
        from scipy.spatial import ConvexHull
        return ConvexHull(xy).volume, len(sel)
    except Exception:
        return float("nan"), len(sel)


print("\n=== 胫骨/腓骨 截面（z 越高越靠近截断顶）===")
for nm in ("tibia_cort__hex8", "tibia_spon__tet4", "tibia_marrow__tet4",
           "fibula_cort__hex8", "fibula_spon__tet4"):
    if nm not in DOMS:
        continue
    nodes = np.unique([n for c in DOMS[nm][1] for n in c])
    idx = [pos[n] for n in nodes if n in pos]
    P = X[idx]
    zmax = P[:, 2].max()
    print(f"\n  {nm}: {len(idx)} 节点  z {P[:,2].min():.2f}..{zmax:.2f}")
    for frac in (1.0, 0.97, 0.95, 0.9):
        zlo = zmax * frac if zmax < 0 else zmax * frac
        # 用绝对高度：顶面 zmax，往下 2 / 5 / 10 mm
        pass
    for dz in (1.0, 3.0, 6.0, 12.0):
        a, n = cross_area_at_z(DOMS[nm], zmax - dz)
        print(f"     顶面下 {dz:>5.1f} mm：凸包面积 {a:>8.1f} mm²（{n} 节点）")

print("\n=== 对照：模型报的加载面 ===")
print("  [模型输出] 顶面投影面积：胫骨 47.8 mm²，腓骨 12.7 mm²")
print("  真实胫骨远端截面参考：≈ 250~450 mm²")
print("  ⇒ 若凸包面积 >> 47.8，说明 BC 只作用在截面一小片")
