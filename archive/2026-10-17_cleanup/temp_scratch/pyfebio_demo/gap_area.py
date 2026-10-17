"""无定向假设的间隙测量：talus 外面稠密采样 → 到 tibia 外面的**无符号**距离，按面积加权。"""
import re
from pathlib import Path

import numpy as np

from closest_vec import closest_vec

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio")


def load(tag):
    txt = (W / f"_{tag}.feb").read_text(encoding="latin-1")
    XYZ = {}
    for m in re.finditer(
            r'<node id="(\d+)">\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)', txt):
        XYZ[int(m.group(1))] = np.array([float(m.group(2)), float(m.group(3)),
                                         float(m.group(4))])
    sf = {}
    for m in re.finditer(r'<Surface name="([^"]+)"[^>]*>(.*?)</Surface>', txt, re.S):
        sf[m.group(1)] = [tuple(int(y) for y in s.split(",")) for s in
                          re.findall(r"<[a-z0-9]+ id=[^>]*>\s*([\d,]+)\s*</", m.group(2))]
    return XYZ, sf


def tri_area(T):
    return 0.5 * np.linalg.norm(np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]), axis=1)


K = 10
g = np.linspace(0, 1, K + 1)
tris_bc = []
for i in range(K):
    for j in range(K - i):
        a = (i / K, j / K)
        b = ((i + 1) / K, j / K)
        c = (i / K, (j + 1) / K)
        d = ((i + 1) / K, (j + 1) / K)
        if i + j + 1 <= K:
            tris_bc.append((a, b, c))
        if i + j + 2 <= K:
            tris_bc.append((b, d, c))
tris_bc = np.array(tris_bc)


def subdivide(T):
    B0 = (tris_bc[:, 0, 0][:, None, None] * T[None, :, 0, :]
          + tris_bc[:, 0, 1][:, None, None] * T[None, :, 1, :]
          + (1 - tris_bc[:, 0].sum(1))[:, None, None] * T[None, :, 2, :])
    B1 = (tris_bc[:, 1, 0][:, None, None] * T[None, :, 0, :]
          + tris_bc[:, 1, 1][:, None, None] * T[None, :, 1, :]
          + (1 - tris_bc[:, 1].sum(1))[:, None, None] * T[None, :, 2, :])
    B2 = (tris_bc[:, 2, 0][:, None, None] * T[None, :, 0, :]
          + tris_bc[:, 2, 1][:, None, None] * T[None, :, 1, :]
          + (1 - tris_bc[:, 2].sum(1))[:, None, None] * T[None, :, 2, :])
    return np.stack([B0, B1, B2], 2).reshape(-1, 3, 3)


for tag, lbl in (("B", "基线 t=0"), ("E", "m2 t=0")):
    fp = W / f"_{tag}.feb"
    if not fp.exists():
        continue
    XYZ, sf = load(tag)
    Tti = np.array([np.array([XYZ[i] for i in f]) for f in sf["surf_tibiotalar_tibia"]])
    Tta = np.array([np.array([XYZ[i] for i in f]) for f in sf["surf_tibiotalar_talus"]])
    sub = subdivide(Tta)
    ar = tri_area(sub)
    cen = sub.mean(1)
    d = closest_vec(cen, Tti)[0]
    Ati, Ata = tri_area(Tti).sum(), tri_area(Tta).sum()
    print("\n" + "=" * 76)
    print(f"{tag} {lbl}: tibia 面 {len(Tti)} 三角 {Ati:.0f} mm² / "
          f"talus 面 {len(Tta)} 三角 {Ata:.0f} mm²")
    print("=" * 76)
    print(f"  采样 {len(cen):,} 点，面积 {ar.sum():.0f} mm²")
    print(f"  面积加权间隙（无符号）：")
    for th in (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0):
        A = ar[d < th].sum()
        print(f"    gap < {th:<5.2f} mm : {A:>8.1f} mm² ({A/ar.sum()*100:5.1f}%)")
    print(f"  统计：min {d.min():.4f}  1% {np.percentile(d,1):.4f}  "
          f"中位 {np.median(d):.4f}  max {d.max():.4f} mm")
    print(f"  ⇒ 参照 Anderson 生理接触面积 ≈ 300 mm²")
