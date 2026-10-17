"""真实接触面积：把 talus 关节面稠密采样，算到 tibia 关节面的有符号距离，按面积加权。
关键判据：gap<0.05mm 的**面积**是多少？≈300mm² 则正常（Anderson 生理值）。"""
import re
from pathlib import Path

import numpy as np

from closest_vec import closest_vec

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio")
FN = W / "_B.feb"
txt = FN.read_text(encoding="latin-1")
XYZ = {}
for m in re.finditer(
        r'<node id="(\d+)">\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)', txt):
    XYZ[int(m.group(1))] = np.array([float(m.group(2)), float(m.group(3)),
                                     float(m.group(4))])
sf = {}
for m in re.finditer(r'<Surface name="([^"]+)"[^>]*>(.*?)</Surface>', txt, re.S):
    sf[m.group(1)] = [tuple(int(y) for y in s.split(",")) for s in
                      re.findall(r"<[a-z0-9]+ id=[^>]*>\s*([\d,]+)\s*</", m.group(2))]
els = {}
for m in re.finditer(
        r'<Elements\b[^>]*type="(tet4|hex8|penta6)"[^>]*name="([^"]+)"[^>]*>(.*?)</Elements>',
        txt, re.S):
    els[m.group(2)] = {int(y) for s in
                       re.findall(r"<elem\b[^>]*>\s*([\d,]+)\s*</elem>", m.group(3))
                       for y in s.split(",")}


def tri_area(T):
    return 0.5 * np.linalg.norm(np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]), axis=1)


Tti = np.array([np.array([XYZ[i] for i in f]) for f in sf["surf_tibiotalar_tibia"]])
Tta = np.array([np.array([XYZ[i] for i in f]) for f in sf["surf_tibiotalar_talus"]])
tibC = np.array([XYZ[i] for i in els["tibia_cort__hex8"]]).mean(0)
print("=" * 78)
print("① 关节面面积")
print("=" * 78)
Ati, Ata = tri_area(Tti).sum(), tri_area(Tta).sum()
print(f"  tibia 侧 {len(Tti)} 三角 {Ati:.1f} mm²   talus 侧 {len(Tta)} 三角 {Ata:.1f} mm²")
print(f"  ⇒ Anderson 生理接触面积 ≈ 300 mm²（§2.4 E·A/t 反推）")

# 朝向：tibia 三角背离胫骨质心
nr = np.cross(Tti[:, 1] - Tti[:, 0], Tti[:, 2] - Tti[:, 0])
flip = np.einsum("ij,ij->i", nr, Tti.mean(1) - tibC) < 0
Tti = Tti.copy()
Tti[flip] = Tti[flip][:, ::-1]

print("\n" + "=" * 78)
print("② 稠密采样 talus 面 → 到 tibia 面的有符号距离，按面积加权")
print("=" * 78)
K = 10
g = np.linspace(0, 1, K + 1)
bc = np.array([(g[i], g[j]) for i in range(len(g)) for j in range(len(g) - i)])
sub = bc[bc[:, 0] + bc[:, 1] <= 1 + 1e-12]
# 每个子三角（顶点重心坐标）
tris_bc = []
for i in range(len(g) - 1):
    for j in range(len(g) - 1 - i):
        a = (i / K, j / K)
        b = ((i + 1) / K, j / K)
        c = (i / K, (j + 1) / K)
        d = ((i + 1) / K, (j + 1) / K)
        if i + j + 1 <= K:
            tris_bc.append((a, b, c))
        if i + 1 + j + 1 <= K:
            tris_bc.append((b, d, c))
tris_bc = np.array(tris_bc)
print(f"  每个 talus 三角细分为 {len(tris_bc)} 个子三角；"
      f"总样本三角 {len(tris_bc)*len(Tta):,}")

cen = []
area = []
for tb in tris_bc:
    pt = (tb[:, 0][:, None, None] * Tta + tb[:, 1][:, None, None] * Tta[:, 1][None]
          if False else None)
# 向量化构造子三角
B0 = (tris_bc[:, 0, 0][:, None, None] * Tta[None, :, 0, :]
      + tris_bc[:, 0, 1][:, None, None] * Tta[None, :, 1, :]
      + (1 - tris_bc[:, 0].sum(1))[:, None, None] * Tta[None, :, 2, :])
B1 = (tris_bc[:, 1, 0][:, None, None] * Tta[None, :, 0, :]
      + tris_bc[:, 1, 1][:, None, None] * Tta[None, :, 1, :]
      + (1 - tris_bc[:, 1].sum(1))[:, None, None] * Tta[None, :, 2, :])
B2 = (tris_bc[:, 2, 0][:, None, None] * Tta[None, :, 0, :]
      + tris_bc[:, 2, 1][:, None, None] * Tta[None, :, 1, :]
      + (1 - tris_bc[:, 2].sum(1))[:, None, None] * Tta[None, :, 2, :])
sub_tris = np.stack([B0, B1, B2], 2).reshape(-1, 3, 3)
areas = tri_area(sub_tris)
cents = sub_tris.mean(1)
print(f"  实际子三角 {len(sub_tris):,}  总面积 {areas.sum():.1f} mm²")

d, ti, foot, nm = closest_vec(cents, Tti)
s = np.einsum("ij,ij->i", cents - foot, nm)

print(f"\n  面积加权间隙分布（n={len(s):,}，总面积 {areas.sum():.0f} mm²）：")
order = np.argsort(s)
cum = np.cumsum(areas[order]) / areas.sum()
for th in (0.0, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0):
    A = areas[s < th].sum()
    print(f"    gap < {th:<5.2f} mm : {A:>8.1f} mm²  ({A/areas.sum()*100:5.1f}%)")

print(f"\n  间隙统计：min {s.min():+.4f}  1% {np.percentile(s,1):+.4f} "
      f"中位 {np.median(s):+.4f}  max {s.max():+.4f} mm")
neg = areas[s < 0].sum()
print(f"  侵入（gap<0）面积 = {neg:.1f} mm² ({neg/areas.sum()*100:.1f}%)")
