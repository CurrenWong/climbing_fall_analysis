"""P2g：用 .feb 连通性 + h5 位移，逐域定位柔度。"""
import re
from pathlib import Path

import numpy as np
import h5py

R = Path(r"D:\Project\climbing_fall_analysis")
H5 = R / "temp/thums/febio_lig/ankle_lig_WK_stress.h5"
FEB = R / "temp/thums/febio_lig/ankle_lig_WK.feb"
F = 592.602

with h5py.File(H5, "r") as f:
    nd = f["meshes/0/nodes"][:]
    nid = np.asarray(nd["id"], int)
    X = np.column_stack([nd["x"], nd["y"], nd["z"]]).astype(float)
    sids = sorted(f["states"].keys(), key=int)
    D = np.asarray(f["states"][sids[-1]]["node_data"]["displacement"]["all"][:],
                   float)
pos = {int(v): i for i, v in enumerate(nid)}
dz = D[:, 2]

txt = FEB.read_text(encoding="latin-1", errors="replace")
DOMS = {}
for m in re.finditer(r'<Elements\s+type="(\w+)"[^>]*?name="([^"]*)"[^>]*>'
                     r'(.*?)</Elements>', txt, re.S):
    kind, nm, body = m.group(1), m.group(2), m.group(3)
    conn = [[int(x) for x in c.split(",")] for c in
            re.findall(r'<elem[^>]*>\s*([\d,\s]+?)\s*</elem>', body)]
    DOMS[nm] = (kind, [c for c in conn if len(c) >= 3])

print("=" * 108)
print("逐域：z 跨度 / dz 跨度 / 隐含轴向应变 / |u|max")
print("=" * 108)
print(f"{'域':<34}{'类型':>7}{'节点':>7}{'z跨度':>9}{'dz跨度':>10}"
      f"{'dz均值':>10}{'εz(‰)':>10}{'|u|max':>10}")
rows = []
for nm, (kind, conn) in sorted(DOMS.items()):
    nodes = np.unique([n for c in conn for n in c])
    idx = np.array([pos[n] for n in nodes if n in pos])
    if len(idx) < 3:
        continue
    zz, dd = X[idx, 2], dz[idx]
    zspan = zz.max() - zz.min()
    dspan = dd.max() - dd.min()
    eps = dspan / zspan * 1000 if zspan > 1e-9 else float("nan")
    umax = np.linalg.norm(D[idx], axis=1).max()
    rows.append((nm, kind, len(idx), zspan, dspan, dd.mean(), eps, umax))
    print(f"{nm:<34}{kind:>7}{len(idx):>7}{zspan:>9.2f}{dspan:>10.4f}"
          f"{dd.mean():>10.4f}{eps:>10.2f}{umax:>10.4f}")

print("\n" + "=" * 108)
print("按 dz 跨度排序（= 该域贡献的变形量）")
print("=" * 108)
for r in sorted(rows, key=lambda x: -x[4])[:8]:
    print(f"  {r[0]:<34} dz跨度 {r[4]:>8.4f} mm  (z跨度 {r[3]:>7.2f} mm, "
          f"εz {r[6]:>7.2f}‰)")

print("\n" + "=" * 108)
print("👉 柔度在哪")
print("=" * 108)
tot = sorted(rows, key=lambda x: -x[4])
print(f"  总指定位移 1.6000 mm；关节（tibiotalar_*）只闭 0.0142 mm")
print(f"  dz 跨度最大的域：{tot[0][0]} = {tot[0][4]:.4f} mm")
print(f"  第二：{tot[1][0]} = {tot[1][4]:.4f} mm")
