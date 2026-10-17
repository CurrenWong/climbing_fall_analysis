"""骨面↔骨面实际间距 ⇒ 为什么 offset 会产生干涉。"""
import re
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from closest_vec import closest_vec

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio")
txt = (W / "_B.feb").read_text(encoding="latin-1")
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

allbones = set().union(*[s for d, s in els.items() if not d.startswith("tibiotalar")])
cti = els["tibiotalar_tibia__tet4"]
cta = els["tibiotalar_talus__tet4"]
# 软骨的「骨侧」节点 = 与骨共享的
bti = cti & allbones
bta = cta & allbones
# 「外侧」节点 = 剩余
oti = cti - allbones
ota = cta - allbones

print("=" * 78)
print("① 两软骨的骨侧/外侧节点")
print("=" * 78)
for nm, b, o in (("胫骨软骨", bti, oti), ("距骨软骨", bta, ota)):
    print(f"  {nm}: 骨侧 {len(b)}  外侧 {len(o)}")

print("\n" + "=" * 78)
print("② 骨面 ↔ 骨面 间距（关节处）")
print("=" * 78)
B1 = np.array([XYZ[i] for i in sorted(bti)])
B2 = np.array([XYZ[i] for i in sorted(bta)])
d, idx = cKDTree(B2).query(B1)
print(f"  胫骨骨侧 {len(B1)} 节点 ↔ 距骨骨侧 {len(B2)} 节点（节点距）")
print(f"    min {d.min():+.4f}  中位 {np.median(d):+.4f}  max {d.max():+.4f} mm")
for th in (0.1, 0.5, 1.0, 2.0, 3.0):
    print(f"    < {th:<4.1f} mm : {int((d<th).sum()):>3d}/{len(d)}")

print("\n" + "=" * 78)
print("③ 两个 offset 的距离（外侧 ↔ 外侧，即软骨外面之间）")
print("=" * 78)
O1 = np.array([XYZ[i] for i in sorted(oti)])
O2 = np.array([XYZ[i] for i in sorted(ota)])
d2, _ = cKDTree(O2).query(O1)
print(f"    min {d2.min():+.4f}  中位 {np.median(d2):+.4f}  max {d2.max():+.4f} mm")

print("\n" + "=" * 78)
print("④ 一致性检验：外面间距 ≈ 骨面间距 − 两软骨厚度和？")
print("=" * 78)
# 厚度：骨侧节点 → 外侧节点（用最近骨侧节点近似）
dti = cKDTree(B1).query(O1)[0]
dta = cKDTree(B2).query(O2)[0]
print(f"  胫骨软骨厚度（外侧→最近骨侧）：中位 {np.median(dti):.3f}  "
      f"[{dti.min():.3f},{dti.max():.3f}]")
print(f"  距骨软骨厚度：中位 {np.median(dta):.3f} [{dta.min():.3f},{dta.max():.3f}]")
print(f"  厚度和（中位）= {np.median(dti)+np.median(dta):.3f} mm")
print(f"  骨面间距（中位）= {np.median(d):.3f} mm")
print(f"  ⇒ 预期外面间距 = {np.median(d)-np.median(dti)-np.median(dta):+.3f} mm"
      f"   （实测外面间距中位 {np.median(d2):+.3f} mm）")

print("\n" + "=" * 78)
print("⑤ 判据")
print("=" * 78)
print(f"  骨面最小间距 = {d.min():.4f} mm；两软骨厚度和 ≈ "
      f"{np.median(dti)+np.median(dta):.3f} mm")
if d.min() < np.median(dti) + np.median(dta) - 0.5:
    print("  ⇒ ❌ **骨面本来就靠得比「两软骨厚度之和」更近 ⇒ offset 必然造成干涉**")
    print(f"     ⇒ 干涉深度 ≈ {(np.median(dti)+np.median(dta)) - d.min():.3f} mm")
    print("     ⇒ THUMS 的骨面很可能**已含软骨**（或尚未摆到负重位）")
else:
    print("  ⇒ 骨面间距够，offset 合理")
