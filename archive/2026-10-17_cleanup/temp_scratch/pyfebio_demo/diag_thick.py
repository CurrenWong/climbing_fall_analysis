"""几何核对：软骨厚度 vs 骨面间距 ⇒ 两个 offset 到底合不合得上。"""
import re
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio")
txt = (W / "_B.feb").read_text(encoding="latin-1")
els, XYZ, sf = {}, {}, {}
for m in re.finditer(
        r'<Elements\b[^>]*type="(tet4|hex8|penta6)"[^>]*name="([^"]+)"[^>]*>(.*?)</Elements>',
        txt, re.S):
    els[m.group(2)] = {int(y) for s in
                       re.findall(r"<elem\b[^>]*>\s*([\d,]+)\s*</elem>", m.group(3))
                       for y in s.split(",")}
for m in re.finditer(
        r'<node id="(\d+)">\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)', txt):
    XYZ[int(m.group(1))] = np.array([float(m.group(2)), float(m.group(3)),
                                     float(m.group(4))])
for m in re.finditer(r'<Surface name="([^"]+)"[^>]*>(.*?)</Surface>', txt, re.S):
    sf[m.group(1)] = {int(y) for s in
                      re.findall(r"<[a-z0-9]+ id=[^>]*>\s*([\d,]+)\s*</", m.group(2))
                      for y in s.split(",")}

bones = set().union(*[s for d, s in els.items() if not d.startswith("tibiotalar")])

print("=" * 76)
print("① 各软骨的厚度场（offset 面节点 → 最近骨节点）")
print("=" * 76)
for side, cart, outer in (("胫骨", "tibiotalar_tibia__tet4", "surf_tibiotalar_tibia"),
                          ("距骨", "tibiotalar_talus__tet4", "surf_tibiotalar_talus")):
    cn = els[cart]
    share = cn & bones
    B = np.array([XYZ[i] for i in sorted(share)])
    O = np.array([XYZ[i] for i in sorted(sf[outer])])
    d, _ = cKDTree(B).query(O)
    print(f"  {side}: 软骨 {len(cn)} 节点（共享骨 {len(share)}），"
          f"offset 面 {len(O)} 节点")
    print(f"      厚度（offset → 骨面节点距）：min {d.min():.4f} 中位 "
          f"{np.median(d):.4f} p95 {np.percentile(d,95):.4f} max {d.max():.4f} mm")

print("\n" + "=" * 76)
print("② 骨面之间最小距离（关节处）⇒ 与两软骨厚度之和比较")
print("=" * 76)
AX = np.array([0.8153, 0.135, 0.5631])
AX = AX / np.linalg.norm(AX)
tb = np.array([XYZ[i] for i in sorted(els["tibia_cort__hex8"])])
ta = np.array([XYZ[i] for i in sorted(els["talus_cort__hex8"])])
# 只看关节区（轴位置最远端 15mm）
pr = tb @ AX
tbj = tb[pr > np.percentile(pr, 92)]
prt = ta @ AX
taj = ta[prt < np.percentile(prt, 12)]
d, idx = cKDTree(taj).query(tbj)
print(f"  胫骨关节区 {len(tbj)} 节点 / 距骨关节区 {len(taj)} 节点")
print(f"  骨面→骨面距离：min {d.min():.4f} p5 {np.percentile(d,5):.4f} "
      f"中位 {np.median(d):.4f} mm")
print(f"  ⇒ 若两软骨厚度之和 ≈ 1.04+? 则两软骨外面应在该位置接触")

print("\n" + "=" * 76)
print("③ 两软骨外面（真实间隙）在场里的分布（t=0）")
print("=" * 76)
P = np.array([XYZ[i] for i in sorted(sf["surf_tibiotalar_tibia"])])
Q = np.array([XYZ[i] for i in sorted(sf["surf_tibiotalar_talus"])])
d2, _ = cKDTree(Q).query(P)
print(f"  外面-外面（节点距，上界）：min {d2.min():.4f} 中位 {np.median(d2):.4f} "
      f"max {d2.max():.4f} mm")
print(f"  <0.1mm 的 {int((d2<0.1).sum())}/{len(P)}；<1.0mm 的 {int((d2<1.0).sum())}")
# 关键：外面间隙 是否 ≈ 骨间隙 − 两厚度
print("\n④ 结论判据")
