"""无假设裁决：一个软骨的外表面节点是否落在另一个软骨的四面体内部？（点-在-tet）"""
import re
from pathlib import Path

import numpy as np

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio")
txt = (W / "_B.feb").read_text(encoding="latin-1")
XYZ = {}
for m in re.finditer(
        r'<node id="(\d+)">\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)', txt):
    XYZ[int(m.group(1))] = np.array([float(m.group(2)), float(m.group(3)),
                                     float(m.group(4))])
els = {}
for m in re.finditer(
        r'<Elements\b[^>]*type="(tet4|hex8|penta6)"[^>]*name="([^"]+)"[^>]*>(.*?)</Elements>',
        txt, re.S):
    els[m.group(2)] = [tuple(int(y) for y in s.split(",")) for s in
                       re.findall(r"<elem\b[^>]*>\s*([\d,]+)\s*</elem>", m.group(3))]
sf = {}
for m in re.finditer(r'<Surface name="([^"]+)"[^>]*>(.*?)</Surface>', txt, re.S):
    sf[m.group(1)] = {int(y) for s in
                      re.findall(r"<[a-z0-9]+ id=[^>]*>\s*([\d,]+)\s*</", m.group(2))
                      for y in s.split(",")}
allbones = set().union(*[set(x for e in els[d] for x in e) for d in els
                         if not d.startswith("tibiotalar")])

cti = els["tibiotalar_tibia__tet4"]
cta = els["tibiotalar_talus__tet4"]
oti = sorted(sf["surf_tibiotalar_tibia"])
ota = sorted(sf["surf_tibiotalar_talus"])


def mindet(tet, P):
    """tet:(4,3); P:(n,3) -> 行列式符号 (n,)，>0 表示在同侧"""
    a, b, c, d = tet
    out = []
    for tri, opp in (((a, b, c), d), ((a, d, b), c), ((a, c, d), b), ((b, d, c), a)):
        nvec = np.cross(tri[1] - tri[0], tri[2] - tri[0])
        out.append(np.einsum("ij,j->i", P - tri[0], nvec)
                   * np.dot(opp - tri[0], nvec))
    return np.stack(out, 1)      # (n,4)


# 自检
tet0 = np.array([[0., 0, 0], [1., 0, 0], [0., 1, 0], [0., 0, 1]])
for p, exp in (([0.1, 0.1, 0.1], True), ([0.5, 0.5, 0.5], False),
               ([0.2, 0.2, 0.6], False), ([0.05, 0.05, 0.05], True)):
    s = mindet(tet0, np.array([p], float))[0]
    inside = bool((s > 0).all())
    print(f"  自检 {str(p):18s} 期望 inside={exp} 实得 {inside}  "
          f"[{'OK' if inside == exp else 'XX'}]")

print("\n" + "=" * 74)
print("① talus 软骨外表面节点 是否在 tibia 软骨 tet 内")
print("=" * 74)
Pta = np.array([XYZ[i] for i in ota])
T = np.array([np.array([XYZ[i] for i in e]) for e in cti])
TT = np.transpose(T, (0, 2, 1))
inside = np.zeros(len(Pta), bool)
depth = np.zeros(len(Pta))
for k in range(len(T)):
    mn = np.linalg.norm(T[k].min(0) - Pta, axis=1)
    mx = np.linalg.norm(T[k].max(0) - Pta, axis=1)
    near = mn < 6
    if not near.any():
        continue
    s = mindet(T[k], Pta[near])
    ins = (s > 0).all(1)
    if ins.any():
        inside[near] |= ins
print(f"  talus 外侧节点 {len(Pta)}，落在 tibia 软骨 tets 内的 {int(inside.sum())}")

print("\n" + "=" * 74)
print("② tibia 软骨外表面节点 是否在 talus 软骨 tet 内")
print("=" * 74)
Pti = np.array([XYZ[i] for i in oti])
T2 = np.array([np.array([XYZ[i] for i in e]) for e in cta])
inside2 = np.zeros(len(Pti), bool)
for k in range(len(T2)):
    mn = np.linalg.norm(T2[k].min(0) - Pti, axis=1)
    if (mn < 6).sum() == 0:
        continue
    near = mn < 6
    s = mindet(T2[k], Pti[near])
    ins = (s > 0).all(1)
    if ins.any():
        inside2[near] |= ins
print(f"  tibia 外侧节点 {len(Pti)}，落在 talus 软骨 tets 内的 {int(inside2.sum())}")

print("\n" + "=" * 74)
print("③ 交叉检验：外侧节点到「对侧软骨外面」的最近距离")
print("=" * 74)
# 外侧 ↔ 外侧
d, _ = np.array([0]), 0
from scipy.spatial import cKDTree
d1 = cKDTree(Pta).query(Pti)[0]
print(f"  tibia 外侧 → talus 外侧（节点距）：min {d1.min():.4f} "
      f"中位 {np.median(d1):.4f} max {d1.max():.4f} mm")

print("\n" + "=" * 74)
print("④ 结论")
print("=" * 74)
n_in = int(inside.sum()) + int(inside2.sum())
print(f"  互穿节点总数 = {n_in}")
print(f"  {'⇒ ❌ 存在真实干涉（点-在-tet 确认）' if n_in > 0 else '⇒ ✅ 无干涉（与 cartilage.json 一致）'}")
