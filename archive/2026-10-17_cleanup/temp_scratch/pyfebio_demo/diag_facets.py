"""检查接触面 facet：面积分布 / 节点重合 / 重复面，并交叉验证 FEBio 的 invalid 警告。"""
import re
from pathlib import Path
import numpy as np

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio")
feb = (W / "ankle_tibiotalar.feb").read_text(encoding="latin-1")

log = (W / "ankle_tibiotalar.log").read_text(encoding="latin-1", errors="replace")
warn = dict(re.findall(r'surface "([^"]+)" has (\d+) invalid facets', log))
print("FEBio 报告的 invalid facet 数:", warn)

node = {}
for m in re.finditer(r"<node id=\"(\d+)\">([-0-9.eE+,]+)</node>", feb):
    node[int(m.group(1))] = np.array([float(x) for x in m.group(2).split(",")])
print(f"节点 {len(node):,}")

for name in ("surf_subtalar_calcaneus", "surf_subtalar_talus",
             "surf_tibiotalar_tibia", "surf_tibiotalar_talus"):
    m = re.search(rf'<Surface name="{name}">(.*?)</Surface>', feb, re.S)
    if not m:
        continue
    tri = re.findall(r"<tri3 id=\"\d+\">([0-9,]+)</tri3>", m.group(1))
    quads = re.findall(r"<quad4 id=\"\d+\">([0-9,]+)</quad4>", m.group(1))
    ids = [[int(x) for x in t.split(",")] for t in tri + quads]
    areas, dupnode, zero = [], 0, 0
    for f in ids:
        if len(set(f)) != len(f):
            dupnode += 1
            continue
        p = [node[i] for i in f]
        a = 0.0
        for i in range(len(p)):
            a += float(np.linalg.norm(np.cross(p[i], p[(i + 1) % len(p)]))) * 0.5
        areas.append(a)
        if a < 1e-8:
            zero += 1
    a = np.array(areas)
    keys = [tuple(sorted(f)) for f in ids]
    dupf = len(keys) - len(set(keys))
    print(f"\n{name}: tri3 {len(tri)} / quad4 {len(quads)}"
          f"  (FEBio invalid={warn.get(name, 0)})")
    print(f"  面积 min {a.min():.6e} / p1 {np.percentile(a,1):.4e} "
          f"/ 中位 {np.median(a):.4e} mm^2")
    print(f"  面积 <1e-6 的 {int((a<1e-6).sum())}  <1e-4 的 {int((a<1e-4).sum())}"
          f"  节点重合 {dupnode}  重复面 {dupf}")
    # 边：只被 1 个面使用的边 = 开放边界；>2 = 非流形
    from collections import Counter
    ec = Counter()
    for f in ids:
        for i in range(len(f)):
            e = tuple(sorted((f[i], f[(i + 1) % len(f)])))
            ec[e] += 1
    cnt = Counter(ec.values())
    print(f"  边使用次数分布 {dict(sorted(cnt.items()))}"
          f"（1=开放边界 {cnt.get(1,0)}，2=正常 {cnt.get(2,0)}，>2=非流形）")
