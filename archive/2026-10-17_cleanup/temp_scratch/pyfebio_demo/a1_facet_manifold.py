"""A1 面片合法性预检（只测不改）：2-流形性 / 面积 / 重复面 / FEBio invalid 警告。
正对照 = tibiotalar 的 surf_/tiesurf_（已在跑、能收敛）。"""
import re
from collections import Counter
from pathlib import Path

import numpy as np

W = Path(r"D:\Project\climbing_fall_analysis\temp/thums/febio_lig")
feb = (W / "ankle_lig_WK.feb").read_text(encoding="latin-1")
log = (W / "febio_WK.log").read_text(encoding="latin-1", errors="replace")
warn = dict(re.findall(r'surface "([^"]+)" has (\d+) invalid facets', log))
print("FEBio invalid facet 报告:", warn)

node = {}
for m in re.finditer(r'<node id="(\d+)">([-0-9.eE+,]+)</node>', feb):
    node[int(m.group(1))] = np.array([float(x) for x in m.group(2).split(",")])
print(f"节点 {len(node):,}")

NAMES = ["surf_tibiotalar_tibia", "surf_tibiotalar_talus",
         "surf_subtalar_talus", "surf_subtalar_calcaneus",
         "tiesurf_tibiotalar_tibia", "tiesurf_tibiotalar_talus",
         "tiesurf_subtalar_talus", "tiesurf_subtalar_calcaneus"]
for name in NAMES:
    m = re.search(rf'<Surface name="{name}">(.*?)</Surface>', feb, re.S)
    if not m:
        print(f"\n{name}: ❌ .feb 里没有该 Surface")
        continue
    tri = re.findall(r'<tri3 id="\d+">([0-9,]+)</tri3>', m.group(1))
    quad = re.findall(r'<quad4 id="\d+">([0-9,]+)</quad4>', m.group(1))
    ids = [[int(x) for x in t.split(",")] for t in tri + quad]
    areas, dupnode = [], 0
    for f in ids:
        if len(set(f)) != len(f):
            dupnode += 1
            continue
        p = [node[i] for i in f]
        a = 0.0
        for i in range(len(p)):
            a += float(np.linalg.norm(np.cross(p[i], p[(i + 1) % len(p)]))) * 0.5
        areas.append(a)
    a = np.array(areas) if areas else np.array([0.0])
    keys = [tuple(sorted(f)) for f in ids]
    dupf = len(keys) - len(set(keys))
    ec = Counter()
    for f in ids:
        for i in range(len(f)):
            ec[tuple(sorted((f[i], f[(i + 1) % len(f)])))] += 1
    c = Counter(ec.values())
    print(f"\n{name}: 面 {len(ids)}（tri3 {len(tri)} / quad4 {len(quad)}）"
          f"  FEBio invalid={warn.get(name, '—')}")
    print(f"  面积 min {a.min():.3e} / 中位 {np.median(a):.3e} / "
          f"max {a.max():.3e} mm²   节点重合 {dupnode}  重复面 {dupf}")
    print(f"  边使用次数 {dict(sorted(c.items()))}  "
          f"（1=开放边界 {c.get(1,0)}｜2=正常 {c.get(2,0)}｜>2=非流形 "
          f"{sum(v for k, v in c.items() if k > 2)}）")
