"""核对 tie 的 primary/secondary 到底引用了哪些节点（防"选错骨/选错区域"）。"""
import re
from pathlib import Path

import numpy as np

feb = Path(r"D:\Project\climbing_fall_analysis\temp/thums/febio_lig/ankle_lig_WK.feb")
t = feb.read_text(encoding="latin-1")
xyz = {}
for m in re.finditer(r'<node id="(\d+)">([-0-9.eE+,]+)</node>', t):
    xyz[int(m.group(1))] = np.array([float(x) for x in m.group(2).split(",")])

# 各域的节点集合（用于判断 bonesurf 的节点属于哪块骨）
dom_nodes = {}
for m in re.finditer(r'<Elements type="(\w+)" name="([^"]+)">(.*?)</Elements>',
                     t, re.S):
    name, body = m.group(2), m.group(3)
    base = name.split("__")[0]
    for e in re.finditer(r'<elem id="\d+">([0-9,]+)</elem>', body):
        dom_nodes.setdefault(base, set()).update(
            int(x) for x in e.group(1).split(","))

for name in ("tiesurf_subtalar_calcaneus", "bonesurf_subtalar_calcaneus",
             "tiesurf_subtalar_talus", "bonesurf_subtalar_talus",
             "tiesurf_tibiotalar_talus", "bonesurf_tibiotalar_talus"):
    m = re.search(rf'<Surface name="{name}">(.*?)</Surface>', t, re.S)
    if not m:
        print(f"{name}: 不存在")
        continue
    ids = [int(x) for x in re.findall(r'[0-9]+',
                                      m.group(1).replace("<tri3 id=", "")
                                      .replace(">", ","))]
    ids = sorted({int(x) for grp in re.findall(r'<(?:tri3|quad4) id="\d+">([0-9,]+)</(?:tri3|quad4)>',
                                               m.group(1)) for x in grp.split(",")})
    P = np.array([xyz[i] for i in ids])
    owner = {d: len(set(ids) & s) for d, s in dom_nodes.items()}
    owner = {d: n for d, n in owner.items() if n}
    print(f"{name}: 节点 {len(ids)}  质心 {np.round(P.mean(0),1)}  "
          f"bbox {np.round(P.min(0),1)}~{np.round(P.max(0),1)}")
    print(f"    这些节点属于哪些域: {owner}")

# tie pair 定义
for m in re.finditer(r'<SurfacePair name="([^"]+)">(.*?)</SurfacePair>', t, re.S):
    print("SurfacePair", m.group(1), re.findall(r"<(\w+)>(\w+)</\1>", m.group(2)))
