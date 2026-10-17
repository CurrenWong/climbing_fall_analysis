"""只读：统计各 part 的「面匹配」情况——决定能否对 hex / tet 各自独立加密。

对每个单元集合，把每个单元的每个面取出来（按节点集合规范化），
统计「同一面出现 2 次」= 面匹配（共形）；出现 1 次 = 边界面。
再分别统计：hex-hex 面匹配数、tet-tet 面匹配数、hex-tet 混合面数。
"""
from pathlib import Path
from collections import defaultdict
import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
d = np.load(ROOT / "temp/thums/thums_lowerlimb.npz", allow_pickle=True)

HEX_FACES = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
             (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
WEDGE_FACES = [(0, 2, 1), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)]
TET_FACES = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]

KEYS = [k for k in d.files if k.endswith("_solid")]
print(f"solid parts: {len(KEYS)}")

face_owner = defaultdict(list)   # 规范面 -> [(part, 类型)]
face_kind = defaultdict(set)     # 规范面 -> {3,4} 面节点数

for k in KEYS:
    for row in d[k]:
        uq = sorted({int(x) for x in row if x > 0})
        if len(uq) == 8:
            tmpl, et = HEX_FACES, "hex8"
        elif len(uq) == 6:
            tmpl, et = WEDGE_FACES, "wedge6"
        elif len(uq) == 4:
            tmpl, et = TET_FACES, "tet4"
        else:
            continue
        for f in tmpl:
            if max(f) >= len(uq):
                continue
            fs = tuple(sorted(uq[z] for z in f))
            face_owner[fs].append((k, et))
            face_kind[fs].add(len(fs))

cnt = defaultdict(int)
for fs, owners in face_owner.items():
    n = len(owners)
    if n == 1:
        cnt["边界面(1次)"] += 1
    elif n == 2:
        types = {o[1] for o in owners}
        if len(types) == 1:
            cnt[f"面匹配 {types.pop()}-{next(iter({o[1] for o in owners}))}"] += 1
        else:
            cnt["面匹配 **混合** " + "/".join(sorted(types))] += 1
    else:
        cnt[f"⚠️出现{n}次"] += 1

print("\n=== 面统计（全 solid）===")
for k in sorted(cnt):
    print(f"  {k:40s} {cnt[k]:>8,}")

# 只统计踝部 4 块骨
ANK = ("81000700_solid", "81000600_solid", "81000601_solid",
       "81001100_solid", "81001000_solid",
       "81001300_solid", "81001200_solid",
       "81000900_solid", "81000800_solid", "81000801_solid")
face_owner2 = defaultdict(list)
for k in ANK:
    for row in d[k]:
        uq = sorted({int(x) for x in row if x > 0})
        if len(uq) == 8:
            tmpl, et = HEX_FACES, "hex8"
        elif len(uq) == 6:
            tmpl, et = WEDGE_FACES, "wedge6"
        elif len(uq) == 4:
            tmpl, et = TET_FACES, "tet4"
        else:
            continue
        for f in tmpl:
            if max(f) >= len(uq):
                continue
            fs = tuple(sorted(uq[z] for z in f))
            face_owner2[fs].append((k, et))
cnt2 = defaultdict(int)
for fs, owners in face_owner2.items():
    n = len(owners)
    if n == 1:
        cnt2["边界面"] += 1
    elif n == 2:
        types = sorted({o[1] for o in owners})
        cnt2["面匹配 " + "/".join(types)] += 1
    else:
        cnt2[f"⚠️出现{n}次"] += 1
print("\n=== 面统计（仅踝部 10 个 solid）===")
for k in sorted(cnt2):
    print(f"  {k:40s} {cnt2[k]:>8,}")
