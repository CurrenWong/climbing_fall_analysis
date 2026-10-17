"""查：bone_faces 里为什么有 30/93 个在 elem_blocks 里找不到对应 facet。

也从侧面确认 tile 的 secondary 面到底涵盖多少骨关节面。
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
THUMS = ROOT / "temp/thums"
sys.path.insert(0, str(ROOT / "scripts/ankle_fe"))

cg = np.load(THUMS / "cartilage_mesh_t0d.npz", allow_pickle=True)
db = np.load(THUMS / "thums_lowerlimb.npz")
nid_all, xyz_all = db["node_ids"], db["node_xyz"]
idx_all = {int(v): i for i, v in enumerate(nid_all)}

BONE_PARTS = {
    "tibia_cort": ("81000700_solid",),
    "tibia_spon": ("81000600_solid", "81000601_solid"),
    "fibula_cort": ("81000900_solid",),
    "talus_cort": ("81001100_solid",),
    "talus_spon": ("81001000_solid",),
    "calcaneus_cort": ("81001300_solid",),
}
# 直接从 npz 取骨单元（不做 Jacobian 过滤，只看 facet 存不存在）
TET_FACES = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
HEX_FACES = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
             (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
WEDGE_FACES = [(0, 2, 1), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)]

keys = {k for k in db.files if k.endswith("_solid")}
facets = {}
sizes = {}
for k in keys:
    for row in db[k]:
        v = [int(x) for x in row if int(x) != 0]
        if len(v) == 8:
            tmpl = HEX_FACES
        elif len(v) == 6:
            tmpl = WEDGE_FACES
        elif len(v) == 4:
            tmpl = TET_FACES
        else:
            continue
        for f in tmpl:
            fv = tuple(sorted(v[z] for z in f))
            facets.setdefault(fv, len(f))
            sizes[len(v)] = sizes.get(len(v), 0) + 1
print("npz 骨单元按节点数:", sizes)
print(f"npz 骨 facet 总数 {len(facets):,}；"
      f"其中四边形 {sum(1 for s in facets.values() if s == 4):,} / "
      f"三角 {sum(1 for s in facets.values() if s == 3):,}")

for j in ("tibiotalar_tibia", "tibiotalar_talus"):
    bf = np.asarray(cg[f"{j}_bone_faces"], dtype=np.int64)
    ok = miss3 = miss4 = 0
    widths = {}
    for row in bf:
        v = [int(x) for x in row if int(x) != 0]
        widths[len(v)] = widths.get(len(v), 0) + 1
        k = tuple(sorted(v))
        if k in facets:
            ok += 1
        elif len(v) == 4:
            miss4 += 1
        else:
            miss3 += 1
    print(f"\n[{j}] bone_faces {len(bf)}  宽度 {widths}")
    print(f"   在 npz 骨 facet 里找到 {ok} / 四边形缺 {miss4} / 三角缺 {miss3}")
    # 这些缺的，是不是本来就该是 hex8 的四边形而 npz 里有？
    for row in bf[:0]:
        pass
