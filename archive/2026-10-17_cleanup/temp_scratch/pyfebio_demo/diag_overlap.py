"""检查两层软骨外表面初始是否过盈。
判据：对 src 侧每个 offset 点，算它到 dst 侧**骨面**的距离 d_bone，
      再取最近 dst offset 节点处的软骨厚度 t_local；
      d_bone < t_local  =>  src 已越过 dst 的 offset 面 => 戳入。
"""
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

P = Path(r"D:\Project\climbing_fall_analysis")
db = np.load(P / "temp/thums/thums_lowerlimb.npz")
nid, xyz = db["node_ids"], db["node_xyz"]
idx = {int(v): i for i, v in enumerate(nid)}
cg = np.load(P / "temp/thums/cartilage_mesh.npz")
TET_FACES = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
HEX_T = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6),
         (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]
WEDGE_T = [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]
PAIRS = [("tibiotalar_tibia", "tibiotalar_talus"),
         ("subtalar_talus", "subtalar_calcaneus")]
ALL_J = [x for pr in PAIRS for x in pr]

bone_ids = sorted({int(x) for j in ALL_J for row in cg[f"{j}_elems"]
                   for x in row if int(x) > 0})
bone_pos = {v: i for i, v in enumerate(bone_ids)}
coords = [tuple(xyz[idx[v]]) for v in bone_ids]
off_of = {}
for j in ALL_J:
    for k, p in enumerate(cg[f"{j}_offsets"]):
        off_of[(j, k)] = len(coords)
        coords.append(tuple(p))
coords = np.asarray(coords, float)
OFFSET_IDS = set(off_of.values())


def tets_of(j):
    out = []
    for row in cg[f"{j}_elems"]:
        vals = [int(x) for x in row if int(x) != 0]
        if not vals or vals[0] <= 0:
            continue
        m = len(vals) // 2
        # 骨节点要过 bone_pos 映射成 coords 下标，offset 节点已是 coords 下标
        full = [bone_pos[v] for v in vals[:m]] + \
               [off_of[(j, -x - 1)] for x in vals[m:]]
        for t in (HEX_T if m == 4 else WEDGE_T):
            out.append([full[z] for z in t])
    return out


def outer_faces(j):
    seen, fs = set(), []
    for ids in tets_of(j):
        for f in TET_FACES:
            fv = [ids[z] for z in f]
            if not all(x in OFFSET_IDS for x in fv):
                continue
            k = tuple(sorted(fv))
            if k not in seen:
                seen.add(k)
                fs.append(fv)
    return fs


def bone_of_offset(j):
    m = {}
    for row in cg[f"{j}_elems"]:
        vals = [int(x) for x in row if int(x) != 0]
        if not vals or vals[0] <= 0:
            continue
        k = len(vals) // 2
        for b, o in zip(vals[:k], [off_of[(j, -x - 1)] for x in vals[k:]]):
            m[o] = bone_pos[b]
    return m


print(f"{'接触对':36s} {'A面':>5s} {'B面':>5s} | "
      f"{'A戳入B':>8s} {'最深mm':>8s} | {'B戳入A':>8s} {'最深mm':>8s}")
for ja, jb in PAIRS:
    fa, fb = outer_faces(ja), outer_faces(jb)
    res = {}
    for src, dst, key in ((fa, fb, "A2B"), (fb, fa, "B2A")):
        dmap = bone_of_offset(jb if key == "A2B" else ja)
        dst_off = np.array([dmap[f[0]] for f in dst])
        dst_bone = coords[[dmap[f[0]] - 1 for f in dst]]
        t_dst = np.linalg.norm(coords[[f[0] - 1 for f in dst]] - dst_bone, axis=1)
        pts = coords[[x - 1 for f in src for x in f]]
        d_bone, _ = cKDTree(dst_bone).query(pts)
        _, near = cKDTree(dst_off).query(pts)
        over = t_dst[near] - d_bone
        pen = over > 1e-9
        res[key] = (int(pen.sum()), float(over[pen].max()) if pen.any() else 0.0)
    print(f"{ja+' <-> '+jb:36s} {len(fa):>5,} {len(fb):>5,} | "
          f"{res['A2B'][0]:>8,} {res['A2B'][1]:>8.4f} | "
          f"{res['B2A'][0]:>8,} {res['B2A'][1]:>8.4f}")
print("\n戳入=初始过盈；penalty 接触在 t=0 就产生大反力 => 第一个载荷步压翻单元。")
