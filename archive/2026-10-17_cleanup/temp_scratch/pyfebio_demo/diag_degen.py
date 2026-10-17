"""查我的 cart_patch 输出：退化底面面（零面积 / 重复节点）+ 重复节点。"""
import numpy as np
from pathlib import Path
from collections import Counter

THUMS = Path(r"D:\Project\climbing_fall_analysis\temp\thums")
src = np.load(THUMS / "cartilage_mesh.npz", allow_pickle=True)
bpos = {int(v): src["bone_xyz"][i] for i, v in enumerate(src["bone_node_ids"])}
MIN_AREA = 1e-8

for tag, fn in (("基线 t0d", "cartilage_mesh_t0d.npz"),
                ("我的 m=1", "_pm1.npz"),
                ("我的 m=2", "cartilage_mesh_patch.npz")):
    p = THUMS / fn
    if not p.exists():
        print(f"{tag}: 缺 {fn}")
        continue
    npz = np.load(p, allow_pickle=True)
    print(f"\n===== {tag}  ({fn}) =====")
    for j in ["tibiotalar_tibia", "tibiotalar_talus"]:
        if f"{j}_elems" not in npz.files:
            continue
        nd = npz[f"{j}_nodes"]
        isoff = npz[f"{j}_is_offset"] if f"{j}_is_offset" in npz.files else None
        E = npz[f"{j}_elems"]

        def P_of(r):
            return nd[-r - 1] if r < 0 else bpos[r]

        # 重复节点：位置完全相同的节点对
        key = [tuple(np.round(q, 9)) for q in nd]
        cnt = Counter(key)
        ndup = sum(v - 1 for v in cnt.values() if v > 1)
        # 底面三角（3 节点全非 offset）的退化统计
        noff = 0 if isoff is None else int((~isoff).sum())
        ndeg_area = ndeg_dup = ntri = 0
        minA = 1e9
        for row in E:
            v = [int(x) for x in row if int(x) != 0]
            m = len(v) // 2
            base = v[:m]
            if isoff is not None and not all(not bool(isoff[-r - 1]) for r in base):
                continue
            for k in range(1, m - 1):
                tri = [base[0], base[k], base[k + 1]]
                ntri += 1
                if len(set(tri)) != len(tri):
                    ndeg_dup += 1
                    continue
                q = np.array([P_of(r) for r in tri])
                a = 0.5 * float(np.linalg.norm(np.cross(q[1] - q[0], q[2] - q[0])))
                minA = min(minA, a)
                if a <= MIN_AREA:
                    ndeg_area += 1
        print(f"  {j:20s} 节点 {len(nd):5d} 重复 {ndup:4d} | 底面三角 {ntri:5d} "
              f"零面积 {ndeg_area:3d} 重复节点 {ndeg_dup:3d} 最小面积 {minA:.3e}")
