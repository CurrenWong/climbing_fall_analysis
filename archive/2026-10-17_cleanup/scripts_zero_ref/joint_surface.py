"""B1/B2: 提取 THUMS 骨的真实外表面，并测量关节间隙分布。

为什么必须"测量"而不是"设定"软骨厚度：
  软骨填充的是骨与骨之间的间隙。间隙多大，软骨层就必须多厚 ——
  凭空取 2 mm 可能让两层软骨互相穿透（或悬空）。
  吴恺用 Geomagic "加厚" 生成软骨，本质也是贴合间隙。

外边界提取方法：
  把 CORT ∪ SPON 的所有单元面收集起来，**只出现 1 次的面 = 并集外表面**
  （CORT 与 SPON 之间的界面面会出现 2 次，被自动排除）。
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

NPZ = Path(r"D:\Project\climbing_fall_analysis\temp\thums\thums_lowerlimb.npz")
OUT = Path(r"D:\Project\climbing_fall_analysis\temp\thums")

# 单元类型 -> 面（用局部节点序号；退化面会自动坍缩）
TET_FACES = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
HEX_FACES = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
             (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
WEDGE_FACES = [(0, 2, 1), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)]

GROUPS = {
    "tibia":     ["81000700_solid", "81000600_solid", "81000601_solid"],
    "fibula":    ["81000900_solid", "81000800_solid", "81000801_solid"],
    "talus":     ["81001100_solid", "81001000_solid"],
    "calcaneus": ["81001300_solid", "81001200_solid"],
}

d = np.load(NPZ)


def union_boundary(keys: list[str]):
    """返回 (外表面面片列表 [每项是节点 index 元组], 该组用到的节点 index 集合)"""
    faces: Counter = Counter()
    used: set[int] = set()
    for k in keys:
        arr = d[k]
        for row in arr:
            nodes = [int(x) for x in row]
            uniq = sorted(set(nodes))
            used.update(uniq)
            if len(uniq) == 4:
                tmpl = TET_FACES
            elif len(uniq) == 8:
                tmpl = HEX_FACES
            elif len(uniq) == 6:
                tmpl = WEDGE_FACES
            else:
                continue
            # 把局部模板映射到"去重后"的节点序列上
            for f in tmpl:
                try:
                    face = tuple(sorted(uniq[z] for z in f))
                except IndexError:
                    continue
                if len(set(face)) >= 3:
                    faces[face] += 1
    outer = [f for f, c in faces.items() if c == 1]
    return outer, used


print("=== B1 并集外表面提取 ===")
surf_nodes: dict[str, set[int]] = {}
for g, keys in GROUPS.items():
    outer, used = union_boundary(keys)
    surf_nodes[g] = set().union(*[set(f) for f in outer]) if outer else set()
    tris = sum(1 for f in outer if len(f) == 3)
    quads = sum(1 for f in outer if len(f) == 4)
    print(f"  {g:10s} 单元 {len(keys)} 组 | 外表面面片 {len(outer):>6,} "
          f"(tri {tris:>6,} / quad {quads:>6,}) | 外表面节点 {len(surf_nodes[g]):>6,} "
          f"| 部件总节点 {len(used):>6,}")

# ── 节点坐标查表 ──
nid = d["node_ids"]
xyz = d["node_xyz"]
idx = {int(v): i for i, v in enumerate(nid)}
P = {g: xyz[np.asarray(sorted(s), dtype=np.int64).astype(int)] if False else
        xyz[np.fromiter((idx[int(v)] for v in s), dtype=np.int64, count=len(s))]
     for g, s in surf_nodes.items()}
ID = {g: np.fromiter(sorted(s), dtype=np.int64, count=len(s)) for g, s in surf_nodes.items()}

print("\n=== B2 关节间隙测量（外表面节点间最近距离）===")
PAIRS = [("tibia", "talus", "踝关节（胫骨下关节面）"),
         ("fibula", "talus", "下胫腓/外踝-距骨"),
         ("talus", "calcaneus", "距下关节"),
         ("tibia", "fibula", "下胫腓联合")]

results = {}
for a, b, label in PAIRS:
    A, B = P[a], P[b]
    tree = cKDTree(B)
    dist, j = tree.query(A, k=1)
    # 双向取最小值更稳
    treeA = cKDTree(A)
    dist2, _ = treeA.query(B, k=1)
    both = np.concatenate([dist, dist2])
    qs = np.percentile(both, [0, 1, 5, 25, 50, 75, 95, 100])
    n_close = int((both < 5.0).sum())
    results[f"{a}__{b}"] = {
        "label": label, "n_a": len(A), "n_b": len(B),
        "min_mm": float(both.min()),
        "p1": float(qs[1]), "p5": float(qs[2]), "p25": float(qs[3]),
        "median": float(qs[4]), "p75": float(qs[5]), "p95": float(qs[6]),
        "n_within_5mm": n_close,
    }
    print(f"\n  {label}  ({a} {len(A)} 面节点 <-> {b} {len(B)} 面节点)")
    print(f"    最近 {both.min():7.3f} mm | p1 {qs[1]:7.3f} | p5 {qs[2]:7.3f} | "
          f"中位 {qs[4]:7.3f} | p95 {qs[6]:7.3f} | max {qs[7]:8.2f}")
    print(f"    <5 mm 的节点对数: {n_close:,}")

# ── 关节面区域规模（间隙 < 阈值）──
print("\n=== 候选关节面规模（对某一侧骨，其表面节点距对侧 < 阈值 的数量）===")
print(f"{'pair':28s} " + " ".join(f"{f'<{t}mm':>10s}" for t in (1.5, 2, 3, 5)))
region = {}
for a, b, label in PAIRS:
    A, B = P[a], P[b]
    tree = cKDTree(B)
    dist, _ = tree.query(A, k=1)
    row = []
    for t in (1.5, 2, 3, 5):
        n = int((dist < t).sum())
        row.append(n)
    region[f"{a}__{b}"] = {"a": a, "b": b, "thresholds_mm": [1.5, 2, 3, 5],
                           "counts_on_a": row, "n_a": len(A)}
    print(f"{label:28s} " + " ".join(f"{n:>10,}" for n in row))

(OUT / "joint_gaps.json").write_text(
    json.dumps({"pairs": results, "regions": region}, ensure_ascii=False, indent=2),
    encoding="utf-8")
print(f"\n写出 {OUT/'joint_gaps.json'}")
