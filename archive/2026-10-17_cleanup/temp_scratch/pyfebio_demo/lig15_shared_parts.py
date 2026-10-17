"""诊断：每条韧带与哪些部件共用节点？（找出真正的附着路径）"""
import numpy as np
from pathlib import Path

T = Path(r"D:\Project\climbing_fall_analysis\temp\thums")
z = np.load(T / "thums_lowerlimb.npz", allow_pickle=True)

# 所有非韧带/非节点的部件 -> 其节点集合
part_nodes = {}
for k in z.files:
    if k in ("node_ids", "node_xyz"):
        continue
    part_nodes[k] = {int(v) for v in z[k].ravel()}

LIG = [(81100701, "PTFL"), (81100801, "CFL 跟腓"), (81100901, "ATFL"),
       (81101001, "PTCL"), (81101101, "下胫腓"), (81101201, "内侧1"),
       (81101202, "内侧2 胫跟"), (81101203, "内侧3 胫舟"),
       (81101204, "内侧4"), (81101401, "跖长")]

for pid, nm in LIG:
    k = f"{pid:08d}_shell"
    mine = part_nodes[k]
    print(f"\n### {pid} {nm}  ({len(mine)} 节点)")
    hits = []
    for pk, pn in part_nodes.items():
        if pk == k or pk.startswith(f"{pid:08d}"):
            continue
        inter = mine & pn
        if inter:
            hits.append((len(inter), pk))
    hits.sort(reverse=True)
    if not hits:
        print("    (与任何部件都不共用节点！)")
    for n, pk in hits[:8]:
        print(f"    共用 {n:3d} 节点 <- {pk}")
