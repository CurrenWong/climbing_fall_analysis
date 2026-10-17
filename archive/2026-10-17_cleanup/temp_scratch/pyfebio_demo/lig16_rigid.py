"""关键诊断：韧带壳的边界节点是否全部有骨连接？（刚体模态风险）
同时做 2-流形检查（LRN-019 同款）。"""
from collections import Counter
import numpy as np
from pathlib import Path

T = Path(r"D:\Project\climbing_fall_analysis\temp\thums")
z = np.load(T / "thums_lowerlimb.npz", allow_pickle=True)
xyz = z["node_xyz"]
i_of = {int(v): i for i, v in enumerate(z["node_ids"])}

BONE_PIDS = [81000700, 81000600, 81000900, 81000800, 81001100, 81001000,
             81001300, 81001200, 81001400, 81001401]
bone_nodes = set()
for pid in BONE_PIDS:
    for kind in ("solid", "shell"):
        k = f"{pid:08d}_{kind}"
        if k in z:
            bone_nodes |= {int(v) for v in z[k].ravel()}

LIG = [(81100701, "PTFL 距腓后"), (81100801, "CFL 跟腓"),
       (81100901, "ATFL 距腓前"), (81101001, "PTCL 距跟后"),
       (81101101, "下胫腓"), (81101201, "内侧1"),
       (81101202, "内侧2 胫跟"), (81101203, "内侧3 胫舟"),
       (81101204, "内侧4"), (81101401, "跖长")]

print(f"{'pid':>9} {'名称':14s} {'单元':>4} {'节点':>4} {'边':>4} {'边界节点':>8} "
      f"{'边界已连':>8} {'悬空边界':>8} {'流形':>6}")
print("-" * 92)
susp_all = []
for pid, nm in LIG:
    e = z[f"{pid:08d}_shell"]
    nodes = {int(v) for v in e.ravel()}
    # 每单元 4 条边 -> 统计边使用次数
    edge = Counter()
    for q in e:
        n = [int(v) for v in q]
        for a, b in ((0, 1), (1, 2), (2, 3), (3, 0)):
            edge[tuple(sorted((n[a], n[b])))] += 1
    bnd = {n for (a, b), c in edge.items() if c == 1 for n in (a, b)}
    bnd_conn = bnd & bone_nodes
    susp = bnd - bone_nodes
    counts = Counter(edge.values())
    manifold = "OK" if set(counts) == {2} else f"非2流形{dict(counts)}"
    susp_all += [(pid, nm, n) for n in susp]
    print(f"{pid:>9} {nm:14s} {len(e):>4} {len(nodes):>4} {len(edge):>4} "
          f"{len(bnd):>8} {len(bnd_conn):>8} {len(susp):>8} {manifold:>6}")

print(f"\n\n=== 悬空边界节点合计 {len(susp_all)} ===")
from collections import defaultdict
by = defaultdict(list)
for pid, nm, n in susp_all:
    by[(pid, nm)].append(n)
for (pid, nm), ns in sorted(by.items()):
    P = xyz[[i_of[n] for n in ns]]
    print(f"  {pid} {nm}: {len(ns)} 个 -> {ns}")
    if len(ns) <= 12:
        for n in ns:
            p = xyz[i_of[n]]
            d = np.min(np.linalg.norm(xyz[[i_of[b] for b in bone_nodes]] - p,
                                      axis=1))
            print(f"      node {n} @({p[0]:7.2f},{p[1]:7.2f},{p[2]:7.2f}) "
                  f"离最近骨节点 {d:.3f} mm")
