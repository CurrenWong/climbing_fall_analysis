"""诊断：① V1 的收敛数据 ② tibia_cort hex8 是否含「坐标重合节点」导致退化"""
import subprocess
from pathlib import Path
import numpy as np

P = Path(r"D:\Project\climbing_fall_analysis")
FEBIO = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")
W = P / "temp/thums/febio/bisect"

r = subprocess.run([str(FEBIO), "-i", str(W / "v1_hex8_only.feb")], cwd=str(W),
                   capture_output=True, text=True, timeout=900)
out = (r.stdout or "") + (r.stderr or "")
L = out.split("\n")
print("=== V1 的前两组收敛数据 ===")
for i, l in enumerate(L):
    if "convergence norms" in l:
        for j in range(i, min(i + 6, len(L))):
            print("  " + L[j].rstrip())
        print("  ---")
        break
print("\n=== V1 是否报告负 Jacobian ===")
print("  ", [l.strip() for l in L if "jacobian" in l.lower()][:5] or "无")

# ② 单元内最小节点间距
db = np.load(P / "temp/thums/thums_lowerlimb.npz")
nid, xyz = db["node_ids"], db["node_xyz"]
idx = {int(v): i for i, v in enumerate(nid)}
print("\n=== 单元内最小节点间距（判断退化/重合节点）===")
for key, label in (("81000700_solid", "tibia_cort"), ("81000600_solid", "tibia_end_spon"),
                   ("81001100_solid", "talus_cort"), ("81001000_solid", "talus_spon")):
    mins = []
    n_coincide = 0
    for row in db[key]:
        conn = []
        for x in row:
            xi = int(x)
            if xi not in conn:
                conn.append(xi)
        pts = xyz[[idx[v] for v in conn]]
        if len(pts) < 2:
            continue
        d = np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=-1)
        np.fill_diagonal(d, np.inf)
        m = d.min()
        mins.append(m)
        if m < 1e-6:
            n_coincide += 1
    mins = np.array(mins)
    print(f"  {label:16s} n={len(mins):>6,}  最小间距 min={mins.min():.4f} "
          f"p1={np.percentile(mins,1):.4f} 中位={np.median(mins):.4f}  "
          f"重合(<1e-6) {n_coincide}")

# ③ 全局是否有重合节点（不同 id 同坐标）
print("\n=== 全局重合节点检查（不同 id 同坐标）===")
from scipy.spatial import cKDTree
t = cKDTree(xyz)
pairs = t.query_pairs(1e-6, output_type="ndarray")
print(f"  {len(xyz):,} 个节点中，坐标重合的节点对: {len(pairs):,}")

# ④ 单元是否「非流形」：每个节点被几个单元共用
print("\n=== 节点连通性（应无孤立/自由漂浮分量）===")
import collections
adj = collections.defaultdict(set)
for key in ("81000700_solid", "81000600_solid", "81000601_solid"):
    for row in db[key]:
        conn = []
        for x in row:
            xi = int(x)
            if xi not in conn:
                conn.append(xi)
        for i, a in enumerate(conn):
            for b in conn[i + 1:]:
                adj[a].add(b)
                adj[b].add(a)
# 连通分量
seen, comps = set(), []
for s in adj:
    if s in seen:
        continue
    stack, comp = [s], 0
    seen.add(s)
    while stack:
        u = stack.pop()
        comp += 1
        for v in adj[u]:
            if v not in seen:
                seen.add(v)
                stack.append(v)
    comps.append(comp)
comps.sort(reverse=True)
print(f"  胫骨(全部3部件) 节点 {len(adj):,}，连通分量 {len(comps)} 个，"
      f"最大 {comps[0]:,}，其余 {comps[1:6]}{'...' if len(comps)>6 else ''}")
