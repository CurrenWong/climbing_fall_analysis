"""只读：探查 thums_lowerlimb.npz 结构，为均匀 hex8 细分做准备。"""
from pathlib import Path
import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
d = np.load(ROOT / "temp/thums/thums_lowerlimb.npz", allow_pickle=True)
print("=== keys ===")
for k in d.files:
    a = d[k]
    print(f"  {k:28s} shape={str(a.shape):>16s} dtype={a.dtype}")

print("\n=== 单元连接度统计（按 part）===")
ANKLE = ("81000700_solid", "81000600_solid", "81000601_solid",
         "81000900_solid", "81000800_solid", "81000801_solid",
         "81001100_solid", "81001000_solid",
         "81001300_solid", "81001200_solid",
         "81001400_solid", "81001500_solid")
tot_e = tot_uniq = 0
for k in d.files:
    if not k.endswith(("_solid", "_shell")):
        continue
    a = d[k]
    if a.ndim != 2:
        continue
    cols = a.shape[1]
    # 每行去重后的唯一节点数
    uq = np.array([len(set(int(x) for x in row if x > 0)) for row in a])
    vals, cnts = np.unique(uq, return_counts=True)
    tag = "★踝" if k in ANKLE else ""
    print(f"  {k:24s} n={a.shape[0]:>6,} cols={cols}  "
          f"unique节点数分布 {dict(zip(vals.tolist(), cnts.tolist()))} {tag}")
    if k in ANKLE:
        tot_e += a.shape[0]
print(f"\n踝部单元合计 {tot_e:,}")

# 全局节点使用情况：踝部 part 用到的节点范围
print("\n=== 全局 ===")
print(f"  node_ids {d['node_ids'].shape}, node_xyz {d['node_xyz'].shape}")
print(f"  node_ids min/max = {d['node_ids'].min()}/{d['node_ids'].max()}")
# 检查是否存在退化 hex（唯一节点 <8）
