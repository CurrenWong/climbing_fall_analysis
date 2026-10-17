"""CARTFIX 诊断 1：cart_m1 网格结构 + 抬回机制数据验证。

回答：
  A. cart_m1/cartilage_mesh.npz 有哪些键？是 patch 产物还是 cartilage.py 源？
  B. 顶面位置 = 底面 + 逐节点 disp(中位数) 的重心插值？逐单元验证。
  C. 每个楔形的：bottom/top 三角法向、扭角、面积、厚度、逐顶点 disp。
"""
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\Project\climbing_fall_analysis")
D = ROOT / "temp/thums/cart_m1"
d = np.load(D / "cartilage_mesh.npz", allow_pickle=True)
print("=== npz 键 ===")
for k in d.files:
    a = d[k]
    print(f"  {k:34s} shape={str(a.shape):16s} dtype={a.dtype}")

JOINTS = ["tibiotalar_tibia", "tibiotalar_talus",
          "subtalar_talus", "subtalar_calcaneus"]
for j in JOINTS:
    if f"{j}_elems" not in d.files:
        continue
    E = np.asarray(d[f"{j}_elems"])
    N = np.asarray(d[f"{j}_nodes"], float) if f"{j}_nodes" in d.files else None
    print(f"\n=== {j} ===")
    print(f"  elems shape={E.shape}  uniq counts per row: "
          f"{sorted(set((np.asarray(E)!=0).sum(axis=1).tolist()))}")
    if N is not None:
        print(f"  nodes={N.shape}  is_offset true={(np.asarray(d[f'{j}_is_offset'],bool)).sum()}")
    print(f"  另外的键: {[k for k in d.files if k.startswith(j)]}")
