"""只读探针：摸清 cartilage.json 的结构 + 各关节面/节点规模。"""
import json
from pathlib import Path

import numpy as np

P = Path(r"D:\Project\climbing_fall_analysis\temp\thums\cartilage.json")
d = json.loads(P.read_text(encoding="utf-8"))

print("顶层 key:", list(d.keys()))
for k, v in d.items():
    if isinstance(v, dict):
        print(f"  {k}: dict, 子 key = {list(v.keys())[:6]}")
    elif isinstance(v, list):
        print(f"  {k}: list len={len(v)}  样例={str(v[0])[:110] if v else '-'}")
    else:
        print(f"  {k}: {v!r}")

# 逐关节
joints = d.get("cart") or d.get("joints") or {}
print("\n--- 逐关节 ---")
for j, rec in (joints.items() if isinstance(joints, dict) else []):
    print(f"\n[{j}] keys={list(rec.keys())}")
    for k, v in rec.items():
        if isinstance(v, list):
            print(f"   {k}: n={len(v)}  样例={str(v[0])[:130] if v else '-'}")
        else:
            print(f"   {k}: {str(v)[:100]}")

for extra in ("offset_pts", "OFFSET_PTS", "offsets"):
    if extra in d:
        o = d[extra]
        print(f"\n{extra}: type={type(o).__name__}", end=" ")
        if isinstance(o, dict):
            print("keys=", list(o.keys()))
            for k, v in o.items():
                print(f"   {k}: {np.asarray(v).shape}")
        else:
            print("len=", len(o))
