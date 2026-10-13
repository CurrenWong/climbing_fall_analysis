# -*- coding: utf-8 -*-
"""A2 漏判诊断：股骨/胸椎为何「常值1」（低高度即骨折 ⇒ 无过渡 ⇒ 不能是阶跃）？

只读冻结件，打印：
  - 股骨（长骨：ends/shaft 双阈值）在参考载荷下的 σ 与阈值之比
  - 胸椎（紧凑骨：单一阈值）同上
  - 二者 σ(h) 的最低值（若 h=1 就超阈 ⇒ 常值1）
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "opensim_fe"))
import paper_alignment as pa  # noqa: E402

FM = json.loads((ROOT / "results" / "opensim_fe" / "fracture_matrix_bc2.json").read_text(encoding="utf-8"))
YI = json.loads((ROOT / "results" / "opensim_fe" / "s4_alignment_bc2.json").read_text(encoding="utf-8"))
det = YI["detail_yi"]

print("=== 强度阈值（来自 REGIONS 的 _sc） ===")
for cn, bone, joint, subs in pa.REGIONS:
    if cn in ("股骨", "胸椎", "腓骨", "骨盆", "颅骨"):
        for sub in subs:
            print(f"  {cn:<4} {bone:<14} 子区域 {sub[0]:<10} σ_c = {pa._sc(sub[1]):>7.2f} MPa")
print()

for cn in ("股骨", "胸椎"):
    d = det[cn]
    print(f"=== {cn} ===")
    for k, v in d.items():
        if isinstance(v, (int, float)):
            print(f"  {k} = {v}")
    print()

print("=== 依据矩阵：h=1 是否已骨折（常值1 的判据） ===")
heights = FM["heights_m"]
for cn in ("足部", "胫骨", "腓骨", "股骨", "腰椎", "胸椎", "颈椎", "颅骨"):
    d = det[cn]
    m = d.get("matrix_by_part") or d.get("mask")
    if m is None:
        print(f"  {cn:<4} (无 mask 字段, keys={list(d)[:6]})")
        continue
    arr = np.asarray(m, int)
    print(f"  {cn:<4} h=1:{arr[0]}  h=5:{arr[min(4,len(arr)-1)]}  "
          f"首次骨折 h = {heights[int(np.argmax(arr))] if arr.any() else 'NEVER'}")
