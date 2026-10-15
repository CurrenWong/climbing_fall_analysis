# -*- coding: utf-8 -*-
"""A1 三个锚定部位的 σ 余量（零 FE）—— 要推动首折需多少倍？

目标：把 胫骨 h* 5→6+（#1），足部 h* 13→16-18（#2）。
含义：σ(h) 与阈值 σ_c 的比值，以及在目标高度处还差多少倍。
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "opensim_fe"))

fm = json.loads((ROOT / "results" / "opensim_fe" /
                 "fracture_matrix_bc2.json").read_text(encoding="utf-8"))
heights = np.asarray(fm["heights_m"], dtype=float)
parts = {p["bone"]: p for p in fm["parts"]}
sigs = fm["sigma_vm_mpa"]["median"]

TARGET = {"calcaneus_r": ("足部", 17.0), "tibia_r": ("胫骨", 7.0),
          "parietal_r": ("颅骨", 7.0)}

print("=" * 94)
print("A1 锚定部位的 σ 余量（现役 BC2 · median · mm-N-MPa）")
print("=" * 94)
print("  %-12s %-6s %8s %10s %12s %12s %14s"
      % ("bone", "部位", "σ_c", "首折h", "σ@首折", "σ/σ_c", "目标h 需降"))
for bone, (cn, h_tgt) in TARGET.items():
    row = parts.get(bone)
    if row is None:
        print("  %-12s 缺" % bone)
        continue
    sig = np.asarray(sigs.get(bone, []), dtype=float)
    if sig.size == 0:
        print("  %-12s 无 σ" % bone)
        continue
    sc = float(row.get("sigma_c_mpa") or row.get("sigma_c") or 0)
    i_first = next((i for i, s in enumerate(sig) if s > sc), None)
    h_first = heights[i_first] if i_first is not None else None
    s_first = sig[i_first] if i_first is not None else None
    i_t = int(np.where(heights == h_tgt)[0][0])
    s_t = sig[i_t]
    need = s_t / sc if sc else float("nan")
    print("  %-12s %-6s %8.1f %10s %12s %12s %14s"
          % (bone, cn, sc, ("%.0f m" % h_first) if h_first else ">50",
             ("%.2f" % s_first) if s_first else "—",
             ("%.3f" % (s_first / sc)) if s_first else "—",
             "σ(%.0fm)=%.1f ⇒ 需降 %.2f×" % (h_tgt, s_t, need)))

print()
print("  （「需降 N×」= 该部位在目标高度处的 σ / σ_c；<1 表示该高度本不应折，")
print("    >1 表示要推迟到该高度，σ 必须整体降 N 倍。）")
print()
print("=" * 94)
print("逐高度 σ 剖面（\u00d7σ_c）")
print("=" * 94)
hs = [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 13, 14, 15, 16, 17, 18, 20]
print("  %-6s %10s %10s %10s" % ("h(m)", "胫骨", "足部", "颅骨"))
for h in hs:
    i = int(np.where(heights == h)[0][0])
    cells = []
    for bone in ("tibia_r", "calcaneus_r", "parietal_r"):
        row = parts[bone]
        sc = float(row.get("sigma_c_mpa") or 0)
        v = np.asarray(sigs[bone], dtype=float)[i]
        cells.append("%10.3f" % (v / sc if sc else float("nan")))
    flag = ""
    if h in (6, 7, 16, 17, 18):
        flag = "  ← 目标窗"
    print("  %-6.0f %s%s" % (h, " ".join(cells), flag))
