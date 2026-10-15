# -*- coding: utf-8 -*-
"""足部（跟骨）的 gauge 口径逐高度首折 —— 用它能否把 A1 的 #2 推到 16-18 m？

数据源：results/opensim_fe/s2_thums_height_sweep_anat.json（**方向 C 解剖跟骨，
gauge R=4mm，已存在的扫描结果**）。零 FE。
对照：9 骨矩阵 median 口径首折 = 13 m。
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "results" / "opensim_fe" / "s2_thums_height_sweep_anat.json"
THRESH = 150.0

d = json.loads(SRC.read_text(encoding="utf-8"))
print("=" * 92)
print("s2_thums_height_sweep_anat.json 结构")
print("=" * 92)
if isinstance(d, list):
    print("  顶层是 list, n =", len(d))
    rows = d
else:
    print("  top keys:", list(d.keys())[:12])
    rows = None
    for k in ("rows", "sweep", "results", "cases"):
        if k in d and isinstance(d[k], list):
            rows = d[k]
            print("  rows 键 =", k, " n =", len(rows))
            break
if not rows:
    print(json.dumps(d, ensure_ascii=False, indent=1)[:2500])
    raise SystemExit(0)
print("  首行字段:", list(rows[0].keys()))

print()
print("=" * 92)
print("gauge 口径逐高度（阈值 %.0f MPa）" % THRESH)
print("=" * 92)
key_h = next((k for k in rows[0] if k.lower() in ("h_m", "h", "height_m", "height")), None)
key_g = next((k for k in rows[0] if "gauge" in k.lower() and "max" in k.lower()), None)
key_r = next((k for k in rows[0] if k.lower() in ("max", "raw_max", "sigma_max")), None)
key_a = next((k for k in rows[0] if k.lower() in ("act", "a", "activity")), None)
print("  用列：h='%s'  gauge='%s'  raw='%s'  act='%s'" % (key_h, key_g, key_r, key_a))
print("  %-6s %-6s %12s %10s %10s %s" % ("h(m)", "act", "gauge_max", "×阈", "raw_max", ">阈?"))
first = None
for r in rows:
    if key_a and float(r.get(key_a, 0.0)) != 0.0:
        continue
    h = float(r[key_h])
    g = float(r[key_g])
    if first is None and g > THRESH:
        first = h
    if h <= 20 or h in (25.0, 33.0, 50.0):
        print("  %-6.0f %-6s %12.2f %10.3f %10.2f %s"
              % (h, r.get(key_a, ""), g, g / THRESH,
                 float(r[key_r]) if key_r else float("nan"),
                 "★ 折" if g > THRESH else ""))
print("\n  ⇒ 方向 C 解剖跟骨 gauge 口径首折 h* = %s" % (first if first else ">扫描上界"))
print("  ⇒ 对照 9 骨矩阵 median 口径首折 = 13 m；A1 的 #2 需要 16–18 m")
