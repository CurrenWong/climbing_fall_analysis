# -*- coding: utf-8 -*-
"""修好的 BC 下，胸椎 T6 的首折高度 = ?（零 FE，用 σ∝F 线性外推）

依据
----
V1 实测（注册修正后，能力门 3/4 PASS）：
  h=10  F=19675.9 N  σ_median = 209.6 MPa
  h=30  F=38257.3 N  σ_median = 416.3 MPa
⇒ σ/F = 2.096/196.759 = 0.010652  vs  4.163/382.573 = 0.010881
   两者差 **2.1%** ⇒ **线性**（严格线弹性，本仓库既有口径：σ(h)=σ_ref·load(h)/load_ref）
⇒ k = 0.010766 MPa/N（均值）

然后：σ(h) = k × load_T6(h)，load_T6(h) = lumbar(h) × LOAD_SHARE["T6"]
阈值 σ_c = 150 MPa（[Y25] 材料强度，皮质/压缩）
求首个 σ>150 的高度。

再把这个 h* 喂给官方的 A1 判定，看 #2 边界与 Δ₂。
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "opensim_fe"))

import paper_alignment as pa            # noqa: E402
import s4_alignment_bc2 as s4           # noqa: E402

SIGMA_C = 150.0
SHARE_T6 = 0.6104

d = json.loads((ROOT / "results" / "opensim_fe" /
                "fracture_matrix_bc2.json").read_text(encoding="utf-8"))
heights = np.asarray(d["heights_m"], dtype=float)
lumbar = np.asarray(d["loads_n"]["lumbar"], dtype=float)
load_t6 = lumbar * SHARE_T6

# σ ∝ F 的标定（V1 实测两点）
k1 = 209.626 / 19675.9
k2 = 416.309 / 38257.3
k = 0.5 * (k1 + k2)
print("=" * 88)
print("① σ ∝ F 线性标定（V1 注册修正后的实测两点）")
print("=" * 88)
print("   k1 = 209.626/19675.9  = %.6f MPa/N" % k1)
print("   k2 = 416.309/38257.3  = %.6f MPa/N" % k2)
print("   两者偏差 %.2f%%  ⇒ %s" % (100 * abs(k2 / k1 - 1),
                                     "线性成立" if abs(k2 / k1 - 1) < 0.05 else "非线性"))
print("   采用 k = %.6f MPa/N" % k)

print()
print("=" * 88)
print("② 胸椎 T6：σ(h) 与首折高度（σ_c = %.0f MPa）" % SIGMA_C)
print("=" * 88)
print("   %-5s %12s %11s %10s %s" % ("h(m)", "lumbar(N)", "load_T6(N)", "σ_med(MPa)", ">σ_c?"))
first = None
for i, h in enumerate(heights):
    F = load_t6[i]
    sig = k * F
    hit = sig > SIGMA_C
    if hit and first is None:
        first = h
    if h <= 12 or h in (15.0, 16.0, 18.0, 20.0, 30.0, 50.0) or (hit and h <= 8):
        print("   %-5.0f %12.1f %11.1f %10.1f %s"
              % (h, lumbar[i], F, sig, "★ 折" if hit else ""))

if first is None:
    print("\n   ⇒ 50 m 内不折（常值0）")
else:
    print("\n   ⇒ **胸椎首折 h* = %.0f m**  ⇒ 行为 = **%s**"
          % (first, "阶跃" if first > 1 else "常值1"))

# 反推：要落在 A1 窗口 [16,18] 需要多少载荷
need_F = SIGMA_C / k
h_need = np.interp(need_F, load_t6, heights)
print("   不折条件：σ<%.0f ⇒ load_T6 < %.0f N ⇒ lumbar < %.0f N ⇒ h < **%.1f m**"
      % (SIGMA_C, need_F, need_F / SHARE_T6, h_need))
print("   ⇒ 现模型在 h=%.0f m 就折；A1 要求 16 m 以下不折 ⇒ 载荷需再降 **%.1f×**"
      % (first or 0, load_t6[15] / need_F))

# ------------------------------------------------------------------ A1
print()
print("=" * 88)
print("③ 把首折 h*=%s 喂给官方 A1 判定（只改胸椎行，其余取 BC2 正对照）" % first)
print("=" * 88)
yi = s4.route_yi("median")
names = pa.REGION_NAMES
rows = {cn: np.asarray(yi["matrix_by_part"][cn], dtype=int) for cn in names}


def step_at(h):
    return np.array([1 if (i + 1) >= h else 0 for i in range(50)], dtype=int)


def a1_of(over):
    r2 = {cn: over.get(cn, rows[cn]) for cn in names}
    M = np.array([[int(r2[cn][i]) for cn in names] for i in range(50)], dtype=int)
    clus = pa.jaccard_clusters(M, 5)
    ph = s4.three_phases_seq(clus["labels"], M)
    return s4.a1_verdict(ph), ph


for tag, over in (("BC2 正对照（胸椎 常值1）", {}),
                  ("修好 BC：胸椎 阶跃@%d" % first if first else "—",
                   {"胸椎": step_at(first)} if first else {})):
    v, ph = a1_of(over)
    print("  %-34s 边界=%-14s Δ=%-14s -> %s"
          % (tag, v["ours_boundaries_m"], v["abs_delta"], v["verdict"]))
    print("      相段: %s" % [p["range_m"] for p in ph])

out = {"k_mpa_per_n": k, "k1": k1, "k2": k2,
       "sigma_c_mpa": SIGMA_C, "share_t6": SHARE_T6,
       "first_fracture_m": first,
       "behaviour": "阶跃" if (first and first > 1) else ("常值1" if first == 1 else "常值0"),
       "load_t6_need_n": need_F, "lumbar_need_n": need_F / SHARE_T6,
       "h_need_m": float(h_need),
       "load_reduction_needed_x": float(load_t6[15] / need_F),
       "sigma_at_h10": float(k * load_t6[9]),
       "sigma_at_h16": float(k * load_t6[15])}
(ROOT / "results" / "opensim_fe" / "P2_T6_THRESHOLD_FIXEDBC.json").write_text(
    json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
print("\n产物: results/opensim_fe/P2_T6_THRESHOLD_FIXEDBC.json")
