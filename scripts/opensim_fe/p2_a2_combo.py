# -*- coding: utf-8 -*-
"""组合验证：混 BC（股骨→P5 接触，已入账）＋ 注册修正（T6）⇒ A2 到几？

四个口径并列，全部走官方 route_yi / a1_verdict / a2_verdict：
  ① BC2 纯正对照
  ② 混 BC（股骨修，已入账 → 现役 A2 4/5）
  ③ 注册修正（T6 修）
  ④ **混 BC ＋ 注册修正**（股骨 + 胸椎 都修）  ← 真正的候选
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "opensim_fe"))

import s4_alignment_bc2 as s4           # noqa: E402

P5_JSON = ROOT / "temp" / "b1_official" / "fracture_matrix_bc2_p5femur.json"
P5_DIR = ROOT / "temp" / "b1_official" / "bc2_p5femur"
BC2_JSON = ROOT / "results" / "opensim_fe" / "fracture_matrix_bc2.json"
COMBO = ROOT / "temp" / "opensim_fe" / "fm_bc2_p5femur_t6regfix.json"

# --- 造 ④：在 ② 的基础上换掉 T6 的 σ 剖面 ---
combo = json.loads(P5_JSON.read_text(encoding="utf-8"))
d0 = json.loads(BC2_JSON.read_text(encoding="utf-8"))
heights = np.asarray(d0["heights_m"], dtype=float)
lumbar = np.asarray(d0["loads_n"]["lumbar"], dtype=float)
k = 0.010768
sig_T6 = [float(k * x * 0.6104) for x in lumbar]
i_t6 = next(i for i, x in enumerate(combo["parts"]) if x["bone"] == "T6")
for stat in combo["sigma_vm_mpa"]:
    combo["sigma_vm_mpa"][stat]["T6"] = list(sig_T6)
combo["fracture_matrix_9x50"]["median"][i_t6] = [1 if s > 150.0 else 0 for s in sig_T6]
COMBO.write_text(json.dumps(combo, ensure_ascii=False, indent=1), encoding="utf-8")

_oj, _od = s4.BC2_JSON, s4.BC2_DIR


def run(tag, j, dr):
    s4.BC2_JSON, s4.BC2_DIR = j, dr
    try:
        yi = s4.route_yi("median")
        a1 = s4.a1_verdict(yi["phases"])
        a2 = s4.a2_verdict(yi["logistic"])
        beh = {kk: (vv.get("status") if isinstance(vv, dict) else vv)
               for kk, vv in yi["logistic"].items()}
    finally:
        s4.BC2_JSON, s4.BC2_DIR = _oj, _od
    print("\n  【%s】" % tag)
    print("    A2: n_hit=%d/%d  ⇒ **%s**" % (a2["n_hit"], a2["n_target"], a2["verdict"]))
    print("    A1: 边界=%-14s Δ=%-14s ⇒ %s"
          % (a1["ours_boundaries_m"], a1["abs_delta"], a1["verdict"]))
    print("    行为: 股骨=%-20s 胸椎=%s" % (beh.get("股骨"), beh.get("胸椎")))
    miss = [k2 for k2, v2 in a2.items() if k2.endswith("_miss")]
    if miss:
        print("    目标逐个: %s"
              % json.dumps({k2: v2 for k2, v2 in a2.items()
                            if k2 not in ("n_hit", "n_target", "verdict")},
                           ensure_ascii=False)[:300])
    return {"a1": a1, "a2": a2, "behaviour": beh}


print("=" * 92)
print("四口径官方判定对照")
print("=" * 92)
res = {}
res["① BC2 纯正对照"] = run("① BC2 纯正对照（股骨/胸椎 都未修）",
                        BC2_JSON, ROOT / "temp" / "opensim_fe" / "fracture_matrix_bc2")
res["② 混 BC（股骨修，现役）"] = run("② 混 BC：股骨→P5 接触（已入账，现役 A2 4/5）",
                             P5_JSON, P5_DIR)
res["④ 混 BC + 注册修正"] = run("④ **混 BC ＋ T6 注册修正**（股骨 + 胸椎 都修）",
                           COMBO, P5_DIR)

print()
print("=" * 92)
print("结论")
print("=" * 92)
b, m, c = (res["① BC2 纯正对照"], res["② 混 BC（股骨修，现役）"],
           res["④ 混 BC + 注册修正"])
print("  A2 计数：BC2 %d/5(%s) → 混BC %d/5(%s) → **混BC+T6修 %d/5(%s)**"
      % (b["a2"]["n_hit"], b["a2"]["verdict"], m["a2"]["n_hit"], m["a2"]["verdict"],
         c["a2"]["n_hit"], c["a2"]["verdict"]))
print("  ⇒ 现役 A2 = %s（混 BC，%d/5，**恰在阈值上无余量**）"
      % (m["a2"]["verdict"], m["a2"]["n_hit"]))
print("  ⇒ 加 T6 注册修正后 = **%s（%d/5）** ⇒ 余量 %+d"
      % (c["a2"]["verdict"], c["a2"]["n_hit"], c["a2"]["n_hit"] - 4))
print("  A1：BC2 %s / 混BC %s / 组合 %s"
      % (b["a1"]["verdict"], m["a1"]["verdict"], c["a1"]["verdict"]))
(ROOT / "results" / "opensim_fe" / "P2_A2_COMBO_T6REGFIX.json").write_text(
    json.dumps(res, ensure_ascii=False, indent=1, default=str), encoding="utf-8")
print("\n产物: results/opensim_fe/P2_A2_COMBO_T6REGFIX.json")
