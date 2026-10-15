# -*- coding: utf-8 -*-
"""足部（跟骨）的首折高度对**统计量口径**的敏感性（零 FE）。

若 `max` 给很早、而 `median` 给很晚 ⇒ 首折由**少数峰值单元**决定 = 载荷带伪影签名
⇒ 修掉伪影即可把首折推到目标窗。
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
fm = json.loads((ROOT / "results" / "opensim_fe" /
                 "fracture_matrix_bc2.json").read_text(encoding="utf-8"))
heights = np.asarray(fm["heights_m"], dtype=float)
parts = {p["bone"]: p for p in fm["parts"]}

print("=" * 92)
print("足部 calcaneus_r：各统计量下的首折高度（σ_c = 150 MPa）")
print("=" * 92)
sc = float(parts["calcaneus_r"].get("sigma_c_mpa") or 150.0)
print("  %-8s %10s %10s %12s %12s" % ("stat", "首折h", "σ@首折", "σ(13m)", "σ(17m)"))
for stat, d in fm["sigma_vm_mpa"].items():
    sig = np.asarray(d.get("calcaneus_r", []), dtype=float)
    if sig.size == 0:
        continue
    i1 = next((i for i, s in enumerate(sig) if s > sc), None)
    i13 = int(np.where(heights == 13)[0][0])
    i17 = int(np.where(heights == 17)[0][0])
    print("  %-8s %10s %10s %12.2f %12.2f"
          % (stat, ("%.0f m" % heights[i1]) if i1 is not None else ">50",
             ("%.1f" % sig[i1]) if i1 is not None else "—", sig[i13], sig[i17]))

print()
print("  注：矩阵行用的是哪个统计量？看现役骨折行与各统计量的首折是否吻合。")
i_cal = next(i for i, p in enumerate(fm["parts"]) if p["bone"] == "calcaneus_r")
row = fm["fracture_matrix_9x50"]["median"][i_cal]
fh_row = next((heights[i] for i, v in enumerate(row) if v), None)
print("     现役骨折行（median 矩阵）首折 = %s" % fh_row)

print()
print("=" * 92)
print("各部位：max 口径 vs median 口径的首折对比（伪影放大倍数）")
print("=" * 92)
print("  %-14s %10s %10s %10s" % ("bone", "max首折", "median首折", "max/med"))
for bone in ("calcaneus_r", "tibia_r", "femur_r", "fibula_r", "T6", "L3",
             "C5", "parietal_r"):
    rowp = parts.get(bone)
    if rowp is None:
        continue
    s = float(rowp.get("sigma_c_mpa") or 0)
    out = []
    for stat in ("max", "median"):
        sig = np.asarray(fm["sigma_vm_mpa"][stat].get(bone, []), dtype=float)
        i1 = next((i for i, v in enumerate(sig) if v > s), None) if sig.size else None
        out.append(heights[i1] if i1 is not None else 50.0)
    print("  %-14s %10s %10s %10s" % (bone, "%.0f m" % out[0], "%.0f m" % out[1],
                                      ("%.1f×" % (out[0] / out[1])) if out[1] else "—"))
