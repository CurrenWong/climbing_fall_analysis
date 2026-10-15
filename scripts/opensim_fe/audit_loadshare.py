# -*- coding: utf-8 -*-
"""审计 LOAD_SHARE（中轴骨）：来源核对 + 独立交叉校验。

① 复算 risk_1d_loadshare 的 mass-above 比值，看是否 = 0.6104
② **独立校验**：力学上 F_level(h) = m_above(level)/m_total × GRF(h)（刚体上体），
   所以 F_lumbar(h)/GRF(h) 应 ≈ m_above(L5)/m_total ≈ 0.45–0.50。
   若实测比值远离该值 ⇒ 载荷链有问题。
"""
from __future__ import annotations

import io
import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "opensim_fe"))

import risk_1d_loadshare as r1          # noqa: E402

print("=" * 92)
print("① 复算 mass-above 比值（risk_1d_loadshare 的模型）")
print("=" * 92)
print("  节段质量分数（% BM，Pearsall/Ivancic/Dempster）:")
print("    M_REF(经过 L5/S1) = %.2f %% BM  = %.2f kg" % (r1.M_REF_PCT, r1.M_REF_KG))
print("    ARMS_PCT = %.1f" % r1.ARMS_PCT)
print()
print("  %-12s %6s %8s %10s %10s %12s" %
      ("bone", "level", "arms?", "m_above%", "m_above_kg", "RATIO(m/m_ref)"))
for b in r1.AXIAL_BONES:
    lv, arms = r1.LEVEL_OF[b]
    print("  %-12s %6s %8s %10.2f %10.2f %12.4f"
          % (b, lv, arms, r1.MASS_ABOVE[b], r1.MASS_ABOVE_KG[b], r1.RATIO[b]))

print()
print("  ⚠️ 对照 p2_head_contact.LOAD_SHARE：")
import p2_head_contact as p2            # noqa: E402
for b in r1.AXIAL_BONES:
    ls = p2.LOAD_SHARE[b]
    rr = r1.RATIO[b]
    flag = "✅ 一致" if abs(ls - rr) < 5e-4 else "❌ **不一致（差 %.4f）**" % (ls - rr)
    print("     %-12s LOAD_SHARE=%.4f   mass-above ratio=%.4f   %s"
          % (b, ls, rr, flag))

print()
print("=" * 92)
print("② 独立交叉校验：F_lumbar(h) 应 ≈ m_above(L5)/m_total × GRF(h)")
print("=" * 92)
# 找载荷链缓存
cands = []
for base in (ROOT / "temp", ROOT / "results"):
    for dp, dn, fn in os.walk(base):
        for f in fn:
            if f.endswith(".npz") and ("load" in f.lower() or "chain" in f.lower()):
                cands.append(Path(dp) / f)
print("  候选 npz:", [str(c.relative_to(ROOT)) for c in cands][:10])

fm = json.loads((ROOT / "results" / "opensim_fe" /
                 "fracture_matrix_bc2.json").read_text(encoding="utf-8"))
heights = np.asarray(fm["heights_m"], dtype=float)

for c in cands:
    try:
        z = np.load(c, allow_pickle=True)
    except Exception as e:                       # noqa: BLE001
        print("   [skip] %s: %s" % (c.name, e))
        continue
    if "grf_peak_n" not in z.files or "lumbar__vert" not in z.files:
        continue
    grf = np.asarray(z["grf_peak_n"], dtype=float)
    lum = np.asarray(z["lumbar__vert"], dtype=float)
    h = np.asarray(z["heights"], dtype=float)
    print("\n  ★ 找到载荷链缓存: %s" % c.relative_to(ROOT))
    print("  %-6s %12s %12s %10s" % ("h(m)", "GRF_peak(N)", "F_lumbar(N)", "比值"))
    for hh in (1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 50.0):
        i = int(np.where(h == hh)[0][0])
        if not np.isfinite(grf[i]):
            print("  %-6.0f %12s %12.1f %10s" % (hh, "n/a", lum[i], "—"))
            continue
        print("  %-6.0f %12.1f %12.1f %10.3f" % (hh, grf[i], lum[i], lum[i] / grf[i]))
    ok = np.isfinite(grf)
    print("\n  F_lumbar/GRF：中位 %.4f  范围 [%.4f, %.4f]"
          % (np.nanmedian(lum[ok] / grf[ok]), np.nanmin(lum[ok] / grf[ok]),
             np.nanmax(lum[ok] / grf[ok])))
    print("  ⇒ 力学预期 = m_above(L5)/m_total。THUMS lumbar 子树实测质量 34.24 kg / "
          "总质量 75.337 kg = %.4f" % (34.24 / 75.337))
    print("     而 mass-above 表的 M_REF = %.2f %% BM = %.4f"
          % (r1.M_REF_PCT, r1.M_REF_PCT / 100.0))
    break
