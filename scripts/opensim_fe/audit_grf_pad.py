# -*- coding: utf-8 -*-
"""决定性实测：on_pad True/False 的峰值 GRF 对比（论文场景 = 软垫）。"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from climbing.coupling import opensim_grf as og      # noqa: E402
from climbing import pad as padmod                   # noqa: E402

MASS = 75.337
BW = MASS * 9.81
print("=" * 96)
print("峰值 GRF：硬地面（on_pad=False，现役）vs 软垫（on_pad=True，论文场景）")
print("质量 %.3f kg  BW = %.1f N  ｜ 默认姿势 = %s" % (MASS, BW, og.DEFAULT_POSTURE))
print("可用姿势 = %s" % sorted(padmod.POSTURES))
print("=" * 96)
print("  %-6s %14s %9s %14s %9s %8s" %
      ("h(m)", "硬地面(N)", "×BW", "软垫(N)", "×BW", "降幅"))
rows = []
for h in (1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 50.0):
    got = {}
    for tag, on_pad in (("rigid", False), ("pad", True)):
        try:
            r = og.ground_reaction(height_m=h, mass_kg=MASS, on_pad=on_pad)
            f = np.asarray(r.f_total_n if hasattr(r, "f_total_n")
                           else r.meta.get("f_total_n"), dtype=float)
            got[tag] = float(np.nanmax(np.abs(f))) if len(f) else float("nan")
        except Exception as e:                        # noqa: BLE001
            got[tag] = float("nan")
            print("   [%s h=%.0f] %s: %s" % (tag, h, type(e).__name__, str(e)[:70]))
    a, b = got.get("rigid", np.nan), got.get("pad", np.nan)
    rows.append((h, a, b))
    print("  %-6.0f %14.1f %9.1f %14.1f %9.1f %8s"
          % (h, a, a / BW, b, b / BW,
             ("%.2f×" % (a / b)) if (np.isfinite(a) and np.isfinite(b) and b) else "—"))

print()
print("文献锚点（同行评审落体实验，屈腿柔性落地）：")
print("   0.3 m/65 kg -> 1910 N ≈ 3.0 ×BW ｜ 0.5 m/65 kg -> 2655 N ≈ 4.2 ×BW")
print("   0.5 m/70 kg -> 2710 N ≈ 3.9 ×BW")
print()
print("读法：若软垫档的 ×BW 落到文献量级（几倍 BW），而硬地面档在几十倍 BW，")
print("      则「现役载荷链用错了场景」成立，且这正是 T6 σ 过高与 A1 窗口不可达的源头。")
