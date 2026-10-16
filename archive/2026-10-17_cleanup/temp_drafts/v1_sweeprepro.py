"""复现 `V1_WHY_20KN_AND_GAP_LOCUS.md` §2 的 k 扫描表（h=1 m 刚性地面）。
用来定「正对照」的真实口径：基座是 `feet-first-stiff` 还是 `controlled-drop`？
期望 F_flex 列：2e4⇒4.26 / 4e4⇒5.50 / 8e4⇒8.09 / 2.2e5⇒15.95 / 4e5⇒20.05 / 1.1e6⇒27.43 kN
"""
import sys
from dataclasses import replace
from pathlib import Path

ROOT = Path(r'D:\Project\climbing_fall_analysis')
sys.path.insert(0, str(ROOT / "src"))
import climbing.pad as P  # noqa: E402

BW = 80.0 * 9.80665
KS = [2.0e4, 4.0e4, 8.0e4, 2.2e5, 4.0e5, 1.1e6]
WANT = {2.0e4: 4.26, 4.0e4: 5.50, 8.0e4: 8.09, 2.2e5: 15.95, 4.0e5: 20.05, 1.1e6: 27.43}

for base_key in ("feet-first-stiff", "controlled-drop"):
    base = P.POSTURES[base_key]
    print(f"\n=== 基座 = {base_key} "
          f"(m_low_frac={base.m_low_frac}, x_flex_max={base.x_flex_max}, "
          f"collapse_ref_m={base.collapse_ref_m}) ===")
    print(f"  {'k':>9} {'F_flex(kN)':>11} {'kBW':>7} {'want':>7} {'Δ%':>8} "
          f"{'F_pad(kN)':>10} {'行程(mm)':>9}")
    key = f"__sweep_{base_key}"
    for k in KS:
        P.POSTURES[key] = replace(base, name="sweep", k_flex=k)
        r = P.simulate_boulder_fall(height_m=1.0, mass_kg=80.0, on_pad=False, posture=key)
        ff = float(r.meta["f_flex_peak_n"])
        tol_ = abs(ff / 1e3 - WANT[k]) / WANT[k] * 100
        # 行程 = F/k
        print(f"  {k:>9.1e} {ff/1e3:>11.2f} {ff/BW:>7.2f} {WANT[k]:>7.2f} {tol_:>7.1f}% "
              f"{r.peak_force_n/1e3:>10.2f} {ff/k*1e3:>9.1f}")
    del P.POSTURES[key]
