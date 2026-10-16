"""正对照排查：`feet-first-stiff` @1 m 的 F_flex 期望 27.43 kN（`V1_WHY_20KN_AND_GAP_LOCUS.md:46`）。
试不同求解器容差，找出口径。"""
import inspect
import sys
from pathlib import Path

ROOT = Path(r'D:\Project\climbing_fall_analysis')
sys.path.insert(0, str(ROOT / "src"))
import climbing.pad as P  # noqa: E402

print("simulate_boulder_fall 签名：")
print(" ", inspect.signature(P.simulate_boulder_fall))
print()

CASES = [
    ("defaults", {}),
    ("tight 1e-4/1e-8/1e-11", dict(max_ode_step=1e-4, rtol=1e-8, atol=1e-11)),
    ("tight, no max_step", dict(rtol=1e-8, atol=1e-11)),
]

print(f"{'口径':<26} {'F_pad(kN)':>10} {'F_flex(kN)':>11} {'F_pad/BW':>9} {'F_flex/BW':>10} "
      f"{'bottomed':>9}")
BW = 80.0 * 9.80665
for lbl, kw in CASES:
    try:
        r = P.simulate_boulder_fall(height_m=1.0, mass_kg=80.0, on_pad=False,
                                   posture="feet-first-stiff", **kw)
        fp, ff = float(r.peak_force_n), float(r.meta["f_flex_peak_n"])
        print(f"{lbl:<26} {fp/1e3:>10.2f} {ff/1e3:>11.2f} {fp/BW:>9.2f} {ff/BW:>10.2f} "
              f"{str(bool(r.bottomed_out)):>9}")
    except TypeError as e:
        print(f"{lbl:<26} TypeError: {e}")
print()
print("期望（同文件 :46）：F_flex = 27.43 kN = 35.0 BW（k=1.1e6 档）")
print()
print("F_pad 与 k 无关的旁证（三档锚定 @0.42 m，defaults）：")
for nm in ("anchored-soft", "anchored-self", "anchored-stiff"):
    r = P.simulate_boulder_fall(height_m=0.42, mass_kg=80.0, on_pad=False, posture=nm)
    print(f"  {nm:<16} k={P.POSTURES[nm].k_flex/1e3:>6.1f} kN/m  "
          f"F_pad={r.peak_force_n/1e3:>6.2f} kN  F_flex={r.meta['f_flex_peak_n']/1e3:>6.2f} kN")
