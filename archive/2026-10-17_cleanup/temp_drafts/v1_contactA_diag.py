"""诊断：软接触档高高度处的 `F_flex` 非单调是什么机制（bottomed_out？行程饱和？）。"""
import sys
from pathlib import Path

ROOT = Path(r'D:\Project\climbing_fall_analysis')
sys.path.insert(0, str(ROOT / "src"))
import climbing.pad as P  # noqa: E402

post = P.POSTURES["anchored-self"]
print(f"anchored-self: x_flex_max={post.x_flex_max} m, collapse_ref_m={post.collapse_ref_m}, "
      f"k_flex={post.k_flex} N/m")
print(f"k_stop 止挡 = 3.0e5 N/m（pad.py 内固定）\n")

orig = P.RIGID_CONTACT["k_n_per_m"]
try:
    for kn in (2.1e6, 1.6e6):
        P.RIGID_CONTACT["k_n_per_m"] = kn
        print(f"=== k_n = {kn:.1e} ===")
        print(f"  {'h(m)':>6} {'F_flex(kN)':>11} {'F_pad(kN)':>11} {'max_comp(mm)':>13} "
              f"{'k_eff(kN/m)':>12} {'bottomed':>9} {'absorbed%':>10} {'resid(J)':>10}")
        for h in (0.5, 1.0, 2.0, 5.0, 9.0):
            r = P.simulate_boulder_fall(height_m=h, mass_kg=80.0, on_pad=False,
                                        posture="anchored-self")
            comp = float(r.max_compression_m)
            print(f"  {h:>6.1f} {r.meta['f_flex_peak_n']/1e3:>11.2f} "
                  f"{r.peak_force_n/1e3:>11.2f} {comp*1e3:>13.2f} "
                  f"{r.peak_force_n/max(comp,1e-12)/1e3:>12.1f} "
                  f"{str(bool(r.bottomed_out)):>9} "
                  f"{float(getattr(r,'absorbed_fraction',float('nan')))*100:>9.1f}% "
                  f"{float(getattr(r,'energy_residual_j',float('nan'))):>10.1f}")
        print()
finally:
    P.RIGID_CONTACT["k_n_per_m"] = orig
