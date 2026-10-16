"""能量账红旗定位：`absorbed_fraction > 100%` 的三项拆解。只读 + 临时门。"""
import sys
from pathlib import Path

ROOT = Path(r'D:\Project\climbing_fall_analysis')
sys.path.insert(0, str(ROOT / "src"))
import climbing.pad as P  # noqa: E402

o = P.RIGID_CONTACT["k_n_per_m"]
print(f"{'k_n':>9} {'h(m)':>5} {'mgh_J':>9} {'impact_J':>10} {'into_pad':>10} "
      f"{'into_flex':>10} {'sum':>10} {'sum/impact':>11} {'resid_J':>9} {'bottomed':>9}")
try:
    for kn in (None, 2.1e6, 1.6e6):
        P.RIGID_CONTACT["k_n_per_m"] = kn
        for h in (0.5, 1.0, 2.0, 5.0):
            r = P.simulate_boulder_fall(height_m=h, mass_kg=80.0, on_pad=False,
                                        posture="anchored-self")
            s = r.energy_into_pad_j + r.energy_into_flex_j
            tag = "None" if kn is None else f"{kn:.1e}"
            print(f"{tag:>9} {h:>5.1f} {80*9.80665*h:>9.1f} {r.impact_energy_j:>10.1f} "
                  f"{r.energy_into_pad_j:>10.1f} {r.energy_into_flex_j:>10.1f} "
                  f"{s:>10.1f} {s/r.impact_energy_j:>10.3f} {r.energy_residual_j:>9.1f} "
                  f"{str(bool(r.bottomed_out)):>9}")
        print()
finally:
    P.RIGID_CONTACT["k_n_per_m"] = o
    print(f"(已恢复 {o!r})")
