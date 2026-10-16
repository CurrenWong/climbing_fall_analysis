"""★ 关键自检：A 的判定档（@0.42 m）是否也在能量越界区？若是 ⇒ A 判定作废。
同时对照 LEG_K_ANCHOR（门关）与各档 3 个姿势。只读 + 临时门。"""
import sys
from pathlib import Path

ROOT = Path(r'D:\Project\climbing_fall_analysis')
sys.path.insert(0, str(ROOT / "src"))
import climbing.pad as P  # noqa: E402

LAB = ("soft", "self", "stiff")
SET = [("None(LEG_K_ANCHOR 口径)", None), ("3.0e6", 3.0e6), ("2.1e6(A-J1 判定档)", 2.1e6),
       ("1.6e6", 1.6e6)]

o = P.RIGID_CONTACT["k_n_per_m"]
print("@0.42 m — 能量账自检（absorbed = (into_pad+into_flex)/impact）")
print(f"{'接触设定':<24} {'档':<6} {'impact_J':>9} {'into_pad':>9} {'into_flex':>10} "
      f"{'absorbed':>9} {'resid_J':>9} {'判定':>8}")
bad = []
try:
    for tag, kn in SET:
        P.RIGID_CONTACT["k_n_per_m"] = kn
        for lab in LAB:
            r = P.simulate_boulder_fall(height_m=0.42, mass_kg=80.0, on_pad=False,
                                        posture=f"anchored-{lab}")
            af = float(r.absorbed_fraction)
            flag = "OK" if af <= 1.0 else "**越界**"
            if af > 1.0:
                bad.append((tag, lab, af))
            print(f"{tag:<24} {lab:<6} {r.impact_energy_j:>9.1f} "
                  f"{r.energy_into_pad_j:>9.1f} {r.energy_into_flex_j:>10.1f} "
                  f"{af:>9.3f} {r.energy_residual_j:>9.1f} {flag:>8}")
finally:
    P.RIGID_CONTACT["k_n_per_m"] = o
    print(f"\n(已恢复 {o!r})")

print()
if bad:
    print("!! 越界格：")
    for t, l, a in bad:
        print(f"   {t} / {l}: absorbed={a:.3f}")
    print("\n!! 若 A-J1 的判定档（2.1e6）在其中 ⇒ A 的判定必须撤回重跑")
else:
    print("=> @0.42 m（A 的判定高度）无能量越界")
