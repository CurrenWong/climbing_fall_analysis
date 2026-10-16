"""A 件之后的力判据首折 h*（附加诊断报告，**不改 A 的冻结判定**）。
网格与阈值与 FP_ANCHOR 完全同口径：h ∈ {0.2,0.5,1,2,5,7,9} m，F_flex 线性内插到 6.2 kN。
只读 + 临时门（finally 恢复）。"""
import sys
from pathlib import Path

ROOT = Path(r'D:\Project\climbing_fall_analysis')
sys.path.insert(0, str(ROOT / "src"))
import climbing.pad as P  # noqa: E402

BW = 80.0 * 9.80665
F_THR = 6200.0
GRID = [0.20, 0.50, 1.00, 2.00, 5.00, 7.00, 9.00]
LAB = ("soft", "self", "stiff")
SETTINGS = [("None(刚性 548kN/m)", None), ("3.0e6(177kN/m)", 3.0e6),
            ("2.1e6(128kN/m)", 2.1e6), ("1.6e6(100kN/m)", 1.6e6)]


def hstar_of(vals):
    for i in range(len(GRID) - 1):
        f0, f1 = vals[i], vals[i + 1]
        if f0 <= F_THR <= f1:
            return GRID[i] + (F_THR - f0) / (f1 - f0) * (GRID[i + 1] - GRID[i])
    return None


orig = P.RIGID_CONTACT["k_n_per_m"]
out = {}
try:
    for tag, kn in SETTINGS:
        P.RIGID_CONTACT["k_n_per_m"] = kn
        print(f"\n=== {tag} ===")
        print(f"  {'档':<6} " + " ".join(f"{h:>7.2f}m" for h in GRID) + f" {'h*(m)':>8}")
        for lab in LAB:
            ff, fp = [], []
            for h in GRID:
                r = P.simulate_boulder_fall(height_m=h, mass_kg=80.0, on_pad=False,
                                            posture=f"anchored-{lab}")
                ff.append(float(r.meta["f_flex_peak_n"]))
                fp.append(float(r.peak_force_n))
            hs = hstar_of(ff)
            out.setdefault(lab, {})[tag] = {"h_star_m": hs, "f_flex_n": ff, "f_pad_n": fp}
            print(f"  {lab:<6} " + " ".join(f"{v/1e3:>7.2f} " for v in ff)
                  + f" {('%.3f' % hs) if hs else '>9':>8}")
finally:
    P.RIGID_CONTACT["k_n_per_m"] = orig
    print(f"\n(已恢复 RIGID_CONTACT = {orig!r})")

print("\n=== 首折 h* 汇总（力判据 F_flex → 6.2 kN）===\n")
print(f"  {'接触设定':<22} " + " ".join(f"{l:>9}" for l in LAB))
for tag, _ in SETTINGS:
    print(f"  {tag:<22} " + " ".join(
        f"{(out[l][tag]['h_star_m'] if out[l][tag]['h_star_m'] else float('nan')):>8.3f}m"
        for l in LAB))
print("\n  参照：流行病学 1.1–4 m；临床 ≤5 ft(1.52 m) 亦可致跟骨骨折；[Y25] 靶 7~9 m")
