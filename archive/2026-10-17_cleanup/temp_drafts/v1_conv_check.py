"""★ 最后一个未决项：软接触下 `F_flex(h)` 的非单调是真物理还是积分精度？

对非单调格（anchored-self @ k_n=2.1e6，h=2 与 5 m）细化解算器设置，看是否收敛。
只读 + 临时门。"""
import sys
from pathlib import Path

ROOT = Path(r'D:\Project\climbing_fall_analysis')
sys.path.insert(0, str(ROOT / "src"))
import climbing.pad as P  # noqa: E402

CASES = [
    ("default 1e-8/1e-11/1e-4", dict()),
    ("mid     1e-9/1e-12/3e-5", dict(rtol=1e-9, atol=1e-12, max_ode_step=3e-5)),
    ("tight   1e-10/1e-13/1e-5", dict(rtol=1e-10, atol=1e-13, max_ode_step=1e-5)),
]
HEIGHTS = (2.0, 5.0, 9.0)

o = P.RIGID_CONTACT["k_n_per_m"]
try:
    for kn in (None, 2.1e6):
        P.RIGID_CONTACT["k_n_per_m"] = kn
        tag = "None" if kn is None else f"{kn:.1e}"
        print(f"\n=== k_n = {tag} · anchored-self ===")
        print(f"  {'解算器':<26} " + " ".join(f"{h:>10.1f}m" for h in HEIGHTS)
              + f" {'单调?':>7}")
        prev = None
        for lbl, kw in CASES:
            vals = []
            for h in HEIGHTS:
                r = P.simulate_boulder_fall(height_m=h, mass_kg=80.0, on_pad=False,
                                            posture="anchored-self", **kw)
                vals.append(float(r.meta["f_flex_peak_n"]))
            mono = all(vals[i] < vals[i + 1] for i in range(len(vals) - 1))
            print(f"  {lbl:<26} " + " ".join(f"{v/1e3:>10.3f}" for v in vals)
                  + f" {('单调' if mono else '非单调'):>7}")
            prev = vals
finally:
    P.RIGID_CONTACT["k_n_per_m"] = o
    print(f"\n(已恢复 {o!r})")
print("\n判据：三个解算器设置给出同一组值 ⇒ 收敛 ⇒ 非单调是真物理；否则是数值。")
