"""跑 `scripts/pad_energy_check.py` 的有效性审计：
 ① A 的 9 个格（@0.42 m，3 接触设定 × 3 档）；② 非单调格（anchored-self，h=2/5/9）；
 ③ 顺带把 `f_flex` 峰值拆成「弹簧项 / 阻尼项」（为定性那个非单调）。
只读 + 临时门（finally 恢复）。"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(r'D:\Project\climbing_fall_analysis')
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
import climbing.pad as P                                    # noqa: E402
from scripts.pad_energy_check import audit_energy, elastic_energy_flex  # noqa: E402

G, LEG0, C_FLEX, K_STOP = 9.80665, 0.85, 400.0, 3.0e5
LAB = ("soft", "self", "stiff")
BW = 80.0 * G


def decompose(r, post, mass=80.0):
    """把 `f_flex` 峰值拆成弹簧项与阻尼项（复现 pad.py 的 flex_force 口径）。"""
    t = np.asarray(r.t_s, float)
    zf, zt = np.asarray(r.z_foot_m, float), np.asarray(r.z_torso_m, float)
    xf = np.maximum(0.0, LEG0 - (zf - zt))
    rate = np.maximum(0.0, np.gradient(xf, t))
    f_spr = np.where(xf <= post.x_flex_max,
                     post.k_flex * xf,
                     post.k_flex * post.x_flex_max + K_STOP * (xf - post.x_flex_max))
    f_dmp = C_FLEX * rate
    return xf, f_spr, f_dmp


o = P.RIGID_CONTACT["k_n_per_m"]
print("=" * 96)
print("① 有效性审计（沿轨迹 E_mech 非增）—— A 的判定格 @0.42 m")
print("=" * 96)
print(f"  {'设定':>9} {'档':<6} {'KE(窗首)J':>10} {'violation(J)':>13} {'容差(J)':>9} "
      f"{'max xf(m)':>10} {'max comp(mm)':>13} {'单调':>6} {'判定':>8}")
try:
    for kn in (None, 3.0e6, 2.1e6, 1.6e6):
        P.RIGID_CONTACT["k_n_per_m"] = kn
        for lab in LAB:
            r = P.simulate_boulder_fall(height_m=0.42, mass_kg=80.0, on_pad=False,
                                        posture=f"anchored-{lab}")
            a = audit_energy(r, P.POSTURES[f"anchored-{lab}"], 80.0)
            tag = "None" if kn is None else f"{kn:.1e}"
            print(f"  {tag:>9} {lab:<6} {a['ke0_j']:>10.1f} {a['violation_j']:>13.2f} "
                  f"{a['tolerance_j']:>9.2f} {a['max_xf_m']:>10.4f} "
                  f"{a['max_comp_m']*1e3:>13.2f} "
                  f"{('是' if a['monotone'] else '否'):>6} "
                  f"{('PASS' if a['pass'] else 'FAIL'):>8}")

    print()
    print("=" * 96)
    print("② 非单调格（anchored-self）+ ③ f_flex 峰值分解")
    print("=" * 96)
    for kn in (None, 2.1e6):
        tag = "None" if kn is None else f"{kn:.1e}"
        P.RIGID_CONTACT["k_n_per_m"] = kn
        print(f"\n  --- k_n = {tag} · anchored-self ---")
        print(f"  {'h(m)':>5} {'F_flex(kN)':>11} {'弹簧峰(kN)':>11} {'阻尼峰(kN)':>11} "
              f"{'max xf(m)':>10} {'max rate(m/s)':>14} {'violation(J)':>13} {'审计':>6}")
        for h in (0.5, 1.0, 2.0, 5.0, 9.0):
            r = P.simulate_boulder_fall(height_m=h, mass_kg=80.0, on_pad=False,
                                        posture="anchored-self")
            post = P.POSTURES["anchored-self"]
            xf, f_spr, f_dmp = decompose(r, post)
            a = audit_energy(r, post, 80.0)
            print(f"  {h:>5.1f} {r.meta['f_flex_peak_n']/1e3:>11.3f} "
                  f"{f_spr.max()/1e3:>11.3f} {f_dmp.max()/1e3:>11.3f} "
                  f"{xf.max():>10.4f} {np.max(np.gradient(xf, np.asarray(r.t_s, float))):>14.3f} "
                  f"{a['violation_j']:>13.2f} {('PASS' if a['pass'] else 'FAIL'):>6}")
finally:
    P.RIGID_CONTACT["k_n_per_m"] = o
    print(f"\n(已恢复 {o!r})")
print("\n审计判据：`violation_j ≤ 1%·mgh` 且 `E_mech` 单调 ⇒ 数值物理有效。")
