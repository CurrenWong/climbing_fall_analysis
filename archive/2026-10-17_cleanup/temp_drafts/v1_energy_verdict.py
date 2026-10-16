"""★ 决定性判据：`absorbed_fraction>1` 是「物理缺陷」还是「口径错配」？

几何（pad.py:632-640）：xf = leg_ext0 − (z_foot − z_torso)，leg_ext0 = 0.85，
          y0 = [z_foot=0, vz=0, z_torso=−0.85, vz=0]，z 向下为正。
⇒ 任一时刻 PE 释放量 = g·[ m_low·z_foot + m_up·(z_torso + 0.85) ]
                   = g·[ m·z_foot + m_up·xf ]        （代入 z_torso = xf − 0.85 + z_foot）
而 impact_energy_j = m·g·h 只在**首次接触**时刻成立。

判据（必要充分，无需速度）：
  末态弹性储能  E_el = energy_into_pad_j + energy_into_flex_j
  必须 ≤ PE 释放量（E_el + 阻尼耗散 + 末态动能 = PE 释放量）
  ⇒ 若 E_el > PE_释放 ⇒ **真缺陷**；若 E_el ≤ PE_释放 ⇒ **口径错配，物理没问题**。
另做两条内部一致性：① E_flex ≟ ½k·xf²；② ∫F_pad·d(comp) ≟ energy_into_pad_j。
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(r'D:\Project\climbing_fall_analysis')
sys.path.insert(0, str(ROOT / "src"))
import climbing.pad as P  # noqa: E402

G = 9.80665
LEG0 = 0.85
LAB = ("soft", "self", "stiff")

o = P.RIGID_CONTACT["k_n_per_m"]
print(f"{'设定':>10} {'档':<6} {'xf_end(m)':>10} {'z_foot_end':>11} {'comp(mm)':>9} "
      f"{'PE释放(J)':>10} {'E_el(J)':>9} {'E_el/PE':>8} {'½kxf²(J)':>10} "
      f"{'E_flex(J)':>10} {'∫Fdx(J)':>9} {'E_pad(J)':>9} {'判定':>10}")
try:
    for kn in (None, 2.1e6):
        P.RIGID_CONTACT["k_n_per_m"] = kn
        for lab in LAB:
            post = P.POSTURES[f"anchored-{lab}"]
            r = P.simulate_boulder_fall(height_m=0.42, mass_kg=80.0, on_pad=False,
                                        posture=f"anchored-{lab}")
            zf, zt = np.asarray(r.z_foot_m), np.asarray(r.z_torso_m)
            zf_e, zt_e = float(zf[-1]), float(zt[-1])
            comp = np.maximum(0.0, zf - 0.42)
            xf_e = LEG0 - (zf_e - zt_e)
            m_low, m_up = post.split_mass(80.0)
            pe = G * (80.0 * zf_e + m_up * xf_e)
            e_el = r.energy_into_pad_j + r.energy_into_flex_j
            # 内部一致性
            kxf2 = 0.5 * post.k_flex * max(xf_e, 0.0) ** 2
            int_fdx = float(np.trapezoid(r.pad_force_n, comp)) if hasattr(np, "trapezoid") \
                else float(np.trapz(r.pad_force_n, comp))
            verdict = "真缺陷" if e_el > pe * 1.001 else "口径错配 ⭕"
            tag = "None" if kn is None else f"{kn:.1e}"
            print(f"{tag:>10} {lab:<6} {xf_e:>10.4f} {zf_e:>11.4f} "
                  f"{max(comp)*1e3:>9.2f} {pe:>10.1f} {e_el:>9.1f} {e_el/pe:>8.3f} "
                  f"{kxf2:>10.1f} {r.energy_into_flex_j:>10.1f} {int_fdx:>9.1f} "
                  f"{r.energy_into_pad_j:>9.1f} {verdict:>10}")
        print()
finally:
    P.RIGID_CONTACT["k_n_per_m"] = o
    print(f"(已恢复 {o!r})")

print("\n注：impact_energy_j（=mgh=329.5 J）只在首次接触时刻成立；")
print("    `absorbed_fraction` 的分母用它、分子是**末态**弹性储能 ⇒ 两者不同时刻。")
