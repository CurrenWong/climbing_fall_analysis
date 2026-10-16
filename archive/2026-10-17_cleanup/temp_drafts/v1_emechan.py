"""★★★ 判定性检验：沿轨迹直接算 `E_mech(t) = KE + PE_grav + E_elastic` 是否单调不增。

参照（E_mech(0) = 0）：
  PE_grav(t) = −[ m_low·g·z_foot(t) + m_up·g·(z_torso(t) + leg_ext0) ]
  KE(t)      = ½·m_low·vzf² + ½·m_up·vzt²，   vz = ∫az dt（az = −g·accel_*_g，仅用模型报的量）
  E_elastic  = ∫F_pad dx + [½k·x²  或  ½k·xm² + k·xm·(x−xm) + ½k_stop·(x−xm)²]
判据：`E_mech(t) ≤ 0`（阻尼只耗散）——**违反 ⇒ 能量被凭空创造（真缺陷）**。
只用模型自报量，不依赖它的 `absorbed_fraction` 口径。
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(r'D:\Project\climbing_fall_analysis')
sys.path.insert(0, str(ROOT / "src"))
import climbing.pad as P  # noqa: E402

G = 9.80665
LEG0 = 0.85
K_STOP = 3.0e5          # pad.py 内固定
H = 0.42
LAB = ("soft", "self", "stiff")


def elastic_flex(x, k, xm):
    x = np.maximum(x, 0.0)
    lo = np.minimum(x, xm)
    e = 0.5 * k * lo ** 2
    over = np.maximum(x - xm, 0.0)
    return e + k * xm * over + 0.5 * K_STOP * over ** 2


o = P.RIGID_CONTACT["k_n_per_m"]
print(f"{'设定':>9} {'档':<6} {'max E_mech(J)':>14} {'@t(s)':>8} {'E_mech(末)(J)':>14} "
      f"{'max xf(m)':>10} {'E_flex报(J)':>12} {'½k·maxxf²(J)':>13} {'判定':>12}")
try:
    for kn in (None, 2.1e6, 1.6e6):
        P.RIGID_CONTACT["k_n_per_m"] = kn
        for lab in LAB:
            post = P.POSTURES[f"anchored-{lab}"]
            r = P.simulate_boulder_fall(height_m=H, mass_kg=80.0, on_pad=False,
                                        posture=f"anchored-{lab}")
            t = np.asarray(r.t_s)
            zf, zt = np.asarray(r.z_foot_m), np.asarray(r.z_torso_m)
            al = np.asarray(r.accel_leg_g) * G       # az_f = −g·accel_leg_g
            at = np.asarray(r.accel_torso_g) * G
            m_low, m_up = post.split_mass(80.0)
            # 速度：∫az dt（v(0)=0）
            vzf = np.concatenate([[0.0], np.cumsum(0.5 * (al[1:] + al[:-1]) * np.diff(t))])
            vzt = np.concatenate([[0.0], np.cumsum(0.5 * (at[1:] + at[:-1]) * np.diff(t))])
            ke = 0.5 * m_low * vzf ** 2 + 0.5 * m_up * vzt ** 2
            pe = -(m_low * G * zf + m_up * G * (zt + LEG0))
            xf = LEG0 - (zf - zt)
            comp = np.maximum(zf - H, 0.0)
            e_flex = elastic_flex(xf, post.k_flex, post.x_flex_max)
            f_pad = np.asarray(r.pad_force_n)
            e_pad = np.concatenate([[0.0], np.cumsum(
                0.5 * (f_pad[1:] + f_pad[:-1]) * np.diff(comp))])
            e_mech = ke + pe + e_flex + e_pad
            i = int(np.argmax(e_mech))
            bad = e_mech[i] > 1.0
            tag = "None" if kn is None else f"{kn:.1e}"
            print(f"{tag:>9} {lab:<6} {e_mech[i]:>14.2f} {t[i]:>8.5f} {e_mech[-1]:>14.2f} "
                  f"{xf.max():>10.4f} {r.energy_into_flex_j:>12.1f} "
                  f"{0.5*post.k_flex*xf.max()**2:>13.1f} "
                  f"{'**创造**' if bad else 'OK':>12}")
        print()
finally:
    P.RIGID_CONTACT["k_n_per_m"] = o
    print(f"(已恢复 {o!r})")
print("\n判据：E_mech(t) 必须 ≤ 0（起点 0，只允许阻尼耗散下降）。")
