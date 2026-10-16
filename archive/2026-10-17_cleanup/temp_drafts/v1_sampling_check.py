"""判定 meta['f_flex_peak_n'] 是否**采样依赖**（伪量）。只读 + 临时门。"""
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

ROOT = Path(r'D:\Project\climbing_fall_analysis')
sys.path.insert(0, str(ROOT / "src"))
import climbing.pad as P  # noqa: E402

MASS = 80.0
L, V_REF, M_ECC, K_TH = 0.425, 1.0, 540.0, 350.0
_T = np.arange(0.0, 135.0 + 1e-9, 2.5)
_A = [(0.0, 1.30), (15.0, 1.79), (45.0, 2.90), (60.0, 3.20), (90.0, 2.78), (135.0, 2.78)]
TAB = (np.radians(_T), PchipInterpolator([a for a, _ in _A], [b * MASS for _, b in _A])(_T))
KEYS = ("k_theta_nm_per_rad", "segment_l_m", "theta_max_rad",
        "m_ecc_max_nm", "ecc_rate_ref_ms", "m_theta_table")
TH_MAX = float(2.0 * np.arccos(1.0 - 0.37 / (2.0 * L)))


def run(h, **kw):
    prev = {k: P.TWO_SEG_LEG.get(k) for k in KEYS}
    c0 = P.RIGID_CONTACT.get("k_n_per_m")
    try:
        P.TWO_SEG_LEG.update({"k_theta_nm_per_rad": K_TH, "segment_l_m": L,
                              "theta_max_rad": TH_MAX, "m_ecc_max_nm": M_ECC,
                              "ecc_rate_ref_ms": V_REF, "m_theta_table": TAB})
        P.RIGID_CONTACT["k_n_per_m"] = 2.1e6
        return P.simulate_boulder_fall(height_m=h, mass_kg=MASS, on_pad=False,
                                       posture="anchored-self", **kw)
    finally:
        for k in KEYS:
            P.TWO_SEG_LEG[k] = prev[k]
        P.RIGID_CONTACT["k_n_per_m"] = c0


print("=" * 84)
print("A. 同配置、只改采样参数 ⇒ f_flex_peak_n 变不变？")
print("=" * 84)
print(f"  {'dt_sample':>10} {'max_ode_step':>13} {'f_flex_peak_N':>14} {'n_t':>7} {'F_pad_N':>10}")
for dt, ms in ((1e-3, None), (1e-4, None), (1e-3, 1e-3), (1e-4, 1e-4), (1e-5, 1e-4)):
    kw = {"dt_sample": dt}
    if ms is not None:
        kw["max_ode_step"] = ms
    r = run(1.0, **kw)
    print(f"  {dt:>10.0e} {str(ms):>13} {r.meta['f_flex_peak_n']:>14.1f} "
          f"{len(r.t_s):>7d} {r.peak_force_n:>10.1f}")
print("  ⇒ 若随采样变 ⇒ **观测量是采样依赖的伪量**，不是物理量。")
