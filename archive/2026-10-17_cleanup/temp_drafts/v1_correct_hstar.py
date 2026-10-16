"""用**正确表**+**已修后处理**重测首折 h*。只读 + 临时门。
表格式：必须是 (TAB_TH_rad, TAB_M) 元组 —— 验证：np.interp(1.0 rad) 应 ≈ 254.7 N·m。
"""
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

ROOT = Path(r'D:\Project\climbing_fall_analysis')
sys.path.insert(0, str(ROOT / "src"))
import climbing.pad as P  # noqa: E402

MASS, L, V_REF, M_ECC, K_TH, THR = 80.0, 0.425, 1.0, 540.0, 350.0, 6200.0
_A = [(0.0, 1.30), (15.0, 1.79), (45.0, 2.90), (60.0, 3.20), (90.0, 2.78), (135.0, 2.78)]
_D = np.arange(0.0, 135.0 + 1e-9, 2.5)
TAB_TH = np.radians(_D)
TAB_M = PchipInterpolator([a for a, _ in _A], [b * MASS for _, b in _A])(_D)
assert abs(float(np.interp(1.0, TAB_TH, TAB_M)) - 254.7) < 2.0, "表格式错！"
print(f"表自检：interp(1.0 rad) = {float(np.interp(1.0, TAB_TH, TAB_M)):.1f} N·m (应 ~254.7) OK")
KEYS = ("k_theta_nm_per_rad", "segment_l_m", "theta_max_rad",
        "m_ecc_max_nm", "ecc_rate_ref_ms", "m_theta_table")

def run(lab, h, dt=1e-4):
    prev = {k: P.TWO_SEG_LEG.get(k) for k in KEYS}
    c0 = P.RIGID_CONTACT.get("k_n_per_m")
    try:
        P.TWO_SEG_LEG.update({"k_theta_nm_per_rad": K_TH, "segment_l_m": L,
                              "theta_max_rad": float(2.0 * np.arccos(
                                  1.0 - {"soft": 0.43, "self": 0.37, "stiff": 0.23}[lab] / (2.0 * L))),
                              "m_ecc_max_nm": M_ECC, "ecc_rate_ref_ms": V_REF,
                              "m_theta_table": (TAB_TH, TAB_M)})
        P.RIGID_CONTACT["k_n_per_m"] = 2.1e6
        return P.simulate_boulder_fall(height_m=h, mass_kg=MASS, on_pad=False,
                                       posture=f"anchored-{lab}", dt_sample=dt)
    finally:
        for k in KEYS:
            P.TWO_SEG_LEG[k] = prev[k]
        P.RIGID_CONTACT["k_n_per_m"] = c0

print()
print(f"  {'lab':>6} " + "".join(f"{h:>8.1f}" for h in (0.5, 1, 2, 3, 5, 7, 9)))
res = {}
for lab in ("soft", "self", "stiff"):
    row = {}
    for h in (0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 9.0):
        r = run(lab, h)
        row[h] = float(r.meta["f_flex_peak_n"])
    res[lab] = row
    print(f"  {lab:>6} " + "".join(f"{row[h]/1e3:>8.3f}" for h in row) + "  kN")
print()
for lab in ("soft", "self", "stiff"):
    row = res[lab]
    hs = sorted(row)
    hsx = None
    for i in range(len(hs) - 1):
        a, b = row[hs[i]], row[hs[i + 1]]
        if a <= THR <= b:
            hsx = hs[i] + (THR - a) * (hs[i + 1] - hs[i]) / (b - a)
            break
    print(f"  {lab:>6} h* = {('%.2f m' % hsx) if hsx else '> 9 m（网格内未交叉）'}")
