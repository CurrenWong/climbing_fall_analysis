"""A/B 对拍：meta['f_flex_peak_n'] 在两个探针与我新脚本间差 35%（7002 vs 5198）。
唯一可疑差异 = 传给 `m_theta_table` 的**格式**（`_TAB` 列表对 vs `(TAB_TH, TAB_M)` 元组）。
只读 + 临时门。"""
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

ROOT = Path(r'D:\Project\climbing_fall_analysis')
sys.path.insert(0, str(ROOT / "src"))
import climbing.pad as P  # noqa: E402

MASS, L, V_REF, M_ECC, K_TH = 80.0, 0.425, 1.0, 540.0, 350.0
_A = [(0.0, 1.30), (15.0, 1.79), (45.0, 2.90), (60.0, 3.20), (90.0, 2.78), (135.0, 2.78)]
_D = np.arange(0.0, 135.0 + 1e-9, 2.5)
_FIT = PchipInterpolator([a for a, _ in _A], [b * MASS for _, b in _A])
# 探针的构造方式（unload_process_probe.py:34-38）
_TAB_PAIRS = [[float(t), float(_FIT(t))] for t in _D]
# 我的构造方式（verify_flex_observable.py）
_TUPLE = (np.radians(_D), _FIT(_D))
KEYS = ("k_theta_nm_per_rad", "segment_l_m", "theta_max_rad",
        "m_ecc_max_nm", "ecc_rate_ref_ms", "m_theta_table")
TH_MAX = float(2.0 * np.arccos(1.0 - 0.37 / (2.0 * L)))


def run(table, h=1.0):
    prev = {k: P.TWO_SEG_LEG.get(k) for k in KEYS}
    c0 = P.RIGID_CONTACT.get("k_n_per_m")
    try:
        P.TWO_SEG_LEG.update({"k_theta_nm_per_rad": K_TH, "segment_l_m": L,
                              "theta_max_rad": TH_MAX, "m_ecc_max_nm": M_ECC,
                              "ecc_rate_ref_ms": V_REF, "m_theta_table": table})
        P.RIGID_CONTACT["k_n_per_m"] = 2.1e6
        return P.simulate_boulder_fall(height_m=h, mass_kg=MASS, on_pad=False,
                                       posture="anchored-self")
    finally:
        for k in KEYS:
            P.TWO_SEG_LEG[k] = prev[k]
        P.RIGID_CONTACT["k_n_per_m"] = c0


print("=" * 90)
print("A/B：table 格式对 meta['f_flex_peak_n'] 的影响（其余完全相同；h=1.0 m）")
print("=" * 90)
# 先看两种格式在 np.interp 下的实际读数
for nm, tab in (("_TAB 列表对（探针）", _TAB_PAIRS), ("(TAB_TH,TAB_M) 元组（我）", _TUPLE)):
    xp, fp = np.asarray(tab[0], float), np.asarray(tab[1], float)
    print(f"  {nm:<24} xp[:3]={xp[:3]}  fp[:3]={fp[:3]}")
    for th in (0.3, 1.0):
        print(f"      interp(θ={th}) = {float(np.interp(th, xp, fp)):9.3f} N·m")
print()
print(f"  {'config':<26} {'f_flex_N':>10} {'F_pad_N':>10} {'n_t':>7} {'window_s':>22}")
for nm, tab in (("_TAB 列表对（探针）", _TAB_PAIRS), ("(TAB_TH,TAB_M) 元组（我）", _TUPLE)):
    r = run(tab)
    print(f"  {nm:<26} {r.meta['f_flex_peak_n']:>10.1f} {r.peak_force_n:>10.1f} "
          f"{len(r.t_s):>7d} {str(r.meta.get('impact_window_s')):>22}")
print()
print("  ⇒ 若两者差 ~35% ⇒ **探针传的表格式是错的**（旧结果 7002 无效，我新脚本的 5198 才对）")
