"""为什么 D1 的力比 C1 高？—— 同格对比「弹性项 / 率项 / 率(rate)」，C1 vs D1。
（C1 = 常数 k_θ·θ；D1 = 力-长表 M(θ)。其余完全相同。）
只读 + 临时门（finally 恢复）。
"""
import sys
from math import degrees
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

ROOT = Path(r'D:\Project\climbing_fall_analysis')
sys.path.insert(0, str(ROOT / "src"))
import climbing.pad as P  # noqa: E402

G, LEG0, L, V_REF, M_ECC, R_REF, RAMP, PAD_KN = 9.80665, 0.85, 0.425, 1.0, 540.0, 0.040, 0.349, 2.1e6
K_TH = 350.0
COM = {"soft": 0.43, "self": 0.37, "stiff": 0.23}
ANCH = [(0.0, 1.30), (15.0, 1.79), (45.0, 2.90), (60.0, 3.20), (90.0, 2.78), (135.0, 2.78)]
_D = np.arange(0.0, 135.0 + 1e-9, 2.5)
TAB = [[float(t), float(PchipInterpolator([a for a, _ in ANCH], [b * 80 for _, b in ANCH])(t))]
       for t in _D]
TAB_TH, TAB_M = np.radians([t for t, _ in TAB]), np.array([m for _, m in TAB])
KEYS = ("k_theta_nm_per_rad", "segment_l_m", "theta_max_rad", "m_ecc_max_nm",
        "ecc_rate_ref_ms", "m_theta_table")


def th_of(s):
    return 2.0 * np.arccos(np.clip(1.0 - np.asarray(s, float) / (2.0 * L), -1.0, 1.0))


def run(lab, h, use_table):
    prev = {k: P.TWO_SEG_LEG.get(k) for k in KEYS}
    c0 = P.RIGID_CONTACT.get("k_n_per_m")
    try:
        P.TWO_SEG_LEG.update({"k_theta_nm_per_rad": K_TH, "segment_l_m": L,
                              "theta_max_rad": float(th_of(COM[lab])),
                              "m_ecc_max_nm": M_ECC, "ecc_rate_ref_ms": V_REF,
                              "m_theta_table": TAB if use_table else None})
        P.RIGID_CONTACT["k_n_per_m"] = PAD_KN
        r = P.simulate_boulder_fall(height_m=h, mass_kg=80.0, on_pad=False,
                                    posture=f"anchored-{lab}")
    finally:
        for k in KEYS:
            P.TWO_SEG_LEG[k] = prev[k]
        P.RIGID_CONTACT["k_n_per_m"] = c0
    t = np.asarray(r.t_s, float)
    zf, zt = np.asarray(r.z_foot_m, float), np.asarray(r.z_torso_m, float)
    t0, t1 = r.meta["impact_window_s"]
    i0 = int(np.searchsorted(t, t0)); i1 = max(int(np.searchsorted(t, t1)), i0 + 2)
    s = np.maximum(0.0, LEG0 - (zf - zt))
    th = th_of(s)
    sh = np.sin(0.5 * np.maximum(th, 1e-9))
    rate = np.maximum(0.0, np.gradient(s, t))
    r_arm = R_REF * np.minimum(1.0, th / RAMP)
    m_el = np.interp(th, TAB_TH, TAB_M) if use_table else K_TH * th
    if use_table:
        m_el = m_el * np.minimum(1.0, th / RAMP)
    el = np.where(sh > 1e-9, m_el / (L * sh), 0.0)
    cap = np.where(sh > 1e-9, (M_ECC * (r_arm / R_REF)) / (L * sh), 0.0)
    rt = cap * rate / (rate + V_REF)
    return {"F": float(r.meta["f_flex_peak_n"]), "el": float(el[i0:i1].max()),
            "rt": float(rt[i0:i1].max()), "rate": float(rate[i0:i1].max()),
            "th": degrees(float(th[i0:i1].max())), "s": float(s[i0:i1].max())}


for lab in ("stiff", "soft"):
    h = 2.0
    a = run(lab, h, False)
    b = run(lab, h, True)
    print(f"=== {lab} @ {h} m ===")
    print(f"  {'':>16} {'F_flex(kN)':>11} {'弹性(kN)':>9} {'率项(kN)':>9} "
          f"{'max rate(m/s)':>14} {'θ@峰°':>7} {'s@峰(m)':>8}")
    for tag, d in (("C1 常数 k_θ·θ", a), ("D1 力-长 M(θ)", b)):
        print(f"  {tag:>16} {d['F']/1e3:>11.3f} {d['el']/1e3:>9.3f} {d['rt']/1e3:>9.3f} "
              f"{d['rate']:>14.3f} {d['th']:>7.1f} {d['s']:>8.4f}")
    print(f"  ⇒ ΔF = {(b['F']-a['F'])/1e3:+.3f} kN ；Δ弹性 = {(b['el']-a['el'])/1e3:+.3f} "
          f"；Δ率项 = {(b['rt']-a['rt'])/1e3:+.3f} ；Δrate = {b['rate']-a['rate']:+.3f}")
    print()
print("(门已恢复)")
