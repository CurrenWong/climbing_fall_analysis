"""核查：修复后的 f_flex 读数 vs 力台外部锚 + 力的饱和性 + 两神读数对比。
只读 + 临时门（finally 恢复）；不改 pad.py。
"""
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

ROOT = Path(r'D:\Project\climbing_fall_analysis')
sys.path.insert(0, str(ROOT / "src"))
import climbing.pad as P  # noqa: E402

MASS, G = 80.0, 9.80665
L, V_REF, M_ECC, R_REF, K_TH = 0.425, 1.0, 540.0, 0.040, 350.0
ANCH = [(0.0, 1.30), (15.0, 1.79), (45.0, 2.90), (60.0, 3.20), (90.0, 2.78), (135.0, 2.78)]
_D = np.arange(0.0, 135.0 + 1e-9, 2.5)
TAB_TH = np.radians(_D)
TAB_M = PchipInterpolator([a for a, _ in ANCH], [b * MASS for _, b in ANCH])(_D)
KEYS = ("k_theta_nm_per_rad", "segment_l_m", "theta_max_rad",
        "m_ecc_max_nm", "ecc_rate_ref_ms", "m_theta_table")


def th_of(s):
    return 2.0 * np.arccos(np.clip(1.0 - np.asarray(s, float) / (2.0 * L), -1.0, 1.0))


def run(lab, h):
    prev = {k: P.TWO_SEG_LEG.get(k) for k in KEYS}
    c0 = P.RIGID_CONTACT.get("k_n_per_m")
    try:
        P.TWO_SEG_LEG.update({"k_theta_nm_per_rad": K_TH, "segment_l_m": L,
                              "theta_max_rad": float(th_of(0.37)),
                              "m_ecc_max_nm": M_ECC, "ecc_rate_ref_ms": V_REF,
                              "m_theta_table": (TAB_TH, TAB_M)})
        P.RIGID_CONTACT["k_n_per_m"] = 2.1e6
        return P.simulate_boulder_fall(height_m=h, mass_kg=MASS, on_pad=False,
                                       posture=f"anchored-{lab}")
    finally:
        for k in KEYS:
            P.TWO_SEG_LEG[k] = prev[k]
        P.RIGID_CONTACT["k_n_per_m"] = c0


print("=" * 78)
print("A. 修复后的 f_flex 读数（meta = rhs 同一个 flex_force）")
print("=" * 78)
print(f"  {'h_m':>5} {'f_flex_N':>10} {'n_samples':>10}")
vals = {}
for h in (0.42, 0.5, 1.0, 2.0, 5.0):
    r = run("self", h)
    vals[h] = float(r.meta["f_flex_peak_n"])
    print(f"  {h:>5.2f} {vals[h]:>10.1f} {len(r.t_s):>10d}")

print()
print("=" * 78)
print("B. ★ 外部锚对照：力台 GRFmax（PMC9465758, 0.42 m 落箱, n=24 运动员）")
print("=" * 78)
anchor = 5.36 * MASS * G
print(f"  力台锚 GRFmax = 5.36 BW = {anchor:.0f} N  （stiff 档；soft/self 2.76/3.73 BW）")
f042 = vals[0.42]
print(f"  模型 @0.42 m f_flex = {f042:.0f} N")
print(f"  => 比值 = {f042 / anchor:.2f}x  （f_flex 是腿内力，GRF ≥ f_flex ⇒ 模型高估至少这个倍数）")
for bw in (2.76, 3.73):
    print(f"     若换 {bw} BW 锚（{bw*MASS*G:.0f} N）：比值 {f042/(bw*MASS*G):.2f}x")

print()
print("=" * 78)
print("C. 力的饱和性：f_ecc_max = M_ach(θ)/(l·sin(θ/2)) 几乎不随 θ 变（r 归零因子抵消）")
print("=" * 78)
print(f"  {'θ_deg':>7} {'r_arm_mm':>9} {'M_ach_Nm':>9} {'f_ecc_max_N':>12} {'rate/(rate+1)@3.1':>18}")
for th in (0.05, 0.15, 0.3, 0.6, 1.0, 1.5):
    r_arm = R_REF * min(1.0, th / 0.349)
    m_ach = M_ECC * (r_arm / R_REF)
    sh = float(np.sin(0.5 * th))
    fe = m_ach / (L * sh) if sh > 1e-9 else 0.0
    print(f"  {np.degrees(th):>7.1f} {r_arm*1e3:>9.2f} {m_ach:>9.1f} {fe:>12.0f} "
          f"{3.13/(3.13+1.0):>18.3f}")
f_cap_min = M_ECC * (2.0 / 0.349) / (L * 1.0) * min(1.0, 1.0)
print(f"  => f_ecc_max 在 θ→0 极限 ≈ {M_ECC*(R_REF/0.349)*(2.0)/(L*R_REF):.0f} N"
      f"  ⇒ 率项上限 × 0.76 ≈ {M_ECC*(R_REF/0.349)*(2.0)/(L*R_REF)*3.13/4.13:.0f} N")
print("  ⇒ 曲线 0.5→5 m 只从 6421 升到 7818 N（+22%）= **接近饱和的平台**")
print("  ⇒ 用固定阈 6200 N 切一条平台 ⇒ h* 对阈值/参数**病态敏感**（这才是『又回去了』的真因）")
