"""SECOND_CONTACT 探针 —— 「二次接触（臀/背）」三问（Q1 解析可达性 / Q2 能量-行程账 / Q3 姿态对比）。

只读 + 临时门（finally 恢复）；**不改 pad.py**。
预注册：`docs/methods/SECOND_CONTACT_预注册.md`（判据与预测均已写死）
产物：`results/second_contact_probe.json`
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

ROOT = Path(r'D:\Project\climbing_fall_analysis')
sys.path.insert(0, str(ROOT / "src"))
import climbing.pad as P  # noqa: E402

# ---- 冻结口径（预注册 §2）------------------------------------------------
MASS, G = 80.0, 9.80665
L, LEG0, V_REF, M_ECC, R_REF, RAMP = 0.425, 0.85, 1.0, 540.0, 0.040, 0.349
K_TH, THRESH_N, PAD_KN = 350.0, 6200.0, 2.1e6
COM_DEPTH = {"soft": 0.23, "self": 0.37, "stiff": 0.43}
H_BODY = 1.75
H_ISCHIUM = (0.46 * H_BODY, 0.50 * H_BODY)      # 坐骨结节站立高度
ANCH = [(0.0, 1.30), (15.0, 1.79), (45.0, 2.90), (60.0, 3.20), (90.0, 2.78), (135.0, 2.78)]
_D = np.arange(0.0, 135.0 + 1e-9, 2.5)
TAB_TH, TAB_M = np.radians(_D), PchipInterpolator(
    [a for a, _ in ANCH], [b * MASS for _, b in ANCH])(_D)
KEYS = ("k_theta_nm_per_rad", "segment_l_m", "theta_max_rad",
        "m_ecc_max_nm", "ecc_rate_ref_ms", "m_theta_table")
GRID = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
POSTURES = [("anchored-self", "self"), ("controlled-drop", "self"), ("tuck-roll", "self")]


def th_of(s):
    return 2.0 * np.arccos(np.clip(1.0 - np.asarray(s, float) / (2.0 * L), -1.0, 1.0))


def run(posture, lab, h):
    prev = {k: P.TWO_SEG_LEG.get(k) for k in KEYS}
    c0 = P.RIGID_CONTACT.get("k_n_per_m")
    try:
        P.TWO_SEG_LEG.update({"k_theta_nm_per_rad": K_TH, "segment_l_m": L,
                              "theta_max_rad": float(th_of(COM_DEPTH[lab])),
                              "m_ecc_max_nm": M_ECC, "ecc_rate_ref_ms": V_REF,
                              "m_theta_table": (TAB_TH, TAB_M)})
        P.RIGID_CONTACT["k_n_per_m"] = PAD_KN
        r = P.simulate_boulder_fall(height_m=h, mass_kg=MASS, on_pad=False,
                                    posture=posture)
    finally:
        for k in KEYS:
            P.TWO_SEG_LEG[k] = prev[k]
        P.RIGID_CONTACT["k_n_per_m"] = c0
    return r


def pad_peak(r):
    for nm in ("f_pad_n", "f_pad", "pad_force_n", "f_contact_n"):
        v = getattr(r, nm, None)
        if v is not None:
            return float(np.max(np.abs(np.asarray(v, float))))
    return float("nan")


def cross_h(hs, fs, thr):
    """线性插值求 f(h)=thr 的 h。"""
    hs, fs = np.asarray(hs, float), np.asarray(fs, float)
    idx = np.where(np.diff(np.sign(fs - thr)) != 0)[0]
    if not len(idx):
        return float("nan")
    i = int(idx[0])
    return float(hs[i] + (thr - fs[i]) * (hs[i + 1] - hs[i]) / (fs[i + 1] - fs[i]))


print("=" * 100)
print("Q1 几何可达性（解析）：1D 垂直落地里，臀/背要碰地还得下降多少？腿给得起吗？")
print("=" * 100)
print(f"  坐骨结节站立高度 = 0.46~0.50 × H({H_BODY} m) = {H_ISCHIUM[0]:.3f} ~ {H_ISCHIUM[1]:.3f} m"
      f"  （与模型 leg_ext0={LEG0} m 同量级 ⇒ 与躯干同乘骨盆下降）")
print(f"  {'行程锚':>10} {'腿给(m)':>9} {'需下降(m)':>10} {'缺口(m)':>9}")
q1 = {}
for lab, d in COM_DEPTH.items():
    lo, hi = min(x - d for x in H_ISCHIUM), max(x - d for x in H_ISCHIUM)
    q1[lab] = {"available_m": d, "gap_lo_m": lo, "gap_hi_m": hi}
    print(f"  {lab:>10} {d:>9.3f} {H_ISCHIUM[0]:.2f}~{H_ISCHIUM[1]:.2f} "
          f"{lo:>9.3f}~{hi:.3f}")
gmin = min(v["gap_lo_m"] for v in q1.values())
gmax = max(v["gap_hi_m"] for v in q1.values())
print(f"  => 缺口 {gmin:.3f}~{gmax:.3f} m（= 腿行程的 {gmin/0.43:.1f}~{gmax/0.23:.1f} 倍）")
print(f"  => J1（缺口>0.20 m ⇒ 1D 不可达）：{'PASS' if gmin > 0.20 else 'FAIL'}\n")

print("=" * 100)
print("Q2 能量-行程账（解析，恒力理想界）：在 F<=6.2 kN 下吸收 mgh 需要多少行程？")
print("=" * 100)
print(f"  s_req(h) = m·g·h / F_thr = {MASS*G:.1f}·h / {THRESH_N:.0f}")
q2 = {}
for lab, d in COM_DEPTH.items():
    h_an = THRESH_N * d / (MASS * G)
    q2[lab] = {"s_avail_m": d, "h_analytic_m": h_an}
print(f"  {'行程锚':>10} {'腿给(m)':>9} {'解析 h*(m)':>12}")
for lab, d in COM_DEPTH.items():
    print(f"  {lab:>10} {d:>9.3f} {q2[lab]['h_analytic_m']:>12.3f}")
for h in (4.0, 7.0, 9.0):
    print(f"  s_req({h:g} m) = {MASS*G*h/THRESH_N:.3f} m"
          f"  ⇒ 相对腿给 0.43 m 的缺口 {MASS*G*h/THRESH_N-0.43:+.3f} m")
print()

print("=" * 100)
print("Q3 姿态对比（实跑）：tuck-roll vs controlled-drop（同锚 self=0.37 m，D1 力-长链）")
print("=" * 100)
res = {}
for posture, lab in POSTURES:
    rows = []
    for h in GRID:
        r = run(posture, lab, h)
        rows.append({"h_m": h, "f_flex_n": float(r.meta["f_flex_peak_n"]),
                     "f_pad_n": pad_peak(r)})
    res[posture] = rows
hs = [r["h_m"] for r in res[POSTURES[0][0]]]
print(f"  {'姿态':>16} {'':>5}" + "".join(f"{h:>9.1f}" for h in hs) + f"{'h*(m)':>9}")
for posture, _ in POSTURES:
    ff = [r["f_flex_n"] / 1e3 for r in res[posture]]
    hstar = cross_h(hs, [r["f_flex_n"] for r in res[posture]], THRESH_N)
    print(f"  {posture:>16} {'F_flex':>5}" + "".join(f"{x:>9.3f}" for x in ff)
          + f"{hstar:>9.3f}")
    res[posture + "__hstar"] = hstar
print()
for posture, _ in POSTURES:
    fp = [r["f_pad_n"] / 1e3 for r in res[posture]]
    print(f"  {posture:>16} {'F_pad':>5}" + "".join(f"{x:>9.2f}" for x in fp))

ha, hb = res["anchored-self__hstar"], res["tuck-roll__hstar"]
hc = res["controlled-drop__hstar"]
dh = (hb - hc) / hc
print(f"\n  h*: anchored-self {ha:.3f} | controlled-drop {hc:.3f} | tuck-roll {hb:.3f}")
print(f"  => 姿态交换 Δh*/h* = {dh*100:+.1f}%   J2（|Δ|<=10%）："
      f"{'PASS' if abs(dh) <= 0.10 else 'FAIL'}")
fp_d = [r["f_pad_n"] for r in res["controlled-drop"]]
fp_t = [r["f_pad_n"] for r in res["tuck-roll"]]
j3 = all(t > d for t, d in zip(fp_t, fp_d))
print(f"  => F_pad 逐格比较（tuck-roll > controlled-drop）：{'全部成立' if j3 else '不全部成立'}"
      f"   J3：{'PASS' if j3 else 'FAIL'}")
h_an = q2["self"]["h_analytic_m"]
print(f"  => J4：解析 h*({h_an:.3f}) > ODE h*({ha:.3f}) ⇒ {h_an > ha}；"
      f"比值 {h_an/ha:.2f}× (<2.5) ⇒ {'PASS' if h_an > ha and h_an/ha < 2.5 else 'FAIL'}")

out = {"q1_reachability": q1, "q2_energy_ledger": q2, "q3_postures": res,
       "threshold_n": THRESH_N, "com_depth_m": COM_DEPTH,
       "h_ischi_m": list(H_ISCHIUM), "grid_m": GRID}
(ROOT / "results" / "second_contact_probe.json").write_text(
    json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
print("\nOUT = results/second_contact_probe.json  （门已恢复）")
