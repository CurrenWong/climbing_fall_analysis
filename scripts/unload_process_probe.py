"""D1 探针：力-长关系 M(θ)（滞回 h=0 = 纯弹性，明确的保守极限）
O1（主）= 弹性项 / 总力 —— 当前非主导（~1.65-2.0 kN vs 实测 4.2 kN）；D1-J1 要求 >0.50。
锚：M(θ) 由 Sci Rep 16:8582 + PMC8554160（等长 MVIC 膝伸肌，N·m/kg × 80 kg）：
  θ=15°→1.79、45°→2.9、60°→3.2、90°→2.78 N·m/kg
  ⚠️ θ<15° 与 θ>90° 无锚 ⇒ **常数外推**（np.interp 夹持），显式标注。
用法: .venv\\Scripts\\python.exe scripts/unload_process_probe.py
"""
from __future__ import annotations

import json
import pathlib
import sys
from math import degrees

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from climbing.pad import RIGID_CONTACT, TWO_SEG_LEG, simulate_boulder_fall  # noqa: E402

G, LEG0 = 9.80665, 0.85
BW = 80.0 * G
F_THR = 6200.0
H_GRID = (0.50, 1.00, 2.00, 5.00, 9.00)
LABS = ("soft", "self", "stiff")
COM = {"soft": 0.43, "self": 0.37, "stiff": 0.23}
L, K_TH, V_REF, M_ECC, R_REF, RAMP, PAD_KN = 0.425, 350.0, 1.0, 540.0, 0.040, 0.349, 2.1e6
S_HARD = 0.3624                                    # 深蹲膝屈曲 110° 代入 s(θ)
NL_PER_KG = [(15.0, 1.79), (45.0, 2.90), (60.0, 3.20), (90.0, 2.78)]   # 锚点
TABLE = [[t, v * 80.0] for t, v in NL_PER_KG]
TAB_TH, TAB_M = np.radians([t for t, _ in TABLE]), np.array([m for _, m in TABLE])
UP = ROOT / "results" / "fp_anchor_probe.json"
OUT = ROOT / "results" / "unload_process_probe.json"
KEYS = ("k_theta_nm_per_rad", "segment_l_m", "theta_max_rad", "m_ecc_max_nm",
        "ecc_rate_ref_ms", "m_theta_table")


def th_of(s):
    return 2.0 * np.arccos(np.clip(1.0 - np.asarray(s, float) / (2.0 * L), -1.0, 1.0))


def run(lab, h, on_d1, kn=PAD_KN):
    prev = {k: TWO_SEG_LEG.get(k) for k in KEYS}
    c0 = RIGID_CONTACT.get("k_n_per_m")
    try:
        for k in KEYS:
            TWO_SEG_LEG[k] = None
        if on_d1:
            TWO_SEG_LEG.update({"k_theta_nm_per_rad": K_TH, "segment_l_m": L,
                                "theta_max_rad": float(th_of(COM[lab])),
                                "m_ecc_max_nm": M_ECC, "ecc_rate_ref_ms": V_REF,
                                "m_theta_table": TABLE})
        RIGID_CONTACT["k_n_per_m"] = kn
        r = simulate_boulder_fall(height_m=h, mass_kg=80.0, on_pad=False,
                                  posture=f"anchored-{lab}")
    finally:
        for k in KEYS:
            TWO_SEG_LEG[k] = prev[k]
        RIGID_CONTACT["k_n_per_m"] = c0
    t = np.asarray(r.t_s, float)
    zf, zt = np.asarray(r.z_foot_m, float), np.asarray(r.z_torso_m, float)
    t0, t1 = r.meta["impact_window_s"]
    i0 = int(np.searchsorted(t, t0)); i1 = max(int(np.searchsorted(t, t1)), i0 + 2)
    s = np.maximum(0.0, LEG0 - (zf - zt))
    th = th_of(s)
    sh = np.sin(0.5 * np.maximum(th, 1e-9))
    rate = np.maximum(0.0, np.gradient(s, t))
    r_arm = R_REF * np.minimum(1.0, th / RAMP)
    el = np.where(sh > 1e-9, np.interp(th, TAB_TH, TAB_M) / (L * sh), 0.0)
    cap = np.where(sh > 1e-9, (M_ECC * (r_arm / R_REF)) / (L * sh), 0.0)
    rt = cap * rate / (rate + V_REF)
    tot = el + rt
    e_max, t_max = float(el[i0:i1].max()), float(tot[i0:i1].max())
    return {"h_m": h, "f_flex_n": float(r.meta["f_flex_peak_n"]),
            "elastic_peak_n": e_max, "rate_peak_n": float(rt[i0:i1].max()),
            "share_elastic": e_max / max(t_max, 1e-9),
            "s_peak_m": float(s[i0:i1].max()),
            "theta_peak_deg": degrees(float(th[i0:i1].max())),
            "over_s_hard": bool(s[i0:i1].max() > S_HARD)}


def hstar(rows):
    for i in range(len(rows) - 1):
        a, b = rows[i], rows[i + 1]
        if a["over_s_hard"] or b["over_s_hard"]:
            continue
        if a["f_flex_n"] <= F_THR <= b["f_flex_n"]:
            return a["h_m"] + (F_THR - a["f_flex_n"]) / (b["f_flex_n"] - a["f_flex_n"]) * \
                (b["h_m"] - a["h_m"])
    return None


print("=" * 100)
print("D1-J0 POSITIVE CONTROL（门关 + 刚性口径 ⇒ 逐位复现 fp_anchor_probe.json）")
print("=" * 100)
up = json.loads(UP.read_text(encoding="utf-8"))
ok = True
for lab in LABS:
    w = up["main_rows"][lab]
    g = run(lab, 0.42, on_d1=False, kn=None)
    good = abs(g["f_flex_n"] - w["f_flex_n"]) < 1e-6
    ok &= good
    print(f"  {lab:<6} want {w['f_flex_n']:.6f} | got {g['f_flex_n']:.6f}  "
          f"{'PASS' if good else 'FAIL'}")
if not ok:
    print("!! INVALID")
    OUT.write_text(json.dumps({"positive_control": False}), encoding="utf-8")
    raise SystemExit(1)
print("=> POSITIVE CONTROL PASS")
print(f"   M(θ) 表（N·m）: " + " ".join(f"{t:.0f}°={m:.1f}" for t, m in TABLE)
      + "  ⚠️ θ<15°/θ>90° 常数外推（无锚）")
print(f"   s_hard = {S_HARD} m（深蹲膝屈曲 110°）")

print()
print("=" * 100)
print("O1-O4  D1 ｜ 实测跟垫口径 ｜ **O1 = 弹性项/总力**（D1-J1 要求 >0.50）")
print("=" * 100)
print(f"  {'档':<6} " + " ".join(f"{h:>7.2f}m" for h in H_GRID)
      + f" {'h*(m)':>9} {'单调':>5} {'s@峰':>7} {'θ@峰°':>7} {'弹性kN':>7} {'份额':>6}")
grid, verdict = {}, {}
for lab in LABS:
    rows = [run(lab, h, on_d1=True) for h in H_GRID]
    grid[lab] = rows
    hs = hstar(rows)
    mono = all(rows[i]["f_flex_n"] < rows[i + 1]["f_flex_n"] for i in range(len(rows) - 1))
    verdict[lab] = {"h_star_m": hs, "monotone": mono,
                    "share_elastic_at_hstar": max(
                        (r["share_elastic"] for r in rows if hs and r["h_m"] <= hs), default=0.0),
                    "share_elastic_max": max(r["share_elastic"] for r in rows),
                    "over_s_hard_at": [r["h_m"] for r in rows if r["over_s_hard"]]}
    print(f"  {lab:<6} " + " ".join(f"{r['f_flex_n']/1e3:>7.2f} " for r in rows)
          + f" {('%.3f' % hs) if hs else '不可定':>9} {('是' if mono else '否'):>5} "
          f" {rows[-1]['s_peak_m']:>7.4f} {rows[-1]['theta_peak_deg']:>7.1f} "
          f" {rows[-1]['elastic_peak_n']/1e3:>7.2f} {rows[-1]['share_elastic']:>6.3f}")

print()
print("=" * 100)
print("判定（预注册 §5）")
print("=" * 100)
j1 = {l: verdict[l]["share_elastic_at_hstar"] > 0.50 for l in LABS}
j2 = {l: verdict[l]["monotone"] for l in LABS}
j3 = {l: (verdict[l]["h_star_m"] is not None
          and verdict[l]["share_elastic_at_hstar"] > 0) for l in LABS}
for l in LABS:
    v = verdict[l]
    print(f"  {l:<6} J1 弹性份额@h*={v['share_elastic_at_hstar']:.3f}>{0.50}→{j1[l]} "
          f"（全网格最大 {v['share_elastic_max']:.3f}）  J2 单调={v['monotone']}  "
          f"J3 h*={('%.3f' % v['h_star_m']) if v['h_star_m'] else '不可定'}")
print(f"  => D1-J1={'PASS' if all(j1.values()) else 'FAIL'}  "
      f"D1-J2={'PASS' if all(j2.values()) else 'FAIL'}  "
      f"D1-J3={'PASS' if all(j3.values()) else 'FAIL'}")
print("  ⛔ J1 若 FAIL ⇒ 如实记「力-长不足以让几何项主导 ⇒ 缺主动控制/滞回」，不得回调 M(θ)")

OUT.write_text(json.dumps({"positive_control": True, "m_theta_table_nm": TABLE,
                           "note_extrap": "θ<15° 与 θ>90° 无锚 ⇒ 常数外推",
                           "s_hard_m": S_HARD, "m_ecc_nm": M_ECC, "pad_kn": PAD_KN,
                           "hysteresis": 0.0, "grid": grid, "verdict": verdict,
                           "D1_J1": j1, "D1_J2": j2, "D1_J3": j3},
                          ensure_ascii=False, indent=2), encoding="utf-8")
print(f"\nOUT = {OUT}")
