"""C2 探针（被动路径 / 有效上限）—— `docs/methods/P3_骨性软组织接触_预注册.md`（§2bis 更正后）

判据（更正版）：
  C2-J0 门关 ⇒ 逐位复现 fp_anchor_probe.json
  C2-J1 实测跟垫口径 F_flex(h) 单调
  C2-J2 ★ **h* 所在高度的 s@峰 < s_hard** ⇒ h* 不依赖自由 k_stop/s_hard（可辩护）
  C2-J3 三档皆可定 h*
配置 = C1（B2a 两段腿 + 饱和离心律，M_ecc = 540 N·m 双腿口径）
s_hard = 0.49 / 0.56 m（人体测量学：深蹲 COM 下降 0.28–0.32·H，H=1.75）

用法: .venv\\Scripts\\python.exe scripts/passive_path_probe.py
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
L, K_TH, V_REF, M_ECC = 0.425, 350.0, 1.0, 540.0      # 双腿口径（C1 §7.3 建议）
PAD_KN = 2.1e6
S_HARD = (0.49, 0.56)
UP = ROOT / "results" / "fp_anchor_probe.json"
OUT = ROOT / "results" / "passive_path_probe.json"
KEYS = ("k_theta_nm_per_rad", "segment_l_m", "theta_max_rad", "m_ecc_max_nm",
        "ecc_rate_ref_ms")


def th_of(s):
    s = np.asarray(s, dtype=float)
    return 2.0 * np.arccos(np.clip(1.0 - s / (2.0 * L), -1.0, 1.0))


def run(lab, h, on_c2, kn=PAD_KN):
    """`kn=None` ⇒ 刚性口径（= `LEG_K_ANCHOR`，正对照必须用这个）。"""
    prev = {k: TWO_SEG_LEG.get(k) for k in KEYS}
    c0 = RIGID_CONTACT.get("k_n_per_m")
    try:
        for k in KEYS:
            TWO_SEG_LEG[k] = None
        if on_c2:
            TWO_SEG_LEG.update({"k_theta_nm_per_rad": K_TH, "segment_l_m": L,
                                "theta_max_rad": float(th_of(COM[lab])),
                                "m_ecc_max_nm": M_ECC, "ecc_rate_ref_ms": V_REF})
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
    s_win = float(s[i0:i1].max())
    th_max = float(th_of(COM[lab]))
    return {"h_m": h, "f_flex_n": float(r.meta["f_flex_peak_n"]),
            "s_peak_m": s_win, "theta_peak_deg": degrees(float(th_of(s_win))),
            "theta_max_deg": degrees(th_max),
            "exhausted": bool(th_of(s_win) > th_max),
            "s_over_hard": {str(x): bool(s_win > x) for x in S_HARD}}


def hstar(rows):
    for i in range(len(rows) - 1):
        a, b = rows[i], rows[i + 1]
        if a["exhausted"] or b["exhausted"]:
            continue
        if a["f_flex_n"] <= F_THR <= b["f_flex_n"]:
            return a["h_m"] + (F_THR - a["f_flex_n"]) / (b["f_flex_n"] - a["f_flex_n"]) * \
                (b["h_m"] - a["h_m"])
    return None


print("=" * 100)
print("C2-J0 POSITIVE CONTROL（门关 ⇒ 逐位复现 fp_anchor_probe.json）")
print("=" * 100)
up = json.loads(UP.read_text(encoding="utf-8"))
ok = True
for lab in LABS:
    w = up["main_rows"][lab]
    g = run(lab, 0.42, on_c2=False, kn=None)   # ★ 正对照必须用刚性口径（k_n=None）
    good = abs(g["f_flex_n"] - w["f_flex_n"]) < 1e-6
    ok &= good
    print(f"  {lab:<6} want {w['f_flex_n']:.6f} | got {g['f_flex_n']:.6f}  "
          f"{'PASS' if good else 'FAIL'}")
if not ok:
    print("!! INVALID"); OUT.write_text(json.dumps({"positive_control": False}), encoding="utf-8")
    raise SystemExit(1)
print("=> POSITIVE CONTROL PASS")
print(f"   配置 = C1（B2a + 饱和离心律，M_ecc = {M_ECC} N·m 双腿）｜s_hard = {S_HARD} m")

print()
print("=" * 100)
print("O1-O4  C2 ｜ 实测跟垫口径")
print("=" * 100)
print(f"  {'档':<6} " + " ".join(f"{h:>7.2f}m" for h in H_GRID)
      + f" {'h*(m)':>9} {'单调':>5} {'s@峰(m)':>8} {'θ@峰°':>7} {'θ_max°':>7} "
      f"{'耗尽':>5} {'>0.49':>6} {'>0.56':>6}")
grid, verdict = {}, {}
for lab in LABS:
    rows = [run(lab, h, on_c2=True) for h in H_GRID]
    grid[lab] = rows
    hs = hstar(rows)
    mono = all(rows[i]["f_flex_n"] < rows[i + 1]["f_flex_n"] for i in range(len(rows) - 1))
    # h* 落在哪一对格之间 ⇒ 用该段的 s@峰 上界判 C2-J2
    s_at_hstar = max(r["s_peak_m"] for r in rows if hs is not None and r["h_m"] <= max(hs, 1e-9))
    verdict[lab] = {"h_star_m": hs, "monotone": mono,
                    "s_peak_at_or_below_hstar_m": s_at_hstar,
                    "s_peak_max_m": max(r["s_peak_m"] for r in rows),
                    "exhausted_at": [r["h_m"] for r in rows if r["exhausted"]]}
    print(f"  {lab:<6} " + " ".join(f"{r['f_flex_n']/1e3:>7.2f} " for r in rows)
          + f" {('%.3f' % hs) if hs else '不可定':>9} {('是' if mono else '否'):>5} "
          f" {rows[-1]['s_peak_m']:>8.4f} {rows[-1]['theta_peak_deg']:>7.1f} "
          f" {rows[-1]['theta_max_deg']:>7.1f} "
          f" {('是' if rows[-1]['exhausted'] else '否'):>5} "
          f" {('是' if rows[-1]['s_over_hard']['0.49'] else '否'):>6} "
          f" {('是' if rows[-1]['s_over_hard']['0.56'] else '否'):>6}")

print()
print("=" * 100)
print("判定（预注册 §5，更正版 C2-J2）")
print("=" * 100)
j1 = {l: verdict[l]["monotone"] for l in LABS}
j2 = {l: (verdict[l]["s_peak_at_or_below_hstar_m"] < S_HARD[0]) for l in LABS}
j3 = {l: verdict[l]["h_star_m"] is not None for l in LABS}
for l in LABS:
    v = verdict[l]
    print(f"  {l:<6} J1单调={v['monotone']}  J2 s@(≤h*)={v['s_peak_at_or_below_hstar_m']:.4f}"
          f"<{S_HARD[0]}→{j2[l]}  J3 h*={('%.3f' % v['h_star_m']) if v['h_star_m'] else '不可定'}")
# J2 需至少一档到 h* 之上（否则判据无检验力），故同时报「h* 之上是否触发 s_hard」
print(f"  => C2-J1={'PASS' if all(j1.values()) else 'FAIL'}  "
      f"C2-J2（h* 附近 s<s_hard）={'PASS' if all(j2.values()) else 'FAIL'}  "
      f"C2-J3={'PASS' if all(j3.values()) else 'FAIL'}")
print("  ⚠️ J2 需注意检验力：若**所有**格都 s<s_hard，则该判据在实测范围内无检验力（须如实标注）")
print(f"     各档 s@峰 最大值："
      + "  ".join(f"{l}={verdict[l]['s_peak_max_m']:.4f}" for l in LABS)
      + f"  vs s_hard={S_HARD}")

OUT.write_text(json.dumps({"positive_control": True, "m_ecc_nm": M_ECC, "s_hard_m": list(S_HARD),
                           "config": "C1 (B2a + saturated eccentric)", "pad_kn": PAD_KN,
                           "h_grid_m": list(H_GRID), "grid": grid, "verdict": verdict,
                           "C2_J1": j1, "C2_J2": j2, "C2_J3": j3},
                          ensure_ascii=False, indent=2), encoding="utf-8")
print(f"\nOUT = {OUT}")
