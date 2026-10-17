"""D1 探针（修正版）：① M(θ) 用 **PCHIP 单调三次拟合**锚点 ⇒ 密集表 ⇒ 导数连续（解刚）；
② **多进程并行**跑；③ 网格缩到 {0.5,1,2,5}（h* 落在 1~5 之间，够用）。
滞回 h = 0（纯弹性）。锚：Sci Rep 16:8582 + PMC8554160（等长 MVIC 膝伸肌，N·m/kg × 80 kg）。
⚠️ 0° 与 >90° 无锚 ⇒ 边缘取**带端点常数**（KEP 外推），显式标注。
用法: .venv\\Scripts\\python.exe scripts/leg1d/unload_process_probe.py
"""
from __future__ import annotations

import json
import os
import pathlib
import sys
from concurrent.futures import ProcessPoolExecutor
from math import degrees

import numpy as np
from scipy.interpolate import PchipInterpolator

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

G, LEG0 = 9.80665, 0.85
BW = 80.0 * G
F_THR = 6200.0
H_GRID = (0.50, 1.00, 2.00, 5.00)
LABS = ("soft", "self", "stiff")
COM = {"soft": 0.43, "self": 0.37, "stiff": 0.23}
L, K_TH, V_REF, M_ECC, R_REF, RAMP, PAD_KN = 0.425, 350.0, 1.0, 540.0, 0.040, 0.349, 2.1e6
S_HARD = 0.3624
#: 锚点（θ_deg, N·m/kg）；0° 取「近全伸展」带下端 1.3，135° 取「屈曲位」2.78（两端【无锚·常数外推】）
ANCH_NL_PER_KG = [(0.0, 1.30), (15.0, 1.79), (45.0, 2.90), (60.0, 3.20),
                  (90.0, 2.78), (135.0, 2.78)]
#: PCHIP ⇒ 密集表（每 2.5°）—— 导数连续，解 ODE 的刚性
_D = np.arange(0.0, 135.0 + 1e-9, 2.5)
_TAB = [[float(t), float(PchipInterpolator([a for a, _ in ANCH_NL_PER_KG],
                                           [b * 80.0 for _, b in ANCH_NL_PER_KG])(t))]
        for t in _D]
TAB_TH, TAB_M = np.radians([t for t, _ in _TAB]), np.array([m for _, m in _TAB])
UP = ROOT / "results" / "fp_anchor_probe.json"
OUT = ROOT / "results" / "unload_process_probe.json"
KEYS = ("k_theta_nm_per_rad", "segment_l_m", "theta_max_rad", "m_ecc_max_nm",
        "ecc_rate_ref_ms", "m_theta_table")


def _th_of(s):
    return 2.0 * np.arccos(np.clip(1.0 - np.asarray(s, float) / (2.0 * L), -1.0, 1.0))


def _one(job) -> dict:
    """一个格子的计算（供并行）。"""
    from climbing.pad import RIGID_CONTACT, TWO_SEG_LEG, simulate_boulder_fall
    lab, h, on_d1, kn = job
    prev = {k: TWO_SEG_LEG.get(k) for k in KEYS}
    c0 = RIGID_CONTACT.get("k_n_per_m")
    try:
        for k in KEYS:
            TWO_SEG_LEG[k] = None
        if on_d1:
            TWO_SEG_LEG.update({"k_theta_nm_per_rad": K_TH, "segment_l_m": L,
                                "theta_max_rad": float(_th_of(COM[lab])),
                                "m_ecc_max_nm": M_ECC, "ecc_rate_ref_ms": V_REF,
                                "m_theta_table": _TAB})
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
    th = _th_of(s)
    sh = np.sin(0.5 * np.maximum(th, 1e-9))
    rate = np.maximum(0.0, np.gradient(s, t))
    r_arm = R_REF * np.minimum(1.0, th / RAMP)
    # ★ 与 pad.py 对齐：弹性项也要乘 r(θ) 归零因子（否则 θ→0 发散；见 commit a319b36）
    el = np.where(sh > 1e-9,
                  (np.interp(th, TAB_TH, TAB_M) * np.minimum(1.0, th / RAMP)) / (L * sh),
                  0.0)
    cap = np.where(sh > 1e-9, (M_ECC * (r_arm / R_REF)) / (L * sh), 0.0)
    rt = cap * rate / (rate + V_REF)
    tot = el + rt
    em, tm = float(el[i0:i1].max()), float(tot[i0:i1].max())
    return {"label": lab, "h_m": h, "on_d1": on_d1,
            "f_flex_n": float(r.meta["f_flex_peak_n"]),
            "elastic_peak_n": em, "rate_peak_n": float(rt[i0:i1].max()),
            "share_elastic": em / max(tm, 1e-9),
            "s_peak_m": float(s[i0:i1].max()),
            "theta_peak_deg": degrees(float(th[i0:i1].max())),
            "over_s_hard": bool(s[i0:i1].max() > S_HARD)}


def _hstar(rows):
    rs = sorted(rows, key=lambda r: r["h_m"])
    for i in range(len(rs) - 1):
        a, b = rs[i], rs[i + 1]
        if a["over_s_hard"] or b["over_s_hard"]:
            continue
        if a["f_flex_n"] <= F_THR <= b["f_flex_n"]:
            return a["h_m"] + (F_THR - a["f_flex_n"]) / (b["f_flex_n"] - a["f_flex_n"]) * \
                (b["h_m"] - a["h_m"])
    return None


def main() -> int:
    nw = max(1, min(8, (os.cpu_count() or 4)))
    jobs = ([(l, 0.42, False, None) for l in LABS]
            + [(l, h, True, PAD_KN) for l in LABS for h in H_GRID])
    print(f"并行跑 {len(jobs)} 格（{nw} 进程）…", flush=True)
    with ProcessPoolExecutor(max_workers=nw) as ex:
        res = list(ex.map(_one, jobs))
    print("跑完。", flush=True)

    print("=" * 96)
    print("D1-J0 POSITIVE CONTROL（门关 + 刚性口径 ⇒ 逐位复现 fp_anchor_probe.json）")
    print("=" * 96)
    up = json.loads(UP.read_text(encoding="utf-8"))
    ok = True
    for r in res:
        if r["on_d1"]:
            continue
        w = up["main_rows"][r["label"]]
        good = abs(r["f_flex_n"] - w["f_flex_n"]) < 1e-6
        ok &= good
        print(f"  {r['label']:<6} want {w['f_flex_n']:.6f} | got {r['f_flex_n']:.6f}  "
              f"{'PASS' if good else 'FAIL'}")
    if not ok:
        print("!! INVALID")
        OUT.write_text(json.dumps({"positive_control": False}), encoding="utf-8")
        return 1
    print("=> POSITIVE CONTROL PASS")
    print(f"   M(θ) = PCHIP 拟合锚点 ⇒ 密集表（2.5° 步，0–135°）；"
          f"两端【无锚·常数外推】｜s_hard = {S_HARD}")

    print()
    print("=" * 96)
    print("O1-O4  D1 ｜ 实测跟垫口径 ｜ **O1 = 弹性项/总力**（D1-J1 要求 >0.50）")
    print("=" * 96)
    print(f"  {'档':<6} " + " ".join(f"{h:>7.2f}m" for h in H_GRID)
          + f" {'h*(m)':>9} {'单调':>5} {'s@峰':>7} {'θ@峰°':>7} {'弹性kN':>7} {'份额':>6}")
    grid, verdict = {}, {}
    for lab in LABS:
        rows = [r for r in res if r["on_d1"] and r["label"] == lab]
        grid[lab] = rows
        hs = _hstar(rows)
        mono = all(rows[i]["f_flex_n"] < rows[i + 1]["f_flex_n"] for i in range(len(rows) - 1))
        sh_at = max((r["share_elastic"] for r in rows if hs and r["h_m"] <= hs), default=0.0)
        verdict[lab] = {"h_star_m": hs, "monotone": mono, "share_elastic_at_hstar": sh_at,
                        "share_elastic_max": max(r["share_elastic"] for r in rows),
                        "over_s_hard_at": [r["h_m"] for r in rows if r["over_s_hard"]]}
        last = max(rows, key=lambda r: r["h_m"])
        print(f"  {lab:<6} " + " ".join(f"{r['f_flex_n']/1e3:>7.2f} " for r in rows)
              + f" {('%.3f' % hs) if hs else '不可定':>9} {('是' if mono else '否'):>5} "
              f" {last['s_peak_m']:>7.4f} {last['theta_peak_deg']:>7.1f} "
              f" {last['elastic_peak_n']/1e3:>7.2f} {last['share_elastic']:>6.3f}")

    print()
    print("=" * 96)
    print("判定（预注册 §5）")
    print("=" * 96)
    j1 = {l: verdict[l]["share_elastic_at_hstar"] > 0.50 for l in LABS}
    j2 = {l: verdict[l]["monotone"] for l in LABS}
    j3 = {l: (verdict[l]["h_star_m"] is not None) for l in LABS}
    for l in LABS:
        v = verdict[l]
        print(f"  {l:<6} J1 弹性份额@h*={v['share_elastic_at_hstar']:.3f}（全网格max "
              f"{v['share_elastic_max']:.3f}）→{j1[l]}  J2 单调={v['monotone']}  "
              f"J3 h*={('%.3f' % v['h_star_m']) if v['h_star_m'] else '不可定'}")
    print(f"  => D1-J1={'PASS' if all(j1.values()) else 'FAIL'}  "
          f"D1-J2={'PASS' if all(j2.values()) else 'FAIL'}  "
          f"D1-J3={'PASS' if all(j3.values()) else 'FAIL'}")
    print("  ⛔ J1 若 FAIL ⇒ 如实记「力-长不足以让几何项主导 ⇒ 缺主动控制/滞回」，不得回调 M(θ)")

    OUT.write_text(json.dumps({"positive_control": True,
                               "m_theta_table_nm": _TAB, "method": "PCHIP dense",
                               "extrap_note": "θ<0° 与 >135° 无锚·常数外推；0°取1.3、135°取2.78 N·m/kg",
                               "s_hard_m": S_HARD, "m_ecc_nm": M_ECC, "pad_kn": PAD_KN,
                               "hysteresis": 0.0, "h_grid_m": list(H_GRID),
                               "grid": grid, "verdict": verdict,
                               "D1_J1": j1, "D1_J2": j2, "D1_J3": j3},
                              ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\nOUT = {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
