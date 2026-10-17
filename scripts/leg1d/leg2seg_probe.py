"""B2a 探针（两段腿等效单关节）—— `docs/methods/P3_两段腿等效关节_预注册.md`

设计定理：两段腿 ≡ 非线性（硬化）弹簧
    θ(s) = 2·arccos(1 − s/(2l))    F(s) = k_θ·θ / (l·sin(θ/2))
⇒ 关节准静态时与单自由度框架等价（无需新增自由度）。

正对照（B2a-J0，不过则整份作废）：门关时必须逐位复现 `results/fp_anchor_probe.json`。

判据（§5，跑前冻结）
    J1 实测跟垫口径下 `F_flex(h)` 在未耗尽的格上**严格单调递增**（线性模型在此非单调）
    J2 实测跟垫口径下**三档都能在未耗尽范围内定出 `h*`**（B1 下只有 soft 可定）
    J3 @0.42 m `F_flex`/BW ÷ 实测靶 ∈ [0.70,1.30]（**只报不判**）
    J5 ⛔ 禁调 `k_θ` 使 J1/J2 变绿

用法::

    .venv\\Scripts\\python.exe scripts/leg1d/leg2seg_probe.py
"""
from __future__ import annotations

import json
import pathlib
import sys
from math import acos, degrees, sin

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from climbing.pad import POSTURES, RIGID_CONTACT, TWO_SEG_LEG, simulate_boulder_fall  # noqa: E402

G, LEG_EXT0 = 9.80665, 0.85
MASS_KG = 80.0
BW = MASS_KG * G
F_THR_N = 6200.0
#: 跑前网格选择：**弃用 0.20 m** —— 该格在实测跟垫口径下每格需 **~186 s**（软化非线性律导致
#: 长时振铃、`settled` 不触发 ⇒ 积分走满 `t_max`），三档合计 ~9.5 min，性价比低且不在 `h*` 附近。
H_GRID = (0.50, 1.00, 2.00, 5.00, 7.00, 9.00)
LABELS = ("soft", "self", "stiff")
COM_DEPTH_M = {"soft": 0.43, "self": 0.37, "stiff": 0.23}
L_SEG = 0.425                      # 2l = 0.85 = leg_ext0
K_THETA = 350.0                    # N·m/rad（0.08 N·m/kg/deg × 76 kg；Shamaei 2013 ~400）
FP_TARGET_BW = {"soft": 3.38, "self": 3.70, "stiff": 5.36}
BAND = (0.70, 1.30)
UPSTREAM = ROOT / "results" / "fp_anchor_probe.json"
OUT = ROOT / "results" / "leg2seg_probe.json"


def theta_max_of(dl: float, l: float = L_SEG) -> float:
    return 2.0 * acos(1.0 - dl / (2.0 * l))


def theta_of(s: float, l: float = L_SEG) -> float:
    return 2.0 * acos(min(1.0, max(-1.0, 1.0 - s / (2.0 * l))))


def _run(label: str, h: float, kn, use_b2a: bool) -> dict:
    keys = ("k_theta_nm_per_rad", "segment_l_m", "theta_max_rad")
    prev = {k: TWO_SEG_LEG.get(k) for k in keys}
    c0 = RIGID_CONTACT.get("k_n_per_m")
    try:
        if use_b2a:
            TWO_SEG_LEG["k_theta_nm_per_rad"] = K_THETA
            TWO_SEG_LEG["segment_l_m"] = L_SEG
            TWO_SEG_LEG["theta_max_rad"] = theta_max_of(COM_DEPTH_M[label])
        else:
            for k in keys:
                TWO_SEG_LEG[k] = None
        RIGID_CONTACT["k_n_per_m"] = kn
        r = simulate_boulder_fall(height_m=h, mass_kg=MASS_KG, on_pad=False,
                                  posture=f"anchored-{label}")
    finally:
        for k in keys:
            TWO_SEG_LEG[k] = prev[k]
        RIGID_CONTACT["k_n_per_m"] = c0
    t = np.asarray(r.t_s, float)
    zf, zt = np.asarray(r.z_foot_m, float), np.asarray(r.z_torso_m, float)
    t0, t1 = r.meta["impact_window_s"]
    i0 = int(np.searchsorted(t, t0))
    i1 = max(int(np.searchsorted(t, t1)), i0 + 2)
    s = np.maximum(0.0, LEG_EXT0 - (zf - zt))
    th = 2.0 * np.arccos(np.clip(1.0 - s / (2.0 * L_SEG), -1.0, 1.0))
    th_win = float(th[i0:i1].max())
    th_max = theta_max_of(COM_DEPTH_M[label]) if use_b2a else None
    ff = float(r.meta["f_flex_peak_n"])
    return {
        "label": label, "h_m": h, "k_n": kn, "use_b2a": use_b2a,
        "f_pad_n": float(r.peak_force_n), "f_flex_n": ff,
        "f_pad_bw": float(r.peak_force_n) / BW, "f_flex_bw": ff / BW,
        "theta_peak_rad": th_win, "theta_peak_deg": degrees(th_win),
        "theta_max_deg": (degrees(th_max) if th_max else None),
        "exhausted": bool(th_max is not None and th_win > th_max),
    }


def _hstar(rows: list[dict]) -> float | None:
    """只在**相邻且两端都未耗尽**的格之间内插（B1 纪律，`PAD_GOTCHAS` §6d）。"""
    for i in range(len(rows) - 1):
        a, b = rows[i], rows[i + 1]
        if a["exhausted"] or b["exhausted"]:
            continue
        if a["f_flex_n"] <= F_THR_N <= b["f_flex_n"]:
            h0, h1 = a["h_m"], b["h_m"]
            f0, f1 = a["f_flex_n"], b["f_flex_n"]
            return h0 + (F_THR_N - f0) / (f1 - f0) * (h1 - h0)
    return None


def positive_control() -> tuple[bool, list[dict]]:
    up = json.loads(UPSTREAM.read_text(encoding="utf-8"))
    rows, ok = [], True
    for lab in LABELS:
        w = up["main_rows"][lab]
        g = _run(lab, 0.42, None, use_b2a=False)
        good = (abs(g["f_pad_n"] - w["f_pad_n"]) < 1e-6
                and abs(g["f_flex_n"] - w["f_flex_n"]) < 1e-6)
        ok &= good
        rows.append({"label": lab, "pass": good,
                     "want": {"f_pad_n": w["f_pad_n"], "f_flex_n": w["f_flex_n"]},
                     "got": {"f_pad_n": g["f_pad_n"], "f_flex_n": g["f_flex_n"]}})
    return ok, rows


def main() -> int:
    print("=" * 94)
    print("B2a-J0  POSITIVE CONTROL（门关 ⇒ 必须逐位复现 fp_anchor_probe.json）")
    print("=" * 94)
    pcok, pcrows = positive_control()
    for r in pcrows:
        print(f"  {r['label']:<6} want {r['want']['f_pad_n']:.6f}/{r['want']['f_flex_n']:.6f}"
              f" | got {r['got']['f_pad_n']:.6f}/{r['got']['f_flex_n']:.6f}"
              f"  {'PASS' if r['pass'] else 'FAIL'}")
    if not pcok:
        print("!! POSITIVE CONTROL FAILED -- INVALID")
        OUT.write_text(json.dumps({"positive_control": False}, ensure_ascii=False, indent=2),
                       encoding="utf-8")
        return 1
    print("=> POSITIVE CONTROL PASS")
    print(f"   l={L_SEG} m（2l=0.85=leg_ext0）｜k_θ={K_THETA} N·m/rad｜"
          f"θ_max = " + " / ".join(f"{l} {degrees(theta_max_of(COM_DEPTH_M[l])):.1f}°"
                                   for l in LABELS))

    grid, verdict = {}, {}
    for tag, kn in (("rigid", None), ("pad", 2.1e6)):
        print()
        print("=" * 94)
        print(f"B2a 两段腿 ｜ 接触口径 = {tag}")
        print("=" * 94)
        print(f"  {'档':<6} " + " ".join(f"{h:>7.2f}m" for h in H_GRID)
              + f" {'h*(m)':>9} {'θ@峰°':>8} {'θ_max°':>8} {'耗尽':>6}")
        for lab in LABELS:
            rows = []
            for h in H_GRID:
                print(f"    [run] {tag}/{lab}/h={h:.2f}", flush=True)
                rows.append(_run(lab, h, kn, use_b2a=True))
            grid[f"{tag}:{lab}"] = rows
            hs = _hstar(rows)
            mono = all(rows[i]["f_flex_n"] < rows[i + 1]["f_flex_n"]
                       for i in range(len(rows) - 1))
            first_ex = next((r["h_m"] for r in rows if r["exhausted"]), None)
            verdict.setdefault(tag, {})[lab] = {
                "h_star_m": hs, "monotone": mono, "first_exhausted_h_m": first_ex,
                "exhausted_at": [r["h_m"] for r in rows if r["exhausted"]],
            }
            print(f"  {lab:<6} " + " ".join(f"{r['f_flex_n']/1e3:>7.2f} " for r in rows)
                  + f" {('%.3f' % hs) if hs else '不可定':>9} "
                  + f" {rows[-1]['theta_peak_deg']:>8.1f} {rows[-1]['theta_max_deg']:>8.1f} "
                  + f" {('@%.2f' % first_ex) if first_ex else '从未':>6}"
                  + ("" if mono else "  ⚠非单调"))

    print()
    print("=" * 94)
    print("判定（预注册 §5，如实二值）")
    print("=" * 94)
    j1 = {l: verdict["pad"][l]["monotone"] for l in LABELS}
    for l in LABELS:
        print(f"  B2a-J1 pad/{l:<6} 单调 = {verdict['pad'][l]['monotone']}  "
              f"{'PASS' if j1[l] else 'FAIL'}")
    j2 = {l: verdict["pad"][l]["h_star_m"] is not None for l in LABELS}
    for l in LABELS:
        v = verdict["pad"][l]
        print(f"  B2a-J2 pad/{l:<6} h* = "
              f"{('%.3f m' % v['h_star_m']) if v['h_star_m'] else '不可定'}  "
              f"{'PASS' if j2[l] else 'FAIL'}")
    print(f"    B2a-J1(三档) = {'PASS' if all(j1.values()) else 'FAIL'} ｜ "
          f"B2a-J2(三档) = {'PASS' if all(j2.values()) else 'FAIL'}")
    print("  B2a-J3 @0.42 m（只报不判）：")
    up = json.loads(UPSTREAM.read_text(encoding="utf-8"))
    for l in LABELS:
        d = _run(l, 0.42, None, use_b2a=True)
        ratio = d["f_flex_bw"] / FP_TARGET_BW[l]
        print(f"      {l:<6} linear {up['main_rows'][l]['f_flex_bw']:.2f} BW -> "
              f"2seg {d['f_flex_bw']:.2f} BW  ratio {ratio:.3f} "
              f"(band {BAND[0]}-{BAND[1]})")

    payload = {
        "positive_control": True, "positive_control_rows": pcrows,
        "l_seg_m": L_SEG, "k_theta_nm_per_rad": K_THETA,
        "com_depth_m": COM_DEPTH_M,
        "theta_max_deg": {l: degrees(theta_max_of(COM_DEPTH_M[l])) for l in LABELS},
        "mass_kg": MASS_KG, "f_threshold_n": F_THR_N, "h_grid_m": list(H_GRID),
        "fp_target_bw": FP_TARGET_BW, "band": list(BAND),
        "grid": grid, "verdict": verdict,
        "B2a_J1": j1, "B2a_J1_all": all(j1.values()),
        "B2a_J2": j2, "B2a_J2_all": all(j2.values()),
    }
    OUT.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print()
    print(f"OUT = {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
