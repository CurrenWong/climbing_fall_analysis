"""C1 探针（肌肉离心吸能）—— `docs/methods/P3_肌肉离心吸能_预注册.md`
把无上限粘性 `c_flex·rate` 换成 Hill 型饱和律 `f = f_ecc_max·rate/(rate+v_ref)`，
`f_ecc_max = M_ecc/(l·sin(θ/2))`（M_ecc = 270 N·m，Ruas 实测）。

正对照 C1-J0：门关 ⇒ 逐位复现 fp_anchor_probe.json。
用法: .venv\\Scripts\\python.exe scripts/ecc_muscle_probe.py
"""
from __future__ import annotations

import json
import pathlib
import sys
from math import acos, degrees, sin

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from climbing.pad import RIGID_CONTACT, TWO_SEG_LEG, simulate_boulder_fall  # noqa: E402

G, LEG_EXT0, C_FLEX = 9.80665, 0.85, 400.0
MASS, BW = 80.0, 80.0 * 9.80665
F_THR = 6200.0
H_GRID = (0.50, 1.00, 2.00, 5.00, 7.00, 9.00)
LABS = ("soft", "self", "stiff")
COM = {"soft": 0.43, "self": 0.37, "stiff": 0.23}
L, K_TH, M_ECC, V_REF = 0.425, 350.0, 270.0, 1.0
PAD_KN = 2.1e6
UP = ROOT / "results" / "fp_anchor_probe.json"
OUT = ROOT / "results" / "ecc_muscle_probe.json"
KEYS = ("k_theta_nm_per_rad", "segment_l_m", "theta_max_rad", "m_ecc_max_nm",
        "ecc_rate_ref_ms")


def th_of(s):
    """向量化版本（`s` 可为标量或数组）。"""
    s = np.asarray(s, dtype=float)
    return 2.0 * np.arccos(np.clip(1.0 - s / (2.0 * L), -1.0, 1.0))


def run(lab, h, kn, on_c1, v_ref=V_REF):
    prev = {k: TWO_SEG_LEG.get(k) for k in KEYS}
    c0 = RIGID_CONTACT.get("k_n_per_m")
    try:
        for k in KEYS:
            TWO_SEG_LEG[k] = None
        if on_c1:
            TWO_SEG_LEG.update({"k_theta_nm_per_rad": K_TH, "segment_l_m": L,
                                "theta_max_rad": th_of(COM[lab]),
                                "m_ecc_max_nm": M_ECC, "ecc_rate_ref_ms": v_ref})
        RIGID_CONTACT["k_n_per_m"] = kn
        r = simulate_boulder_fall(height_m=h, mass_kg=MASS, on_pad=False,
                                  posture=f"anchored-{lab}")
    finally:
        for k in KEYS:
            TWO_SEG_LEG[k] = prev[k]
        RIGID_CONTACT["k_n_per_m"] = c0
    t = np.asarray(r.t_s, float)
    zf, zt = np.asarray(r.z_foot_m, float), np.asarray(r.z_torso_m, float)
    t0, t1 = r.meta["impact_window_s"]
    i0 = int(np.searchsorted(t, t0))
    i1 = max(int(np.searchsorted(t, t1)), i0 + 2)
    s = np.maximum(0.0, LEG_EXT0 - (zf - zt))
    th = th_of(s)
    th_w = float(th[i0:i1].max())
    rate = np.maximum(0.0, np.gradient(s, t))
    f_cap = M_ECC / (L * np.sin(0.5 * np.maximum(th, 1e-9)))
    f_mus = f_cap * rate / (rate + v_ref)
    j = i0 + int(np.argmax(f_mus[i0:i1]))
    ff = float(r.meta["f_flex_peak_n"])
    return {"label": lab, "h_m": h, "f_flex_n": ff, "f_flex_bw": ff / BW,
            "theta_peak_deg": degrees(th_w), "theta_max_deg": degrees(th_of(COM[lab])),
            "exhausted": bool(th_w > th_of(COM[lab])),
            "f_muscle_peak_n": float(f_mus[i0:i1].max()),
            "f_ecc_cap_at_peak_n": float(f_cap[j]),
            "saturation": float(f_mus[i0:i1].max() / max(f_cap[j], 1e-9))}


def hstar(rows):
    for i in range(len(rows) - 1):
        a, b = rows[i], rows[i + 1]
        if a["exhausted"] or b["exhausted"]:
            continue
        if a["f_flex_n"] <= F_THR <= b["f_flex_n"]:
            return a["h_m"] + (F_THR - a["f_flex_n"]) / (b["f_flex_n"] - a["f_flex_n"]) * \
                (b["h_m"] - a["h_m"])
    return None


print("=" * 96)
print("C1-J0 POSITIVE CONTROL（门关 ⇒ 逐位复现 fp_anchor_probe.json）")
print("=" * 96)
up = json.loads(UP.read_text(encoding="utf-8"))
pcok = True
for lab in LABS:
    w = up["main_rows"][lab]
    g = run(lab, 0.42, None, on_c1=False)
    good = abs(g["f_flex_n"] - w["f_flex_n"]) < 1e-6
    pcok &= good
    print(f"  {lab:<6} want {w['f_flex_n']:.6f} | got {g['f_flex_n']:.6f}  "
          f"{'PASS' if good else 'FAIL'}")
if not pcok:
    print("!! POSITIVE CONTROL FAILED -- INVALID")
    OUT.write_text(json.dumps({"positive_control": False}, ensure_ascii=False, indent=2),
                   encoding="utf-8")
    raise SystemExit(1)
print("=> POSITIVE CONTROL PASS")
print(f"   l={L} | k_θ={K_TH} | M_ecc={M_ECC} N·m | v_ref={V_REF} m/s | 实测跟垫 k_n={PAD_KN:.1e}")

print()
print("=" * 96)
print("O1-O4  C1（B2a + 饱和离心律）｜ 实测跟垫口径")
print("=" * 96)
print(f"  {'档':<6} " + " ".join(f"{h:>7.2f}m" for h in H_GRID)
      + f" {'h*(m)':>9} {'单调':>5} {'θ@峰°':>7} {'θ_max°':>7} {'饱和':>6} {'上限kN':>7}")
grid, verdict = {}, {}
for lab in LABS:
    rows = []
    for h in H_GRID:
        print(f"    [run] {lab}/h={h:.2f}", flush=True)
        rows.append(run(lab, h, PAD_KN, on_c1=True))
    grid[lab] = rows
    hs = hstar(rows)
    mono = all(rows[i]["f_flex_n"] < rows[i + 1]["f_flex_n"] for i in range(len(rows) - 1))
    verdict[lab] = {"h_star_m": hs, "monotone": mono,
                    "sat_hi": max(r["saturation"] for r in rows),
                    "exhausted_at": [r["h_m"] for r in rows if r["exhausted"]]}
    print(f"  {lab:<6} " + " ".join(f"{r['f_flex_n']/1e3:>7.2f} " for r in rows)
          + f" {('%.3f' % hs) if hs else '不可定':>9} {('是' if mono else '否'):>5} "
          f" {rows[-1]['theta_peak_deg']:>7.1f} {rows[-1]['theta_max_deg']:>7.1f} "
          f" {verdict[lab]['sat_hi']:>6.2f} {rows[-1]['f_ecc_cap_at_peak_n']/1e3:>7.2f}")

print()
print("=" * 96)
print("判定（预注册 §5）")
print("=" * 96)
j1 = {l: verdict[l]["monotone"] for l in LABS}
j2 = {l: verdict[l]["h_star_m"] is not None for l in LABS}
j3 = {l: verdict[l]["sat_hi"] > 0.95 for l in LABS}
for l in LABS:
    v = verdict[l]
    print(f"  {l:<6} J1单调={v['monotone']}  J2 h*={('%.3f' % v['h_star_m']) if v['h_star_m'] else '不可定'}"
          f"  J3饱和={v['sat_hi']:.2f}")
print(f"  => C1-J1(三档)={'PASS' if all(j1.values()) else 'FAIL'}  "
      f"C1-J2={'PASS' if all(j2.values()) else 'FAIL'}  "
      f"C1-J3(深高度饱和)={'PASS' if any(j3.values()) else 'FAIL'}")

print()
print("=== v_ref 敏感性（弱锚，§9-1）===")
for vr in (0.5, 2.0):
    row = [run("self", h, PAD_KN, on_c1=True, v_ref=vr) for h in (1.0, 2.0, 5.0)]
    print(f"  v_ref={vr}: self F_flex = "
          + " / ".join(f"{r['f_flex_n']/1e3:.2f}" for r in row) + " kN")

OUT.write_text(json.dumps({"positive_control": True, "l_seg_m": L, "k_theta": K_TH,
                           "m_ecc_nm": M_ECC, "v_ref_ms": V_REF, "pad_kn": PAD_KN,
                           "grid": grid, "verdict": verdict,
                           "C1_J1": j1, "C1_J1_all": all(j1.values()),
                           "C1_J2": j2, "C1_J2_all": all(j2.values()),
                           "C1_J3": j3}, ensure_ascii=False, indent=2), encoding="utf-8")
print(f"\nOUT = {OUT}")
