"""力台锚定探针（LEG_K_ANCHOR 预注册 v2）—— `docs/methods/P3_腿刚度力台锚定_预注册.md`

为什么
------
`pad.py` 三档落地姿势的 `k_flex`（`controlled-drop` 2.2e5、`feet-first-stiff` 1.1e6 N/m）
是**未标定的自由参数**。外部锚定值来自 `PMC9465758`（Ward 2019 口径
``Stiffness = GRFmax / Depth``；该表 `GRFmax` 以 **N·kg⁻¹** 报 ⇒ 表值量纲实为
**N·kg⁻¹·m⁻¹**，须 ×76 kg 才是物理刚度）⇒ soft/self/stiff = **6.3 / 8.2 / 18.5 kN/m**。

正对照（positive control，不过则**整份结果作废**）
--------------------------------------------------
同一姿势 `feet-first-stiff` @1 m 刚性地面，**两种调用必须给两个已知值**：
  * 传**键** `"feet-first-stiff"` ⇒ `F_pad` **37.46** / `F_flex` **31.10** kN（**含 L2 门**：k×2.0）
  * 传 `replace(...)` **对象**         ⇒ `F_pad` **33.65** / `F_flex` **27.43** kN（**无门**）
这同时证明：① 代码路径可复现；② `pad.py` 的 **L2 刚性门**存在且**按姿势*键*判定**
（传对象静默绕门，差 **13.4%**）。锚定档的键是 `anchored-*` ⇒ **门不生效**（本件要的口径）。

判据（§5，跑前冻结）
--------------------
J1a `F_flex/BW` ÷ 实测靶 ∈ [0.70,1.30]；J1b `F_pad/BW` ÷ 实测靶 同带；J1 = J1a AND J1b；
J2 `h*` ≤ 1.52 m（弱外部判据）；J3 单调 `k↑ ⇒ h*↓`。
⛔ 禁止以「使 V1 进 6–9 m」为任何选择标准。

用法::

    .venv\\Scripts\\python.exe scripts/1d/fp_anchor_probe.py
"""
from __future__ import annotations

import json
import pathlib
import sys
from dataclasses import replace

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from climbing.pad import POSTURES, simulate_boulder_fall  # noqa: E402

MASS_KG = 80.0
BW = MASS_KG * 9.80665
F_THR_N = 6200.0                     # 尸体跟骨阈（Yoganandan），继承不新设
H_MAIN = 0.42                        # 力台落箱高度
H_GRID = (0.20, 0.50, 1.00, 2.00, 5.00, 7.00, 9.00)
LABELS = ("soft", "self", "stiff")
#: 该表 `GRFmax (N/kg)` ÷ 9.80665 = BW 倍数
FP_TARGET_BW = {"soft": 3.38, "self": 3.70, "stiff": 5.36}
BAND = (0.70, 1.30)
H_STAR_MAX_M = 1.52

#: J0 正对照（逐位）
POSCTRL = {
    "by_key": {"posture": "feet-first-stiff", "F_pad_kN": 37.46, "F_flex_kN": 31.10},
    "by_object": {"posture": "feet-first-stiff", "F_pad_kN": 33.65, "F_flex_kN": 27.43},
}
POSCTRL_TOL_KN = 0.05

OUT = ROOT / "results" / "fp_anchor_probe.json"


def _run(posture, h: float) -> dict:
    r = simulate_boulder_fall(height_m=h, mass_kg=MASS_KG, on_pad=False, posture=posture)
    return {
        "h_m": h,
        "k_flex_n_per_m": float(POSTURES[posture].k_flex) if isinstance(posture, str)
        else float(posture.k_flex),
        "f_pad_n": float(r.peak_force_n),
        "f_flex_n": float(r.meta["f_flex_peak_n"]),
        "f_pad_bw": float(r.peak_force_n) / BW,
        "f_flex_bw": float(r.meta["f_flex_peak_n"]) / BW,
        "bottomed_out": bool(r.bottomed_out),
    }


def positive_control() -> tuple[bool, list[dict]]:
    """J0：同一姿势两种调用的双正对照。"""
    rows = []
    ok = True
    # 传键（含门）
    r = simulate_boulder_fall(height_m=1.0, mass_kg=MASS_KG, on_pad=False,
                             posture="feet-first-stiff")
    got = (float(r.peak_force_n) / 1e3, float(r.meta["f_flex_peak_n"]) / 1e3)
    want = (POSCTRL["by_key"]["F_pad_kN"], POSCTRL["by_key"]["F_flex_kN"])
    good = all(abs(a - b) <= POSCTRL_TOL_KN for a, b in zip(got, want))
    ok &= good
    rows.append({"call": "by_key(gated)", "want": want, "got": got, "pass": good})
    # 传对象（无门）
    obj = replace(POSTURES["feet-first-stiff"], name="__posctrl")
    r = simulate_boulder_fall(height_m=1.0, mass_kg=MASS_KG, on_pad=False, posture=obj)
    got = (float(r.peak_force_n) / 1e3, float(r.meta["f_flex_peak_n"]) / 1e3)
    want = (POSCTRL["by_object"]["F_pad_kN"], POSCTRL["by_object"]["F_flex_kN"])
    good = all(abs(a - b) <= POSCTRL_TOL_KN for a, b in zip(got, want))
    ok &= good
    rows.append({"call": "by_object(ungated)", "want": want, "got": got, "pass": good})
    return ok, rows


def main() -> int:
    print("=" * 78)
    print("J0  POSITIVE CONTROL  (feet-first-stiff @1 m: key vs object)")
    print("=" * 78)
    pcok, pcrows = positive_control()
    for r_ in pcrows:
        print(f"  {r_['call']:<22} want=({r_['want'][0]:.2f}, {r_['want'][1]:.2f})  "
              f"got=({r_['got'][0]:.2f}, {r_['got'][1]:.2f})  "
              f"{'PASS' if r_['pass'] else 'FAIL'}")
    if not pcok:
        print("!! POSITIVE CONTROL FAILED -- results INVALID, do not use")
        OUT.write_text(json.dumps({"positive_control": False, "rows": pcrows},
                                  ensure_ascii=False, indent=2), encoding="utf-8")
        return 1
    print("=> POSITIVE CONTROL PASS")
    print(f"   (门的存在性证据：键 vs 对象 @F_flex "
          f"{pcrows[0]['got'][1]:.2f} vs {pcrows[1]['got'][1]:.2f} kN ⇒ "
          f"{pcrows[0]['got'][1] / pcrows[1]['got'][1] - 1:+.1%})")

    print()
    print("=" * 78)
    print(f"O1/O2  @{H_MAIN} m 硬地面（对照该表 GRFmax，BW 归一化）")
    print("=" * 78)
    print(f"  {'档':<7} {'k(kN/m)':>8} {'F_pad(kN)':>10} {'F_pad/BW':>9} "
          f"{'F_flex(kN)':>11} {'F_flex/BW':>10} {'实测靶':>7} {'J1a ratio':>10} "
          f"{'J1b ratio':>10}")
    main_rows = {}
    for lab in LABELS:
        key = f"anchored-{lab}"
        row = _run(key, H_MAIN)
        tgt = FP_TARGET_BW[lab]
        r_a = row["f_flex_bw"] / tgt
        r_b = row["f_pad_bw"] / tgt
        main_rows[lab] = {**row, "fp_target_bw": tgt, "J1a_ratio": r_a, "J1b_ratio": r_b}
        print(f"  {lab:<7} {row['k_flex_n_per_m']/1e3:>8.1f} {row['f_pad_n']/1e3:>10.2f} "
              f"{row['f_pad_bw']:>9.2f} {row['f_flex_n']/1e3:>11.2f} "
              f"{row['f_flex_bw']:>10.2f} {tgt:>7.2f} {r_a:>10.3f} {r_b:>10.3f}")

    print()
    print("=" * 78)
    print("O3/J2/J3  力判据首折 h*（F_flex 折线内插到 6.2 kN）")
    print("=" * 78)
    print(f"  {'档':<7} {'k(kN/m)':>8} " + " ".join(f"{h:>7.2f}m" for h in H_GRID)
          + f" {'h*(m)':>8}")
    grid = {}
    hstar = {}
    for lab in LABELS:
        key = f"anchored-{lab}"
        vals = []
        for h in H_GRID:
            row = _run(key, h)
            vals.append(row["f_flex_n"])
            grid.setdefault(lab, []).append({"h_m": h, **row})
        hs = None
        for i in range(len(H_GRID) - 1):
            f0, f1 = vals[i], vals[i + 1]
            if f0 <= F_THR_N <= f1:
                hs = H_GRID[i] + (F_THR_N - f0) / (f1 - f0) * (H_GRID[i + 1] - H_GRID[i])
                break
        hstar[lab] = hs
        print(f"  {lab:<7} {grid[lab][0]['k_flex_n_per_m']/1e3:>8.1f} "
              + " ".join(f"{v/1e3:>7.2f} " for v in vals)
              + f" {('%.3f' % hs) if hs else '>9':>8}")

    # ---- 判定 ---------------------------------------------------------
    j1a = {lab: (BAND[0] <= main_rows[lab]["J1a_ratio"] <= BAND[1]) for lab in LABELS}
    j1b = {lab: (BAND[0] <= main_rows[lab]["J1b_ratio"] <= BAND[1]) for lab in LABELS}
    j1 = all(j1a.values()) and all(j1b.values())
    j2 = all((hstar[lab] is not None and hstar[lab] <= H_STAR_MAX_M) for lab in LABELS)
    j3 = (hstar["soft"] is not None and hstar["self"] is not None
          and hstar["stiff"] is not None
          and hstar["soft"] > hstar["self"] > hstar["stiff"])

    print()
    print("=" * 78)
    print("判定（预注册 §5，如实二值）")
    print("=" * 78)
    for lab in LABELS:
        print(f"  J1a {lab:<6} ratio={main_rows[lab]['J1a_ratio']:.3f}  "
              f"{'PASS' if j1a[lab] else 'FAIL'}")
    for lab in LABELS:
        print(f"  J1b {lab:<6} ratio={main_rows[lab]['J1b_ratio']:.3f}  "
              f"{'PASS' if j1b[lab] else 'FAIL'}")
    print(f"  J1  = J1a AND J1b                              "
          f"{'PASS' if j1 else 'FAIL'}")
    for lab in LABELS:
        hs = hstar[lab]
        print(f"  J2  {lab:<6} h*={('%.3f' % hs) if hs else '>9'} m  "
              f"{'PASS' if (hs and hs <= H_STAR_MAX_M) else 'FAIL'}")
    print(f"  J3  单调 k↑⇒h*↓                                {'PASS' if j3 else 'FAIL'}")

    payload = {
        "positive_control": True,
        "positive_control_rows": pcrows,
        "gate_note": ("pad.py L2 刚性门按姿势*键*判定；传 Posture 对象静默绕门。"
                      "锚定档键 anchored-* ⇒ 门不生效。"),
        "mass_kg": MASS_KG, "bw_n": BW, "f_threshold_n": F_THR_N,
        "h_main_m": H_MAIN, "h_grid_m": list(H_GRID),
        "fp_target_bw": FP_TARGET_BW, "band": list(BAND), "h_star_max_m": H_STAR_MAX_M,
        "anchored_k_n_per_m": {lab: float(POSTURES[f"anchored-{lab}"].k_flex)
                               for lab in LABELS},
        "main_rows": main_rows, "grid": grid, "h_star_m": hstar,
        "verdict": {"J1a": j1a, "J1b": j1b, "J1": j1, "J2": j2, "J3": j3},
    }
    OUT.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print()
    print(f"OUT = {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
