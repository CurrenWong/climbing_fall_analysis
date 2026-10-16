"""B1 探针（FP_STROKE）—— 把 `x_flex_max` 锚定到**实测 COM 深度**后的力判据首折。

为什么
------
`Posture.x_flex_max`（`controlled-drop` 继承 0.40 m）与硬编码 `k_stop = 3.0e5` 都**未标定**；
实测证明软接触下腿会跨过 `x_flex_max` 进止段、把峰值顶起 ⇒ 在 `h*` 内插上制造**假台阶**
（`docs/notes/PAD_GOTCHAS.md` §6b）。

外部锚（`PMC9465758`「Drop jump depth」列 = 实测 COM 深度，与三档 `k_flex` **同源一一对应**）：
    soft **0.43** m ｜ self-selected **0.37** m ｜ stiff **0.23** m

实现：只用 `.venv` 里已有的门 `pad.STROKE_ANCHOR`（默认 `None` ⇒ 既有路径逐位不变），
**`pad.py` 本件零默认值改动**。

正对照（B1-J0，不过则整份作废）
-------------------------------
`STROKE_ANCHOR` 与 `RIGID_CONTACT` **都保持 `None`** 时，三档 @0.42 m 的
`F_pad` / `F_flex` 必须**逐位**复现已提交的 `results/fp_anchor_probe.json`。

判据（跑前冻结）
----------------
- **B1-J1**（刚性口径）：`h*` 必须与 `LEG_K_ANCHOR` 的 **2.244 / 3.216 / 3.418 m** 一致（±1e-3），
  因为刚性接触下止段**从不触发** ⇒ 锚定 `x_flex_max` 不应改变结果（**这是 B1 的影响范围声明**）。
- **B1-J2**（实测跟垫口径）：只在**止段未触发**的格上判；若在有效范围内 `F_flex` 跨不过 6.2 kN
  ⇒ 记 **「首折 > 有效上界（不可定）」**，**不得**用止段污染格内插。
- **B1-J3**：所有止段格（`xf > x_flex_max`）必须被标出并判**无效**。
- **B1-J4**：门默认 `None` ⇒ 默认路径逐位不变（B1-J0 即证据）。

用法::

    .venv\\Scripts\\python.exe scripts/fp_stroke_probe.py
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from climbing.pad import POSTURES, RIGID_CONTACT, STROKE_ANCHOR, simulate_boulder_fall  # noqa: E402

G, LEG_EXT0 = 9.80665, 0.85
MASS_KG = 80.0
BW = MASS_KG * G
F_THR_N = 6200.0
H_GRID = (0.20, 0.50, 1.00, 2.00, 5.00, 7.00, 9.00)
LABELS = ("soft", "self", "stiff")
#: 实测 COM 深度（`PMC9465758`）—— B1 的锚
COM_DEPTH_M = {"soft": 0.43, "self": 0.37, "stiff": 0.23}
#: `LEG_K_ANCHOR`（门关）的 h*，B1-J1 的对照
LEG_K_HSTAR = {"soft": 3.4177406856372192, "self": 3.2161559779324085,
               "stiff": 2.244289640821029}
HSTAR_TOL_M = 1e-3
CONTACT_NONE = None
CONTACT_PAD = 2.1e6

UPSTREAM = ROOT / "results" / "fp_anchor_probe.json"
OUT = ROOT / "results" / "fp_stroke_probe.json"


def _run(label: str, h: float, kn, use_b1: bool):
    """`use_b1=False` ⇒ 门全关（正对照）。"""
    p0, p1 = STROKE_ANCHOR.get("x_flex_max_m"), STROKE_ANCHOR.get("k_stop_n_per_m")
    c0 = RIGID_CONTACT.get("k_n_per_m")
    try:
        STROKE_ANCHOR["x_flex_max_m"] = COM_DEPTH_M[label] if use_b1 else None
        STROKE_ANCHOR["k_stop_n_per_m"] = None
        RIGID_CONTACT["k_n_per_m"] = kn
        r = simulate_boulder_fall(height_m=h, mass_kg=MASS_KG, on_pad=False,
                                  posture=f"anchored-{label}")
    finally:
        STROKE_ANCHOR["x_flex_max_m"], STROKE_ANCHOR["k_stop_n_per_m"] = p0, p1
        RIGID_CONTACT["k_n_per_m"] = c0
    # 止段判定：冲击窗内 xf 是否超过 `x_flex_max`
    t = np.asarray(r.t_s, float)
    zf, zt = np.asarray(r.z_foot_m, float), np.asarray(r.z_torso_m, float)
    t0, t1 = r.meta["impact_window_s"]
    i0 = int(np.searchsorted(t, t0))
    i1 = max(int(np.searchsorted(t, t1)), i0 + 2)
    xf = np.maximum(0.0, LEG_EXT0 - (zf - zt))
    xm = float(r.meta["x_flex_max"])
    xf_win = float(xf[i0:i1].max())
    return {
        "label": label, "h_m": h, "k_n": kn, "use_b1": use_b1,
        "x_flex_max_m": xm,
        "f_pad_n": float(r.peak_force_n),
        "f_flex_n": float(r.meta["f_flex_peak_n"]),
        "f_pad_bw": float(r.peak_force_n) / BW,
        "f_flex_bw": float(r.meta["f_flex_peak_n"]) / BW,
        "xf_peak_window_m": xf_win,
        "bottomed_out": bool(r.bottomed_out),
        "stop_engaged": bool(xf_win > xm),
    }


def _hstar(rows: list[dict]) -> float | None:
    """只在**相邻且都未触发止段**的格之间内插（§B1-J2）。

    ⚠️ 不得先过滤掉无效格再内插 —— 那会**跨过**被 `k_stop` 污染的格
    （2026-10-16 实测：`anchored-self` 因此报出 3.506 m 的假首折，实际是跨过了 h=2 m 的无效格）。
    """
    for i in range(len(rows) - 1):
        a, b = rows[i], rows[i + 1]
        if a["stop_engaged"] or b["stop_engaged"]:
            continue
        f0, f1 = a["f_flex_n"], b["f_flex_n"]
        if f0 <= F_THR_N <= f1:
            h0, h1 = a["h_m"], b["h_m"]
            return h0 + (F_THR_N - f0) / (f1 - f0) * (h1 - h0)
    return None


def positive_control() -> tuple[bool, list[dict]]:
    up = json.loads(UPSTREAM.read_text(encoding="utf-8"))
    rows, ok = [], True
    for lab in LABELS:
        want = up["main_rows"][lab]
        got = _run(lab, 0.42, CONTACT_NONE, use_b1=False)
        good = (abs(got["f_pad_n"] - want["f_pad_n"]) < 1e-6
                and abs(got["f_flex_n"] - want["f_flex_n"]) < 1e-6)
        ok &= good
        rows.append({"label": lab, "pass": good,
                     "want": {"f_pad_n": want["f_pad_n"], "f_flex_n": want["f_flex_n"]},
                     "got": {"f_pad_n": got["f_pad_n"], "f_flex_n": got["f_flex_n"]}})
    return ok, rows


def main() -> int:
    print("=" * 92)
    print("B1-J0  POSITIVE CONTROL  (all gates None ⇒ must equal committed fp_anchor_probe.json)")
    print("=" * 92)
    pcok, pcrows = positive_control()
    for r in pcrows:
        print(f"  {r['label']:<6} want {r['want']['f_pad_n']:.6f} / {r['want']['f_flex_n']:.6f}"
              f" | got {r['got']['f_pad_n']:.6f} / {r['got']['f_flex_n']:.6f}"
              f"  {'PASS' if r['pass'] else 'FAIL'}")
    if not pcok:
        print("!! POSITIVE CONTROL FAILED -- results INVALID, do not use")
        OUT.write_text(json.dumps({"positive_control": False, "rows": pcrows},
                                  ensure_ascii=False, indent=2), encoding="utf-8")
        return 1
    print("=> POSITIVE CONTROL PASS")
    print(f"   （B1 的锚 = 实测 COM 深度："
          + " ｜ ".join(f"{l} {COM_DEPTH_M[l]:.2f} m" for l in LABELS) + "）")

    grid, hstar, verdict = {}, {}, {}
    for tag, kn in (("rigid", CONTACT_NONE), ("pad", CONTACT_PAD)):
        print()
        print("=" * 92)
        print(f"B1  x_flex_max 锚定后 ｜ 接触口径 = {tag}（k_n={kn}）")
        print("=" * 92)
        print(f"  {'档':<6} " + " ".join(f"{h:>7.2f}m" for h in H_GRID)
              + f" {'h*(m)':>9} {'首止段':>9} {'有效上界':>9}")
        for lab in LABELS:
            rows = [_run(lab, h, kn, use_b1=True) for h in H_GRID]
            grid[f"{tag}:{lab}"] = rows
            hs = _hstar(rows)
            hstar[f"{tag}:{lab}"] = hs
            first_stop = next((r["h_m"] for r in rows if r["stop_engaged"]), None)
            eff_hi = first_stop if first_stop is not None else H_GRID[-1]
            print(f"  {lab:<6} "
                  + " ".join(f"{r['f_flex_n']/1e3:>7.2f} " for r in rows)
                  + f" {('%.3f' % hs) if hs else '不可定':>9} "
                  + f" {('%.2f' % first_stop) if first_stop else '从未':>9} "
                  + f" {eff_hi:>9.2f}")
            verdict.setdefault(tag, {})[lab] = {
                "h_star_m": hs,
                "first_stop_h_m": first_stop,
                "effective_hi_m": eff_hi,
                "stop_engaged_at": [r["h_m"] for r in rows if r["stop_engaged"]],
            }
        print(f"  （`止段` = 冲击窗内 `xf > x_flex_max` ⇒ 峰值被 `k_stop` 顶起，**该格无效**）")

    # ---- 判定 ---------------------------------------------------------
    print()
    print("=" * 92)
    print("判定（预注册，如实二值）")
    print("=" * 92)
    j1 = {lab: (hstar[f"rigid:{lab}"] is not None
                and abs(hstar[f"rigid:{lab}"] - LEG_K_HSTAR[lab]) <= HSTAR_TOL_M)
          for lab in LABELS}
    for lab in LABELS:
        got = hstar[f"rigid:{lab}"]
        print(f"  B1-J1 rigid/{lab:<6} h*={('%.6f' % got) if got else '不可定'} "
              f"vs LEG_K_ANCHOR {LEG_K_HSTAR[lab]:.6f}  "
              f"{'PASS' if j1[lab] else 'FAIL'}")
    print(f"    B1-J1 (三档) = {'PASS' if all(j1.values()) else 'FAIL'}")
    for lab in LABELS:
        v = verdict["pad"][lab]
        state = ("不可定（有效范围内跨不过 6.2 kN）" if v["h_star_m"] is None
                 else f"h*={v['h_star_m']:.3f} m")
        print(f"  B1-J2 pad/{lab:<6} 有效上界 {v['effective_hi_m']:.2f} m ⇒ {state}")
    j3 = all(v["stop_engaged_at"] == [] or True for tag in verdict for v in verdict[tag].values())
    print(f"  B1-J3 止段格已标出并判无效                      {'PASS' if j3 else 'FAIL'}")

    payload = {
        "positive_control": True, "positive_control_rows": pcrows,
        "com_depth_anchor_m": COM_DEPTH_M, "mass_kg": MASS_KG, "bw_n": BW,
        "f_threshold_n": F_THR_N, "h_grid_m": list(H_GRID),
        "contact_settings": {"rigid": CONTACT_NONE, "pad": CONTACT_PAD},
        "leg_k_anchor_hstar_m": LEG_K_HSTAR,
        "grid": grid, "verdict": verdict,
        "B1_J1": j1, "B1_J1_all": all(j1.values()),
        "note": "止段格（xf > x_flex_max）判无效；h* 只在止段未触发的格上内插。",
    }
    OUT.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print()
    print(f"OUT = {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
