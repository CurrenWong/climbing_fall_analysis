"""P5 前置侦察 · Step 2 —— 垫刚度连续扫描（α 倍泡沫应力）。

为什么
------
Step 1 只有两个工况端点（发表【垫】 vs 发表【刚】）。两点差不能证明
"工况"是**单调**起作用的 —— 可能是本构形状差异造成的端点偶然。
本脚本在**垫家族内部**连续加硬：CrashPad 的 (s_plateau, s_dense,
bottoming_modulus) 三者同乘 α，α=1 即发表基线，α→∞ ⇒ 接近刚性。

正对照（positive control，不过则**整份结果作废**）
--------------------------------------------------
α=1 必须逐位复现 results/rigid_vtarget_probe.json 里【垫】的
    h=2: k200=37.0 / k2500=31.4
    h=3: k200=33.8 / k2500=25.6
    h=5: k200=31.2 / k2500=20.4

判据
----
若 t_torso / tpeak 随 α **单调下降**并逼近刚性端（10.8 / 10.2 / 9.4）
⇒ 工况效应是连续、可由"垫子加硬"实现的 ⇒ **P5 的垫子物理路线有真实收益**。
若早段饱和 / 非单调 ⇒ 端点差是本构形状差，P5 加硬垫子**买不到**时刻。

用法::

    .venv\\Scripts\\python.exe temp/p5_cond/pad_stiffness_sweep.py
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from climbing.head_branch import (  # noqa: E402
    HeadBranchParams,
    simulate_head_branch,
)
from climbing.pad import CrashPad, simulate_boulder_fall  # noqa: E402

# ---- 与 scripts/rigid_vtarget_probe.py 逐字同口径 ------------------------
SOLVER_PAD = dict(max_ode_step=1e-4, rtol=1e-8, atol=1e-11)   # 【垫】紧容差
HEAD_POSTURE = "feet-first-stiff"
HEAD_MASS_KG = 80.0
HEIGHTS = [2.0, 3.0, 5.0]
K_NECKS = (200.0, 2500.0)

ALPHAS = (1.0, 3.0, 10.0, 30.0, 100.0)

# 已发表【垫】端点 = 正对照
POSITIVE_CONTROL = {
    2.0: {200: 37.0, 2500: 31.4},
    3.0: {200: 33.8, 2500: 25.6},
    5.0: {200: 31.2, 2500: 20.4},
}
# 已发表【刚】端点 = 单调终点参考（k2500 列 = V3 判据）
RIGID_REF = {2.0: 10.8, 3.0: 10.2, 5.0: 9.4}

V3_THR_MS = 12.5
OUT = ROOT / "results" / "p5_pad_stiffness_sweep.json"


def scaled_pad(alpha: float) -> CrashPad:
    """α=1 逐字等于 CrashPad()；只缩放应力三件套，几何与形状参数不动。"""
    base = CrashPad()
    return CrashPad(
        thickness_m=base.thickness_m,
        s_plateau_pa=base.s_plateau_pa * alpha,
        s_dense_pa=base.s_dense_pa * alpha,
        xi_plateau=base.xi_plateau,
        n_dense=base.n_dense,
        xi_dense_max=base.xi_dense_max,
        bottoming_modulus_pa=base.bottoming_modulus_pa * alpha,
        linear_modulus_pa=base.linear_modulus_pa,
    )


def peak_rel_ms(r, arr) -> float:
    """冲击窗内 |arr| 首峰的 (t - t0)，ms —— 与探针同款窗口守卫。"""
    t = np.asarray(r.t_s, dtype=float)
    a = np.asarray(arr, dtype=float)
    t0, t1 = r.meta["impact_window_s"]
    m = (t >= t0) & (t <= t1)
    i = int(np.argmax(np.abs(a[m])))
    return round(float(t[m][i] - t0) * 1e3, 2)


def one(h: float, alpha: float) -> dict:
    r = simulate_boulder_fall(
        height_m=h, mass_kg=HEAD_MASS_KG, posture=HEAD_POSTURE,
        on_pad=True, pad=scaled_pad(alpha), **SOLVER_PAD,
    )
    row = {
        "h_m": h, "alpha": alpha,
        "window_ms": round(float(r.meta["impact_window_s"][1]
                                 - r.meta["impact_window_s"][0]) * 1e3, 2),
        "tF_ms": peak_rel_ms(r, r.pad_force_n),
        "t_torso_ms": peak_rel_ms(r, r.accel_torso_g),
        "peak_pad_n": float(np.max(np.abs(r.pad_force_n))),
    }
    for kn in K_NECKS:
        v = simulate_head_branch(r, HeadBranchParams(k_neck_n_per_mm=kn))
        row[f"tpeak_ms_k{int(kn)}"] = round(v.peak_time_s * 1e3, 2)
    return row


def main() -> int:
    rows = []
    for h in HEIGHTS:
        for a in ALPHAS:
            print(f"[run] h={h:>4.1f} m  alpha={a:>6.1f}  ...", flush=True)
            rows.append(one(h, a))

    # ---- 正对照 -------------------------------------------------------
    print()
    print("=" * 74)
    print("POSITIVE CONTROL  (alpha=1 must equal published pad column)")
    print("=" * 74)
    ok = True
    for h in HEIGHTS:
        got = next(r for r in rows if r["h_m"] == h and r["alpha"] == 1.0)
        for kn, want in POSITIVE_CONTROL[h].items():
            have = got[f"tpeak_ms_k{kn}"]
            good = abs(have - want) < 0.05
            ok &= good
            print(f"  h={h:>4.1f} k={kn:>4}  want={want:>6.1f}  "
                  f"got={have:>6.2f}  {'PASS' if good else 'FAIL'}")
    if not ok:
        print("!! POSITIVE CONTROL FAILED -- results INVALID, do not use")
        OUT.write_text(json.dumps({"rows": rows, "positive_control": False},
                                  ensure_ascii=False, indent=2), encoding="utf-8")
        return 1
    print("=> POSITIVE CONTROL PASS")

    # ---- 扫描表 -------------------------------------------------------
    print()
    print("=" * 74)
    print("PAD STIFFNESS SWEEP  (alpha x foam stress; rigid ref in last col)")
    print("=" * 74)
    print(f"{'h(m)':>5} {'alpha':>7} {'window':>7} {'t_torso':>8} "
          f"{'tF':>6} {'tpk200':>7} {'tpk2500':>8} {'peakPadN':>9} {'V3@2500':>8}")
    print("-" * 74)
    for h in HEIGHTS:
        ref = RIGID_REF[h]
        for a in ALPHAS:
            r = next(x for x in rows if x["h_m"] == h and x["alpha"] == a)
            tag = "PASS" if r["tpeak_ms_k2500"] <= V3_THR_MS else "fail"
            print(f"{r['h_m']:>5.1f} {r['alpha']:>7.1f} {r['window_ms']:>7.2f} "
                  f"{r['t_torso_ms']:>8.2f} {r['tF_ms']:>6.2f} "
                  f"{r['tpeak_ms_k200']:>7.2f} {r['tpeak_ms_k2500']:>8.2f} "
                  f"{r['peak_pad_n']:>9.1f} {tag:>8}")
        print(f"{'':>5} {'rigid':>7} {'':>7} {'':>8} {'':>6} "
              f"{'':>7} {ref:>8.2f} {'':>9} {'(ref)':>8}")
        print("-" * 74)

    # ---- 单调性 -------------------------------------------------------
    print()
    print("MONOTONICITY (t_torso and tpk2500 must be non-increasing in alpha):")
    mono = {}
    for h in HEIGHTS:
        seq_t = [next(x for x in rows if x["h_m"] == h and x["alpha"] == a)
                 ["t_torso_ms"] for a in ALPHAS]
        seq_k = [next(x for x in rows if x["h_m"] == h and x["alpha"] == a)
                 ["tpeak_ms_k2500"] for a in ALPHAS]
        inc_t = any(seq_t[i] > seq_t[i - 1] + 1e-9 for i in range(1, len(seq_t)))
        inc_k = any(seq_k[i] > seq_k[i - 1] + 1e-9 for i in range(1, len(seq_k)))
        mono[f"{h:g}"] = {
            "t_torso_ms": seq_t, "tpeak_ms_k2500": seq_k,
            "t_torso_non_increasing": not inc_t,
            "tpeak_non_increasing": not inc_k,
            "rigid_ref_ms": RIGID_REF[h],
            "alpha100_gap_to_rigid_ms": round(seq_k[-1] - RIGID_REF[h], 2),
        }
        print(f"  h={h:>4.1f}  t_torso={'MONO' if not inc_t else 'NOT'} {seq_t}"
              f"  | tpk2500={'MONO' if not inc_k else 'NOT'} {seq_k}"
              f"  | rigid ref {RIGID_REF[h]:.1f}")

    payload = {
        "rows": rows, "positive_control": True, "mono": mono,
        "alpha_note": "s_plateau/s_dense/bottoming x alpha; geometry fixed",
        "solver": SOLVER_PAD, "posture": HEAD_POSTURE, "mass_kg": HEAD_MASS_KG,
        "k_necks": list(K_NECKS), "alphas": list(ALPHAS),
        "v3_threshold_ms": V3_THR_MS,
        "rigid_ref_ms": {str(h): RIGID_REF[h] for h in HEIGHTS},
    }
    OUT.write_text(json.dumps(payload, ensure_ascii=False, indent=2),
                   encoding="utf-8")
    print()
    print(f"OUT = {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
