"""P5 前置侦察 · Step 2b —— 补 20 m（V8b 判据高度）的垫刚度扫描。

Step 2 只跑了 V3 的高度 [2,3,5]。V8b = t_torso @20 m <= 5.0 ms，
所以 20 m 必须补。同时把 Step 2 的判定边界写清楚：
  * V3 临界 α（刚好过 12.5 ms 的最小 α）
  * V8b @20 m 能不能靠加硬垫子过 5.0 ms
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from climbing.head_branch import HeadBranchParams, simulate_head_branch  # noqa: E402
from climbing.pad import CrashPad, simulate_boulder_fall  # noqa: E402

SOLVER_PAD = dict(max_ode_step=1e-4, rtol=1e-8, atol=1e-11)
HEAD_POSTURE = "feet-first-stiff"
HEAD_MASS_KG = 80.0
H = 20.0
ALPHAS = (1.0, 3.0, 10.0, 30.0, 100.0, 300.0)
K_NECKS = (200.0, 2500.0)
V8B_THR_MS = 5.0

# 已发表端点（正对照：α=1 必须复现 pad@20 的 29.4 / 26.6 / 14.6）
POS_CONTROL = {200: 26.6, 2500: 14.6}
# rigid 端（同 JSON）：t_torso=14.0, t_moment=21.4, head@k2500=8.6
RIGID_REF = {"t_torso_ms": 14.0, "tpeak_ms_k2500": 8.6}

OUT = ROOT / "results" / "p5_pad_stiffness_sweep_h20.json"


def scaled_pad(alpha: float) -> CrashPad:
    b = CrashPad()
    return CrashPad(
        thickness_m=b.thickness_m,
        s_plateau_pa=b.s_plateau_pa * alpha,
        s_dense_pa=b.s_dense_pa * alpha,
        xi_plateau=b.xi_plateau,
        n_dense=b.n_dense,
        xi_dense_max=b.xi_dense_max,
        bottoming_modulus_pa=b.bottoming_modulus_pa * alpha,
        linear_modulus_pa=b.linear_modulus_pa,
    )


def peak_rel_ms(r, arr) -> float:
    t = np.asarray(r.t_s, dtype=float)
    a = np.asarray(arr, dtype=float)
    t0, t1 = r.meta["impact_window_s"]
    m = (t >= t0) & (t <= t1)
    i = int(np.argmax(np.abs(a[m])))
    return round(float(t[m][i] - t0) * 1e3, 2)


def main() -> int:
    rows = []
    for a in ALPHAS:
        print(f"[run] h={H} m  alpha={a:>6.1f}  ...", flush=True)
        r = simulate_boulder_fall(
            height_m=H, mass_kg=HEAD_MASS_KG, posture=HEAD_POSTURE,
            on_pad=True, pad=scaled_pad(a), **SOLVER_PAD)
        w = r.meta["impact_window_s"]
        row = {"h_m": H, "alpha": a,
               "window_ms": round(float(w[1] - w[0]) * 1e3, 2),
               "tF_ms": peak_rel_ms(r, r.pad_force_n),
               "t_torso_ms": peak_rel_ms(r, r.accel_torso_g),
               "peak_pad_n": float(np.max(np.abs(r.pad_force_n)))}
        for kn in K_NECKS:
            v = simulate_head_branch(r, HeadBranchParams(k_neck_n_per_mm=kn))
            row[f"tpeak_ms_k{int(kn)}"] = round(v.peak_time_s * 1e3, 2)
        rows.append(row)

    # 正对照：alpha=1 复现已发表【垫】@20 m
    print()
    print("=" * 72)
    print(f"POSITIVE CONTROL  h={H} m, alpha=1 vs published pad column")
    print("=" * 72)
    a1 = rows[0]
    ok = True
    for kn, want in POS_CONTROL.items():
        have = a1[f"tpeak_ms_k{kn}"]
        good = abs(have - want) < 0.05
        ok &= good
        print(f"  head k={kn:>4}  want={want:>6.1f}  got={have:>6.2f}  "
              f"{'PASS' if good else 'FAIL'}")
    # t_torso 参考值（同 JSON pad@20 t_torso=29.4）
    good_tt = abs(a1["t_torso_ms"] - 29.4) < 0.05
    ok &= good_tt
    print(f"  t_torso     want=  29.4  got={a1['t_torso_ms']:>6.2f}  "
          f"{'PASS' if good_tt else 'FAIL'}")
    if not ok:
        print("!! POSITIVE CONTROL FAILED -- INVALID")
        return 1
    print("=> POSITIVE CONTROL PASS")

    print()
    print("=" * 72)
    print(f"PAD STIFFNESS SWEEP @ {H:.0f} m   (V8b threshold = "
          f"{V8B_THR_MS} ms)")
    print("=" * 72)
    print(f"{'alpha':>7} {'window':>7} {'t_torso':>8} {'tF':>6} "
          f"{'tpk200':>7} {'tpk2500':>8} {'peakPadN':>10} "
          f"{'V8b':>6} {'V3@2500':>8}")
    print("-" * 72)
    for r in rows:
        v8b = "PASS" if r["t_torso_ms"] <= V8B_THR_MS else "fail"
        v3 = "PASS" if r["tpeak_ms_k2500"] <= 12.5 else "fail"
        print(f"{r['alpha']:>7.1f} {r['window_ms']:>7.2f} "
              f"{r['t_torso_ms']:>8.2f} {r['tF_ms']:>6.2f} "
              f"{r['tpeak_ms_k200']:>7.2f} {r['tpeak_ms_k2500']:>8.2f} "
              f"{r['peak_pad_n']:>10.1f} {v8b:>6} {v3:>8}")
    print(f"{'rigid':>7} {'':>7} {RIGID_REF['t_torso_ms']:>8.2f} {'':>6} "
          f"{'':>7} {RIGID_REF['tpeak_ms_k2500']:>8.2f}")
    print("-" * 72)

    seq = [r["t_torso_ms"] for r in rows]
    mono = all(seq[i] <= seq[i - 1] + 1e-9 for i in range(1, len(seq)))
    print()
    print(f"t_torso monotone in alpha : {'YES' if mono else 'NO'}  {seq}")
    print(f"min t_torso               : {min(seq):.2f} ms "
          f"(threshold {V8B_THR_MS})")
    print(f"V8b reachable by stiffening pad : "
          f"{'YES' if min(seq) <= V8B_THR_MS else 'NO'}")

    OUT.write_text(json.dumps({
        "rows": rows, "positive_control": True,
        "v8b_threshold_ms": V8B_THR_MS, "rigid_ref": RIGID_REF,
        "t_torso_monotone": mono,
        "v8b_reachable": bool(min(seq) <= V8B_THR_MS),
    }, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"OUT = {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
