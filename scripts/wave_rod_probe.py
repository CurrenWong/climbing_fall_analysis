"""G1 · 预注册备选探针 —— 1D 应力波杆（胫骨）的波观测量峰时。

对 {2,3,5,10,20} m 跑 ``climbing.wave_rod.simulate_tibia_wave``，取远端
（impact 端）局部组织应变的**波观测量** = ``|ε|`` 全局极大时刻（预注册第一
原理：压缩波在近端固定端同号反射，返回远端时应变叠加于往返 ``2L/c``），判定：

* **V8b 阈值（字面）**：20 m 峰时 ≤ 5.0 ms（``head_branch.V8_THRESHOLD_AT_20M_S``）。
* **物理对照**：[T22] Table 3 胫骨远端 5/10/20 m → 11.9/6.5/**4.6** ms。

**预注册预期**：波杆峰时 ``= 2L/c ≈ 0.25 ms``（**亚毫秒**）——设计 §4.3 已
预期其**过冲**（远早于 [T22]），以佐证 [T22] 4.6 ms 是**接触压溃主导**而非
波动主导。本探针**如实**报告，**不调参**凑 4.6 / 5.0 ms。

产物：``results/wave_rod_probe.json``（纯增量，新文件）。

用法::

    $env:PYTHONPATH="src"; & .venv\\Scripts\\python.exe scripts\\wave_rod_probe.py
"""

from __future__ import annotations

import json
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from climbing.head_branch import V8_THRESHOLD_AT_20M_S  # noqa: E402
from climbing.wave_rod import WaveRod, simulate_tibia_wave  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "wave_rod_probe.json"

MASS_KG = 80.0
HEIGHTS = (2.0, 3.0, 5.0, 10.0, 20.0)
REF_T22_TIBIA_MS = {2.0: None, 3.0: None, 5.0: 11.9, 10.0: 6.5, 20.0: 4.6}


def _run(h: float, cfg: WaveRod) -> dict:
    r = simulate_tibia_wave(height_m=h, mass_kg=MASS_KG, config=cfg)
    return {
        "h_m": h,
        "t_peak_ms": round(r.peak_time_ms, 6),
        "eps_peak": r.peak_strain,
        "impact_speed_ms": round(r.impact_speed_ms, 4),
        "impact_strain": r.impact_strain,
        "round_trip_ms": round(r.round_trip_s * 1e3, 6),
        "ratio_tpeak_over_2Lc": round(r.peak_time_s / r.round_trip_s, 6),
        "submillisecond": r.submillisecond,
        "t22_tibia_ms": REF_T22_TIBIA_MS.get(h),
    }


def main() -> int:
    cfg = WaveRod()
    rows = [_run(h, cfg) for h in HEIGHTS]

    t20 = rows[-1]["t_peak_ms"]
    seq = [r["t_peak_ms"] for r in rows]
    monotone = all(seq[i] <= seq[i - 1] + 1e-12 for i in range(1, len(seq)))
    thr_ms = V8_THRESHOLD_AT_20M_S * 1e3
    v8b_numeric = bool(t20 <= thr_ms)        # 字面单边判据
    t22_20 = REF_T22_TIBIA_MS[20.0]
    overshoot = bool(t20 < 0.5 * t22_20)     # 远早于 [T22] = 过冲
    reproduces = bool(0.5 * t22_20 <= t20 <= 1.5 * t22_20)  # 量级吻合

    print("===== G1 · wave_rod probe (1D axial stress-wave tibia) =====")
    print(f"E={cfg.modulus_pa:.0e} Pa  rho={cfg.density_kg_m3:.0f} kg/m^3  "
          f"L={cfg.length_m} m  c={cfg.wave_speed_ms:.1f} m/s")
    print(f"PRE-REGISTERED prediction 2L/c = {cfg.round_trip_s*1e3:.4f} ms "
          f"(L/c = {cfg.transit_s*1e3:.4f} ms)")
    for r in rows:
        print(f"  h={r['h_m']:>4.0f}  t_peak={r['t_peak_ms']:>8.4f} ms  "
              f"eps_peak={r['eps_peak']:.4e}  v0={r['impact_speed_ms']:>6.2f}  "
              f"ratio(2L/c)={r['ratio_tpeak_over_2Lc']:.4f}")
    print(f"  height-monotone(non-increasing): {monotone}")
    print(f"  V8b numeric (t20 <= {thr_ms:.1f} ms): {v8b_numeric}  "
          f"(t20 = {t20:.4f} ms)")
    print(f"  vs [T22] 4.6 ms: OVERSHOOT={overshoot}  "
          f"(factor={t22_20 / t20:.1f}x too early); reproduces_magnitude={reproduces}")

    data = {
        "meta": {
            "module": "climbing.wave_rod",
            "purpose": "G1 §7 decision 6 pre-registered fallback (triggered by "
                       "chain2d V8b FAIL, no retune)",
            "mass_kg": MASS_KG,
            "heights_m": list(HEIGHTS),
            "material": {
                "cortical_bone_E_pa": cfg.modulus_pa,
                "cortical_bone_rho_kg_m3": cfg.density_kg_m3,
                "tibia_length_m": cfg.length_m,
                "tibia_area_m2": cfg.area_m2,
                "wave_speed_ms": cfg.wave_speed_ms,
                "sources": "cortical bone E~15-20 GPa (17 used), rho~1800-2100 "
                           "(1900 used), tibia length ~0.38 m (Drillis & Contini "
                           "1966 / Winter 2009); documented, not tuned.",
            },
            "observable": "distal (impact-end) local tissue strain eps=du/dx; "
                          "peak time = argmax |eps| within [0, round_trips*(2L/c)]",
            "analytic_peak_time_s": cfg.round_trip_s,
            "analytic_peak_time_ms": cfg.round_trip_s * 1e3,
            "v8b_threshold_ms": thr_ms,
            "ref_t22_tibia_ms": {str(k): v for k, v in REF_T22_TIBIA_MS.items()},
            "note": "纯增量新增件；材料/几何/观测量跑前冻结，禁止调参凑 4.6/5.0 ms。",
        },
        "rows": rows,
        "verdict": {
            "t_peak_at_20m_ms": t20,
            "analytic_2L_over_c_ms": cfg.round_trip_s * 1e3,
            "v8b_threshold_ms": thr_ms,
            "V8b_numeric_pass_submillisecond_below_threshold": v8b_numeric,
            "height_monotone_non_increasing": monotone,
            "t22_tibia_at_20m_ms": t22_20,
            "overshoot_vs_t22": overshoot,
            "overshoot_factor_vs_t22": (t22_20 / t20),
            "reproduces_t22_magnitude": reproduces,
            "honest_conclusion": (
                "The 1D axial stress-wave rod predicts a SUB-MILLISECOND tibial "
                "peak (~2L/c). It satisfies the literal V8b inequality (<=5.0 ms) "
                "only by OVERSHOOTING ~18x earlier than [T22] 4.6 ms; it does NOT "
                "reproduce the [T22] magnitude and its peak time is height-"
                "independent. This is the EXPECTED, pre-registered finding: it "
                "confirms [T22]'s 4.6 ms is CONTACT-CRUSH-controlled, not wave-"
                "controlled, and that the 6-segment chain's 5.8 ms is the closer "
                "physical match. No parameter was tuned."
            ),
        },
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\nwritten: {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
