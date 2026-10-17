"""G1 · V8b 新观测量探针 —— 6 段平面软关节链的局部组织峰时。

对 {2,3,5,10,20} m 刚性地面工况跑 ``simulate_planar_chain_fall``，取
**胫骨段**的局部组织峰 ``(t_peak, ε_peak)``（``segment_local_peak``，
= ``ε=(F/A)/E`` 在 ``meta["impact_window_s"]`` 内首显著局部极大，复用
``head_branch._first_significant_peak``），判定：

* **V8a**：峰时随高度**单调不增**（[T22] Table 3 胫骨列：11.9/6.5/4.6 ms）。
* **V8b**：20 m 时 **≤ 5.0 ms**（``head_branch.V8_THRESHOLD_AT_20M_S``）。

口径（跑前冻结，禁止回解凑靶）
* 姿势 ``feet-first-stiff``、质量 80 kg、``on_pad=False``（刚地面）——
  与既有 ``scripts/paper/rigid_vtarget_probe.py`` 的 V8 口径一致。
* 观测量先验冻结为**胫骨段**（索引 1）；逐段峰全报（分布）但不改主判据。

产物：``results/chain_vtarget_probe.json``（纯增量，新文件）。

用法::

    $env:PYTHONPATH="src"; & .venv\\Scripts\\python.exe scripts\\chain_vtarget_probe.py
"""

from __future__ import annotations

import json
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "src"))

from climbing.chain2d import (  # noqa: E402
    SEGMENTS_6,
    PlanarChainConfig,
    segment_local_peak,
    simulate_planar_chain_fall,
)
from climbing.head_branch import V8_THRESHOLD_AT_20M_S  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "chain_vtarget_probe.json"

POSTURE = "feet-first-stiff"
MASS_KG = 80.0
HEIGHTS = (2.0, 3.0, 5.0, 10.0, 20.0)
TIBIA_INDEX = 1   # SEGMENTS_6[1] = "小腿-胫骨"

# 对照：1D 双质点链的 t_torso@20（既有 V8b 口径下的现状值）。
REF_1D_TTORSO_AT_20M_MS = 14.0
REF_T22_MS = {2.0: None, 3.0: None, 5.0: 11.9, 10.0: 6.5, 20.0: 4.6}


def _run(h: float) -> dict:
    cfg = PlanarChainConfig(mass_kg=MASS_KG)
    r = simulate_planar_chain_fall(
        height_m=h, mass_kg=MASS_KG, posture=POSTURE, on_pad=False, config=cfg,
    )
    t0, t1 = r.meta["impact_window_s"]
    t_tibia, eps_tibia = segment_local_peak(r, TIBIA_INDEX)
    per_seg = {}
    for i, s in enumerate(SEGMENTS_6):
        tp, val = segment_local_peak(r, i)
        per_seg[s.name] = {"t_peak_ms": round(tp * 1e3, 2), "eps_peak": val}
    return {
        "h_m": h,
        "t_tibia_ms": round(t_tibia * 1e3, 3),
        "eps_tibia": eps_tibia,
        "window_ms": round((t1 - t0) * 1e3, 2),
        "impact_speed_ms": round(r.impact_speed_ms, 3),
        "peak_force_kn": round(r.peak_force_kn, 2),
        "max_compression_mm": round(r.max_compression_m * 1e3, 2),
        "monotone_through_collapse": bool(
            np.all(np.diff(r.contact_force_n.sum(axis=1)[
                : int(np.argmax(r.contact_force_n.sum(axis=1))) + 1
            ]) >= -1e-6)
        ),
        "per_segment_local_peak": per_seg,
        "t22_tibia_ms": REF_T22_MS.get(h),
    }


def main() -> int:
    rows = [_run(h) for h in HEIGHTS]

    seq = [r["t_tibia_ms"] for r in rows]
    v8a = all(seq[i] <= seq[i - 1] + 1e-9 for i in range(1, len(seq)))
    t20 = seq[-1]
    v8b_threshold_ms = V8_THRESHOLD_AT_20M_S * 1e3
    v8b = bool(t20 <= v8b_threshold_ms)

    print("===== G1 · chain2d V8b probe =====")
    print(f"posture={POSTURE} mass={MASS_KG} on_pad=False")
    for r in rows:
        print(
            f"  h={r['h_m']:>4.0f}  t_tibia={r['t_tibia_ms']:>7.3f} ms  "
            f"eps={r['eps_tibia']:.4e}  win={r['window_ms']:>6.2f} ms  "
            f"Fpk={r['peak_force_kn']:>6.2f} kN  "
            f"monotone={r['monotone_through_collapse']}"
        )
    print(f"  V8a (monotone non-increasing): {v8a}")
    print(f"  V8b @20m = {t20:.3f} ms (threshold ≤ {v8b_threshold_ms:.1f} ms): {v8b}")

    data = {
        "meta": {
            "module": "climbing.chain2d",
            "posture": POSTURE,
            "mass_kg": MASS_KG,
            "on_pad": False,
            "heights_m": list(HEIGHTS),
            "observable": "eps_tibia = (F_ankle/A_tibia)/E_tibia; first-significant "
                          "local peak inside meta['impact_window_s'] "
                          "(head_branch._first_significant_peak)",
            "tibia_segment": SEGMENTS_6[TIBIA_INDEX].name,
            "v8b_threshold_ms": v8b_threshold_ms,
            "ref_1d_ttorso_at_20m_ms": REF_1D_TTORSO_AT_20M_MS,
            "ref_t22_tibia_ms": {str(k): v for k, v in REF_T22_MS.items()},
            "note": "纯增量新增件；观测量跑前冻结（胫骨段首显著局部峰），"
                    "禁止事后换口径。",
        },
        "rows": rows,
        "verdict": {
            "V8a_monotone_non_increasing": v8a,
            "V8b_tibia_at_20m_ms": t20,
            "V8b_threshold_ms": v8b_threshold_ms,
            "V8b_pass": v8b,
            "improvement_vs_1d_ttorso": (
                REF_1D_TTORSO_AT_20M_MS / t20 if t20 > 0 else None
            ),
        },
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\nwritten: {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
