"""G4 · 压缩相本构探针 —— V8a/V8b + 回弹对比 + 能量账（vs G1 基线）。

对 {2,3,5,10,20} m 刚性地面工况同时跑 G1 基线（``chain2d.py``）与 G4 压溃
（``chain2d_crush.py``），取**胫骨段**局部组织峰（``segment_local_peak`` 同
口径 = ``scripts/paper/chain_vtarget_probe.py``），判定：

* **V8a**：峰时随高度单调不增（crush 路径）。
* **V8b**：20 m 时 ≤ 5.0 ms（``head_branch.V8_THRESHOLD_AT_20M_S``）。
* **回弹对比**：扩展观察窗（t0→t0+250ms）的 ``head_rise_m / foot_lift_mm``，
  基线与 crush 各算一份，量化 G4 的「整身回弹消除」声明。注：模块
  ``rebound_metrics`` 的窗口 = 首次接触到首次分离，回弹发生在分离后，
  窗口内测不到，故探针用扩展观察窗（基线与 crush 同口径）。
* **能量账**：``crush_dissipated_energy_j`` vs 输入能量（``m·g·h``）。

口径（跑前冻结，禁止回解凑靶）
* 姿势 ``feet-first-stiff``、质量 80 kg、``on_pad=False``（刚地面）——
  与既有 ``scripts/paper/chain_vtarget_probe.py`` 口径一致。
* 观测量先验冻结为**胫骨段**（索引 1）；逐段峰全报（分布）但不改主判据。

产物：``results/chain_crush_probe.json``（纯增量，新文件）。

用法::

    $env:PYTHONPATH="src"; & .venv\\Scripts\\python.exe scripts\\chain_crush_probe.py [h1 h2 ...]
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
from climbing.chain2d_crush import (  # noqa: E402
    CrushChainConfig,
    crush_dissipated_energy_j,
    rebound_metrics,
    simulate_crush_chain_fall,
)
from climbing.head_branch import V8_THRESHOLD_AT_20M_S  # noqa: E402

REBOUND_HORIZON_S = 0.25   # 观察窗 = 窗口起点 + 250 ms（覆盖整链回弹振荡）


def _extended_rebound(result: "PlanarChainResult", horizon: float = REBOUND_HORIZON_S) -> dict:
    """扩展观察窗回弹度量：t0 → t0+horizon 内头节点最大上升 + 足节点最大离地。

    注：``rebound_metrics`` 的窗口 = 首次接触→首次分离，回弹发生在分离**之后**，
    窗口内测不到；故此处用扩展观察窗（基线与 crush 同口径）量化整身回弹。
    """
    t0, _ = result.meta["impact_window_s"]
    t = result.t_s
    m = (t >= t0) & (t <= t0 + horizon)
    if not m.any():
        return {"head_rise_m": 0.0, "foot_lift_mm": 0.0}
    z_head = result.node_z[m, -1]
    z_foot = result.node_z[m, 0]
    head_rise = max(0.0, float(z_head[0] - z_head.min()))
    foot_lift = max(0.0, float(-z_foot.min())) * 1e3
    return {"head_rise_m": head_rise, "foot_lift_mm": foot_lift}

ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "chain_crush_probe.json"

POSTURE = "feet-first-stiff"
MASS_KG = 80.0
HEIGHTS = (2.0, 3.0, 5.0, 10.0, 20.0)
TIBIA_INDEX = 1   # SEGMENTS_6[1] = "小腿-胫骨"

REF_T22_MS = {2.0: None, 3.0: None, 5.0: 11.9, 10.0: 6.5, 20.0: 4.6}


def _run(h: float) -> dict:
    cfg_base = PlanarChainConfig(mass_kg=MASS_KG)
    cfg_crush = CrushChainConfig(mass_kg=MASS_KG)
    r_base = simulate_planar_chain_fall(
        height_m=h, mass_kg=MASS_KG, posture=POSTURE, on_pad=False, config=cfg_base,
    )
    r_crush = simulate_crush_chain_fall(
        height_m=h, mass_kg=MASS_KG, posture=POSTURE, on_pad=False, config=cfg_crush,
    )

    def row(label: str, r) -> dict:
        tp, eps = segment_local_peak(r, TIBIA_INDEX)
        rb = _extended_rebound(r)
        return {
            "label": label,
            "t_tibia_ms": round(tp * 1e3, 3),
            "eps_tibia": eps,
            "window_ms": round((r.meta["impact_window_s"][1]
                                - r.meta["impact_window_s"][0]) * 1e3, 2),
            "impact_speed_ms": round(r.impact_speed_ms, 3),
            "peak_force_kn": round(r.peak_force_n / 1e3, 2),
            "max_compression_mm": round(r.max_compression_m * 1e3, 2),
            "head_rise_m": rb["head_rise_m"],
            "foot_lift_mm": rb["foot_lift_mm"],
            "per_segment_local_peak": {
                s.name: [round(tp_i * 1e3, 2), val_i]
                for s, (tp_i, val_i) in zip(
                    SEGMENTS_6, [segment_local_peak(r, i) for i in range(len(SEGMENTS_6))]
                )
            },
        }

    base = row("baseline", r_base)
    crush = row("crush", r_crush)
    crush["dissipated_j"] = crush_dissipated_energy_j(r_crush)
    crush["impact_energy_j"] = MASS_KG * cfg_crush.g * h
    return {
        "h_m": h,
        "t22_tibia_ms": REF_T22_MS.get(h),
        "baseline": base,
        "crush": crush,
    }


def main() -> int:
    heights = [float(x) for x in sys.argv[1:]] if len(sys.argv) > 1 else list(HEIGHTS)
    rows = [_run(h) for h in heights]

    seq = [r["crush"]["t_tibia_ms"] for r in rows]
    v8a = all(seq[i] <= seq[i - 1] + 1e-9 for i in range(1, len(seq)))
    t20 = seq[-1]
    v8b_threshold_ms = V8_THRESHOLD_AT_20M_S * 1e3
    v8b = bool(t20 <= v8b_threshold_ms)

    print("===== G4 · chain2d_crush probe (baseline vs crush) =====")
    print(f"posture={POSTURE} mass={MASS_KG} on_pad=False")
    for r in rows:
        b, c = r["baseline"], r["crush"]
        print(
            f"  h={r['h_m']:>4.0f}  t_tibia base={b['t_tibia_ms']:>7.3f} / "
            f"crush={c['t_tibia_ms']:>7.3f} ms  eps_crush={c['eps_tibia']:.4e}  "
            f"Fpk={c['peak_force_kn']:>6.2f}kN  comp={c['max_compression_mm']:>6.2f}mm  "
            f"head_rise base={b['head_rise_m']*1e3:>7.1f} / crush={c['head_rise_m']*1e3:>7.1f} mm  "
            f"foot_lift base={b['foot_lift_mm']:>6.1f} / crush={c['foot_lift_mm']:>6.1f} mm  "
            f"E_diss={c.get('dissipated_j', 0)/1e3:>6.1f}kJ / E_in={c.get('impact_energy_j', 0)/1e3:>6.1f}kJ"
        )
    print(f"  V8a (crush monotone non-increasing): {v8a}")
    print(f"  V8b @20m = {t20:.3f} ms (threshold ≤ {v8b_threshold_ms:.1f} ms): {v8b}")

    data = {
        "meta": {
            "module": "climbing.chain2d_crush (G4 压缩相本构) vs climbing.chain2d (G1 基线)",
            "posture": POSTURE,
            "mass_kg": MASS_KG,
            "on_pad": False,
            "heights_m": list(heights),
            "observable": "eps_tibia = (F_ankle/A_tibia)/E_tibia; first-significant "
                          "local peak inside meta['impact_window_s'] "
                          "(head_branch._first_significant_peak, 同 chain_vtarget_probe)",
            "tibia_segment": SEGMENTS_6[TIBIA_INDEX].name,
            "v8b_threshold_ms": v8b_threshold_ms,
            "ref_t22_tibia_ms": {str(k): v for k, v in REF_T22_MS.items()},
            "note": "纯增量新增件；观测量跑前冻结（胫骨段首显著局部峰），"
                    "禁止事后换口径。回弹 = rebound_metrics 同口径。",
        },
        "rows": rows,
        "verdict": {
            "V8a_monotone_non_increasing": v8a,
            "V8b_tibia_at_20m_ms": t20,
            "V8b_threshold_ms": v8b_threshold_ms,
            "V8b_pass": v8b,
        },
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\nwritten: {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
