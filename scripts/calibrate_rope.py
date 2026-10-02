"""UIAA 101 条件下的绳索本构标定。

标定目标
--------
UIAA 101 动态绳坠落试验：80 kg 质量、坠落系数 FF = 1.77、动态单绳、
标准保护系统，**峰值力 < 12 kN**；公开文献中该条件下的实测峰值
典型落在 **8 - 9 kN**，伸长率约 **30-40%**。

⚠️ 这两个目标**在幂律本构下互斥**（详见 RopeModel 的 docstring）：
标定网格里峰值力随 ``t_ref`` 单调上升、伸长率随之单调下降，要 35% 伸长
就得把峰值压到 7 kN 附近。根因是真实动态绳有明显的**趾部**（护套与
编织层先被拉紧，曲线开头很软，纱线芯才咬合），幂律从零点就硬化，
画不出这个拐点。

所以取舍是：**优先保峰值力**（12 kN 判据直接作用在它上面，
低估冲击力是安全分析里更危险的偏置方向），牺牲伸长率。
正解是换分段（趾部 + 硬化）本构，属于 Phase A-2。

目标函数用相对误差的**均方根**，不能用相加 —— 相加会让"峰值偏低"与
"伸长偏高"互相抵消，出现总误差接近 0 的假优解。
"""

from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd

from climbing.rope import RopeModel, simulate_rope_fall

TARGET_PEAK_N = 8_500.0
TARGET_ELONG_PCT = 35.0
UIAA_LIMIT_N = 12_000.0
# 目标函数权重，见 probe() 里的说明
W_PEAK, W_ELONG = 0.7, 0.3


def probe(t_ref: float, exponent: float, eps_ref: float = 0.30, c: float = 100.0) -> dict:
    rope = RopeModel(t_ref_n=t_ref, eps_ref=eps_ref, exponent=exponent, c_damp=c)
    r = simulate_rope_fall(
        mass_kg=80.0, rope_length_m=2.0, fall_factor=1.77,
        t_hold_n=1e9, device="calibration", rope=rope,
        max_ode_step=2e-4, rtol=1e-8, atol=1e-10,
        # 标定只需要首次冲击的峰值。默认 t_impact_max=3.0 会把回弹振荡
        # 一起积分进来，网格要跑十几分钟；1.2 s 已经覆盖 t_peak≈0.6-0.8 s
        # 及其后的峰值窗口。峰值统计本身也已限定在首次冲击窗内
        # (t_peak_window_s)，回弹不会污染结果。
        t_impact_max=1.2,
    )
    # 相对误差的**加权**均方根。不能用相加 —— 那会让"峰值偏低"与
    # "伸长偏高"互相抵消，出现总误差接近 0 的假优解。
    #
    # 权重体现模块 docstring 里写明的取舍：峰值力是安全关键量
    # （UIAA 12 kN 判据直接作用在它上面），低估冲击力是安全分析里
    # 更危险的偏置方向，所以给峰值 0.7、伸长 0.3。
    e1 = (r.peak_tension_kn * 1e3 - TARGET_PEAK_N) / TARGET_PEAK_N
    e2 = (r.max_elongation_pct - TARGET_ELONG_PCT) / TARGET_ELONG_PCT
    return {
        "t_ref_N": t_ref,
        "exponent": exponent,
        "eps_ref": eps_ref,
        "c_damp": c,
        "peak_kN": r.peak_tension_kn,
        "peak_g": r.peak_accel_g,
        "stroke_m": r.max_elongation,
        "elong_pct": r.max_elongation_pct,
        "v_taut_ms": float(np.sqrt(2 * 9.80665 * r.free_fall_m)),
        "err_peak": e1,
        "err_elong": e2,
        "err_rss": float(np.hypot(W_PEAK * e1, W_ELONG * e2)),
    }


def main() -> None:
    out_dir = pathlib.Path(__file__).resolve().parents[1] / "results"
    out_dir.mkdir(exist_ok=True)

    # 网格搜索。范围覆盖最终写入默认值的那一带
    # (t_ref=10000, exponent=2.30, c_damp=100)。
    t_refs = np.arange(6_000, 12_001, 500.0)
    exponents = np.arange(2.10, 2.41, 0.10)
    c_damps = [100.0, 200.0]
    rows = []
    for n in exponents:
        for t in t_refs:
            for c in c_damps:
                rows.append(probe(float(t), float(n), c=c))
    grid = pd.DataFrame(rows)
    best = grid.loc[grid["err_rss"].idxmin()]
    grid.sort_values("err_rss").head(20).to_csv(out_dir / "calibration_rope.csv", index=False)

    print("=== UIAA 条件 (80 kg, FF=1.77, L0=2 m) 标定网格 top-10 ===")
    print(grid.sort_values("err_rss").head(10).to_string(index=False))

    print(f"\n最优: t_ref={best.t_ref_N:.0f} N, exponent={best.exponent:.2f}, "
          f"c_damp={best.c_damp:.0f}")
    print(f"  -> 峰值 {best.peak_kN:.2f} kN ({best.peak_g:.1f} g), "
          f"冲击位移 {best.stroke_m:.3f} m ({best.elong_pct:.1f}%), "
          f"绳绷紧时 v={best.v_taut_ms:.2f} m/s")

    # UIAA 合格性检查：整个 FF 范围内峰值都 < 12 kN，且对 FF 单调
    print("\n=== 标定后模型在 FF 范围内的表现 (ATC 锁死, 80 kg, L0=2 m) ===")
    rope = RopeModel(t_ref_n=float(best.t_ref_N), eps_ref=float(best.eps_ref),
                     exponent=float(best.exponent), c_damp=float(best.c_damp))
    scan = []
    for ff in [0.3, 0.5, 0.75, 1.0, 1.25, 1.5, 1.77, 2.0]:
        r = simulate_rope_fall(mass_kg=80, rope_length_m=2.0, fall_factor=ff,
                               t_hold_n=1e9, device="locked", rope=rope,
                               t_impact_max=1.2)
        scan.append({"FF": ff, "free_fall_m": r.free_fall_m,
                     "peak_kN": r.peak_tension_kn, "peak_g": r.peak_accel_g,
                     "drop_m": r.max_drop_m, "elong_pct": r.max_elongation_pct,
                     "UIAA_ok": r.peak_tension_n < UIAA_LIMIT_N})
    df = pd.DataFrame(scan)
    print(df.to_string(index=False))
    df.to_csv(out_dir / "calibration_ff_scan.csv", index=False)

    mono = df.loc[df["FF"] > 1.0, "peak_kN"].is_monotonic_increasing
    print(f"\n峰值对 FF(>1) 单调递增: {'是' if mono else '否 —— 回归了'}")
    print(f"\n# 建议写入 RopeModel 默认值：")
    print(f"RopeModel(t_ref_n={best.t_ref_N:.0f}, eps_ref={best.eps_ref:.2f}, "
          f"exponent={best.exponent:.2f}, c_damp={best.c_damp:.0f})")


if __name__ == "__main__":
    main()
