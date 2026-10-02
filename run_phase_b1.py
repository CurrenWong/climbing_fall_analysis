"""
run_phase_b1.py — 抱石坠落摔伤分析 第一阶段 主入口

功能：
  1. 运行单场景详细仿真（默认 80kg @ 3m，标准垫）
  2. 能量守恒自检
  3. 多高度峰值对比（2/3/4/5m）
  4. 软垫力-位移曲线
  5. 高度×质量 峰值力热力图（标准垫 + 厚垫）
  6. 汇总报告（控制台 + 保存 CSV/JSON）

用法：
  python run_phase_b1.py            # 默认全套
  python run_phase_b1.py --quick    # 仅单场景 + 时程图
"""

import argparse
import csv
import json
"""⚠️ 已废弃（2026-10-02）—— 请改用 ``scripts/phase_b2_analysis.py``

这是最初的"质点 + 双线性软垫弹簧"模型（Phase B-1 初稿），已被
``src/climbing/pad.py`` 完全取代。当前实现加入了原方案初稿漏掉的三件事：

1. **渐进接触面积**：脚先落只有 ~0.06 m²，塌下去才展开到臀/背，
   力-位移曲线是**时变**的，不是单一线性弹簧。
2. **屈曲缓冲**：软垫之外，肌肉/关节的离心收缩才是第二个吸能器。
   软垫 + 屈曲的**串联行程**才决定峰值力。
3. **头/躯干减速不同步**：同样的垫反力，作用在小质量上加速度远大。

原模型的结论（"厚软垫显著降峰值"）在当前模型下依然成立，但所有绝对
数值都不可用。保留此文件仅作历史参考。
"""

import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(_ROOT))

from src.boulder_fall_1d import (
    PAD_PRESETS, SAFE_THRESHOLDS, SimConfig, assess, energy_check,
    scan_heights, scan_matrix, simulate,
)
from src.visualize_b1 import (
    OUT, plot_height_comparison, plot_matrix_heatmap, plot_pad_force_curve,
    plot_single_timeseries,
)


def print_result(res, tag=""):
    a = assess(res)
    chk = energy_check(res.cfg, res)
    print(f"\n[{tag}] {res.cfg.mass:.0f}kg @ {res.cfg.height:.1f}m  "
          f"(pad={getattr(res, 'pad_name', 'default')})")
    print(f"  撞击速度      : {res.impact_velocity:.2f} m/s")
    print(f"  撞击动能      : {res.kinetic_energy:.0f} J")
    print(f"  峰值反力      : {a['peak_force_kN']:.1f} kN   "
          f"(阈值 {SAFE_THRESHOLDS['peak_force_kN']:.0f})  "
          f"{'✓' if a['checks']['peak_force'] else '✗'}")
    print(f"  峰值减速度    : {a['peak_accel_g']:.1f} g    "
          f"(阈值 {SAFE_THRESHOLDS['peak_accel_g']:.0f})  "
          f"{'✓' if a['checks']['peak_accel'] else '✗'}")
    print(f"  最大压缩      : {a['max_compression_cm']:.1f} cm  "
          f"(上限 {SAFE_THRESHOLDS['max_compression_cm']:.0f})  "
          f"{'✓' if a['checks']['max_compression'] else '✗'}")
    print(f"  吸能占比      : {a['energy_ratio']*100:.1f}%  "
          f"{'✓' if a['checks']['energy_ratio'] else '✗'}")
    print(f"  压穿到底      : {res.bottomed_out}")
    print(f"  能量自检      : rel_err={chk['rel_err']*100:.3f}%  "
          f"{'PASS' if chk['pass'] else 'FAIL'}")
    print(f"  >>> 风险等级  : {a['risk']}")
    return a, chk


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="仅单场景+时程图")
    args = ap.parse_args()

    print("=" * 64)
    print("  抱石坠落摔伤分析 · Phase B-1 (1D 简化模型)")
    print("=" * 64)

    # 1. 单场景详细仿真
    pad_std = PAD_PRESETS["standard"]
    res = simulate(SimConfig(mass=80.0, height=3.0, pad=pad_std, t_max=2.0))
    a, chk = print_result(res, "单场景")
    f_ts = plot_single_timeseries(res)
    print(f"  [图] 时程: {f_ts}")

    if args.quick:
        return

    # 2. 多高度对比
    print("\n--- 多高度对比 (80kg, 标准垫) ---")
    heights = (2.0, 3.0, 4.0, 5.0)
    results_std = scan_heights(heights=heights, pad=pad_std, t_max=2.0)
    for r in results_std:
        print_result(r, f"h={r.cfg.height:.0f}m")
    f_hc = plot_height_comparison(results_std, pad_name="standard")
    print(f"  [图] 高度对比: {f_hc}")

    # 3. 软垫力-位移曲线
    f_pad = plot_pad_force_curve(pad_std, impact_energy_J=res.kinetic_energy,
                                fname_tag="standard")
    print(f"  [图] 软垫曲线: {f_pad}")

    # 4. 高度×质量 热力图（标准垫 + 厚垫）
    print("\n--- 矩阵扫描 (高度×质量×软垫) ---")
    matrix_std = scan_matrix(heights=heights, masses=(60.0, 80.0, 100.0),
                            pad_presets=("standard",), t_max=2.0)
    f_hm_std = plot_matrix_heatmap(matrix_std, pad_name="standard")
    matrix_thick = scan_matrix(heights=heights, masses=(60.0, 80.0, 100.0),
                              pad_presets=("thick",), t_max=2.0)
    f_hm_thick = plot_matrix_heatmap(matrix_thick, pad_name="thick")
    print(f"  [图] 热力图(标准): {f_hm_std}")
    print(f"  [图] 热力图(厚垫): {f_hm_thick}")

    # 5. 汇总 CSV
    _export_csv(matrix_std + matrix_thick)
    _export_json(res, a, chk)
    print(f"\n输出目录: {OUT}")


def _export_csv(results):
    fname = OUT / "summary_matrix.csv"
    with open(fname, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["pad", "mass_kg", "height_m", "impact_v_ms", "KE_J",
                    "peak_force_kN", "peak_accel_g", "max_comp_cm",
                    "energy_ratio", "bottomed", "risk"])
        for r in results:
            a = assess(r)
            w.writerow([getattr(r, "pad_name", "default"),
                        getattr(r, "mass", r.cfg.mass), r.cfg.height,
                        f"{r.impact_velocity:.2f}", f"{r.kinetic_energy:.0f}",
                        f"{a['peak_force_kN']:.1f}", f"{a['peak_accel_g']:.1f}",
                        f"{a['max_compression_cm']:.1f}", f"{a['energy_ratio']:.3f}",
                        r.bottomed_out, a["risk"]])
    print(f"  [数据] {fname}")


def _export_json(res, a, chk):
    fname = OUT / "single_case.json"
    data = {
        "config": {"mass": res.cfg.mass, "height": res.cfg.height,
                   "pad": {"k1": res.cfg.pad.k1, "k2": res.cfg.pad.k2,
                           "x_switch": res.cfg.pad.x_switch,
                           "max_compression": res.cfg.pad.max_compression}},
        "results": {"impact_velocity": res.impact_velocity,
                    "kinetic_energy": res.kinetic_energy,
                    "peak_force_kN": a["peak_force_kN"],
                    "peak_accel_g": a["peak_accel_g"],
                    "max_compression_cm": a["max_compression_cm"],
                    "energy_ratio": a["energy_ratio"],
                    "bottomed_out": res.bottomed_out, "risk": a["risk"]},
        "energy_check": {"rel_err": chk["rel_err"], "pass": chk["pass"]},
    }
    with open(fname, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    print(f"  [数据] {fname}")


if __name__ == "__main__":
    main()
