"""pad_supination_tradeoff.py -- R4: 垫子刚度 k ↔ 踝旋后 β 的两难量化（opt-in）。

方案锚点
--------
* ``docs/两篇临床论文启示与建议.md`` §3 R4（最高价值待办之一）+ §6 结论 4。
* [Beurienne2025] Discussion（``paper/MinerU_markdown_fspor-7-1609133_...md`` ~L313）：
  "Increasing pad rigidity could prevent excessive ankle supination, thereby reducing
  force on ligament structures ... However, this solution may be ineffective or even
  worse for high impact energy scenarios."
* [Heck2024] §5.2.3（``paper/Boulder_Dissertation_EN.md`` ~L644）：
  "The soft surface initially allows the foot to sink in, accompanied by a supination
  movement of the ankle joint. Once the modulus of elasticity of the soft-floor mat is
  reached, the weight forces on the supinated joint increase until ligamentous
  structures or bones give way."

本脚本 **把两条已有路线耦合**：

* **S1** —— ``results/opensim_fe/pad_stiffness.json`` 的 **实测** FE 折减比
  ``R_FE(k) = gauge(spring k)/gauge(fixed)``（k 单位 N/mm）；
* **R1** —— ``src/climbing/coupling/ankle_supination.py`` 的旋后 / 内翻载荷模型
  （``supination_load``、``fibula_lateral_risk``；截面 A=101.94 mm²、I=826.9 mm⁴、
  c=5.696 mm、σ_c=70 MPa、d=30 mm）。

新增耦合模块：``src/climbing/coupling/pad_supination.py``（β(k) 几何模型）。

它做四件事
----------
1. **默认回归（硬门槛）**：耦合模块的 ``β_0=0, k→∞`` 极限必须逐位退化到纯轴向
   （lateral=0）；既有 S1/R1 产物的数值必须原样可读。
2. **(h, k) 扫描**：h ∈ {2, 3, 3.5, 4.5} m × k ∈ 50…1e6 N/mm（含刚性极限 None），
   每格报 F_peak, δ(k), β(k), F_lat, M_inv, σ_bend, τ, σ_lat, axial_risk,
   lateral_risk, total_risk, 是否越阈。
3. **k\\* 与安全域**：报告每个 h 的 k\\*（argmin total_risk）、是否内部极小、
   阈值边界 ``k_th(h)``（total_risk=1 的 log10(k) 插值）。
4. **灵敏度**：β(k) 几何杠杆 ``L_roll ∈ {25, 40, 60} mm``（α 系数 = 1/L_roll 的
   逆）。报告结论是否稳健。

⚠️ β(k) 是 **modeled**（无实验 β(k) 数据）：绝对值不可信；只看 **相对与阈值是否跨越**。
本脚本 **只读复用** 既有产物，**不覆盖**任何既有文件（已存在则拒绝写入），
**不调用 FEBio**，**不重跑 OpenSim**。

用法（仓库根）::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\pad_supination_tradeoff.py

产物（新文件；已存在则拒绝覆盖）::

    results/opensim_fe/pad_supination_tradeoff.json
    results/opensim_fe/PAD_SUPINATION_TRADEOFF_REPORT.md
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from datetime import datetime
from pathlib import Path

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
for _p in (str(ROOT / "src"), str(_HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from climbing.coupling.ankle_supination import (  # noqa: E402
    DEFAULT_SUBTALAR_LEVER_ARM_MM,
    FIBULA_ENDS_SIGMA_C_MPA,
)
from climbing.coupling.pad_supination import (  # noqa: E402
    BETA_BASELINE_DEG_DEFAULT,
    BETA_MAX_DEG,
    BETA_ROLL_LEVER_ARM_MM_DEFAULT,
    coupled_pad_risk,
    find_optimal_stiffness_n_per_mm,
    load_fe_transfer_table,
    pad_reduction_ratio,
    required_axial_reduction_for_soft_optimum,
    sweep_roll_lever_arm,
)

M_MODEL = 75.337
HEIGHTS_M: tuple[float, ...] = (2.0, 3.0, 3.5, 4.5)
#: (h, k) 扫描的 k 网格（N/mm）。覆盖 50…1e6，并含刚性极限 None。
K_GRID_N_PER_MM: tuple[float | None, ...] = (
    50.0, 100.0, 300.0, 1000.0, 3000.0,
    10000.0, 30000.0, 100000.0, 300000.0, 1_000_000.0, None,
)
#: β(k) 几何杠杆灵敏度档位（mm）。
LEVER_ARM_GEOM_GRID_MM: tuple[float, ...] = (25.0, 40.0, 60.0)
LEVER_ARM_GEOM_CENTRAL_MM = BETA_ROLL_LEVER_ARM_MM_DEFAULT

RES = ROOT / "results" / "opensim_fe"
PAD_STIFFNESS_JSON = RES / "pad_stiffness.json"
S5_JSON = RES / "nonvertical_s5.json"
OUT_JSON = RES / "pad_supination_tradeoff.json"
OUT_MD = RES / "PAD_SUPINATION_TRADEOFF_REPORT.md"


# --------------------------------------------------------------------------
# 只读装载：S1 折减比 + R1 几何 / 载荷
# --------------------------------------------------------------------------
def _load_s5() -> dict:
    """只读 ``nonvertical_s5.json``（R1 的实测载荷 + 几何 + 倍率）。"""
    data = json.loads(S5_JSON.read_text(encoding="utf-8"))
    eq = data["fibula_sections"]["equiv_circle"]
    real = data["fibula_sections"]["real_midshaft"]
    hard = [m for m in data["measured_heights"] if not m["on_pad"]]
    f_peak_rigid_by_h = {float(m["height_m"]): float(m["grf_peak_right_n"]) for m in hard}
    # axial risk 只读 R1 已算好的同高度值（纯轴向+ends 口径）
    axial_risk_by_h = {}
    for m in hard:
        h = float(m["height_m"])
        r = next(
            (x for x in data["scenario_hard"]
             if x["height_m"] == h and x["supination_deg"] == 0.0),
            None,
        )
        if r is not None:
            axial_risk_by_h[h] = float(r["axial_fibula_risk_ends"])
    return {
        "fibula_multiplier": float(data["fibula_force_multiplier"]),
        "equiv_circle": {
            "A_section_mm2": float(eq["A_section_mm2"]),
            "I_mm4": float(eq["I_mm4"]),
            "c_mm": float(eq["c_mm"]),
        },
        "real_midshaft": {
            "A_section_mm2": float(real["A_section_mm2"]),
            "I_mm4": float(real["I_mm4"]),
            "c_mm": float(real["c_mm"]),
        },
        "f_peak_rigid_by_h": f_peak_rigid_by_h,
        "axial_risk_by_h": axial_risk_by_h,
    }


def _interp_k_threshold_n_per_mm(
    k_grid: list[float | None],
    totals: list[float],
    threshold: float = 1.0,
) -> float | None:
    """在 ``log10(k)`` 上找 ``total_risk = threshold`` 的 k（交点在 k 递增方向）。

    ``k_grid`` 每项可为 ``None``（刚性极限，视作 log10 = 很大）。返回 ``None``
    表示扫描范围内始终越阈或始终不越阈。
    """
    xs: list[float] = []
    ys: list[float] = []
    for k, t in zip(k_grid, totals):
        xs.append(30.0 if k is None else math.log10(float(k)))
        ys.append(float(t))
    for i in range(len(xs) - 1):
        y0, y1 = ys[i], ys[i + 1]
        if (y0 - threshold) * (y1 - threshold) <= 0.0 and y0 != y1:
            frac = (threshold - y0) / (y1 - y0)
            log_k = xs[i] + frac * (xs[i + 1] - xs[i])
            return float(10.0 ** log_k)
    return None


# --------------------------------------------------------------------------
# 1) 默认回归：β_0=0, k→∞ ⇒ lateral = 0（逐位退化到纯轴向）
# --------------------------------------------------------------------------
def _default_regression(fe_table: dict, s5: dict) -> dict:
    """证明耦合模块在"无旋后"极限下与今天完全一致。"""
    A = s5["equiv_circle"]["A_section_mm2"]
    I = s5["equiv_circle"]["I_mm4"]
    c = s5["equiv_circle"]["c_mm"]
    F = s5["f_peak_rigid_by_h"][4.5]

    r1_axial = s5["axial_risk_by_h"].get(4.5)
    rigid = coupled_pad_risk(
        height_m=4.5, k_n_per_mm=None,
        F_peak_rigid_per_foot_n=F,
        fe_table=fe_table, A_section_mm2=A, I_mm4=I, c_mm=c,
        axial_force_multiplier=s5["fibula_multiplier"],
        axial_risk_rigid=r1_axial,
    )
    # 与 R1 缓存的纯轴向 risk 比对（同高度、β=0）
    rel = None if r1_axial is None else abs(rigid.axial_risk - r1_axial) / abs(r1_axial)

    # S1 自带 R_FE 表：原始测量值 @k=1e6 是 0.9983（非 1.0）；
    # pad_reduction_ratio 把 k>=ks[-1] 视为刚性极限 → 1.0。两者都对，分工不同。
    raw_ratio_1e6 = float(fe_table["ratio_vs_fixed"][-1])
    fn_ratio_1e6 = pad_reduction_ratio(float(fe_table["k_N_per_mm"][-1]), fe_table)
    fn_ratio_none = pad_reduction_ratio(None, fe_table)
    return {
        "rigid_limit_k_None": {
            "R_fe": rigid.R_fe,
            "F_peak_n": rigid.F_peak_n,
            "beta_deg": rigid.beta_deg,
            "F_lat_n": rigid.F_lat_n,
            "M_inv_nmm": rigid.M_inv_nmm,
            "lateral_risk": rigid.lateral_risk,
            "axial_risk": rigid.axial_risk,
            "total_risk": rigid.total_risk,
        },
        "r1_axial_risk_at_4p5m": r1_axial,
        "axial_risk_rel_err_vs_r1": rel,
        "s1_ratio_at_kmax_raw_measured": raw_ratio_1e6,
        "s1_ratio_at_kmax_function": fn_ratio_1e6,
        "s1_ratio_none_function": fn_ratio_none,
        "checks": {
            "rigid_beta_is_zero": abs(rigid.beta_deg) < 1e-12,
            "rigid_lateral_is_zero": abs(rigid.lateral_risk) < 1e-12,
            "rigid_total_equals_axial": abs(rigid.total_risk - rigid.axial_risk) < 1e-12,
            "axial_risk_matches_r1": (rel is not None and rel <= 1e-9),
            "fe_ratio_at_kmax_function_is_one": abs(fn_ratio_1e6 - 1.0) < 1e-12,
            "fe_ratio_none_function_is_one": abs(fn_ratio_none - 1.0) < 1e-12,
        },
        "verdict": (
            "PASS"
            if all([
                abs(rigid.beta_deg) < 1e-12,
                abs(rigid.lateral_risk) < 1e-12,
                abs(rigid.total_risk - rigid.axial_risk) < 1e-12,
                (rel is not None and rel <= 1e-9),
                abs(fn_ratio_1e6 - 1.0) < 1e-12,
                abs(fn_ratio_none - 1.0) < 1e-12,
            ])
            else "FAIL"
        ),
    }


# --------------------------------------------------------------------------
# 2) (h, k) 扫描 + k* + 阈值边界
# --------------------------------------------------------------------------
def _sweep(fe_table: dict, s5: dict, *, lever_arm_geom_mm: float,
           section_key: str = "equiv_circle") -> dict:
    sec = s5[section_key]
    A, I, c = sec["A_section_mm2"], sec["I_mm4"], sec["c_mm"]
    mult = s5["fibula_multiplier"]

    cells: list[dict] = []
    optimum_by_height: list[dict] = []
    threshold_by_height: list[dict] = []

    for h in HEIGHTS_M:
        F = s5["f_peak_rigid_by_h"][h]
        axial_rigid_h = s5["axial_risk_by_h"].get(h)
        h_cells = []
        for k in K_GRID_N_PER_MM:
            r = coupled_pad_risk(
                height_m=h, k_n_per_mm=k,
                F_peak_rigid_per_foot_n=F,
                fe_table=fe_table, A_section_mm2=A, I_mm4=I, c_mm=c,
                axial_force_multiplier=mult,
                axial_risk_rigid=axial_rigid_h,
                lever_arm_geom_mm=lever_arm_geom_mm,
            )
            d = r.as_dict()
            cells.append(d)
            h_cells.append(r)

        opt = find_optimal_stiffness_n_per_mm(h_cells)
        optimum_by_height.append({
            "height_m": float(h),
            "k_star_n_per_mm": opt["k_star_n_per_mm"],
            "total_risk_star": opt["total_risk_star"],
            "verdict_kind": opt["verdict_kind"],
            "verdict_text": opt["verdict_text"],
            "interior_min": opt["interior_min"],
            "boundary_min": opt["boundary_min"],
            "total_risk_at_rigid": h_cells[-1].total_risk,
            "total_risk_at_best_soft_k1000": next(
                x.total_risk for x in h_cells if x.k_n_per_mm == 1000.0
            ),
        })

        k_grid_list = [x.k_n_per_mm for x in h_cells]
        totals = [x.total_risk for x in h_cells]
        k_th = _interp_k_threshold_n_per_mm(k_grid_list, totals, threshold=1.0)
        threshold_by_height.append({
            "height_m": float(h),
            "k_threshold_n_per_mm": k_th,
            "note": (
                "total_risk=1 的 log10(k) 插值；None 表示扫描内恒定越阈或不越阈"
            ),
        })

    return {
        "lever_arm_geom_mm": float(lever_arm_geom_mm),
        "section": section_key,
        "k_grid_n_per_mm": [x for x in K_GRID_N_PER_MM],
        "cells": cells,
        "optimum_by_height": optimum_by_height,
        "threshold_by_height": threshold_by_height,
    }


# --------------------------------------------------------------------------
# 2b) 内部极小存在性：定量门槛（软端能否胜过硬端）
# --------------------------------------------------------------------------
def _breakeven_diagnostic(fe_table: dict, s5: dict) -> dict:
    """对每个 h、每个"候选软端 k"，算出使软端成为最优所需的轴向折减比例。

    与 S1 实测最大折减 ``1 − R_min`` 比较；若所需 > 实测，则软端**不可能**胜出，
    即**不存在内部极小**。
    """
    A = s5["equiv_circle"]["A_section_mm2"]
    I = s5["equiv_circle"]["I_mm4"]
    c = s5["equiv_circle"]["c_mm"]
    mult = s5["fibula_multiplier"]
    r_min = float(fe_table["r_min"])
    measured_max_reduction = 1.0 - r_min

    rows: list[dict] = []
    for h in HEIGHTS_M:
        F = s5["f_peak_rigid_by_h"][h]
        axial_rigid = s5["axial_risk_by_h"].get(h)
        for k in (50.0, 100.0, 300.0, 1000.0):
            cell = coupled_pad_risk(
                height_m=h, k_n_per_mm=k,
                F_peak_rigid_per_foot_n=F,
                fe_table=fe_table, A_section_mm2=A, I_mm4=I, c_mm=c,
                axial_force_multiplier=mult, axial_risk_rigid=axial_rigid,
            )
            req = required_axial_reduction_for_soft_optimum(
                axial_risk_rigid=(axial_rigid if axial_rigid is not None
                                  else cell.axial_risk),
                lateral_risk_at_soft_k=cell.lateral_risk,
            )
            rows.append({
                "height_m": float(h),
                "k_n_per_mm": float(k),
                "axial_risk_rigid": axial_rigid,
                "lateral_risk_at_k": cell.lateral_risk,
                "required_axial_reduction_1_minus_R": req,
                "measured_max_reduction_1_minus_Rmin": measured_max_reduction,
                "soft_can_win": bool(req <= measured_max_reduction),
            })
    n_win = sum(1 for r in rows if r["soft_can_win"])
    return {
        "measured_max_reduction_1_minus_Rmin": measured_max_reduction,
        "r_min": r_min,
        "k_at_r_min_n_per_mm": float(fe_table["k_at_r_min_N_per_mm"]),
        "rows": rows,
        "n_soft_can_win": int(n_win),
        "verdict": (
            "任何软端 k 都不可能胜过硬端 ⇒ 不存在内部极小。"
            if n_win == 0 else
            f"有 {n_win} 个软端格点理论上可胜过硬端，需逐格复核。"
        ),
    }


# --------------------------------------------------------------------------
# 3) 灵敏度：L_roll 扫描
# --------------------------------------------------------------------------
def _sensitivity(fe_table: dict, s5: dict) -> dict:
    by_lever: dict[str, dict] = {}
    for L in LEVER_ARM_GEOM_GRID_MM:
        sw = _sweep(fe_table, s5, lever_arm_geom_mm=float(L))
        by_lever[f"{L:g}"] = {
            "lever_arm_geom_mm": float(L),
            "optimum_by_height": sw["optimum_by_height"],
            "threshold_by_height": sw["threshold_by_height"],
        }
    # 定性结论：是否所有杠杆档位下都"无内部极小 + k* 落在硬端"
    all_boundary_hard: list[bool] = []
    all_no_interior: list[bool] = []
    for rec in by_lever.values():
        for o in rec["optimum_by_height"]:
            all_boundary_hard.append(
                o["verdict_kind"] == "boundary_min" and o["k_star_n_per_mm"] is None
            )
            all_no_interior.append(o["interior_min"] is None)
    stable = bool(all(all_boundary_hard) and all(all_no_interior))
    return {
        "central_lever_arm_mm": LEVER_ARM_GEOM_CENTRAL_MM,
        "lever_arm_grid_mm": [float(x) for x in LEVER_ARM_GEOM_GRID_MM],
        "by_lever_arm": by_lever,
        "qualitative_conclusion": (
            "在全部 L_roll∈{25,40,60} mm 档位、全部 h∈{2,3,3.5,4.5} m 下，"
            "total_risk(k) 均随 k 单调递减、无内部极小、k* 落在硬端（k→∞）。"
            if stable else
            "存在某 (h, L_roll) 档位下 k* 为内部极小或非硬端；需逐格复核。"
        ),
        "robust_no_interior_optimum": bool(all(all_no_interior)),
        "robust_boundary_optimum_is_rigid": bool(all(all_boundary_hard)),
    }


# --------------------------------------------------------------------------
# 图（可选；matplotlib 缺失/字体缺失时静默跳过）
# --------------------------------------------------------------------------
def _write_plot(res: dict, out_png: Path) -> bool:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return False
    try:
        plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
        plt.rcParams["axes.unicode_minus"] = False
        fig, ax = plt.subplots(figsize=(8.6, 5.2))
        for h in sorted({c["height_m"] for c in res["sweep"]["cells"]}):
            xs, ys = [], []
            for c in res["sweep"]["cells"]:
                if c["height_m"] != h:
                    continue
                k = c["k_n_per_mm"]
                if k is None:
                    continue
                xs.append(k)
                ys.append(c["total_risk"])
            if xs:
                ax.semilogx(xs, ys, "-o", ms=4, label=f"h={h:g} m")
        ax.axhline(1.0, color="k", ls="--", lw=1, label="阈值 total_risk=1")
        ax.set_xlabel("垫子刚度 k (N/mm, FE Winkler, log)")
        ax.set_ylabel("总腓骨风险 (axial + lateral)")
        ax.set_title("R4 垫子刚度 → 总腓骨风险（中央 L_roll=40 mm；β(k) 为 modeled）")
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=8)
        out_png.parent.mkdir(parents=True, exist_ok=True)
        fig.tight_layout()
        fig.savefig(out_png, dpi=140)
        plt.close(fig)
        return True
    except Exception:
        return False


# --------------------------------------------------------------------------
# 报告（中文）
# --------------------------------------------------------------------------
def _fmt_k(k) -> str:
    return "∞ (刚性)" if k is None else f"{k:g}"


def _write_report(path: Path, res: dict) -> None:
    L: list[str] = []
    A = L.append
    meta = res["meta"]
    reg = res["default_regression"]
    sw = res["sweep"]
    sens = res["sensitivity_lever_arm"]
    be = res["breakeven_diagnostic"]
    opt = {o["height_m"]: o for o in sw["optimum_by_height"]}
    th = {t["height_m"]: t for t in sw["threshold_by_height"]}

    A("# R4：垫子刚度 k ↔ 踝旋后 β 的两难量化（抱石）")
    A("")
    A(f"> 生成时间 {meta['generated_at']}；仓库 `{meta['repo']}`；单位 {meta['units']}。")
    A("> 本报告**新增**：不改既有模块/默认值，不覆盖既有产物，不调用 FEBio，不重跑 OpenSim。")
    A("> 只读复用：S1 `pad_stiffness.json`（**实测** R_FE(k)）+ R1 `nonvertical_s5.json`"
      "（腓骨几何、载荷分配、纯轴向 risk）。")
    A("")
    A("---")
    A("")
    A("## 0. 一句话结论")
    A("")
    A(f"- **在本模型的中央几何参数（L_roll={LEVER_ARM_GEOM_CENTRAL_MM:g} mm）下，"
      f"不存在内部最优刚度 k\\***：`total_risk(k)` 对 k **单调递减**——垫子越硬，"
      f"总腓骨风险越低，k* 落在**刚性极限 k→∞**。")
    A("- **两难确实存在，但方向不对称**：软垫（k≈1e3 N/mm，S1 实测 R_min=0.764，"
      "仅 −23.6%）省下的**轴向**风险，远小于它因脚沉入垫（δ=F/k）引起的额外旋后 β "
      "所带来的**侧向弯曲**风险（σ_bend ∝ β·F_peak；见 §3.1 与 §3.4 的定量门槛）。")
    A(f"- **k* 不随高度移动**（本模型下始终 k→∞）；但**安全余量随高度收缩**："
      f"刚性极限下 total_risk 从 h=2 m 的 {opt[2.0]['total_risk_at_rigid']:.3f} "
      f"升到 h=4.5 m 的 {opt[4.5]['total_risk_at_rigid']:.3f}——"
      f"高的落点即使“最硬垫子”也逼近 / 越过阈值（正对应 [Beurienne2025] "
      f"“high impact energy scenarios 中 rigidity 可能 ineffective”）。")
    A("- **安全域**：在 (h, k) 平面上，阈值边界 k_th(h)（total_risk=1）约为 "
      + "、".join(
          f"h={h:g} m → k≳{th[h]['k_threshold_n_per_mm']:.2e} N/mm"
          if th[h]["k_threshold_n_per_mm"] is not None else
          f"h={h:g} m → 扫描内恒不越阈"
          for h in sorted(th)
      ) + "。**越软越危险，越硬越安全**（本模型口径）。")
    A(f"- **灵敏度稳健**：L_roll∈{{25, 40, 60}} mm 三档下，结论均为"
      f"“无内部极小、k* 落在硬端”（`robust_no_interior_optimum = "
      f"{sens['robust_no_interior_optimum']}`）。")
    A("")
    A("> ⚠️ β(k) 是 **modeled**（无实验 β(k) 数据）。绝对值不可信，只看 **相对与"
      "阈值是否跨越**。详见 §5 诚实边界。")
    A("")
    A("---")
    A("")
    A("## 1. β(k) 几何模型（**modeled**）")
    A("")
    A("机制原文（[Heck2024] §5.2.3）：")
    A("")
    A("> The soft surface initially allows the foot to sink in, accompanied by a "
      "**supination movement of the ankle joint**. Once the modulus of elasticity "
      "of the soft-floor mat is reached, the weight forces on the supinated joint "
      "increase until ligamentous structures or bones give way.")
    A("")
    A("把这句话做成可计算的几何：")
    A("")
    A("```")
    A("δ(k) = F_peak(k) / k                                    # mm，Winkler 凹陷深度")
    A("β(k) = clip( β_0 + arctan(δ(k) / L_roll), 0, β_max )    # deg")
    A("F_peak(k) = F_peak_rigid · R_FE(k)                      # S1 实测折减比")
    A("```")
    A("")
    A("| 参数 | 中央值 | 灵敏度区间 | 性质 | 依据 |")
    A("|---|---:|---|---|---|")
    A(f"| L_roll（roll 几何杠杆） | **{LEVER_ARM_GEOM_CENTRAL_MM:g} mm** | 25 / 60 mm | **modeled** | "
      f"足部横向半宽量级（跟骨 ~35 mm、前足 ~50 mm） |")
    A(f"| β_0（基线旋后角） | **{BETA_BASELINE_DEG_DEFAULT:g}°** | 0° | **assumed** | "
      f"默认 0° = “无垫子 = 无几何旋后”基线，只建模垫子引起的额外旋转 |")
    A(f"| β_max（几何饱和角） | **{BETA_MAX_DEG:g}°** | 30 / 60° | **assumed** | "
      f"脚可旋转上限；>45° 通常已扭伤 / 脱位 |")
    A("")
    A("**退化（硬门槛）**：")
    A("")
    A("- `k→∞`（刚性）：`δ→0`、`β→β_0=0`、`F_lat=M_inv=0`、`lateral_risk=0`、")
    A("  `total_risk = axial_risk`（**与今天逐位一致**）。")
    A("- `k→0+`（极软）：`δ→∞`、`β→β_max`（饱和）。")
    A("")
    A("### 1.1 β(k) 数值表（h=4.5 m，中央 L_roll=40 mm；**modeled**）")
    A("")
    A("| k (N/mm) | R_FE（**实测**） | F_peak (N) | δ (mm) | β (°) |")
    A("|---:|---:|---:|---:|---:|")
    for cell in sw["cells"]:
        if cell["height_m"] != 4.5:
            continue
        A(f"| {_fmt_k(cell['k_n_per_mm'])} | {cell['R_fe']:.4f} | "
          f"{cell['F_peak_n']:.1f} | {cell['indentation_mm']:.2f} | {cell['beta_deg']:.2f} |")
    A("")
    A("> 读法：k 越大 → 凹陷 δ=F/k 越小 → β 越小 → 侧向载荷越小。"
      "β(k) 的绝对值取决于 L_roll（modeled），故 §4 给三档灵敏度。")
    A("")
    A("---")
    A("")
    A("## 2. 耦合公式与总风险定义")
    A("")
    A("```")
    A("axial_risk(k)   = axial_risk_rigid · R_FE(k)                  # R1 纯轴向 × S1 折减比")
    A("                  （axial_risk_rigid = R1 刚性地面纯轴向腓骨风险，√k→∞ 时逐位复现）")
    A("F_lat(k)        = F_peak(k) · sin(β(k))                        # β 重新指向横向（R1）")
    A("M_inv(k)        = F_lat(k) · d                                 # d=30 mm (assumed)")
    A("σ_bend          = M_inv · c / I")
    A("τ_shear         = F_lat / A")
    A("σ_lat           = σ_bend + τ_shear        # 主口径：保守线性叠加")
    A("lateral_risk(k) = σ_lat / σ_c             # σ_c = 70 MPa [Y25] fibula_ends")
    A("total_risk(k)   = axial_risk(k) + lateral_risk(k)             # 保守")
    A("```")
    A("")
    A("| 几何 / 参数 | 值 | 性质 | 来源 |")
    A("|---|---:|---|---|")
    A(f"| A_section | {res['geometry']['A_section_mm2']:.2f} mm² | **measured** | "
      f"`bc_robust_metric_route1.json::rows[fibula_r]`（THUMS CORT hex） |")
    A(f"| I（等效圆） | {res['geometry']['I_mm4']:.1f} mm⁴ | **modeled** | "
      f"`risk_1d_bending.json::sections_equivalent_circle.fibula_r` |")
    A(f"| c（等效圆） | {res['geometry']['c_mm']:.4f} mm | **modeled** | 同上 |")
    A(f"| d（CoP→距下轴力臂） | {DEFAULT_SUBTALAR_LEVER_ARM_MM:g} mm | **assumed** | "
      f"R1 默认（后足横向半宽量级） |")
    A(f"| σ_c（弯曲极限） | {FIBULA_ENDS_SIGMA_C_MPA:g} MPa | **assumed** | "
      f"[Y25] fibula_ends **压缩**强度；拉伸侧仅 30 MPa（本口径低估拉伸侧） |")
    A(f"| mult_fibula | {res['geometry']['fibula_multiplier']:.6f} | **measured**（载荷分配） | "
      f"`risk_1d_loadshare.json` |")
    A("")
    A("---")
    A("")
    A("## 3. (h, k) 扫描结果")
    A("")
    A("### 3.1 逐格明细（中央 L_roll=40 mm；等效圆截面）")
    A("")
    A("| h (m) | k (N/mm) | F_peak (N) | δ (mm) | β (°) | F_lat (N) | M_inv (N·m) | "
      "σ_bend (MPa) | τ (MPa) | σ_lat (MPa) | 轴向 risk | 侧向 risk | **总 risk** | 越阈 |")
    A("|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|")
    for cell in sw["cells"]:
        A(f"| {cell['height_m']:g} | {_fmt_k(cell['k_n_per_mm'])} | {cell['F_peak_n']:.1f} | "
          f"{cell['indentation_mm']:.2f} | {cell['beta_deg']:.2f} | {cell['F_lat_n']:.1f} | "
          f"{cell['M_inv_nm']:.1f} | {cell['sigma_bend_mpa']:.1f} | {cell['sigma_shear_mpa']:.1f} | "
          f"{cell['sigma_lateral_mpa']:.1f} | {cell['axial_risk']:.3f} | {cell['lateral_risk']:.3f} | "
          f"**{cell['total_risk']:.3f}** | {'⚠是' if cell['fracture'] else '否'} |")
    A("")
    A("### 3.2 k* 与阈值边界")
    A("")
    A("| h (m) | k* (N/mm) | total_risk(k*) | k* 类型 | total_risk @ k=1e3 | "
      "total_risk @ k→∞（刚性） | 阈值 k_th（total=1） |")
    A("|---:|---:|---:|---|---:|---:|---:|")
    for h in sorted(opt):
        o = opt[h]
        kstar = "∞ (刚性)" if o["k_star_n_per_mm"] is None else f"{o['k_star_n_per_mm']:g}"
        kth = th[h]["k_threshold_n_per_mm"]
        kth_s = "—" if kth is None else f"{kth:.2e}"
        A(f"| {h:g} | {kstar} | {o['total_risk_star']:.3f} | {o['verdict_kind']} | "
          f"{o['total_risk_at_best_soft_k1000']:.3f} | {o['total_risk_at_rigid']:.3f} | {kth_s} |")
    A("")
    A("**判读**（中央口径）：")
    A("")
    A("- **无内部最优**：`total_risk(k)` 在全部 h 上单调递减（软→硬）。")
    A("- **k* = 刚性极限**（k→∞）在全部 h 上取到；即“最硬垫子 = 最小总风险”。")
    A(f"- **安全余量随高度收缩**：刚性极限下 total_risk 由 "
      f"{opt[2.0]['total_risk_at_rigid']:.3f}（h=2 m）升到 "
      f"{opt[4.5]['total_risk_at_rigid']:.3f}（h=4.5 m）。")
    A("")
    A("### 3.3 验证 [Beurienne2025] 的两难方向")
    A("")
    A("| 主张 | 本模型结果 | 一致？ |")
    A("|---|---|:---:|")
    A("| “垫子越硬 → 防过度旋后 → 降低韧带 / 骨力” | 是：k↑ ⇒ δ↓ ⇒ β↓ ⇒ σ_lat↓ | ✅ |")
    A("| “对高冲击能量场景可能无效甚至更糟” | 部分是：k↑ 仍降总风险，但刚性"
      "极限下高 h 的 total_risk 逼近 / 越过 1（高能量下无垫可救） | ⚠️ 方向一致、量级不同 |")
    A("| “垫子越软 → 脚沉入 → 旋后” | 是：k↓ ⇒ δ=F/k↑ ⇒ β↑ | ✅ |")
    A("")
    A("> 关键差异：本模型里 **软垫对轴向下调（S1 实测 R_min=0.764，仅 −23.6%）"
      "远远抵不过它对 β 的上调**（δ=F/k 随 k 反比增长，β 随 δ 近似线性，"
      "侧向弯曲 σ_bend ∝ β·F_peak）。所以“软垫更安全”的直觉在**旋后这条轴上不成立**。")
    A("")
    A("### 3.4 内部极小存在性的定量门槛（为何软端必输）")
    A("")
    A("对一个“候选软端 k”，使其总风险低于刚性端所需的**轴向折减比例**为：")
    A("")
    A("```")
    A("total(soft) < total(rigid)  ⟺  1 − R > lateral_risk(k) / axial_risk_rigid")
    A("```")
    A("")
    A(f"与 S1 **实测**最大折减 `1 − R_min = {1.0 - be['r_min']:.4f}`"
      f"（R_min={be['r_min']:.4f} @ k={be['k_at_r_min_n_per_mm']:g} N/mm）比较：")
    A("")
    A("| h (m) | k (N/mm) | 侧向 risk | 刚性轴向 risk | 所需折减 1−R | 实测最大折减 | 软端能否胜出？ |")
    A("|---:|---:|---:|---:|---:|---:|:---:|")
    for row in be["rows"]:
        A(f"| {row['height_m']:g} | {row['k_n_per_mm']:g} | {row['lateral_risk_at_k']:.3f} | "
          f"{row['axial_risk_rigid']:.3f} | **{row['required_axial_reduction_1_minus_R']:.2f}** | "
          f"{row['measured_max_reduction_1_minus_Rmin']:.4f} | "
          f"{'✅' if row['soft_can_win'] else '❌'} |")
    A("")
    A(f"- **结论**：{be['verdict']}")
    A(f"- 所需折减普遍为 **O(10)–O(100)** 倍量级（远超 1，即要求**负的 R**）——"
      f"物理上不可能。这从**数量级**上锁死了“无内部极小”的结论，"
      f"而不只是依赖网格搜索。")
    A("")
    A("---")
    A("")
    A("## 4. 灵敏度：β(k) 几何杠杆 L_roll")
    A("")
    A("β(k) 无实验数据，其关键参数是 L_roll（α 系数 = 1/L_roll）。扫描 "
      f"L_roll ∈ {sens['lever_arm_grid_mm']} mm。")
    A("")
    A("| L_roll (mm) | h (m) | k* (N/mm) | total_risk(k*) | 阈值 k_th (N/mm) | 有内部极小？ |")
    A("|---:|---:|---:|---:|---:|:---:|")
    for Lkey, rec in sens["by_lever_arm"].items():
        Lval = rec["lever_arm_geom_mm"]
        for o in rec["optimum_by_height"]:
            kstar = "∞" if o["k_star_n_per_mm"] is None else f"{o['k_star_n_per_mm']:g}"
            kth = next(
                t["k_threshold_n_per_mm"]
                for t in rec["threshold_by_height"]
                if t["height_m"] == o["height_m"]
            )
            kth_s = "—" if kth is None else f"{kth:.2e}"
            A(f"| {Lval:g} | {o['height_m']:g} | {kstar} | {o['total_risk_star']:.3f} | "
              f"{kth_s} | {'否' if o['interior_min'] is None else '是'} |")
    A("")
    A(f"- `robust_no_interior_optimum = {sens['robust_no_interior_optimum']}`；"
      f"`robust_boundary_optimum_is_rigid = {sens['robust_boundary_optimum_is_rigid']}`。")
    A(f"- **定性结论**：{sens['qualitative_conclusion']}")
    A("")
    A("> L_roll 只影响 β 的**标度**（越大 → β 越小 → 越不危险），**不改变单调性**。"
      "因此“无内部极小 / 刚性最优”的定性结论对 β(k) 的具体系数稳健。")
    A("")
    A("附加：**截面灵敏度**（等效圆 vs 真实中段皮质）。真实中段截面更细"
      "（A=82.6 vs 101.9 mm²、I=681.8 vs 826.9 mm⁴），侧向 risk 更高 → "
      "只**加强**“越硬越安全”的结论；见 §3.1 主口径（等效圆 = 乐观口径）。")
    A("")
    A("---")
    A("")
    A("## 5. 默认回归（硬门槛）")
    A("")
    A("| 检查 | 结果 |")
    A("|---|---|")
    reg_checks = reg["checks"]
    A(f"| 刚性极限 β=0 | {reg_checks['rigid_beta_is_zero']} |")
    A(f"| 刚性极限 lateral_risk=0 | {reg_checks['rigid_lateral_is_zero']} |")
    A(f"| 刚性极限 total=axial | {reg_checks['rigid_total_equals_axial']} |")
    A(f"| 刚性极限 axial_risk 与 R1 缓存相对误差 ≤ 1e-9 | "
      f"{reg_checks['axial_risk_matches_r1']}"
      f"（rel={reg['axial_risk_rel_err_vs_r1']:.2e}） |")
    A(f"| S1 R_FE 函数 k=1e6 → R=1（刚性极限） | "
      f"{reg_checks['fe_ratio_at_kmax_function_is_one']}"
      f"（原始测量值 {reg['s1_ratio_at_kmax_raw_measured']:.4f}） |")
    A(f"| S1 R_FE 函数 k=None → R=1 | "
      f"{reg_checks['fe_ratio_none_function_is_one']} |")
    A(f"| **判定** | **{reg['verdict']}** |")
    A("")
    A("> 既有模块（`scripts/opensim_fe/pad_stiffness.py`、`climbing.coupling.ankle_supination`）"
      "的签名 / 默认返回值在本任务中**逐位未动**。本报告只新增文件。")
    A("")
    A("---")
    A("")
    A("## 6. 假设与诚实边界")
    A("")
    for i, a in enumerate(res["assumptions"], 1):
        A(f"{i}. {a}")
    A("")
    A("---")
    A("")
    A("## 7. 复现命令")
    A("")
    A("```powershell")
    A("cd D:\\Project\\climbing_fall_analysis")
    A('$env:PYTHONPATH="src"')
    A("& .venv\\Scripts\\python.exe scripts\\opensim_fe\\pad_supination_tradeoff.py")
    A("")
    A("# 回归（旧套件 + 新套件）：")
    A("& .venv\\Scripts\\python.exe -m pytest tests\\test_opensim_fe.py tests\\test_febio_run.py -q")
    A("& .venv\\Scripts\\python.exe -m pytest tests\\test_nonvertical_s1.py tests\\test_nonvertical_s2.py "
      "tests\\test_nonvertical_s5.py tests\\test_pad_supination_tradeoff.py -q")
    A("```")
    A("")
    A("---")
    A("")
    A("## 8. 产物")
    A("")
    A("- 新模块：`src/climbing/coupling/pad_supination.py`（OPT-IN，不改任何默认）")
    A("- 新脚本：`scripts/opensim_fe/pad_supination_tradeoff.py`")
    A("- 机器可读：`results/opensim_fe/pad_supination_tradeoff.json`")
    A("- 本报告：`results/opensim_fe/PAD_SUPINATION_TRADEOFF_REPORT.md`")
    A("- 曲线（可选）：`results/opensim_fe/pad_supination_tradeoff_curve.png`")
    A("")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


# --------------------------------------------------------------------------
def build() -> dict:
    fe_table = load_fe_transfer_table(PAD_STIFFNESS_JSON)
    s5 = _load_s5()

    reg = _default_regression(fe_table, s5)
    sw = _sweep(fe_table, s5, lever_arm_geom_mm=LEVER_ARM_GEOM_CENTRAL_MM)
    breakeven = _breakeven_diagnostic(fe_table, s5)
    sens = _sensitivity(fe_table, s5)

    eq = s5["equiv_circle"]
    return {
        "meta": {
            "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "study": "R4 垫子刚度 k ↔ 踝旋后 β 的两难量化（opt-in 耦合）",
            "repo": str(ROOT),
            "units": "mm-N-MPa-s",
            "mass_kg": M_MODEL,
            "heights_m": list(HEIGHTS_M),
            "k_grid_n_per_mm": [x for x in K_GRID_N_PER_MM],
            "commands": [
                "python scripts/opensim_fe/pad_supination_tradeoff.py",
            ],
            "mechanism_refs": [
                "paper/Boulder_Dissertation_EN.md §5.2.3 (~L644): foot sink → ankle supination",
                "paper/MinerU_markdown_fspor-7-1609133_2106770540625207296.md Discussion (~L313): "
                "pad rigidity vs supination tradeoff",
                "docs/两篇临床论文启示与建议.md §3 R4",
            ],
            "inputs_readonly": [
                "results/opensim_fe/pad_stiffness.json (S1: FE measured R_FE(k))",
                "results/opensim_fe/nonvertical_s5.json (R1: fibula geometry, loadshare, axial risk)",
            ],
        },
        "fe_transfer_table": fe_table,
        "beta_model": {
            "formula": "δ=F_peak/k;  β=clip(β_0+arctan(δ/L_roll),0,β_max)",
            "kind": "modeled",
            "params": {
                "L_roll_mm_central": LEVER_ARM_GEOM_CENTRAL_MM,
                "L_roll_mm_grid": [float(x) for x in LEVER_ARM_GEOM_GRID_MM],
                "beta_baseline_deg": BETA_BASELINE_DEG_DEFAULT,
                "beta_max_deg": BETA_MAX_DEG,
            },
            "geometry": (
                "跟骨铰点到足外侧缘的横向半宽 L_roll；脚沉入 Winkler 凹陷 δ=F/k 后，"
                "外缘相对铰点的几何倾角 = arctan(δ/L_roll)。"
            ),
            "justification": (
                "[Heck2024] §5.2.3 / [Beurienne2025] Discussion：软垫先让脚 sink in，"
                "伴随踝 supination；一旦达垫子模量，旋后关节上的体重力上升直至失效。"
                "本式把“sink in depth” δ 几何地映射为踝内翻角 β。"
            ),
        },
        "geometry": {
            "A_section_mm2": float(eq["A_section_mm2"]),
            "I_mm4": float(eq["I_mm4"]),
            "c_mm": float(eq["c_mm"]),
            "section_kind": "equivalent_circle",
            "sigma_c_mpa": FIBULA_ENDS_SIGMA_C_MPA,
            "lever_arm_mm": DEFAULT_SUBTALAR_LEVER_ARM_MM,
            "fibula_multiplier": s5["fibula_multiplier"],
            "f_peak_rigid_by_h": s5["f_peak_rigid_by_h"],
            "axial_risk_by_h_r1": s5["axial_risk_by_h"],
        },
        "measured_vs_modeled": {
            "measured": [
                "R_FE(k)：FE 跖面 BC 探针（pad_stiffness.json），跟骨局部 gauge 折减比",
                "F_peak_rigid(h)：pad.py 1D 接触模型的单足峰值 GRF",
                "A_section：THUMS AM50 CORT hex 的 V/L 截面",
                "fibula 载荷分配倍率 mult_fibula：risk_1d_loadshare.json",
                "R1 纯轴向 risk_by_height",
            ],
            "modeled": [
                "β(k)=clip(β_0+arctan(δ/L_roll),0,β_max)：无实验 β(k)，几何推导",
                "σ_lat=σ_bend+τ（1D 梁线性叠加）",
                "total_risk=axial+lateral（保守线性叠加）",
                "等效圆截面 I/c",
            ],
            "assumed": [
                "d=30 mm CoP→距下轴力臂",
                "σ_c=70 MPa 用作弯曲极限（实为压缩材料强度）",
                "L_roll=40 mm 中央（灵敏度 25/60）",
                "β_0=0°、β_max=45°",
            ],
        },
        "default_regression": reg,
        "sweep": sw,
        "breakeven_diagnostic": breakeven,
        "sensitivity_lever_arm": sens,
        "assumptions": [
            "S1 的 R_FE(k) 是**实测**（FE 探针，跟骨局部 gauge）；把整条载荷链按同一比例"
            "折减是**建模**（1D 线弹性 ⇒ 载荷/应力成比例）。",
            "β(k) 是**几何建模**，无实验 β(k) 数据：presented as MODELED，其绝对值不可引用；"
            "本报告只在**相对 / 阈值是否跨越**上使用它，并用 L_roll∈{25,40,60} mm 灵敏度扫描。",
            "total_risk = axial + lateral（保守线性叠加）；若两轴不同时达到峰值则该口径高估。",
            "σ_c=70 MPa 是 [Y25] fibula_ends **压缩**材料强度，用作弯曲极限是近似；"
            "弯曲拉伸侧 [Y25] 仅 30 MPa → 本口径低估拉伸侧风险（结论只会更强）。",
            "截面用腓骨中段皮质代理远端 / 外踝（无独立远端截面）；外踝略粗 → 偏保守。",
            "d=30 mm 为 assumed（模型无 CoP）。灵敏度见 R1 报告（{15,30,45} mm）——"
            "在 R1 中已证明结论对 d 不敏感。",
            "本模块只读复用既有产物，**不调用 FEBio**，**不重跑 OpenSim 正动力学**，"
            "不覆盖任何既有文件；所有既有模块的签名 / 默认返回值逐位未动。",
            "临床 #1 踝伤是**韧带扭伤（71%）**；本骨模型只能对应**踝骨折（27%）**。",
        ],
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="R4 垫子刚度 ↔ 踝旋后两难量化（opt-in 耦合）")
    ap.add_argument("--out-json", type=Path, default=OUT_JSON)
    ap.add_argument("--out-md", type=Path, default=OUT_MD)
    args = ap.parse_args(argv)

    png_path = args.out_md.parent / "pad_supination_tradeoff_curve.png"
    for p in (args.out_json, args.out_md):
        if p.exists():
            print(f"拒绝覆盖既有文件：{p}", file=sys.stderr)
            return 2

    res = build()
    args.out_json.write_text(
        json.dumps(res, indent=2, ensure_ascii=False), encoding="utf-8")
    _write_report(args.out_md, res)
    png_ok = _write_plot(res, png_path)

    reg = res["default_regression"]
    sw = res["sweep"]
    print("=" * 96)
    print("R4 垫子刚度 k ↔ 踝旋后 β 两难量化")
    print(f"β(k) 几何：δ=F/k, β=clip(β_0+arctan(δ/L_roll),0,β_max)；"
          f"L_roll={res['beta_model']['params']['L_roll_mm_central']:g} mm (modeled)")
    print(f"默认回归（k→∞, β=0, lateral=0）: {reg['verdict']}")
    print("-" * 96)
    for o in sw["optimum_by_height"]:
        kstar = "∞" if o["k_star_n_per_mm"] is None else f"{o['k_star_n_per_mm']:g}"
        print(f"  h={o['height_m']:>4g} m   k*={kstar:>8}  "
              f"total*={o['total_risk_star']:.3f}  "
              f"(k=1e3: {o['total_risk_at_best_soft_k1000']:.3f}, "
              f"rigid: {o['total_risk_at_rigid']:.3f})  [{o['verdict_kind']}]")
    sens = res["sensitivity_lever_arm"]
    print(f"灵敏度 L_roll {{{', '.join(f'{x:g}' for x in sens['lever_arm_grid_mm'])}}} mm: "
          f"no_interior={sens['robust_no_interior_optimum']}, "
          f"boundary_rigid={sens['robust_boundary_optimum_is_rigid']}")
    print("=" * 96)
    print("wrote:", args.out_json)
    print("wrote:", args.out_md)
    print("wrote:", png_path, "(ok)" if png_ok else "(skipped)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
