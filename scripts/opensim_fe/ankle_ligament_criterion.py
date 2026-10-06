"""R6 probe —— 踝韧带扭伤（sprain）判据（**opt-in**，默认零）。

方案锚点
--------
* ``docs/踝韧带判据方案.md``（R6 方案本体）。
* ``docs/两篇临床论文启示与建议.md`` §2 结构缺口 + §7 待办（韧带判据）。
* [Heck2024]（``paper/Boulder_Dissertation_EN.md``）：踝伤中扭伤 **71.3%** vs 骨折 27.4%。
* [Beurienne2025]（``paper/MinerU_markdown_fspor-7-1609133_...md``）：踝扭伤 #1（28%）。
* [Tochigi2006]（``paper/MinerU_markdown_2006FAIankleligaments_...md``）：正常步态
  ATFL 应变 ≤ 6.2 %、CFL ≤ 2.1 %；"taut" 阈值 1.0 %。

本脚本把 **R1 旋后/内翻载荷模型**（``ankle_supination``，内翻力矩 ``M_inv``）
接到 **R6 韧带判据**（``ankle_ligament``）上，回答："在什么 (h, β) 格点，
ATFL / CFL 会（按哪种判据）扭伤？"

它做五件事
----------
1. **默认回归（硬门槛）**：新模块在 β=0（``M_inv=0``）时全部 risk 为 0；既有
   R1 缓存 ``nonvertical_s5.json`` 的 ``scenario_hard`` 逐格可复现（相对误差 0）；
   既有 ``ankle_supination`` 行为不变（只读复用，无改动）。
2. **场景 (h, β)**：h ∈ {2, 3, 3.5, 4.5} m × β ∈ {0, 5, 10, 15, 30}°，对每格报
   ``M_inv``、``β``、ATFL 力 / 应变，及 ATFL / CFL 是否扭伤（关节级 + 力/应变）。
3. **两个预载口径**：(a) 保守 / 无预载 ``M_fail=21 N·m``；(b) 高承载 / 2 kN 预载
   ``M_fail=77 N·m``（Funk 2002）。并**明确标注**：抱石冲击轴向 ≈ 21 kN，
   预载外推**不确定**。
4. **r_lig 灵敏度**：20 / 22 / 25 mm（**assumed** 力臂）。
5. **临床对照**：我们的韧带扭伤预测 vs [Heck2024] 踝伤 71% 扭伤 / 27% 骨折。

⚠️ 诚实边界：ATFL 力/应变判据把**全部** ``M_inv`` 归给单条韧带（无骨/其他结构
分担），故是 **上界式**；关节级判据是主口径。绝对值不可引用，只看相对与是否越阈。

只读复用既有产物，**不覆盖**既有文件（已存在则拒绝写入），**不调用 FEBio**。

用法（仓库根）::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\ankle_ligament_criterion.py

产物（新文件；已存在则拒绝覆盖）::

    results/opensim_fe/ankle_ligament_criterion.json
    results/opensim_fe/ANKLE_LIGAMENT_CRITERION_REPORT.md
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
for _p in (str(ROOT / "src"), str(_HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from climbing.coupling.ankle_supination import (  # noqa: E402
    frontal_lateral_force_n,
    inversion_moment_nmm,
    supination_load,
)
from climbing.coupling.ankle_ligament import (  # noqa: E402
    ATFL_FAILURE_STRAIN,
    ATFL_INVERSION_FAIL_ANGLE_2KN_DEG,
    ATFL_INVERSION_FAIL_ANGLE_DEG,
    ATFL_INVERSION_FAIL_MOMENT_2KN_PRELOAD_NM,
    ATFL_INVERSION_FAIL_MOMENT_NO_PRELOAD_HIGH_NM,
    ATFL_INVERSION_FAIL_MOMENT_NO_PRELOAD_NM,
    ATFL_PEAK_FORCE_STANCE_N,
    ATFL_PEAK_FORCE_STANCE_STD_N,
    ATFL_RESTING_LENGTH_MM,
    ATFL_STIFFNESS_N_PER_MM,
    ATFL_ULTIMATE_LOAD_N,
    CFL_FAILURE_STRAIN,
    CFL_INVERSION_FAIL_ANGLE_DEG,
    CFL_ULTIMATE_LOAD_N,
    DEFAULT_LIGAMENT_LEVER_ARM_MM,
    R_LIG_SENSITIVITY_MM,
    TOCHIGI_ATFL_NORMAL_STRAIN_MAX,
    TOCHIGI_CFL_NORMAL_STRAIN_MAX,
    TOCHIGI_TAUT_STRAIN,
    ankle_ligament_criterion,
    atfl_sprain_risk,
    cfl_sprain_risk,
    critical_inversion_angle_deg,
    critical_inversion_moment_nmm,
    preload_capacity_interp,
    sprain_risk,
)

M_MODEL = 75.337
HEIGHTS_M: tuple[float, ...] = (2.0, 3.0, 3.5, 4.5)
ANGLES_DEG: tuple[float, ...] = (0.0, 5.0, 10.0, 15.0, 30.0)
LEVER_ARM_PRIMARY_MM = 30.0  # CoP→距下轴 force arm（R1 主口径）

#: 预载口径标签 → 中文说明
PRELOAD_LABELS = {
    "conservative": "保守口径（无预载，M_fail=21 N·m）",
    "elevated": "高承载口径（2 kN 预载，M_fail=77 N·m）",
}

RES = ROOT / "results" / "opensim_fe"
S5_JSON = RES / "nonvertical_s5.json"       # R1 缓存（只读）
OUT_JSON = RES / "ankle_ligament_criterion.json"
OUT_MD = RES / "ANKLE_LIGAMENT_CRITERION_REPORT.md"


# --------------------------------------------------------------------------
# 1) 默认回归（硬门槛）
# --------------------------------------------------------------------------
def _default_regression() -> dict:
    """证明：新模块默认路径零风险 + 既有 R1 缓存可逐格复现（rel err 0）。"""
    checks: dict[str, bool] = {}

    # (a) 新模块默认：M_inv=0 → 全部 risk=0（逐位复现旧行为）
    a = atfl_sprain_risk(0.0)
    c = cfl_sprain_risk(0.0)
    fs_a = sprain_risk(0.0, ligament="ATFL")
    fs_c = sprain_risk(0.0, ligament="CFL")
    checks["atfl_joint_zero_moment_risk_is_zero"] = bool(a.risk == 0.0 and not a.sprain)
    checks["cfl_joint_zero_moment_risk_is_zero"] = bool(c.risk == 0.0 and not c.sprain)
    checks["atfl_force_strain_zero_is_zero"] = bool(
        fs_a["F_lig_n"] == 0.0 and fs_a["strain"] == 0.0 and not fs_a["sprain"]
    )
    checks["cfl_force_strain_zero_is_zero"] = bool(
        fs_c["F_lig_n"] == 0.0 and fs_c["strain"] == 0.0 and not fs_c["sprain"]
    )

    # (b) 既有 ankle_supination 行为不变：默认 supination_load 仍为零
    sl = supination_load((0.0, 21000.0, 0.0))
    checks["supination_load_default_zero"] = bool(sl.is_zero)

    # (c) R1 缓存 scenario_hard 逐格复现（rel err 0）：重算 F_lat / M_inv
    s5 = json.loads(S5_JSON.read_text(encoding="utf-8"))
    max_rel = 0.0
    n_rows = 0
    for r in s5["scenario_hard"]:
        if r["supination_deg"] == 0.0:
            continue
        fvec = (0.0, float(r["F_vert_n"]), 0.0)
        fl = abs(frontal_lateral_force_n(fvec, float(r["supination_deg"])))
        m = inversion_moment_nmm(fl, LEVER_ARM_PRIMARY_MM)
        for got, exp in ((fl, r["F_lat_n"]), (m, r["M_inv_nmm"])):
            denom = abs(exp) if exp != 0.0 else 1.0
            max_rel = max(max_rel, abs(got - exp) / denom)
        n_rows += 1
    checks["r1_cached_scenario_reproduced"] = bool(max_rel <= 1e-12)

    # (d) 逐高度单足峰值 GRF：fresh re-sim（pad.py，纯 Python）vs R1 缓存
    from climbing.coupling.opensim_grf import ground_reaction

    grf_rel = 0.0
    for m in s5["measured_heights"]:
        if m.get("on_pad", True):
            continue
        h = float(m["height_m"])
        g = ground_reaction(height_m=h, mass_kg=M_MODEL)
        f_right = float(np.max(g.f_right_n))
        exp = float(m["grf_peak_right_n"])
        grf_rel = max(grf_rel, abs(f_right - exp) / exp)
    checks["fresh_grf_matches_r1_cached"] = bool(grf_rel <= 1e-9)

    passed = all(checks.values())
    return {
        "checks": checks,
        "n_nonzero_beta_rows_checked": n_rows,
        "scenario_max_rel_err": float(max_rel),
        "grf_max_rel_err": float(grf_rel),
        "verdict": "PASS" if passed else "FAIL",
        "note": (
            "新模块 ankle_ligament 是纯增量（opt-in）；既有 ankle_supination / "
            "opensim_grf 无任何改动，默认路径逐位复现 R1 缓存。"
        ),
    }


# --------------------------------------------------------------------------
# 2) 场景 (h, β)
# --------------------------------------------------------------------------
def _scenario(measured_heights: list[dict]) -> list[dict]:
    """(h, β) × 两预载口径 × ATFL/CFL 双判据。"""
    by_h = {float(m["height_m"]): m for m in measured_heights
            if not m.get("on_pad", True)}
    rows: list[dict] = []
    for h in HEIGHTS_M:
        m = by_h.get(float(h))
        if m is None:
            continue
        f_vert = float(m["grf_peak_right_n"])
        fvec = (0.0, f_vert, 0.0)
        for beta in ANGLES_DEG:
            fl = abs(frontal_lateral_force_n(fvec, float(beta)))
            m_inv = inversion_moment_nmm(fl, LEVER_ARM_PRIMARY_MM)
            for pm in ("conservative", "elevated"):
                atfl = atfl_sprain_risk(m_inv, preload_model=pm)
                cfl = cfl_sprain_risk(m_inv, preload_model=pm)
                fa = sprain_risk(m_inv, ligament="ATFL",
                                 r_lig_mm=DEFAULT_LIGAMENT_LEVER_ARM_MM)
                fc = sprain_risk(m_inv, ligament="CFL",
                                 r_lig_mm=DEFAULT_LIGAMENT_LEVER_ARM_MM)
                rows.append({
                    "height_m": float(h),
                    "supination_deg": float(beta),
                    "preload_model": pm,
                    "F_vert_n": f_vert,
                    "F_lat_n": float(fl),
                    "M_inv_nmm": float(m_inv),
                    "M_inv_nm": float(m_inv) / 1000.0,
                    "atfl_joint_risk": float(atfl.risk),
                    "atfl_joint_sprain": bool(atfl.sprain),
                    "cfl_joint_risk": float(cfl.risk),
                    "cfl_joint_sprain": bool(cfl.sprain),
                    "atfl_force_n": float(fa["F_lig_n"]),
                    "atfl_strain": float(fa["strain"]),
                    "atfl_force_risk": float(fa["risk_force"]),
                    "atfl_strain_risk": float(fa["risk_strain"]),
                    "atfl_fe_sprain": bool(fa["sprain"]),
                    "cfl_force_n": float(fc["F_lig_n"]),
                    "cfl_strain": float(fc["strain"]),
                    "cfl_fe_sprain": bool(fc["sprain"]),
                })
    return rows


# --------------------------------------------------------------------------
# 3) r_lig 灵敏度
# --------------------------------------------------------------------------
def _r_lig_sensitivity(measured_heights: list[dict]) -> dict:
    by_h = {float(m["height_m"]): m for m in measured_heights
            if not m.get("on_pad", True)}
    grid = (DEFAULT_LIGAMENT_LEVER_ARM_MM,) + tuple(
        r for r in R_LIG_SENSITIVITY_MM if abs(r - DEFAULT_LIGAMENT_LEVER_ARM_MM) > 1e-12
    )
    out: dict[str, list[dict]] = {f"{r:g}": [] for r in grid}
    for h in HEIGHTS_M:
        m = by_h.get(float(h))
        if m is None:
            continue
        f_vert = float(m["grf_peak_right_n"])
        fvec = (0.0, f_vert, 0.0)
        for beta in ANGLES_DEG:
            fl = abs(frontal_lateral_force_n(fvec, float(beta)))
            m_inv = inversion_moment_nmm(fl, LEVER_ARM_PRIMARY_MM)
            for r in grid:
                fs = sprain_risk(m_inv, ligament="ATFL", r_lig_mm=float(r))
                out[f"{r:g}"].append({
                    "height_m": float(h),
                    "supination_deg": float(beta),
                    "r_lig_mm": float(r),
                    "M_inv_nm": float(m_inv) / 1000.0,
                    "atfl_force_n": float(fs["F_lig_n"]),
                    "atfl_strain": float(fs["strain"]),
                    "atfl_force_risk": float(fs["risk_force"]),
                    "atfl_strain_risk": float(fs["risk_strain"]),
                    "atfl_fe_sprain": bool(fs["sprain"]),
                })
    return out


# --------------------------------------------------------------------------
# 4) Tochigi 参考带（正常步态 vs 我们的极端载荷）
# --------------------------------------------------------------------------
def _tochigi_band(measured_heights: list[dict]) -> dict:
    """把 [Tochigi2006] 的正常步态应变上限，和本模型在 β=5° 的应变对照。

    **只作参照**（Tochigi 是矢状面步态、600 N 轴向；本模型是冠状面冲击）。
    """
    by_h = {float(m["height_m"]): m for m in measured_heights
            if not m.get("on_pad", True)}
    ref_beta = 5.0
    rows = []
    for h in HEIGHTS_M:
        m = by_h.get(float(h))
        if m is None:
            continue
        f_vert = float(m["grf_peak_right_n"])
        fvec = (0.0, f_vert, 0.0)
        fl = abs(frontal_lateral_force_n(fvec, ref_beta))
        m_inv = inversion_moment_nmm(fl, LEVER_ARM_PRIMARY_MM)
        fs_a = sprain_risk(m_inv, ligament="ATFL")
        rows.append({
            "height_m": float(h),
            "ref_beta_deg": ref_beta,
            "M_inv_nm": m_inv / 1000.0,
            "atfl_strain": float(fs_a["strain"]),
            "atfl_strain_over_tochigi_normal": float(
                fs_a["strain"] / TOCHIGI_ATFL_NORMAL_STRAIN_MAX
            ),
        })
    return {
        "tochigi_atfl_normal_max": TOCHIGI_ATFL_NORMAL_STRAIN_MAX,
        "tochigi_cfl_normal_max": TOCHIGI_CFL_NORMAL_STRAIN_MAX,
        "tochigi_taut_threshold": TOCHIGI_TAUT_STRAIN,
        "reference_beta_deg": ref_beta,
        "rows": rows,
        "caveat": (
            "Tochigi 是矢状面步态 600 N 轴向的**正常**应变；本模型是冠状面冲击的"
            "**极端**载荷。两者不同轴、不同载荷量级 → 只作'非生理参照带'，"
            "不作为失效阈值。"
        ),
    }


# --------------------------------------------------------------------------
# 5) 临床对照
# --------------------------------------------------------------------------
def _clinical_comparison(rows: list[dict]) -> dict:
    cons = [r for r in rows if r["preload_model"] == "conservative"]
    cons_sprain = [r for r in cons if r["atfl_joint_sprain"]]
    return {
        "heck2024_ankle_sprain_share_pct": 71.3,
        "heck2024_ankle_fracture_share_pct": 27.4,
        "beurienne2025_ankle_sprain_rank": "#1（28%）",
        "our_conservative_sprain_cells": len(cons_sprain),
        "our_conservative_total_cells": len(cons),
        "note": (
            "骨模型只覆盖踝骨折 27%；本韧带模块补上踝扭伤 71% 这条主线。"
            "在保守口径下，凡 β≥5° 的格点 ATFL 关节级即越阈 → 定性上与"
            "'扭伤是踝伤 #1' 一致。"
        ),
        "bone_only_axial_reference": {
            "fibula_first_fracture_h": 7.0,
            "note": "R1 纯轴向腓骨首骨折 ≈ 7 m（真实高度全安全）；补旋后后 2–4.5 m 即越阈。",
        },
    }


# --------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="R6 踝韧带扭伤判据 probe（opt-in）")
    ap.add_argument("--out-json", type=Path, default=OUT_JSON)
    ap.add_argument("--out-md", type=Path, default=OUT_MD)
    args = ap.parse_args(argv)

    for p in (args.out_json, args.out_md):
        if p.exists():
            print(f"拒绝覆盖既有文件：{p}", file=sys.stderr)
            return 2

    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("=== R6 踝韧带扭伤（sprain）判据 ===", flush=True)

    print("[1/4] 默认回归（硬门槛）...", flush=True)
    reg = _default_regression()
    print(f"      verdict={reg['verdict']}  scenario_rel_err={reg['scenario_max_rel_err']:.2e}  "
          f"grf_rel_err={reg['grf_max_rel_err']:.2e}", flush=True)

    s5 = json.loads(S5_JSON.read_text(encoding="utf-8"))
    measured = s5["measured_heights"]

    print("[2/4] 场景 (h × β × 预载口径)...", flush=True)
    rows = _scenario(measured)
    print(f"      {len(rows)} 格（{len(HEIGHTS_M)}h × {len(ANGLES_DEG)}β × 2 口径）", flush=True)

    print("[3/4] r_lig 灵敏度 + Tochigi 参照带...", flush=True)
    sens = _r_lig_sensitivity(measured)
    band = _tochigi_band(measured)
    clinical = _clinical_comparison(rows)

    # 组合 API 单格演示（(h=2m, β=15°) 的 M_inv）
    demo_row = next(
        r for r in rows
        if r["preload_model"] == "conservative"
        and r["height_m"] == 2.0 and r["supination_deg"] == 15.0
    )
    api_demo = ankle_ligament_criterion(
        demo_row["M_inv_nmm"], r_lig_mm=DEFAULT_LIGAMENT_LEVER_ARM_MM,
    )

    # 预载承载力曲线
    preload_curve = [
        {"axial_preload_n": p, "M_fail_nm": preload_capacity_interp(p),
         "M_fail_nm_linear_extrap": preload_capacity_interp(p, extrapolation="linear")}
        for p in (0.0, 500.0, 1000.0, 2000.0, 5000.0, 10000.0, 21000.0)
    ]

    # 临界量
    fvec_ref = (0.0, 21000.0, 0.0)
    critical = {
        "joint_moment_thresholds_nm": {
            "conservative_no_preload": (
                critical_inversion_moment_nmm(preload_model="conservative") / 1000.0
            ),
            "conservative_no_preload_high": ATFL_INVERSION_FAIL_MOMENT_NO_PRELOAD_HIGH_NM,
            "elevated_2kn_preload": (
                critical_inversion_moment_nmm(preload_model="elevated") / 1000.0
            ),
        },
        "joint_angle_thresholds_deg": {
            "funk2002_no_preload": ATFL_INVERSION_FAIL_ANGLE_DEG,
            "parenteau1998_no_preload": 34.3,
            "funk2002_2kn_preload": ATFL_INVERSION_FAIL_ANGLE_2KN_DEG,
            "cfl_combined": CFL_INVERSION_FAIL_ANGLE_DEG,
        },
        "our_model_critical_angle_deg_at_21kn": {
            "conservative": critical_inversion_angle_deg(fvec_ref, preload_model="conservative"),
            "elevated": critical_inversion_angle_deg(fvec_ref, preload_model="elevated"),
        },
        "force_thresholds_n": {
            "atfl_ultimate": ATFL_ULTIMATE_LOAD_N,
            "cfl_ultimate": CFL_ULTIMATE_LOAD_N,
            "atfl_stance_peak_bahr1998": ATFL_PEAK_FORCE_STANCE_N,
            "atfl_stance_peak_std": ATFL_PEAK_FORCE_STANCE_STD_N,
        },
        "force_strain_criterion_moment_threshold_nm_at_default_rlig": (
            ATFL_ULTIMATE_LOAD_N * DEFAULT_LIGAMENT_LEVER_ARM_MM / 1000.0
        ),
    }

    result = {
        "meta": {
            "generated_at": ts,
            "study": "R6 踝韧带扭伤（sprain）判据 —— 对接 R1 旋后/内翻载荷",
            "plan_doc": "docs/踝韧带判据方案.md + docs/两篇临床论文启示与建议.md §7",
            "clinic_refs": [
                "paper/Boulder_Dissertation_EN.md ([Heck2024]/[Müller2022]): 踝伤 扭伤 71.3% / 骨折 27.4%",
                "paper/MinerU_markdown_fspor-7-1609133_2106770540625207296.md ([Beurienne2025]): 踝扭伤 #1 28%",
                "paper/MinerU_markdown_2006FAIankleligaments_2106701085136199680.md ([Tochigi2006])",
            ],
            "units": "mm-N-MPa-s",
            "mass_kg": M_MODEL,
            "heights_m": list(HEIGHTS_M),
            "angles_deg": list(ANGLES_DEG),
            "lever_arm_primary_mm": LEVER_ARM_PRIMARY_MM,
            "r_lig_default_mm": DEFAULT_LIGAMENT_LEVER_ARM_MM,
            "r_lig_sensitivity_mm": list(R_LIG_SENSITIVITY_MM),
            "preload_models": PRELOAD_LABELS,
            "commands": [
                "python scripts/opensim_fe/ankle_ligament_criterion.py",
                "python -m pytest tests/test_ankle_ligament_criterion.py",
            ],
            "measured_vs_modeled": {
                "measured": [
                    "地面反力 GRF(t) 由 pad.py 1D 接触模型产出（measured，R1 复用）",
                    "R1 缓存 nonvertical_s5.json 的 F_vert / F_lat / M_inv（只读，逐格复现）",
                    "韧带极限力学参数（ATFL/CFL）：极限载荷、失效应变、刚度、长度、CSA（尸体文献）",
                    "关节级内翻失效力矩 / 角（Funk 2002 / Parenteau 1998，尸体）",
                    "Tochigi 2006 正常步态韧带应变上限（尸体）",
                    "Bahr 1998 步态期 ATFL 峰值力（在体）",
                ],
                "modeled_assumed": [
                    "关节级判据 risk = M_inv / M_fail（主口径）",
                    "预载承载力线性插值 preload_capacity_interp（0→2 kN 锚点）",
                    "力/应变判据 F_lig = M_inv / r_lig；ε ≈ F_lig/(k·L)（1D 线弹性）",
                    "r_lig = 22 mm（**assumed**；灵敏度 20/25 mm）",
                    "把全部 M_inv 归给单条韧带（无骨/其他结构分担）→ 上界式",
                ],
            },
        },
        "default_regression": reg,
        "ligament_params": {
            "ATFL": {
                "ultimate_load_n": ATFL_ULTIMATE_LOAD_N,
                "failure_strain": ATFL_FAILURE_STRAIN,
                "stiffness_n_per_mm": ATFL_STIFFNESS_N_PER_MM,
                "resting_length_mm": ATFL_RESTING_LENGTH_MM,
                "csa_mm2": 13.0,
                "E_mpa": 260.0,
                "stance_peak_force_n": ATFL_PEAK_FORCE_STANCE_N,
                "stance_peak_force_std_n": ATFL_PEAK_FORCE_STANCE_STD_N,
                "sources": "St.Pierre 1983; Chapman 2019; Siegler 1988; Bahr 1998",
            },
            "CFL": {
                "ultimate_load_n": CFL_ULTIMATE_LOAD_N,
                "failure_strain": CFL_FAILURE_STRAIN,
                "stiffness_n_per_mm": 12.0,
                "resting_length_mm": 28.5,
                "csa_mm2": 10.0,
                "E_mpa": 510.0,
                "sources": "St.Pierre 1983; Chapman 2019; Siegler 1988",
            },
        },
        "preload_capacity": {
            "anchors": {
                "no_preload_nm": ATFL_INVERSION_FAIL_MOMENT_NO_PRELOAD_NM,
                "preload_2kn_nm": ATFL_INVERSION_FAIL_MOMENT_2KN_PRELOAD_NM,
            },
            "curve": preload_curve,
            "uncertainty_note": (
                "抱石冲击轴向 ≈ 21 kN 是 2 kN 锚点的 10.5×；文献仅在 0–2 kN 有数据，"
                ">2 kN 的外推**不确定**。clamp 口径给出下界 77 N·m；linear 外推给出上界。"
            ),
        },
        "critical_thresholds": critical,
        "api_demo_single_cell": {
            "note": "ankle_ligament_criterion() 组合 API 在 (h=2m, β=15°) 的单格输出摘要",
            "M_inv_nmm": float(api_demo["m_inv_nmm"]),
            "r_lig_mm": float(api_demo["r_lig_mm"]),
            "atfl_joint_risk": float(api_demo["atfl_joint"].risk),
            "atfl_joint_sprain": bool(api_demo["atfl_joint"].sprain),
            "atfl_F_lig_n": float(api_demo["atfl_force_strain"]["F_lig_n"]),
            "atfl_strain": float(api_demo["atfl_force_strain"]["strain"]),
            "atfl_fe_sprain": bool(api_demo["atfl_force_strain"]["sprain"]),
            "r_lig_sensitivity_keys": sorted(
                float(k) for k in api_demo["r_lig_sensitivity"]
            ),
        },
        "scenario": {
            "grid_heights_m": list(HEIGHTS_M),
            "grid_angles_deg": list(ANGLES_DEG),
            "rows": rows,
        },
        "r_lig_sensitivity": sens,
        "tochigi_reference_band": band,
        "clinical_comparison": clinical,
        "assumptions": [
            "新模块 ankle_ligament 全 opt-in；默认（β=0 / M_inv=0）全部 risk=0，"
            "既有 ankle_supination / opensim_grf / pad_supination 行为逐位不变。",
            "关节级判据（M_inv/M_fail）是**主口径**——避开了 r_lig 力臂假设。",
            "r_lig = 22 mm 是 **assumed**（无实测韧带力臂）；灵敏度 20/25 mm。",
            "力/应变判据把**全部** M_inv 归给单条韧带（无骨/其他韧带分担）→ **上界式**，"
            "其力阈值对应的关节力矩（4.4 N·m）远低于关节级 21 N·m —— 二者差 ≈ 5×，"
            "正说明现实中载荷由多结构分担。",
            "预载承载力对抱石冲击轴向（≈21 kN）的外推**不确定**：同时给保守（无预载）"
            "与高承载（2 kN）两口径，不假装单一阈值。",
            "失效阈值来自尸体文献（Funk/Parenteau/Siegler/Bahr），非本仓库实测；"
            "[Tochigi2006] 只给正常步态应变带，**不是**失效阈值。",
            "临床 71% 是踝**伤**内占比（急诊选择偏倚）；与我们的相对/定性对照即可。",
        ],
    }

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(
        json.dumps(result, indent=2, ensure_ascii=False, default=float), encoding="utf-8"
    )
    _write_report(args.out_md, result, ts)
    print(f"wrote: {args.out_json}", flush=True)
    print(f"wrote: {args.out_md}", flush=True)
    return 0 if reg["verdict"] == "PASS" else 1


# --------------------------------------------------------------------------
def _bool_cn(v: bool) -> str:
    return "是" if v else "否"


def _write_report(path: Path, res: dict, ts: str) -> None:
    reg = res["default_regression"]
    crit = res["critical_thresholds"]
    rows = res["scenario"]["rows"]
    cons = [r for r in rows if r["preload_model"] == "conservative"]
    elev = [r for r in rows if r["preload_model"] == "elevated"]
    sens = res["r_lig_sensitivity"]
    band = res["tochigi_reference_band"]
    clin = res["clinical_comparison"]
    pre = res["preload_capacity"]

    def _sprain_cells(rs, key):
        return sorted({(r["height_m"], r["supination_deg"]) for r in rs if r[key]})

    def _cells_str(cells):
        if not cells:
            return "无"
        return "、".join(f"({h:g}m, {b:g}°)" for h, b in cells)

    L: list[str] = []
    A = L.append
    A("# R6 踝韧带扭伤（sprain）判据报告")
    A("")
    A(f"> 生成时间 {ts}；仓库 `D:\\Project\\climbing_fall_analysis`；单位 mm-N-MPa-s。")
    A("> 方案 `docs/踝韧带判据方案.md` + `docs/两篇临床论文启示与建议.md` §7；"
      "对接 R1 `results/opensim_fe/NONVERTICAL_S5_SUPINATION_REPORT.md`。")
    A("> 参考论文：[Tochigi2006]（`paper/MinerU_markdown_2006FAIankleligaments_...md`）"
      "+ 失效参数文献（Funk 2002 / Parenteau 1998 / St.Pierre 1983 / Siegler 1988 / Bahr 1998）。")
    A("> 本报告**新增**，不改任何模块默认、不覆盖既有产物、**不调用 FEBio**。")
    A("")
    A("---")
    A("")
    A("## 0. 一句话结论")
    A("")
    cons_j = _sprain_cells(cons, "atfl_joint_sprain")
    elev_j = _sprain_cells(elev, "atfl_joint_sprain")
    cons_fe = _sprain_cells(cons, "atfl_fe_sprain")
    A(f"1. **默认回归（硬门槛）: {reg['verdict']}** —— 新模块在 β=0（`M_inv=0`）时"
      f"全部 risk=0；R1 缓存 `scenario_hard` 逐格复现相对误差 "
      f"**{reg['scenario_max_rel_err']:.2e}**；fresh GRF 与 R1 缓存相对误差 "
      f"**{reg['grf_max_rel_err']:.2e}**。")
    A("2. **API 全 opt-in**：新模块 `ankle_ligament` 纯增量；既有模块无改动。")
    A(f"3. **关节级 ATFL 扭伤**（主口径，`M_inv/M_fail`）：保守口径（无预载 21 N·m）下"
      f" **β≥5° 即全部越阈**：{_cells_str(cons_j) if cons_j else '—'}；"
      f"高承载口径（2 kN 预载 77 N·m）下 **β≥10° 越阈**：{_cells_str(elev_j)}。")
    A(f"4. **力/应变 ATFL 扭伤**（`r_lig=22 mm`，上界式）：**β≥5° 全部越阈**："
      f"{_cells_str(cons_fe)}（其阈值关节力矩仅 "
      f"{crit['force_strain_criterion_moment_threshold_nm_at_default_rlig']:.1f} N·m）。")
    A(f"5. **临床对照**：模型在保守口径下 β≥5° 即预测 ATFL 扭伤 → 定性对齐 "
      f"[Heck2024] 踝伤**扭伤 71.3%** / 骨折 27.4%，补上了骨模型（仅 27%）漏掉的#1 主线。")
    A("")
    A("> ⚠️ **两个预载口径都报**：抱石冲击轴向 ≈ 21 kN，远超文献 0–2 kN 范围，"
      "预载外推**不确定**。**不假装单一阈值。** 力/应变判据是**上界式**"
      "（全部 M_inv 归单条韧带）。**绝对值不可引用，只看相对与是否越阈。**")
    A("")
    A("---")
    A("")
    A("## 1. 新 API（全部 opt-in、默认复现旧行为）")
    A("")
    A("| 符号 | 参数（默认） | 单位 | 作用 | 类别 |")
    A("|---|---|---|---|---|")
    A("| `atfl_sprain_risk(M_inv, preload_model='conservative', axial_preload_n=0, m_fail_nm=None)` | 关节级 | N·mm→risk | ATFL `risk=M_inv/M_fail` | 关节级（主口径） |")
    A("| `cfl_sprain_risk(...)` | 同 ATFL | N·mm→risk | CFL（同口径保守下限） | 关节级 |")
    A("| `preload_capacity_interp(P, extrapolation='clamp')` | P (N) | N·m | 轴向预载 P → M_fail | 模型 |")
    A("| `ligament_force_n(M_inv, r_lig_mm)` | r_lig | N | `F_lig = M_inv/r_lig` | 力/应变 |")
    A("| `ligament_strain(F, stiffness, length)` | k, L | – | `ε = F/(k·L)` | 力/应变 |")
    A("| `sprain_risk(M_inv, ligament='ATFL', r_lig_mm=22)` | ligament | – | 力/应变双判据 | 力/应变 |")
    A("| `critical_inversion_moment_nmm(ligament='ATFL', preload_model=...)` | – | N·mm | 首次扭伤的最小 M_inv | 关节级 |")
    A("| `critical_inversion_angle_deg(force_vec, preload_model=...)` | – | deg | 首次扭伤的最小 β | 关节级 |")
    A("| `ankle_ligament_criterion(M_inv, r_lig_mm=22, preload_model=...)` | – | dict | 单格联合判据 | 综合 |")
    A("")
    A("**关键默认**：`preload_model='conservative'`、`r_lig_mm=22`、`M_inv=0` ⇒ 全部"
      " risk=0、`sprain=False`（旧行为）。")
    A("")
    A("---")
    A("")
    A("## 2. 默认回归（硬门槛）")
    A("")
    A("| 检查 | 结果 |")
    A("|---|---|")
    for k, v in reg["checks"].items():
        A(f"| `{k}` | **{_bool_cn(v)}** |")
    A(f"| 非零 β 的 R1 缓存格点数 | {reg['n_nonzero_beta_rows_checked']} |")
    A(f"| 逐格最大相对误差 | **{reg['scenario_max_rel_err']:.2e}** |")
    A(f"| fresh GRF vs 缓存最大相对误差 | **{reg['grf_max_rel_err']:.2e}** |")
    A(f"| **判定** | **{reg['verdict']}** |")
    A("")
    A(f"> {reg['note']}")
    A("")
    A("---")
    A("")
    A("## 3. 韧带参数（measured，文献）")
    A("")
    A("| 参数 | ATFL | CFL | 来源 |")
    A("|---|---:|---:|---|")
    A(f"| 极限载荷 | **{ATFL_ULTIMATE_LOAD_N:g} N**（140–350） | **{CFL_ULTIMATE_LOAD_N:g} N**（300–370） | St.Pierre 1983; Chapman 2019; Siegler 1988 |")
    A(f"| 失效应变 | **{ATFL_FAILURE_STRAIN*100:.0f}%**（11–15） | **{CFL_FAILURE_STRAIN*100:.0f}%**（9–13） | Siegler 1988; Chapman 2019 |")
    A(f"| 割线刚度 | {ATFL_STIFFNESS_N_PER_MM:g}±4 N/mm | 12±4 N/mm | 近期共识 |")
    A(f"| 静息长度 | {ATFL_RESTING_LENGTH_MM:g} mm（18–25） | 28.5 mm（25–32） | Siegler 1988; Golanó |")
    A("| CSA | 13 mm² | 10 mm² | Siegler 1988 |")
    A("| E | 260 MPa | 510 MPa | Siegler 1988 |")
    A(f"| 步态期峰值力 | {ATFL_PEAK_FORCE_STANCE_N:g}±{ATFL_PEAK_FORCE_STANCE_STD_N:g} N | — | Bahr 1998（安全系数 ≈ 2） |")
    A("")
    A("> ⚠️ 用户列表中的 **Woo 1983 是 ACL（膝）**，**不用于踝**；踝韧带初级数值源是 "
      "Siegler 1988。本题眼下的 `r_lig` 是 **assumed**。")
    A("")
    A("---")
    A("")
    A("## 4. 预载依赖的内翻承载力（**两口径**）")
    A("")
    A("`[Funk 2002]` 给出两个锚点：**0 N 预载 → 21 N·m**、**2 kN 预载 → 77 N·m**；"
      "`preload_capacity_interp` 在 `[0, 2 kN]` 线性插值。")
    A("")
    A("| 轴向预载 P (N) | M_fail (N·m)（clamp） | M_fail (N·m)（linear 外推） |")
    A("|---:|---:|---:|")
    for c in pre["curve"]:
        A(f"| {c['axial_preload_n']:,.0f} | {c['M_fail_nm']:.1f} | {c['M_fail_nm_linear_extrap']:.1f} |")
    A("")
    A(f"> **不确定性**：{pre['uncertainty_note']}")
    A("")
    A("**阈值汇总**：")
    A("")
    A("| 口径 | M_fail (N·m) | 失效角 (°) |")
    A("|---|---:|---:|")
    A(f"| 无预载 · Funk 2002 | {crit['joint_moment_thresholds_nm']['conservative_no_preload']:.1f} | {crit['joint_angle_thresholds_deg']['funk2002_no_preload']:.0f} |")
    A(f"| 无预载 · Parenteau 1998 | {crit['joint_moment_thresholds_nm']['conservative_no_preload_high']:.1f} | {crit['joint_angle_thresholds_deg']['parenteau1998_no_preload']:.1f} |")
    A(f"| 2 kN 预载 · Funk 2002 | {crit['joint_moment_thresholds_nm']['elevated_2kn_preload']:.0f} | {crit['joint_angle_thresholds_deg']['funk2002_2kn_preload']:.0f} |")
    A(f"| CFL（组合载荷，晚于 ATFL） | — | {crit['joint_angle_thresholds_deg']['cfl_combined']:.0f} |")
    A("")
    A(f"**本模型在 |F_vert|≈21 kN 下的临界旋后角**（由 GRF 驱动，非尸体测试角）："
      f"保守口径 **{crit['our_model_critical_angle_deg_at_21kn']['conservative']:.2f}°**、"
      f"高承载口径 **{crit['our_model_critical_angle_deg_at_21kn']['elevated']:.2f}°**。")
    A("")
    A("> 两类角不同：Funk/Parenteau 的 33–40° 是**尸体测试的失效角**；上表的 2.6°/9.6° "
      "是**本模型在 21 kN 冲击 GRF 下**首次越阈的输入旋后角。后者更小，因为 21 kN "
      "远大于尸体测试加载。**不可混用。**")
    A("")
    A("---")
    A("")
    A("## 5. 场景：h × β（刚性地面）")
    A("")
    A("F_vert = 单足峰值 GRF（R1 measured）；`F_lat=F_vert·sinβ`；`M_inv=F_lat·d`（d=30 mm）。")
    A("`r_lig=22 mm`。ATFL 力阈值 = 200 N、应变阈值 = 14 %。")
    A("")
    A("### 5.1 保守口径（无预载，M_fail=21 N·m）")
    A("")
    A("| h (m) | β (°) | F_lat (N) | M_inv (N·m) | ATFL risk (关节) | ATFL 扭伤(关节) | ATFL F (N) | ATFL ε | ATFL 扭伤(力/应变) | CFL 扭伤(关节) |")
    A("|---:|---:|---:|---:|---:|:--:|---:|---:|:--:|:--:|")
    for r in cons:
        A(f"| {r['height_m']:g} | {r['supination_deg']:g} | {r['F_lat_n']:,.0f} | "
          f"{r['M_inv_nm']:.1f} | {r['atfl_joint_risk']:.2f} | {_bool_cn(r['atfl_joint_sprain'])} | "
          f"{r['atfl_force_n']:,.0f} | {r['atfl_strain']*100:.1f}% | "
          f"{_bool_cn(r['atfl_fe_sprain'])} | {_bool_cn(r['cfl_joint_sprain'])} |")
    A("")
    A("### 5.2 高承载口径（2 kN 预载，M_fail=77 N·m）")
    A("")
    A("| h (m) | β (°) | F_lat (N) | M_inv (N·m) | ATFL risk (关节) | ATFL 扭伤(关节) | ATFL F (N) | ATFL ε | ATFL 扭伤(力/应变) | CFL 扭伤(关节) |")
    A("|---:|---:|---:|---:|---:|:--:|---:|---:|:--:|:--:|")
    for r in elev:
        A(f"| {r['height_m']:g} | {r['supination_deg']:g} | {r['F_lat_n']:,.0f} | "
          f"{r['M_inv_nm']:.1f} | {r['atfl_joint_risk']:.2f} | {_bool_cn(r['atfl_joint_sprain'])} | "
          f"{r['atfl_force_n']:,.0f} | {r['atfl_strain']*100:.1f}% | "
          f"{_bool_cn(r['atfl_fe_sprain'])} | {_bool_cn(r['cfl_joint_sprain'])} |")
    A("")
    A("> 力/应变判据**不随预载口径变**（它不含预载项）；关节级随口径变。")
    A("")
    A("### 5.3 扭伤格点汇总")
    A("")
    A("| 判据 | 保守口径 M_fail=21 N·m | 高承载口径 M_fail=77 N·m |")
    A("|---|---|---|")
    A(f"| ATFL 关节级 | {_cells_str(cons_j)} | {_cells_str(elev_j)} |")
    A(f"| CFL 关节级 | {_cells_str(_sprain_cells(cons, 'cfl_joint_sprain'))} | {_cells_str(_sprain_cells(elev, 'cfl_joint_sprain'))} |")
    A(f"| ATFL 力/应变 | {_cells_str(cons_fe)} | {_cells_str(_sprain_cells(elev, 'atfl_fe_sprain'))} |")
    A(f"| CFL 力/应变 | {_cells_str(_sprain_cells(cons, 'cfl_fe_sprain'))} | {_cells_str(_sprain_cells(elev, 'cfl_fe_sprain'))} |")
    A("")
    A("> `β=0` 全部格点安全（`M_inv=0`）。**CFL 与 ATFL 关节级同步越阈**（本模块用同一 "
      "M_fail 口径，CFL 是保守下限；现实中 CFL 晚于 ATFL）。")
    A("")
    A("---")
    A("")
    A("## 6. r_lig 灵敏度（**assumed**）")
    A("")
    A(f"r_lig 是韧带到内翻轴的力臂（无实测）。中央 **{DEFAULT_LIGAMENT_LEVER_ARM_MM:g} mm**，"
      f"灵敏度 {list(R_LIG_SENSITIVITY_MM)} mm。")
    A("")
    A("| h (m) | β (°) | r_lig=20: F (N) | r_lig=22: F (N) | r_lig=25: F (N) | 三者是否都扭伤 |")
    A("|---:|---:|---:|---:|---:|:--:|")
    keys = list(sens.keys())
    for i in range(len(sens[keys[0]])):
        r20 = sens["20"][i] if "20" in sens else None
        r22 = sens["22"][i] if "22" in sens else None
        r25 = sens["25"][i] if "25" in sens else None
        allsp = all(x["atfl_fe_sprain"] for x in (r20, r22, r25) if x is not None)
        A(f"| {r22['height_m']:g} | {r22['supination_deg']:g} | "
          f"{r20['atfl_force_n']:,.0f} | {r22['atfl_force_n']:,.0f} | {r25['atfl_force_n']:,.0f} | "
          f"{_bool_cn(allsp)} |")
    A("")
    A("> 三档 r_lig 下，**所有 β≥5° 格点都越阈** → 结论对 r_lig 不敏感"
      "（ATFL 阈值仅 200 N，而 M_inv/r_lig 量级达 10²–10³ N）。r_lig 越小 → F 越大 → 越保守。")
    A("")
    A("---")
    A("")
    A("## 7. [Tochigi2006] 正常步态参照带（**不是失效阈值**）")
    A("")
    A(f"- 正常步态 ATFL 应变 ≤ **{band['tochigi_atfl_normal_max']*100:.1f}%**；"
      f"CFL ≤ **{band['tochigi_cfl_normal_max']*100:.1f}%**；'taut' 阈值 = "
      f"**{band['tochigi_taut_threshold']*100:.1f}%**。")
    A("- 我们把这当成生理参照带：正常步态下 ATFL 几乎不绷紧（Tochigi：11 例中仅 1 例 "
      "在运动弧内 taut）。")
    A("")
    A("| h (m) | β=5° 的 M_inv (N·m) | ATFL ε | ε / Tochigi 正常上限 |")
    A("|---:|---:|---:|---:|")
    for r in band["rows"]:
        A(f"| {r['height_m']:g} | {r['M_inv_nm']:.1f} | {r['atfl_strain']*100:.1f}% | "
          f"**{r['atfl_strain_over_tochigi_normal']:.0f}×** |")
    A("")
    A(f"> {band['caveat']} 结论：冲击载荷下 ATFL 应变远超正常步态带（>10³×），"
      "且远高于 1% taut 阈值 → **是明确的外伤事件**，与 Tochigi "
      "'扭伤 = 极端位置事件'前提一致。")
    A("")
    A("---")
    A("")
    A("## 8. 与临床踝伤构成的对照")
    A("")
    A("| 口径 | 占比 | 我们的模型能否表示 |")
    A("|---|---:|---|")
    A(f"| 踝扭伤（[Heck2024]） | **{clin['heck2024_ankle_sprain_share_pct']:.1f}%** | ✅ 本模块（R6） |")
    A(f"| 踝骨折（[Heck2024]） | {clin['heck2024_ankle_fracture_share_pct']:.1f}% | ✅ R1 骨口径 |")
    A("")
    A(f"- [Beurienne2025]：踝扭伤 **{clin['beurienne2025_ankle_sprain_rank']}**。")
    A(f"- 保守口径下，{clin['our_conservative_sprain_cells']}/{clin['our_conservative_total_cells']} "
      "格点预测 ATFL 扭伤 → 与'扭伤是踝伤 #1'**定性一致**。")
    A(f"- 骨模型（R1）此前只能对应踝骨折 27%；{clin['note']}")
    A("")
    A("---")
    A("")
    A("## 9. 假设（ASSUMPTIONS）与诚实边界")
    A("")
    for i, a in enumerate(res["assumptions"], 1):
        A(f"{i}. {a}")
    A("")
    A("**结构边界**：")
    A("")
    A("- 关节级判据是**主口径**；力/应变判据把全部 M_inv 归单条韧带 → 上界式。")
    A("- r_lig 是 **assumed**；换 r_lig 只改力/应变数值，不改'是否越阈'。")
    A("- 预载承载力外推**不确定**；已给两口径 + clamp/linear 边界。")
    A("- Tochigi 是**矢状面步态**、600 N；我们的是**冠状面冲击**、21 kN —— 不同轴/"
      "量级，只作参照。")
    A("- 本脚本只读复用既有产物，不调用 FEBio，不覆盖既有产物。")
    A("")
    A("---")
    A("")
    A("## 10. 复现命令")
    A("")
    A("```powershell")
    A("cd D:\\Project\\climbing_fall_analysis")
    A('$env:PYTHONPATH="src"')
    A("& .venv\\Scripts\\python.exe scripts\\opensim_fe\\ankle_ligament_criterion.py")
    A("# 只读复用：")
    A("#   results/opensim_fe/nonvertical_s5.json  (R1 旋后/内翻缓存：F_vert/F_lat/M_inv)")
    A("# 回归：")
    A("& .venv\\Scripts\\python.exe -m pytest tests\\test_opensim_fe.py tests\\test_febio_run.py -q")
    A("& .venv\\Scripts\\python.exe -m pytest tests\\test_nonvertical_s1.py tests\\test_nonvertical_s2.py tests\\test_nonvertical_s5.py tests\\test_pad_supination_tradeoff.py tests\\test_ankle_ligament_criterion.py -q")
    A("```")
    A("")
    A("## 11. 产物")
    A("")
    A(f"- `{OUT_MD.relative_to(ROOT)}`（本报告）")
    A(f"- `{OUT_JSON.relative_to(ROOT)}`（机器可读）")
    A(f"- `src/climbing/coupling/ankle_ligament.py`（新模块）")
    A(f"- `scripts/opensim_fe/ankle_ligament_criterion.py`（本脚本）")
    A(f"- `tests/test_ankle_ligament_criterion.py`（新测试）")
    A("")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())