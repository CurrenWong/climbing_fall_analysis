"""R4 —— 垫子刚度 k ↔ 踝旋后 β 两难耦合模型的回归测试。

覆盖新模块 ``climbing.coupling.pad_supination`` 与脚本产物
``results/opensim_fe/pad_supination_tradeoff.json``。

**硬门槛**：耦合模块在"无旋后"极限（β_0=0、k→∞）下必须**逐位退化**到纯轴向 —
``β=0``、``lateral_risk=0``、``total_risk=axial_risk``，且 ``axial_risk`` 与 R1
缓存 ``nonvertical_s5.json`` 的纯轴向值相对误差 ≤ 1e-9。

快测试只读既有 JSON 产物（毫秒级），不重跑 OpenSim / FEBio。
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from climbing.coupling.pad_supination import (
    BETA_BASELINE_DEG_DEFAULT,
    BETA_MAX_DEG,
    BETA_ROLL_LEVER_ARM_MM_DEFAULT,
    CoupledPadResult,
    coupled_pad_risk,
    coupled_pad_risk_at_grid,
    find_optimal_stiffness_n_per_mm,
    indentation_depth_mm,
    load_fe_transfer_table,
    pad_reduction_ratio,
    required_axial_reduction_for_soft_optimum,
    supination_angle_deg,
    sweep_roll_lever_arm,
)

M_MODEL = 75.337
ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results" / "opensim_fe"
PAD_STIFFNESS_JSON = RES / "pad_stiffness.json"
S5_JSON = RES / "nonvertical_s5.json"
TRADEOFF_JSON = RES / "pad_supination_tradeoff.json"

#: 主口径腓骨远端截面（等效圆，来自 risk_1d_bending.json / R1）。
_EQ = dict(A_section_mm2=101.93640134318377, I_mm4=826.8918876962899,
           c_mm=5.696258799381752)
#: h → 单足峰值 GRF（刚性地面，来自 R1 measured_heights）。
#: 与 nonvertical_s5.json 一致；测试里独立硬编码以锁定不变量。
F_RIGID_BY_H = {
    2.0: 21131.84246941172,
    3.0: 22123.30461015003,
    3.5: 22895.856216565313,
    4.5: 25095.73293947136,
}
#: R1 刚性纯轴向腓骨 risk（scenario_hard, β=0）。
AXIAL_RIGID_BY_H = {
    2.0: 0.6590277487199715,
    3.0: 0.6886929667323994,
    3.5: 0.7120194762480256,
    4.5: 0.7789583294729311,
}


def _fe_table() -> dict:
    return load_fe_transfer_table(PAD_STIFFNESS_JSON)


def _cell(h: float, k, *, lever_arm_geom_mm: float = BETA_ROLL_LEVER_ARM_MM_DEFAULT):
    return coupled_pad_risk(
        height_m=h, k_n_per_mm=k,
        F_peak_rigid_per_foot_n=F_RIGID_BY_H[h],
        fe_table=_fe_table(), **_EQ,
        axial_risk_rigid=AXIAL_RIGID_BY_H[h],
        lever_arm_geom_mm=lever_arm_geom_mm,
    )


# --------------------------------------------------------------------------
# S1 折减比复用（measured）
# --------------------------------------------------------------------------
@pytest.mark.skipif(not PAD_STIFFNESS_JSON.is_file(), reason="pad_stiffness.json 未生成")
def test_pad_reduction_ratio_matches_s1_table() -> None:
    fe = _fe_table()
    # 表中每一格必须精确复现（除 k>=max 视为刚性极限 1.0）
    for k, ratio in zip(fe["k_N_per_mm"], fe["ratio_vs_fixed"]):
        k = float(k)
        if k >= float(fe["k_N_per_mm"][-1]):
            assert pad_reduction_ratio(k, fe) == 1.0
        else:
            assert pad_reduction_ratio(k, fe) == pytest.approx(float(ratio), rel=1e-12)


@pytest.mark.skipif(not PAD_STIFFNESS_JSON.is_file(), reason="pad_stiffness.json 未生成")
def test_pad_reduction_ratio_limits() -> None:
    fe = _fe_table()
    assert pad_reduction_ratio(None, fe) == 1.0                 # 刚性极限
    assert pad_reduction_ratio(1e9, fe) == 1.0                  # 远超表 → 刚性
    assert pad_reduction_ratio(1.0, fe) == pytest.approx(      # 低于表 → 平推软端
        float(fe["ratio_vs_fixed"][0]), rel=1e-12)
    # R_min 落在 k=1e3（S1 结论）
    assert fe["k_at_r_min_N_per_mm"] == pytest.approx(1000.0)
    assert pad_reduction_ratio(1000.0, fe) == pytest.approx(fe["r_min"], rel=1e-12)


def test_pad_reduction_ratio_bad_table_raises() -> None:
    with pytest.raises(ValueError):
        pad_reduction_ratio(1000.0, {"k_N_per_mm": [], "ratio_vs_fixed": []})


# --------------------------------------------------------------------------
# β(k) 几何模型（modeled）
# --------------------------------------------------------------------------
def test_indentation_depth_is_force_over_k() -> None:
    assert indentation_depth_mm(20000.0, 1000.0) == pytest.approx(20.0)
    assert indentation_depth_mm(0.0, 1000.0) == 0.0
    assert indentation_depth_mm(20000.0, 0.0) == 0.0     # k<=0 → 不凹陷
    assert indentation_depth_mm(-1.0, 1000.0) == 0.0


def test_supination_angle_decreases_with_k() -> None:
    f = 21000.0
    betas = [supination_angle_deg(f, k) for k in (10.0, 100.0, 1000.0, 10000.0, 1e6)]
    assert all(b2 <= b1 for b1, b2 in zip(betas, betas[1:])), betas


def test_supination_angle_rigid_limit_is_baseline() -> None:
    # k=inf（数值上极大）→ δ=0 → β=β_0
    assert supination_angle_deg(21000.0, float("inf")) == BETA_BASELINE_DEG_DEFAULT
    assert BETA_BASELINE_DEG_DEFAULT == 0.0


def test_supination_angle_saturates_at_max() -> None:
    # 极软 → δ 巨大 → β=β_max
    assert supination_angle_deg(21000.0, 1.0) == pytest.approx(BETA_MAX_DEG)


def test_supination_angle_baseline_offset() -> None:
    b = supination_angle_deg(21000.0, 1000.0, beta_baseline_deg=5.0)
    b0 = supination_angle_deg(21000.0, 1000.0, beta_baseline_deg=0.0)
    assert b == pytest.approx(b0 + 5.0, rel=1e-9)


def test_supination_angle_bad_lever_raises() -> None:
    with pytest.raises(ValueError):
        supination_angle_deg(21000.0, 1000.0, lever_arm_mm=0.0)


def test_supination_angle_monotone_in_lever_arm() -> None:
    # 杠杆越大 → 同 δ 下 β 越小
    b_small = supination_angle_deg(21000.0, 1000.0, lever_arm_mm=25.0)
    b_large = supination_angle_deg(21000.0, 1000.0, lever_arm_mm=60.0)
    assert b_large < b_small


# --------------------------------------------------------------------------
# 耦合：(h, k) → total_risk
# --------------------------------------------------------------------------
def test_coupled_rigid_limit_degenerates_to_axial() -> None:
    r = _cell(4.5, None)
    assert r.R_fe == 1.0
    assert r.beta_deg == 0.0
    assert r.F_lat_n == 0.0
    assert r.M_inv_nmm == 0.0
    assert r.lateral_risk == 0.0
    assert r.total_risk == pytest.approx(r.axial_risk, rel=0, abs=1e-15)
    assert r.fracture is False


def test_coupled_axial_matches_r1_baseline() -> None:
    for h in F_RIGID_BY_H:
        r = _cell(h, None)
        assert r.axial_risk == pytest.approx(AXIAL_RIGID_BY_H[h], rel=1e-12)


def test_coupled_axial_scales_with_reduction_ratio() -> None:
    for k in (100.0, 1000.0, 10000.0):
        r = _cell(4.5, k)
        expected = AXIAL_RIGID_BY_H[4.5] * pad_reduction_ratio(k, _fe_table())
        assert r.axial_risk == pytest.approx(expected, rel=1e-12)


def test_coupled_soft_pad_raises_lateral_risk() -> None:
    # 同一高度：软垫（k=1e3）侧向风险应高于硬垫（k=1e6）
    soft = _cell(4.5, 1000.0)
    hard = _cell(4.5, 1e6)
    assert soft.beta_deg > hard.beta_deg
    assert soft.lateral_risk > hard.lateral_risk
    assert soft.total_risk > hard.total_risk


def test_coupled_total_is_axial_plus_lateral() -> None:
    r = _cell(3.0, 1000.0)
    assert r.total_risk == pytest.approx(r.axial_risk + r.lateral_risk, rel=1e-12)


def test_coupled_fracture_flag_matches_threshold() -> None:
    assert _cell(2.0, 10.0).fracture is True        # 极软 → 越阈
    assert _cell(2.0, None).fracture is False       # 刚性 → 不越阈
    for k in (10.0, 1000.0, 1e6, None):
        r = _cell(4.5, k)
        assert r.fracture == (r.total_risk >= 1.0)


def test_coupled_as_dict_roundtrip() -> None:
    d = _cell(4.5, 1000.0).as_dict()
    for key in ("height_m", "k_n_per_mm", "R_fe", "beta_deg", "F_lat_n",
                "M_inv_nmm", "axial_risk", "lateral_risk", "total_risk",
                "fracture", "A_section_mm2", "sigma_c_mpa"):
        assert key in d
    assert d["M_inv_nm"] == pytest.approx(d["M_inv_nmm"] / 1000.0, rel=1e-12)


def test_coupled_sigma_lateral_is_bend_plus_shear() -> None:
    r = _cell(4.5, 1000.0)
    assert r.sigma_lateral_mpa == pytest.approx(
        r.sigma_bend_mpa + r.sigma_shear_mpa, rel=1e-12)


# --------------------------------------------------------------------------
# k* 寻找
# --------------------------------------------------------------------------
def test_find_optimal_no_interior_returns_rigid_boundary() -> None:
    ks = [50.0, 100.0, 300.0, 1000.0, 3000.0, 10000.0, 1e5, 1e6, None]
    cells = [_cell(4.5, k) for k in ks]
    opt = find_optimal_stiffness_n_per_mm(cells)
    assert opt["verdict_kind"] == "boundary_min"
    assert opt["k_star_n_per_mm"] is None
    assert opt["interior_min"] is None
    assert opt["total_risk_star"] == pytest.approx(cells[-1].total_risk, rel=1e-12)


def test_find_optimal_detects_interior_min() -> None:
    """构造一个合成内部极小（超大 L_roll，数值实验）→ 必须识别为 interior_min。

    这不是物理参数，只用于验证 :func:`find_optimal_stiffness_n_per_mm`
    的**内部极小分支**（真实 L_roll∈{25,40,60} 下无内部极小，见 §4 灵敏度）。
    """
    ks = [100.0, 300.0, 1000.0, 3000.0, 10000.0, 30000.0, 1e5, 1e6, None]
    cells = [_cell(2.0, k, lever_arm_geom_mm=10000.0) for k in ks]
    # 内部极小确实存在（total_risk 在 k=1e3 最低）
    totals = {c.k_n_per_mm: c.total_risk for c in cells}
    assert totals[1000.0] < totals[100.0]
    assert totals[1000.0] < totals[None]
    opt = find_optimal_stiffness_n_per_mm(cells)
    assert opt["verdict_kind"] == "interior_min"
    assert opt["k_star_n_per_mm"] == pytest.approx(1000.0)


def test_find_optimal_empty_raises() -> None:
    with pytest.raises(ValueError):
        find_optimal_stiffness_n_per_mm([])


# --------------------------------------------------------------------------
# 内部极小存在性的定量门槛
# --------------------------------------------------------------------------
def test_required_axial_reduction_formula() -> None:
    # total(soft) < total(rigid) ⟺ 1-R > lat/axial
    assert required_axial_reduction_for_soft_optimum(0.5, 0.1) == pytest.approx(0.2)
    assert required_axial_reduction_for_soft_optimum(0.0, 0.1) == float("inf")


def test_required_axial_reduction_exceeds_measured_for_real_fibula() -> None:
    """真实腓骨下：所需折减 ≫ S1 实测最大折减（1−R_min）⇒ 软端不可能胜出。"""
    fe = _fe_table()
    measured = 1.0 - float(fe["r_min"])
    for h in F_RIGID_BY_H:
        for k in (50.0, 100.0, 1000.0):
            r = _cell(h, k)
            req = required_axial_reduction_for_soft_optimum(
                AXIAL_RIGID_BY_H[h], r.lateral_risk)
            assert req > measured, (h, k, req, measured)


# --------------------------------------------------------------------------
# 灵敏度扫描（L_roll）
# --------------------------------------------------------------------------
def test_sweep_roll_lever_arm_shape() -> None:
    ks = [100.0, 1000.0, 10000.0, 1e6, None]
    out = sweep_roll_lever_arm(
        [2.0, 4.5], ks,
        F_peak_rigid_per_foot_by_h=F_RIGID_BY_H,
        fe_table=_fe_table(), **_EQ,
        axial_risk_rigid_by_h=AXIAL_RIGID_BY_H,
        lever_arm_geom_grid_mm=(25.0, 40.0, 60.0),
    )
    assert set(out.keys()) == {25.0, 40.0, 60.0}
    for L, cells in out.items():
        assert len(cells) == 2 * len(ks)
        assert all(isinstance(c, CoupledPadResult) for c in cells)


def test_sweep_larger_lever_arm_reduces_lateral() -> None:
    ks = [1000.0]
    out = sweep_roll_lever_arm(
        [4.5], ks,
        F_peak_rigid_per_foot_by_h=F_RIGID_BY_H,
        fe_table=_fe_table(), **_EQ,
        axial_risk_rigid_by_h=AXIAL_RIGID_BY_H,
        lever_arm_geom_grid_mm=(25.0, 40.0, 60.0),
    )
    lat25 = out[25.0][0].lateral_risk
    lat40 = out[40.0][0].lateral_risk
    lat60 = out[60.0][0].lateral_risk
    assert lat60 < lat40 < lat25


def test_coupled_pad_risk_at_grid_count() -> None:
    ks = [100.0, 1000.0, None]
    cells = coupled_pad_risk_at_grid(
        [2.0, 3.0], ks,
        F_peak_rigid_per_foot_by_h=F_RIGID_BY_H,
        fe_table=_fe_table(), **_EQ,
        axial_risk_rigid_by_h=AXIAL_RIGID_BY_H,
    )
    assert len(cells) == 2 * len(ks)


# --------------------------------------------------------------------------
# 产物一致性（若已生成脚本产物）
# --------------------------------------------------------------------------
@pytest.mark.skipif(not TRADEOFF_JSON.is_file(), reason="tradeoff.json 未生成")
def test_artifact_default_regression_passes() -> None:
    data = json.loads(TRADEOFF_JSON.read_text(encoding="utf-8"))
    reg = data["default_regression"]
    assert reg["verdict"] == "PASS"
    assert reg["axial_risk_rel_err_vs_r1"] <= 1e-9
    for name, ok in reg["checks"].items():
        assert ok is True, name


@pytest.mark.skipif(not TRADEOFF_JSON.is_file(), reason="tradeoff.json 未生成")
def test_artifact_grid_and_optimum() -> None:
    data = json.loads(TRADEOFF_JSON.read_text(encoding="utf-8"))
    heights = sorted({c["height_m"] for c in data["sweep"]["cells"]})
    assert heights == [2.0, 3.0, 3.5, 4.5]
    # 中央口径：无内部极小、k* 落在刚性极限
    for o in data["sweep"]["optimum_by_height"]:
        assert o["verdict_kind"] == "boundary_min"
        assert o["k_star_n_per_mm"] is None
        assert o["interior_min"] is None


@pytest.mark.skipif(not TRADEOFF_JSON.is_file(), reason="tradeoff.json 未生成")
def test_artifact_sensitivity_robust() -> None:
    data = json.loads(TRADEOFF_JSON.read_text(encoding="utf-8"))
    sens = data["sensitivity_lever_arm"]
    assert sens["robust_no_interior_optimum"] is True
    assert sens["robust_boundary_optimum_is_rigid"] is True
    assert set(sens["by_lever_arm"].keys()) == {"25", "40", "60"}


@pytest.mark.skipif(not TRADEOFF_JSON.is_file(), reason="tradeoff.json 未生成")
def test_artifact_cells_monotone_total_risk() -> None:
    """中央口径：每个 h 的 total_risk 随 k 单调递减（软→硬）。"""
    data = json.loads(TRADEOFF_JSON.read_text(encoding="utf-8"))
    from collections import defaultdict
    by_h: dict[float, list[dict]] = defaultdict(list)
    for c in data["sweep"]["cells"]:
        by_h[c["height_m"]].append(c)
    for h, cells in by_h.items():
        cells = sorted(cells, key=lambda c: 1e30 if c["k_n_per_mm"] is None
                       else c["k_n_per_mm"])
        totals = [c["total_risk"] for c in cells]
        assert all(b <= a + 1e-12 for a, b in zip(totals, totals[1:])), h


@pytest.mark.skipif(not TRADEOFF_JSON.is_file(), reason="tradeoff.json 未生成")
def test_artifact_breakeven_blocks_soft_optimum() -> None:
    data = json.loads(TRADEOFF_JSON.read_text(encoding="utf-8"))
    be = data["breakeven_diagnostic"]
    assert be["n_soft_can_win"] == 0
    measured = float(be["measured_max_reduction_1_minus_Rmin"])
    assert measured == pytest.approx(1.0 - float(data["fe_transfer_table"]["r_min"]), rel=1e-12)


@pytest.mark.skipif(not TRADEOFF_JSON.is_file(), reason="tradeoff.json 未生成")
def test_artifact_beta_monotone_decreasing_in_k() -> None:
    data = json.loads(TRADEOFF_JSON.read_text(encoding="utf-8"))
    from collections import defaultdict
    by_h: dict[float, list[dict]] = defaultdict(list)
    for c in data["sweep"]["cells"]:
        by_h[c["height_m"]].append(c)
    for h, cells in by_h.items():
        cells = sorted(cells, key=lambda c: 1e30 if c["k_n_per_mm"] is None
                       else c["k_n_per_mm"])
        betas = [c["beta_deg"] for c in cells]
        assert all(b <= a + 1e-12 for a, b in zip(betas, betas[1:])), h
