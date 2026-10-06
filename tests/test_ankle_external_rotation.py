"""Phase-S6 外旋（external rotation / SER）—— 踝外旋 / 远端腓骨扭转 + 弯曲风险 回归测试。

**硬门槛（与 Phase-S5 同源）**：外旋角 = 0（默认）时必须**逐位**复现扩展前的
轴向行为。本文件锁定接口契约（无需 OpenSim）：默认零、纯轴向基线、单调性、
阈值跨越，以及与现有 ``ankle_supination`` 模块在 opt-in 上**正交**。

约定与 ``tests/test_nonvertical_s5.py`` 同。
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from climbing.coupling.ankle_rotation import (
    DEFAULT_TIBIAL_LEVER_ARM_MM,
    FIBULA_ENDS_SIGMA_C_MPA,
    critical_external_rotation_angle_deg,
    external_rotation_load,
    fibula_external_rotation_risk,
    fibula_torsional_stress,
    horizontal_torque_nmm,
)
from climbing.coupling.ankle_supination import (
    critical_inversion_angle_deg,
    fibula_lateral_risk,
    supination_load,
)

ROOT = Path(__file__).resolve().parents[1]
S6_JSON = ROOT / "results" / "opensim_fe" / "ankle_external_rotation.json"

#: 主口径腓骨远端截面（等效圆，来自 risk_1d_bending.json），
#: 与 test_nonvertical_s5 共享同一组数（不重写）。
_EQ = dict(i_mm4=826.8918876962899, c_mm=5.696258799381752)
#: A_section 在外旋模型里**用不到**（J=2·I，c=r，无 A），但保留供对照。
_A_SECTION_MM2 = 101.93640134318377

# 典型竖直 GRF：5 m / 75 kg ≈ 21 kN（与 R1 共享同一力）。
_F_VERT_TEST_N = 21000.0


# --------------------------------------------------------------------------
# L0 —— 默认 / 接口契约
# --------------------------------------------------------------------------
def test_default_external_rotation_load_is_zero() -> None:
    """θ=0 ⇒ T_ext=0 ⇒ risk=0 ⇒ **逐位复现**旧轴向行为。"""
    load = external_rotation_load((0.0, _F_VERT_TEST_N, 0.0))
    assert load.external_rotation_deg == 0.0
    assert load.lever_arm_mm == DEFAULT_TIBIAL_LEVER_ARM_MM
    assert load.t_ext_nmm == 0.0
    assert load.t_ext_nm == 0.0
    assert load.is_zero

    r = fibula_external_rotation_risk(load.t_ext_nmm, **_EQ)
    assert r["risk"] == 0.0
    assert r["fracture"] is False
    assert r["sigma_vonmises_mpa"] == 0.0
    assert r["sigma_torsion_mpa"] == 0.0


def test_default_does_not_affect_supination_module() -> None:
    """S6 默认路径不应影响 S5 默认路径（两个 opt-in 维度正交）。"""
    rot = external_rotation_load((0.0, _F_VERT_TEST_N, 0.0))
    sup = supination_load((0.0, _F_VERT_TEST_N, 0.0))
    assert rot.is_zero
    assert sup.is_zero
    assert rot.t_ext_nmm == 0.0 == sup.m_inv_nmm


def test_sigma_c_constant_matches_supination() -> None:
    """与 ``ankle_supination.FIBULA_ENDS_SIGMA_C_MPA`` 同源（70 MPa）。"""
    assert FIBULA_ENDS_SIGMA_C_MPA == 70.0
    # 显式从 supination 模块导入同源常量，避免单向漂移。
    from climbing.coupling.ankle_supination import FIBULA_ENDS_SIGMA_C_MPA as FIB_SUP
    assert FIB_SUP == FIBULA_ENDS_SIGMA_C_MPA


def test_external_rotation_load_bad_vector_raises() -> None:
    with pytest.raises(ValueError):
        external_rotation_load((0.0, 1.0))
    with pytest.raises(ValueError):
        external_rotation_load((0.0, np.nan, 0.0))


def test_external_rotation_load_force_vec_stored() -> None:
    f = (1.5, _F_VERT_TEST_N, -2.5)
    load = external_rotation_load(f, 10.0, 30.0)
    assert load.force_vec_n == (1.5, _F_VERT_TEST_N, -2.5)
    # 力臂 / 角被记录
    assert load.lever_arm_mm == 30.0
    assert load.external_rotation_deg == 10.0


# --------------------------------------------------------------------------
# L1 —— 力矩定义
# --------------------------------------------------------------------------
def test_horizontal_torque_zero_at_theta0() -> None:
    assert horizontal_torque_nmm((0.0, _F_VERT_TEST_N, 0.0), 0.0, 30.0) == 0.0


def test_horizontal_torque_definition() -> None:
    # T_ext = F_y · d · sinθ
    expected = _F_VERT_TEST_N * 30.0 * np.sin(np.deg2rad(15.0))
    got = horizontal_torque_nmm((0.0, _F_VERT_TEST_N, 0.0), 15.0, 30.0)
    assert got == pytest.approx(expected, rel=1e-12)


def test_horizontal_torque_sign_insensitive() -> None:
    """返回值为绝对值（与 ``inversion_moment_nmm`` 同口径）。"""
    a = horizontal_torque_nmm((0.0, _F_VERT_TEST_N, 0.0), 20.0, 30.0)
    b = horizontal_torque_nmm((0.0, _F_VERT_TEST_N, 0.0), -20.0, 30.0)
    assert a == pytest.approx(b, rel=1e-12)


def test_external_rotation_load_t_ext_known() -> None:
    load = external_rotation_load((0.0, _F_VERT_TEST_N, 0.0), 15.0, 30.0)
    expected = _F_VERT_TEST_N * 30.0 * np.sin(np.deg2rad(15.0))
    assert load.t_ext_nmm == pytest.approx(expected, rel=1e-12)
    assert load.t_ext_nm == pytest.approx(expected / 1000.0, rel=1e-12)
    assert not load.is_zero


# --------------------------------------------------------------------------
# L2 —— 扭转 / 弯曲 应力
# --------------------------------------------------------------------------
def test_fibula_torsional_stress_known_values() -> None:
    """τ = T · c / J，J=2·I（等效圆），σ_vm = √3·τ。"""
    load = external_rotation_load((0.0, _F_VERT_TEST_N, 0.0), 15.0, 30.0)
    s = fibula_torsional_stress(load.t_ext_nmm, **_EQ)
    j = 2.0 * _EQ["i_mm4"]
    exp_tau = load.t_ext_nmm * _EQ["c_mm"] / j
    assert s["sigma_torsion_mpa"] == pytest.approx(exp_tau, rel=1e-12)
    assert s["sigma_shear_mpa"] == pytest.approx(exp_tau, rel=1e-12)
    assert s["sigma_vonmises_mpa"] == pytest.approx(np.sqrt(3.0) * exp_tau, rel=1e-12)
    assert s["j_mm4"] == pytest.approx(j, rel=1e-12)


def test_fibula_torsional_stress_zero_path() -> None:
    s = fibula_torsional_stress(0.0, **_EQ)
    assert s["sigma_torsion_mpa"] == 0.0
    assert s["sigma_vonmises_mpa"] == 0.0


def test_fibula_torsional_stress_bad_section_raises() -> None:
    with pytest.raises(ValueError):
        fibula_torsional_stress(1000.0, i_mm4=0.0, c_mm=5.0)
    with pytest.raises(ValueError):
        fibula_torsional_stress(1000.0, i_mm4=826.0, c_mm=-1.0)


# --------------------------------------------------------------------------
# L3 —— 风险 / 阈值跨越
# --------------------------------------------------------------------------
def test_fibula_external_rotation_risk_default_zero() -> None:
    """默认 θ=0 ⇒ risk=0 ⇒ 不越阈 ⇒ **纯轴向基线**。"""
    load = external_rotation_load((0.0, _F_VERT_TEST_N, 0.0))
    r = fibula_external_rotation_risk(load.t_ext_nmm, **_EQ)
    assert r["risk"] == 0.0
    assert r["fracture"] is False
    assert r["sigma_c_mpa"] == FIBULA_ENDS_SIGMA_C_MPA


def test_fibula_external_rotation_risk_crosses_at_theta10() -> None:
    """θ=10°，F_y=21 kN，d=30 mm ⇒ 应**越阈**（按 R1 等口径上界式 1D 筛选）。"""
    load = external_rotation_load((0.0, _F_VERT_TEST_N, 0.0), 10.0, 30.0)
    r = fibula_external_rotation_risk(load.t_ext_nmm, **_EQ)
    assert r["risk"] >= 1.0
    assert r["fracture"] is True


def test_fibula_external_rotation_risk_threshold_scaling() -> None:
    """risk ∝ 1/σ_c：阈值翻倍 ⇒ risk 减半。"""
    load = external_rotation_load((0.0, _F_VERT_TEST_N, 0.0), 15.0, 30.0)
    r70 = fibula_external_rotation_risk(load.t_ext_nmm, sigma_c_mpa=70.0, **_EQ)
    r140 = fibula_external_rotation_risk(load.t_ext_nmm, sigma_c_mpa=140.0, **_EQ)
    assert r70["risk"] == pytest.approx(2.0 * r140["risk"], rel=1e-12)


def test_fibula_external_rotation_risk_bad_sigma_c_raises() -> None:
    with pytest.raises(ValueError):
        fibula_external_rotation_risk(1000.0, sigma_c_mpa=0.0, **_EQ)


# --------------------------------------------------------------------------
# L4 —— 单调性（外旋角 → risk 单调非降）
# --------------------------------------------------------------------------
@pytest.mark.parametrize("theta_a,theta_b", [(0.0, 5.0), (5.0, 15.0), (10.0, 30.0)])
def test_risk_monotone_in_external_rotation(theta_a: float, theta_b: float) -> None:
    a = fibula_external_rotation_risk(
        external_rotation_load((0.0, _F_VERT_TEST_N, 0.0), theta_a, 30.0).t_ext_nmm,
        **_EQ,
    )
    b = fibula_external_rotation_risk(
        external_rotation_load((0.0, _F_VERT_TEST_N, 0.0), theta_b, 30.0).t_ext_nmm,
        **_EQ,
    )
    assert b["risk"] >= a["risk"]


def test_risk_monotone_in_axial_load() -> None:
    """F_y 越大 ⇒ T_ext 越大 ⇒ risk 越大（保持外旋角与力臂固定）。"""
    risks = []
    for fy in (5000.0, 10000.0, 15000.0, 21000.0):
        load = external_rotation_load((0.0, fy, 0.0), 10.0, 30.0)
        r = fibula_external_rotation_risk(load.t_ext_nmm, **_EQ)
        risks.append(r["risk"])
    assert all(b >= a for a, b in zip(risks, risks[1:])), risks


def test_risk_monotone_in_lever_arm() -> None:
    """d 越大 ⇒ T_ext 越大 ⇒ risk 越大。"""
    risks = []
    for d in (10.0, 20.0, 30.0, 45.0):
        load = external_rotation_load((0.0, _F_VERT_TEST_N, 0.0), 10.0, d)
        r = fibula_external_rotation_risk(load.t_ext_nmm, **_EQ)
        risks.append(r["risk"])
    assert all(b >= a for a, b in zip(risks, risks[1:])), risks


# --------------------------------------------------------------------------
# L5 —— 临界外旋角
# --------------------------------------------------------------------------
def test_critical_external_rotation_angle_monotone_bracket() -> None:
    b = critical_external_rotation_angle_deg((0.0, _F_VERT_TEST_N, 0.0), **_EQ)
    assert b is not None and 0.0 < b < 60.0
    # 在临界角处 risk ≈ 1
    load = external_rotation_load((0.0, _F_VERT_TEST_N, 0.0), b, 30.0)
    r = fibula_external_rotation_risk(load.t_ext_nmm, **_EQ)
    assert r["risk"] == pytest.approx(1.0, rel=1e-6)


def test_critical_external_rotation_angle_none_when_below() -> None:
    """极小竖直力 → 60° 内不越阈。"""
    assert critical_external_rotation_angle_deg((0.0, 100.0, 0.0), **_EQ) is None


def test_critical_external_rotation_angle_approaches_zero_at_huge_axial() -> None:
    """极大竖直力 → θ_c → 0°（在 ``θ→0`` 即越阈，二分收敛到极小角）。

    注意：``risk(θ=0)=0``（sin 0 = 0）恒成立，故**严格 θ=0 永不越阈**；
    二分返回的是「risk 首次达到 1」的最小正角。F_y 越大 ⇒ θ_c 越接近 0。
    二分在 [0, 60] 上做 60 次，绝对精度 ≈ 5e-17°；在 1e-5° 量级上相对误差
    约 1e-5，故断言用 ``rel=1e-3``。
    """
    b = critical_external_rotation_angle_deg((0.0, 1.0e9, 0.0), **_EQ)
    assert b is not None
    assert b < 1.0e-3  # 1e9 N ≈ 4.7 万倍冲击 ⇒ θ_c < 1 mdeg
    # 与 21 kN 同力臂下对比：F_y 大 4.76e4× ⇒ sin θ_c 同比例小 ⇒ θ_c ≈ θ_c21k / 4.76e4
    b21k = critical_external_rotation_angle_deg((0.0, 21000.0, 0.0), **_EQ)
    assert b21k is not None
    assert b == pytest.approx(b21k / (1.0e9 / 21000.0), rel=1e-3)


# --------------------------------------------------------------------------
# L6 —— 与 S5 旋后正交（同一 F_y 下两种机制各自的临界角是同一量级，但驱动量不同）
# --------------------------------------------------------------------------
def test_critical_angles_in_same_order_of_magnitude() -> None:
    """同一 F_y=21 kN 下，旋后 θ_c 与外旋 θ_c 都在 0.5–2° 区间 —— 与 R1 同口径
    '上界式 1D 模型对旋转角极敏感'。这不是数值巧合：两者都是 1D 梁 / 圆柱扭转
    在 30 mm 力臂、70 MPa、21 kN 下的同一量级解。"""
    b_sup = critical_inversion_angle_deg(
        (0.0, _F_VERT_TEST_N, 0.0), a_section_mm2=_A_SECTION_MM2, **_EQ,
    )
    b_ext = critical_external_rotation_angle_deg(
        (0.0, _F_VERT_TEST_N, 0.0), **_EQ,
    )
    assert b_sup is not None and 0.0 < b_sup < 5.0
    assert b_ext is not None and 0.0 < b_ext < 5.0


def test_default_path_matches_supination_module_default() -> None:
    """两个模块在默认路径上同时**逐位零** —— 共同支撑 'opt-in 维度正交'。"""
    fvec = (0.0, _F_VERT_TEST_N, 0.0)
    r_sup = fibula_lateral_risk(
        supination_load(fvec).f_lat_n, supination_load(fvec).m_inv_nmm,
        a_section_mm2=_A_SECTION_MM2, **_EQ,
    )
    r_ext = fibula_external_rotation_risk(
        external_rotation_load(fvec).t_ext_nmm, **_EQ,
    )
    assert r_sup["risk"] == 0.0 == r_ext["risk"]
    assert r_sup["fracture"] is False and r_ext["fracture"] is False


# --------------------------------------------------------------------------
# L7 —— 产物一致性（若已生成脚本产物）
# --------------------------------------------------------------------------
@pytest.mark.skipif(not S6_JSON.is_file(), reason="ankle_external_rotation.json 未生成")
def test_artifact_default_rows_have_zero_torque() -> None:
    data = json.loads(S6_JSON.read_text(encoding="utf-8"))
    iface = data["interface_regression"]
    assert iface["external_rotation_deg0_is_zero"] is True
    # 默认（θ=0）行：T_ext=0、σ_vm=0、risk=0、fracture=False
    for r in data["scenario_hard"]:
        if r["external_rotation_deg"] == 0.0:
            assert r["T_ext_nmm"] == 0.0
            assert r["sigma_vonmises_mpa"] == 0.0
            assert r["fibula_external_rotation_risk"] == 0.0
            assert r["fibula_external_rotation_fracture"] is False


@pytest.mark.skipif(not S6_JSON.is_file(), reason="ankle_external_rotation.json 未生成")
def test_artifact_grid_and_threshold_crossing() -> None:
    data = json.loads(S6_JSON.read_text(encoding="utf-8"))
    heights = sorted({r["height_m"] for r in data["scenario_hard"]})
    angles = sorted({r["external_rotation_deg"] for r in data["scenario_hard"]})
    assert heights == [2.0, 3.0, 3.5, 4.5]
    # 角度集合包含 {0, 5, 10, 15, 20, 30}（具体由脚本决定）
    assert 0.0 in angles
    assert any(a >= 5.0 for a in angles)
    # 一旦外旋角≥10°，全部测点应越阈（与 R1 同口径的上界式 1D 敏感性）
    for r in data["scenario_hard"]:
        if r["external_rotation_deg"] >= 10.0:
            assert r["fibula_external_rotation_fracture"] is True, r
