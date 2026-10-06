"""Phase-S5 非垂直落地（踝旋后 / 内翻，冠状面）回归测试。

**硬门槛（方案 §7）**：冠状地面倾角 / 旋后角 = 0（全部默认）时必须**逐位**复现
扩展前的轴向行为。本文件的快测试锁定接口契约（无需 OpenSim）；慢测试
（``RUN_OPEN_FE=1``）直接对比改前缓存 ``results/opensim_fe/s1_h5_opensim.json``。

慢测试放在 ``RUN_OPEN_FE`` 门后，与 ``tests/test_nonvertical_s1.py`` 同一约定。
"""

from __future__ import annotations

import inspect
import json
import os
from pathlib import Path

import numpy as np
import pytest

from climbing.coupling.ankle_supination import (
    DEFAULT_SUBTALAR_LEVER_ARM_MM,
    FIBULA_ENDS_SIGMA_C_MPA,
    critical_inversion_angle_deg,
    fibula_lateral_risk,
    fibula_lateral_stress,
    frontal_lateral_force_n,
    inversion_moment_nmm,
    supination_load,
)
from climbing.coupling.opensim_grf import ground_reaction

M_MODEL = 75.337
ROOT = Path(__file__).resolve().parents[1]
CACHED_AXIAL = ROOT / "results" / "opensim_fe" / "s1_h5_opensim.json"
S5_JSON = ROOT / "results" / "opensim_fe" / "nonvertical_s5.json"

run_open = pytest.mark.skipif(
    os.environ.get("RUN_OPEN_FE") != "1",
    reason="慢（正动力学 ~20 s）；设 RUN_OPEN_FE=1 启用",
)

_K_ARRAYS = ("t_s", "f_total_n", "f_left_n", "f_right_n")

#: 主口径腓骨远端截面（等效圆，来自 risk_1d_bending.json）。
_EQ = dict(a_section_mm2=101.93640134318377, i_mm4=826.8918876962899, c_mm=5.696258799381752)


# --------------------------------------------------------------------------
# L1 —— roll_deg（冠状地面倾角；默认 None = 旧行为）
# --------------------------------------------------------------------------
def test_signature_roll_deg_default_none() -> None:
    sig = inspect.signature(ground_reaction)
    assert sig.parameters["roll_deg"].default is None


def test_roll_deg_zero_bitwise_equals_default() -> None:
    a = ground_reaction(height_m=5.0, mass_kg=M_MODEL)
    b = ground_reaction(height_m=5.0, mass_kg=M_MODEL, roll_deg=0.0)
    assert b.ground_normal == (0.0, 1.0, 0.0)
    for k in _K_ARRAYS:
        assert np.array_equal(getattr(a, k), getattr(b, k)), k


def test_explicit_vertical_normal_bitwise_equals_default() -> None:
    a = ground_reaction(height_m=5.0, mass_kg=M_MODEL)
    c = ground_reaction(height_m=5.0, mass_kg=M_MODEL, ground_normal=(0.0, 1.0, 0.0))
    for k in _K_ARRAYS:
        assert np.array_equal(getattr(a, k), getattr(c, k)), k


def test_roll_deg_coronal_normal_and_vector() -> None:
    r = ground_reaction(height_m=5.0, mass_kg=M_MODEL, roll_deg=20.0)
    assert r.ground_normal == pytest.approx(
        (0.0, np.cos(np.deg2rad(20.0)), np.sin(np.deg2rad(20.0)))
    )
    assert not r.is_vertical
    v = r.force_vec_n("r")
    # 标量大小保留；横向 z = f·sinθ（矢状 x=0）
    assert np.allclose(np.linalg.norm(v, axis=1), r.f_right_n, rtol=1e-12)
    assert np.allclose(v[:, 2], r.f_right_n * np.sin(np.deg2rad(20.0)), rtol=1e-12)
    assert np.allclose(v[:, 0], 0.0)


def test_tilt_plus_roll_compose() -> None:
    x = ground_reaction(height_m=5.0, mass_kg=M_MODEL, tilt_deg=10.0, roll_deg=5.0)
    th, ph = np.deg2rad(10.0), np.deg2rad(5.0)
    n = np.array([np.sin(th), np.cos(th) * np.cos(ph), np.cos(th) * np.sin(ph)])
    n = n / np.linalg.norm(n)
    assert x.ground_normal == pytest.approx(tuple(n))


def test_roll_default_force_vector_pure_vertical() -> None:
    g = ground_reaction(height_m=5.0, mass_kg=M_MODEL)
    v = g.force_vec_n("r")
    assert np.all(v[:, 0] == 0.0) and np.all(v[:, 2] == 0.0)


# --------------------------------------------------------------------------
# 旋后载荷模型（F_lat / M_inv）
# --------------------------------------------------------------------------
def test_frontal_lateral_force_beta0_is_fz() -> None:
    assert frontal_lateral_force_n((10.0, 1000.0, 250.0), 0.0) == pytest.approx(250.0)


def test_frontal_lateral_force_vertical_redirect() -> None:
    f = 21000.0
    assert frontal_lateral_force_n((0.0, f, 0.0), 30.0) == pytest.approx(f * 0.5, rel=1e-12)
    assert frontal_lateral_force_n((0.0, f, 0.0), 0.0) == 0.0


def test_inversion_moment_definition() -> None:
    assert inversion_moment_nmm(500.0, 30.0) == pytest.approx(15000.0)
    assert inversion_moment_nmm(-500.0, 30.0) == pytest.approx(15000.0)


def test_supination_load_default_is_zero() -> None:
    load = supination_load((0.0, 21000.0, 0.0))
    assert load.supination_deg == 0.0
    assert load.lever_arm_mm == DEFAULT_SUBTALAR_LEVER_ARM_MM
    assert load.f_lat_n == 0.0
    assert load.m_inv_nmm == 0.0
    assert load.is_zero


def test_supination_load_beta15_values() -> None:
    load = supination_load((0.0, 21000.0, 0.0), 15.0, 30.0)
    assert load.f_lat_n == pytest.approx(21000.0 * np.sin(np.deg2rad(15.0)), rel=1e-12)
    assert load.m_inv_nmm == pytest.approx(load.f_lat_n * 30.0, rel=1e-12)
    assert load.m_inv_nm == pytest.approx(load.m_inv_nmm / 1000.0, rel=1e-12)
    assert not load.is_zero


def test_supination_load_bad_vector_raises() -> None:
    with pytest.raises(ValueError):
        supination_load((0.0, 1.0))
    with pytest.raises(ValueError):
        supination_load((0.0, np.nan, 0.0))


# --------------------------------------------------------------------------
# 腓骨远端 1D 侧向应力 / 风险
# --------------------------------------------------------------------------
def test_fibula_lateral_stress_known_values() -> None:
    # beta=15°, f=21000 N, d=30 mm, 等效圆截面
    load = supination_load((0.0, 21000.0, 0.0), 15.0, 30.0)
    s = fibula_lateral_stress(load.f_lat_n, load.m_inv_nmm, **_EQ)
    exp_bend = load.m_inv_nmm * _EQ["c_mm"] / _EQ["i_mm4"]
    exp_shear = load.f_lat_n / _EQ["a_section_mm2"]
    assert s["sigma_bend_mpa"] == pytest.approx(exp_bend, rel=1e-12)
    assert s["sigma_shear_mpa"] == pytest.approx(exp_shear, rel=1e-12)
    assert s["sigma_lateral_mpa"] == pytest.approx(exp_bend + exp_shear, rel=1e-12)
    assert s["sigma_vonmises_mpa"] > 0.0


def test_fibula_lateral_risk_zero_at_default() -> None:
    load = supination_load((0.0, 21000.0, 0.0))
    r = fibula_lateral_risk(load.f_lat_n, load.m_inv_nmm, **_EQ)
    assert r["risk"] == 0.0
    assert r["fracture"] is False


def test_fibula_lateral_risk_crosses_at_beta15() -> None:
    load = supination_load((0.0, 21000.0, 0.0), 15.0, 30.0)
    r = fibula_lateral_risk(load.f_lat_n, load.m_inv_nmm, **_EQ)
    assert r["risk"] >= 1.0
    assert r["fracture"] is True
    # 阈值就用 fibula_ends σ_c=70
    assert r["sigma_c_mpa"] == FIBULA_ENDS_SIGMA_C_MPA


def test_fibula_lateral_risk_threshold_scaling() -> None:
    load = supination_load((0.0, 21000.0, 0.0), 15.0, 30.0)
    r70 = fibula_lateral_risk(load.f_lat_n, load.m_inv_nmm, sigma_c_mpa=70.0, **_EQ)
    r140 = fibula_lateral_risk(load.f_lat_n, load.m_inv_nmm, sigma_c_mpa=140.0, **_EQ)
    assert r70["risk"] == pytest.approx(2.0 * r140["risk"], rel=1e-12)


def test_critical_inversion_angle_monotone_bracket() -> None:
    b = critical_inversion_angle_deg((0.0, 21000.0, 0.0), **_EQ)
    assert b is not None and 0.0 < b < 45.0
    # 在临界角处 risk≈1
    r = fibula_lateral_risk(*_load_m_beta(21000.0, b), **_EQ)
    assert r["risk"] == pytest.approx(1.0, rel=1e-6)


def test_critical_inversion_angle_none_when_below() -> None:
    # 极小力 → 45° 内不越阈
    assert critical_inversion_angle_deg((0.0, 100.0, 0.0), **_EQ) is None


def _load_m_beta(f_vert: float, beta_deg: float, arm: float = DEFAULT_SUBTALAR_LEVER_ARM_MM):
    load = supination_load((0.0, f_vert, 0.0), beta_deg, arm)
    return load.f_lat_n, load.m_inv_nmm


# --------------------------------------------------------------------------
# 产物一致性（若已生成脚本产物）：默认格点必须是零横向载荷
# --------------------------------------------------------------------------
@pytest.mark.skipif(not S5_JSON.is_file(), reason="nonvertical_s5.json 未生成")
def test_artifact_default_rows_have_zero_lateral_load() -> None:
    data = json.loads(S5_JSON.read_text(encoding="utf-8"))
    assert data["interface_regression"]["roll_deg_0_bitwise_equals_default"] is True
    assert data["interface_regression"]["supination_beta0_is_zero"] is True
    for r in data["scenario_hard"]:
        if r["supination_deg"] == 0.0:
            assert r["F_lat_n"] == 0.0
            assert r["M_inv_nmm"] == 0.0
            assert r["fibula_lateral_risk"] == 0.0
            assert r["fibula_lateral_fracture"] is False


@pytest.mark.skipif(not S5_JSON.is_file(), reason="nonvertical_s5.json 未生成")
def test_artifact_grid_and_threshold() -> None:
    data = json.loads(S5_JSON.read_text(encoding="utf-8"))
    heights = sorted({r["height_m"] for r in data["scenario_hard"]})
    angles = sorted({r["supination_deg"] for r in data["scenario_hard"]})
    assert heights == [2.0, 3.0, 3.5, 4.5]
    assert angles == [0.0, 15.0, 30.0]
    # 非零旋后在全部真实高度越阈（机制结论）
    for r in data["scenario_hard"]:
        if r["supination_deg"] >= 15.0:
            assert r["fibula_lateral_fracture"] is True


# --------------------------------------------------------------------------
# HARD GATE —— 默认路径数值复现既有轴向结果（慢，默认跳过）
# --------------------------------------------------------------------------
@run_open
def test_dead_drop_default_reproduces_cached_axial() -> None:
    from climbing.coupling.joint_loads import subtalar_reaction
    from climbing.coupling.opensim_fall import run_dead_drop

    assert CACHED_AXIAL.is_file(), f"缺少轴向回归基准 {CACHED_AXIAL}"
    cached = json.loads(CACHED_AXIAL.read_text(encoding="utf-8"))
    grf = ground_reaction(height_m=5.0, mass_kg=M_MODEL)
    assert grf.impulse_ns == pytest.approx(cached["impulse_ns"], rel=1e-9)
    fall = run_dead_drop(grf)
    jr = subtalar_reaction(fall, grf, side="r")
    assert jr.peak_vertical_n == pytest.approx(cached["subtalar_peak_vertical_n"], rel=1e-6)
    assert jr.peak_force_n == pytest.approx(cached["subtalar_peak_force_n"], rel=1e-6)
