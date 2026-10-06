"""Phase-S1 非垂直落地（矢量通路）回归测试。

**硬门槛（方案 §7）**：倾角 / 水平速度 = 0（全部默认）时，必须数值复现
现有轴向结果。本文件的快测试锁定"默认路径逐位不变"的接口契约；慢测试
（``RUN_OPEN_FE=1``）直接对比缓存 ``results/opensim_fe/s1_h5_opensim.json``。

慢测试放在 ``RUN_OPEN_FE`` 门后，不拖慢常规套件 —— 与
``tests/test_opensim_fe.py`` 同一约定。
"""

from __future__ import annotations

import inspect
import json
import os
from pathlib import Path

import numpy as np
import pytest

from climbing.coupling.load_transfer import SUBTALAR_DIR, transfer
from climbing.coupling.opensim_grf import ground_reaction

M_MODEL = 75.337
ROOT = Path(__file__).resolve().parents[1]
CACHED_AXIAL = ROOT / "results" / "opensim_fe" / "s1_h5_opensim.json"

run_open = pytest.mark.skipif(
    os.environ.get("RUN_OPEN_FE") != "1",
    reason="慢（正动力学 ~20 s）；设 RUN_OPEN_FE=1 启用",
)


# --------------------------------------------------------------------------
# L1 —— GRF 方向（默认竖直 = 旧行为）
# --------------------------------------------------------------------------
def test_grf_default_normal_is_vertical() -> None:
    grf = ground_reaction(height_m=5.0, mass_kg=M_MODEL)
    assert grf.ground_normal == (0.0, 1.0, 0.0)
    assert grf.is_vertical


def test_grf_explicit_vertical_bitwise_equals_default() -> None:
    a = ground_reaction(height_m=5.0, mass_kg=M_MODEL)
    b = ground_reaction(height_m=5.0, mass_kg=M_MODEL, ground_normal=(0.0, 1.0, 0.0))
    for k in ("t_s", "f_total_n", "f_left_n", "f_right_n"):
        assert np.array_equal(getattr(a, k), getattr(b, k)), k


def test_grf_tilt_zero_bitwise_equals_default() -> None:
    a = ground_reaction(height_m=5.0, mass_kg=M_MODEL)
    c = ground_reaction(height_m=5.0, mass_kg=M_MODEL, tilt_deg=0.0)
    assert c.ground_normal == (0.0, 1.0, 0.0)
    for k in ("t_s", "f_total_n", "f_left_n", "f_right_n"):
        assert np.array_equal(getattr(a, k), getattr(c, k)), k


def test_grf_default_force_vector_is_pure_vertical() -> None:
    grf = ground_reaction(height_m=5.0, mass_kg=M_MODEL)
    v = grf.force_vec_n("r")
    assert v.shape == (grf.t_s.size, 3)
    assert np.all(v[:, 0] == 0.0) and np.all(v[:, 2] == 0.0)
    assert np.array_equal(v[:, 1], grf.f_right_n)


def test_grf_tilt_vector_magnitude_and_direction() -> None:
    grf = ground_reaction(height_m=5.0, mass_kg=M_MODEL, tilt_deg=20.0)
    assert grf.ground_normal == pytest.approx(
        (np.sin(np.deg2rad(20.0)), np.cos(np.deg2rad(20.0)), 0.0)
    )
    assert np.linalg.norm(grf.ground_normal) == pytest.approx(1.0)
    v = grf.force_vec_n("r")
    # 标量大小保留：|F_vec| == f_scalar；纵向分量 = f·cosθ
    assert np.allclose(np.linalg.norm(v, axis=1), grf.f_right_n, rtol=1e-12)
    assert np.allclose(v[:, 1], grf.f_right_n * np.cos(np.deg2rad(20.0)), rtol=1e-12)


def test_grf_bad_normal_raises() -> None:
    with pytest.raises(ValueError):
        ground_reaction(height_m=5.0, mass_kg=M_MODEL, ground_normal=(0.0, 0.0, 0.0))
    with pytest.raises(ValueError):
        ground_reaction(height_m=5.0, mass_kg=M_MODEL, ground_normal=(np.nan, 1.0, 0.0))


# --------------------------------------------------------------------------
# L3 —— FeLoadSpec 保留完整 3D wrench（默认 = 旧标量口径）
# --------------------------------------------------------------------------
def test_feloadspec_default_is_scalar_dir() -> None:
    spec = transfer(1000.0)
    assert spec.subtalar_n == 1000.0
    assert spec.subtalar_dir == SUBTALAR_DIR
    assert spec.subtalar_force is None
    assert spec.subtalar_moment is None
    # 旧口径等价物
    assert spec.force_vector_n == (0.0, -1000.0, 0.0)
    assert spec.moment_vector_nm == (0.0, 0.0, 0.0)


def test_feloadspec_vector_preserved() -> None:
    f = (300.0, -800.0, 100.0)
    m = (1.0, 2.0, 3.0)
    spec = transfer(1000.0, subtalar_force=f, subtalar_moment=m)
    assert spec.force_vector_n == f
    assert spec.moment_vector_nm == m
    # 标量字段不动（调用方显式给出的纵向值仍保留）
    assert spec.subtalar_n == 1000.0


def test_transfer_default_achilles_unchanged() -> None:
    spec = transfer(1234.0, achilles_activation=0.5)
    assert spec.subtalar_n == 1234.0
    assert spec.achilles_n == pytest.approx(0.5 * 10885.0)


# --------------------------------------------------------------------------
# L2 —— 签名默认值（默认 = 旧行为）
# --------------------------------------------------------------------------
def test_add_foot_force_signature_defaults() -> None:
    from climbing.coupling.opensim_fall import _add_foot_force

    sig = inspect.signature(_add_foot_force)
    assert sig.parameters["f_x"].default is None
    assert sig.parameters["f_z"].default is None
    assert sig.parameters["point"].default is None


def test_run_dead_drop_signature_defaults() -> None:
    from climbing.coupling.opensim_fall import run_dead_drop

    sig = inspect.signature(run_dead_drop)
    assert sig.parameters["vx0_ms"].default == 0.0
    assert sig.parameters["vz0_ms"].default == 0.0


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

    assert jr.peak_vertical_n == pytest.approx(
        cached["subtalar_peak_vertical_n"], rel=1e-6)
    assert jr.peak_force_n == pytest.approx(
        cached["subtalar_peak_force_n"], rel=1e-6)
    assert jr.peak_moment_nm == pytest.approx(
        cached["subtalar_peak_moment_nm"], rel=1e-6)
