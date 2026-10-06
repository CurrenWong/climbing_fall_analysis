"""Phase-S2 非垂直落地（单脚 / 不对称，逐足 GRF）回归测试。

**硬门槛（方案 §7）**：默认 ``split=(0.5, 0.5)``、``foot_forces=None`` 必须**逐位**
复现扩展前行为。本文件的快测试锁定接口契约（无需 OpenSim）；慢测试
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

from climbing.coupling.opensim_grf import (
    GroundReaction,
    ground_reaction,
    single_foot_split,
)

M_MODEL = 75.337
ROOT = Path(__file__).resolve().parents[1]
CACHED_AXIAL = ROOT / "results" / "opensim_fe" / "s1_h5_opensim.json"

run_open = pytest.mark.skipif(
    os.environ.get("RUN_OPEN_FE") != "1",
    reason="慢（正动力学 ~20 s）；设 RUN_OPEN_FE=1 启用",
)

_K_ARRAYS = ("t_s", "f_total_n", "f_left_n", "f_right_n")


def _bitwise_equal(a: GroundReaction, b: GroundReaction) -> bool:
    return all(np.array_equal(getattr(a, k), getattr(b, k)) for k in _K_ARRAYS)


# --------------------------------------------------------------------------
# 接口：签名默认 = 旧行为
# --------------------------------------------------------------------------
def test_signature_foot_forces_default_none() -> None:
    sig = inspect.signature(ground_reaction)
    assert sig.parameters["foot_forces"].default is None
    assert sig.parameters["split"].default == (0.5, 0.5)


def test_single_foot_split_values() -> None:
    assert single_foot_split("l") == (1.0, 0.0)
    assert single_foot_split("r") == (0.0, 1.0)
    with pytest.raises(ValueError):
        single_foot_split("x")


# --------------------------------------------------------------------------
# 默认回归（快）：默认 == 显式 (0.5,0.5) == foot_forces 重构 —— 逐位
# --------------------------------------------------------------------------
def test_default_bitwise_equals_explicit_half_half() -> None:
    a = ground_reaction(height_m=5.0, mass_kg=M_MODEL)
    b = ground_reaction(height_m=5.0, mass_kg=M_MODEL, split=(0.5, 0.5))
    assert _bitwise_equal(a, b)
    assert a.per_foot_forces is False
    assert a.loaded_sides == ("l", "r")


def test_foot_forces_reconstructs_default_bitwise() -> None:
    a = ground_reaction(height_m=5.0, mass_kg=M_MODEL)
    rec = ground_reaction(
        height_m=5.0, mass_kg=M_MODEL,
        foot_forces=((a.t_s, a.f_total_n * 0.5), (a.t_s, a.f_total_n * 0.5)),
    )
    assert _bitwise_equal(a, rec)
    assert rec.per_foot_forces is True


# --------------------------------------------------------------------------
# 单脚（split）
# --------------------------------------------------------------------------
def test_single_foot_right_only_right_loaded() -> None:
    g = ground_reaction(height_m=5.0, mass_kg=M_MODEL, split=single_foot_split("r"))
    assert g.loaded_sides == ("r",)
    assert np.all(g.f_left_n == 0.0)
    assert np.array_equal(g.f_right_n, g.f_total_n)
    assert g.peak_share == (0.0, 1.0)
    assert g.is_vertical  # 单脚仍是竖直方向 → 走旧的逐足施加路径


def test_single_foot_left_only_left_loaded() -> None:
    g = ground_reaction(height_m=5.0, mass_kg=M_MODEL, split=single_foot_split("l"))
    assert g.loaded_sides == ("l",)
    assert np.all(g.f_right_n == 0.0)
    assert np.array_equal(g.f_left_n, g.f_total_n)


# --------------------------------------------------------------------------
# 校验
# --------------------------------------------------------------------------
def test_split_nonnormalized_raises() -> None:
    with pytest.raises(ValueError):
        ground_reaction(height_m=5.0, mass_kg=M_MODEL, split=(0.3, 0.3))


def test_split_negative_raises() -> None:
    with pytest.raises(ValueError):
        ground_reaction(height_m=5.0, mass_kg=M_MODEL, split=(-0.1, 1.1))


def test_foot_forces_must_be_pair() -> None:
    with pytest.raises(ValueError):
        ground_reaction(height_m=5.0, mass_kg=M_MODEL, foot_forces=(0.0,))


def test_foot_forces_bad_shape_raises() -> None:
    with pytest.raises(ValueError):
        ground_reaction(
            height_m=5.0, mass_kg=M_MODEL,
            foot_forces=((np.zeros((2, 2)), np.zeros((2, 2))), 0.0),
        )
    with pytest.raises(ValueError):
        ground_reaction(
            height_m=5.0, mass_kg=M_MODEL,
            foot_forces=((np.array([0.0, 1.0]), np.array([1.0, 2.0, 3.0])), 0.0),
        )


def test_foot_forces_negative_raises() -> None:
    with pytest.raises(ValueError):
        ground_reaction(
            height_m=5.0, mass_kg=M_MODEL,
            foot_forces=(np.array([1.0, -1.0, 1.0]), 0.0),
        )


# --------------------------------------------------------------------------
# 泛化能力：逐足独立时程 / 非对称占比
# --------------------------------------------------------------------------
def test_foot_forces_asymmetric_share_and_total() -> None:
    a = ground_reaction(height_m=5.0, mass_kg=M_MODEL)
    g = ground_reaction(
        height_m=5.0, mass_kg=M_MODEL,
        foot_forces=((a.t_s, a.f_total_n * 0.7), (a.t_s, a.f_total_n * 0.3)),
    )
    assert g.per_foot_forces is True
    assert g.peak_share == pytest.approx((0.7, 0.3))
    assert g.loaded_sides == ("l", "r")
    # 总力逐位守恒（0.7+0.3 的浮点误差内）
    assert np.allclose(g.f_total_n, a.f_total_n, rtol=0.0, atol=1e-6)
    assert g.impulse_ns == pytest.approx(float(np.trapezoid(a.f_total_n, a.t_s)), rel=1e-9)


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
