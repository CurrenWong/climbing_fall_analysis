"""S1 多体段（GRF / FD / 关节反力）回归测试。

快测试只跑 ``pad.py`` 的 1D 求解（毫秒级）；正动力学很慢，放到
``RUN_OPEN_FE=1`` 门后（默认跳过，不拖慢 12 分钟的回归套件）。

不变量对应 ``docs/OpenSim_FE复现方案.md`` §七：T1 冲量、T2 重力基线。
"""

from __future__ import annotations

import os

import numpy as np
import pytest

from climbing.coupling.opensim_grf import (
    DEFAULT_POSTURE,
    IMPULSE_TOL,
    ground_reaction,
)

M_MODEL = 75.337


# --------------------------------------------------------------------------
# S1.1 —— GRF(t)
# --------------------------------------------------------------------------
@pytest.mark.parametrize("h", [1.0, 5.0, 10.0])
def test_t1_impulse_invariant(h: float) -> None:
    """T1：∫F dt / (m·v0) 应在 10% 内（截断到动量归零的耗散冲击）。"""
    grf = ground_reaction(height_m=h, mass_kg=M_MODEL)
    assert grf.impulse_ok, f"h={h}: T1 误差 {grf.impulse_rel_err*100:+.1f}%"


def test_impulse_tol_constant() -> None:
    assert IMPULSE_TOL == 0.10


def test_left_right_split_sums_to_total() -> None:
    grf = ground_reaction(height_m=5.0, mass_kg=M_MODEL)
    assert np.allclose(grf.f_left_n + grf.f_right_n, grf.f_total_n)


def test_window_starts_at_zero_and_monotone() -> None:
    grf = ground_reaction(height_m=5.0, mass_kg=M_MODEL)
    assert grf.t_s[0] == pytest.approx(0.0)
    assert np.all(np.diff(grf.t_s) > 0)


def test_peak_force_monotone_in_height() -> None:
    peaks = [ground_reaction(height_m=h, mass_kg=M_MODEL).peak_total_n
             for h in (1.0, 2.0, 5.0, 10.0)]
    assert all(b > a for a, b in zip(peaks, peaks[1:])), peaks


def test_speed_matches_freefall() -> None:
    from climbing import G

    h = 5.0
    grf = ground_reaction(height_m=h, mass_kg=M_MODEL)
    assert grf.impact_speed_ms == pytest.approx(np.sqrt(2 * G * h), rel=1e-6)


def test_bad_posture_raises() -> None:
    with pytest.raises(KeyError):
        ground_reaction(posture="not-a-posture", mass_kg=M_MODEL)


def test_bad_split_raises() -> None:
    with pytest.raises(ValueError):
        ground_reaction(mass_kg=M_MODEL, split=(0.7, 0.7))


def test_default_posture_is_feet_first() -> None:
    assert DEFAULT_POSTURE == "feet-first-stiff"


# --------------------------------------------------------------------------
# S1.2/S1.3 —— 正动力学 + 反力（慢，默认跳过）
# --------------------------------------------------------------------------
run_open = pytest.mark.skipif(
    os.environ.get("RUN_OPEN_FE") != "1",
    reason="慢（正动力学 ~20 s）；设 RUN_OPEN_FE=1 启用",
)


@run_open
def test_dead_drop_t2_freefall_residual_and_reaction() -> None:
    from climbing.coupling.joint_loads import gravity_baseline, subtalar_reaction
    from climbing.coupling.opensim_fall import LEG_JOINTS, run_dead_drop

    import opensim as os
    from climbing.coupling.opensim_fall import DEFAULT_MODEL

    # T2：刚性腿 + 无外力时，自由落体的远端残差应≈0（远小于冲击峰值）
    m = os.Model(str(DEFAULT_MODEL))
    for c in LEG_JOINTS:
        m.getCoordinateSet().get(c).set_locked(True)
    assert abs(gravity_baseline(m, side="r")) < 5.0

    grf = ground_reaction(height_m=5.0, mass_kg=M_MODEL)
    fall = run_dead_drop(grf, accuracy=1e-3, max_step_s=5e-5)
    jr = subtalar_reaction(fall, grf, side="r")

    # 反力非零、量级与施加的足底力同量级（脚很轻，反力≈外载）
    applied_peak = grf.f_right_n.max()
    assert jr.peak_force_n > 0.5 * applied_peak
    assert jr.peak_force_n < 2.0 * applied_peak
    # 左右对称
    jl = subtalar_reaction(fall, grf, side="l")
    assert jl.peak_vertical_n == pytest.approx(jr.peak_vertical_n, rel=1e-6)
