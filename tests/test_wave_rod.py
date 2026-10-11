"""G1 · 1D 应力波杆（``climbing.wave_rod``）—— 一致性 / 复现性测试。

对应冻结设计 ``docs/methods/G1_多段链设计.md`` §7 决策 6 的预注册备选件。
本件验证的是**模型自洽性**（教材级不变量），不是"凑 4.6 ms"：

1. 轴向波速 ``c = sqrt(E/rho)``；
2. 峰时 ``=`` 往返时间 ``2L/c``（预注册第一原理定义）；
3. 峰时随 ``L`` 线性（``t ∝ L/c``）；
4. 峰时**亚毫秒**（预注册"过冲"预期）且**与高度无关**（波时标不含 v0）；
5. 确定性 / 可复现（两次运行逐位相同）；
6. 无 NaN / Inf；
7. 字段齐备、门控默认关闭。

本文件**只** import 新增的 ``climbing.wave_rod``；不改任何既有文件。
"""

from __future__ import annotations

import numpy as np
import pytest

from climbing import wave_rod as wr
from climbing.wave_rod import (
    WAVE_ROD,
    WaveRod,
    WaveRodResult,
    simulate_tibia_wave,
    wave_observable,
)


# --------------------------------------------------------------------------
# 结构 / 门控
# --------------------------------------------------------------------------
def test_gate_default_none():
    # 门默认关闭（镜像 RIGID_STANDING_CHAIN / PLANAR_CHAIN 范式）。
    assert WAVE_ROD is None


def test_wave_speed_is_sqrt_E_over_rho():
    cfg = WaveRod()
    assert cfg.wave_speed_ms == pytest.approx(
        np.sqrt(cfg.modulus_pa / cfg.density_kg_m3), rel=1e-12
    )
    # 文档值量级：皮质骨 E≈17 GPa、rho≈1900 ⇒ c≈3.0e3 m/s。
    assert cfg.wave_speed_ms == pytest.approx(3000.0, rel=0.05)


def test_result_fields_and_shapes():
    r = simulate_tibia_wave(height_m=20.0, mass_kg=80.0)
    t = r.t_s.size
    n = WaveRod().n_cells
    assert r.strain.shape == (t, n + 1)
    assert r.distal_strain.shape == (t,)
    assert r.x_m.shape == (n + 1,)
    assert r.t_s[0] == 0.0
    assert np.all(np.diff(r.t_s) > 0.0)
    assert r.peak_strain < 0.0          # 压缩为负
    assert r.peak_time_s > 0.0
    assert r.window_s > 0.0
    assert "impact_window_s" not in r.meta  # 本件用自己的 window_s
    assert r.meta["impedance_ratio_m_above_over_rod"] > 10.0  # 固定端判据


# --------------------------------------------------------------------------
# 峰时 = 往返时间 2L/c（预注册第一原理）
# --------------------------------------------------------------------------
def test_peak_time_matches_round_trip():
    cfg = WaveRod()
    r = simulate_tibia_wave(height_m=20.0, mass_kg=80.0, config=cfg)
    analytic = cfg.round_trip_s
    assert analytic == pytest.approx(2.0 * cfg.length_m / cfg.wave_speed_ms, rel=1e-12)
    assert r.analytic_peak_time_s == pytest.approx(analytic, rel=1e-12)
    # 数值峰时 ≈ 2L/c（离散误差以内）。
    assert r.peak_time_s == pytest.approx(analytic, rel=0.05)


def test_peak_time_scales_with_length():
    cfg1 = WaveRod(length_m=0.38)
    cfg2 = WaveRod(length_m=0.76)
    r1 = simulate_tibia_wave(20.0, 80.0, config=cfg1)
    r2 = simulate_tibia_wave(20.0, 80.0, config=cfg2)
    # t ∝ L/c ⇒ 加倍长度 ⇒ 加倍峰时（往返时间）。
    assert r2.round_trip_s == pytest.approx(2.0 * r1.round_trip_s, rel=1e-12)
    assert r2.peak_time_s == pytest.approx(2.0 * r1.peak_time_s, rel=0.02)


# --------------------------------------------------------------------------
# 亚毫秒 + 与高度无关（本件的核心结论）
# --------------------------------------------------------------------------
def test_peak_time_submillisecond():
    r = simulate_tibia_wave(height_m=20.0, mass_kg=80.0)
    assert r.submillisecond
    assert r.peak_time_s < 1.0e-3


def test_peak_time_height_independent():
    # 波速/长度固定 ⇒ 峰时与撞击速度 v0 无关（v0 只定标幅值）。
    t2 = simulate_tibia_wave(height_m=2.0, mass_kg=80.0).peak_time_s
    t20 = simulate_tibia_wave(height_m=20.0, mass_kg=80.0).peak_time_s
    assert t2 == t20
    # 但幅值确实随高度增大（能量口径仍在）。
    e2 = simulate_tibia_wave(height_m=2.0, mass_kg=80.0).impact_strain
    e20 = simulate_tibia_wave(height_m=20.0, mass_kg=80.0).impact_strain
    assert e20 > e2


def test_impact_strain_is_v0_over_c():
    r = simulate_tibia_wave(height_m=20.0, mass_kg=80.0)
    c = WaveRod().wave_speed_ms
    assert r.impact_speed_ms == pytest.approx(np.sqrt(2.0 * 9.80665 * 20.0), rel=1e-9)
    assert r.impact_strain == pytest.approx(r.impact_speed_ms / c, rel=1e-12)


# --------------------------------------------------------------------------
# 确定性 / 复现性 / 数值健全
# --------------------------------------------------------------------------
def test_two_identical_runs_bit_identical():
    cfg = WaveRod()
    r1 = simulate_tibia_wave(20.0, 80.0, config=cfg)
    r2 = simulate_tibia_wave(20.0, 80.0, config=cfg)
    assert np.array_equal(r1.t_s, r2.t_s)
    assert np.array_equal(r1.distal_strain, r2.distal_strain)
    assert np.array_equal(r1.strain, r2.strain)
    assert r1.peak_time_s == r2.peak_time_s
    assert r1.peak_strain == r2.peak_strain


def test_no_nan_or_inf():
    r = simulate_tibia_wave(height_m=20.0, mass_kg=80.0)
    assert np.all(np.isfinite(r.strain))
    assert np.all(np.isfinite(r.distal_strain))
    assert np.isfinite(r.peak_time_s)
    assert np.isfinite(r.peak_strain)


# --------------------------------------------------------------------------
# 波观测量函数 + 门控叠加
# --------------------------------------------------------------------------
def test_wave_observable_matches_result_fields():
    r = simulate_tibia_wave(20.0, 80.0)
    tp, ep = wave_observable(r)
    assert tp == r.peak_time_s
    assert ep == r.peak_strain


def test_gate_overlay_overrides_config(monkeypatch):
    # 门非空 ⇒ 按字段名覆盖（镜像 pad/chain2d）。
    monkeypatch.setattr(wr, "WAVE_ROD", {"length_m": 0.76}, raising=True)
    r = simulate_tibia_wave(20.0, 80.0)
    assert r.meta["length_m"] == pytest.approx(0.76)
    assert r.round_trip_s == pytest.approx(2.0 * 0.76 / r.wave_speed_ms, rel=1e-12)
    # 恢复后回到默认。
    monkeypatch.setattr(wr, "WAVE_ROD", None, raising=True)
    r_def = simulate_tibia_wave(20.0, 80.0)
    assert r_def.meta["length_m"] == pytest.approx(0.38)
