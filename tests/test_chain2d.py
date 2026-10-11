"""G1 · 6 段平面软关节链 —— 不变量 + 复现性测试。

对应冻结设计 ``docs/methods/G1_多段链设计.md`` §6：
1. 旋转→0 收敛回纯压缩（不变量 ①）；
2. 能量守恒（容差在此文档化）；
3. 局部组织峰随高度单调提前（不变量 ④）；
4. 两次相同运行给相同结果（确定性、种子无关）；
5. 结果字段齐备 + 冲击窗口存在。

本文件**只** import 新增的 ``climbing.chain2d``；不改任何既有文件。
"""

from __future__ import annotations

import numpy as np
import pytest

from climbing.chain2d import (
    PLANAR_CHAIN,
    SEGMENTS_6,
    PlanarChainConfig,
    PlanarChainResult,
    _build_ctx,
    _contact_force,
    _kink_forces,
    segment_local_peak,
    simulate_planar_chain_fall,
)
from climbing.pad import POSTURES

#: 能量守恒容差（相对）。**保守配置**（无阻尼、无摩擦）下 Radau 的机械能
#: 漂移。实测 ≈0.00%（bug 修复后）；此处留 2% 余量以免数值平台差异。
ENERGY_TOL = 0.02

#: 探针高度集（冻结设计 §6）。
HEIGHTS = (2.0, 3.0, 5.0, 10.0, 20.0)


def _conservative_cfg(**kw) -> PlanarChainConfig:
    """保守配置：无接触阻尼、无关节阻尼、无摩擦 ⇒ 机械能应守恒。"""
    base = dict(zeta_joint=0.0, zeta_contact=0.0, c_axial_n_s_per_m=0.0, mu_foot=0.0)
    base.update(kw)
    return PlanarChainConfig(mass_kg=80.0, **base)


def _total_energy(r: PlanarChainResult, cfg: PlanarChainConfig) -> np.ndarray:
    """机械能 = KE + 重力势 + 关节弹性势（含止挡）+ 接触弹性势。"""
    post = POSTURES["feet-first-stiff"]
    ctx = _build_ctx(post, False, cfg)
    m = r.meta["node_mass_kg"]
    g = cfg.g
    ke = 0.5 * np.sum(m[None, :] * (r.node_vx**2 + r.node_vz**2), axis=1)
    pe_g = -np.sum(m[None, :] * g * r.node_z, axis=1)
    dx = r.node_x[:, 1:] - r.node_x[:, :-1]
    dz = r.node_z[:, 1:] - r.node_z[:, :-1]
    ell = np.maximum(np.hypot(dx, dz), 1e-9)
    comp = ctx.rest[None, :] - ell
    ext = ell - ctx.rest[None, :]
    pe_j = 0.5 * np.sum(ctx.k_ax[None, :] * (ell - ctx.rest[None, :]) ** 2, axis=1)
    pe_j += np.sum(
        0.5 * ctx.k_comp_stop[None, :]
        * np.maximum(comp - ctx.comp_limit[None, :], 0.0) ** 2, axis=1
    )
    pe_j += np.sum(
        0.5 * ctx.k_ext_stop[None, :]
        * np.maximum(ext - ctx.ext_limit[None, :], 0.0) ** 2, axis=1
    )
    pe_c = np.zeros(r.t_s.size)
    for i in range(r.t_s.size):
        z = r.node_z[i, 0]
        if z > 0.0:
            xs = np.linspace(0.0, z, 80)
            fs = np.array([
                _contact_force(x, post.area_at(x / post.collapse_ref_m), ctx.pad)
                for x in xs
            ])
            pe_c[i] = np.trapezoid(fs, xs)
    return ke + pe_g + pe_j + pe_c


# --------------------------------------------------------------------------
# 结构 / 人体测量学
# --------------------------------------------------------------------------
def test_segments_are_six_and_mass_fractions_sum_to_one():
    assert len(SEGMENTS_6) == 6
    total = sum(s.mass_frac for s in SEGMENTS_6)
    assert total == pytest.approx(1.0, abs=1e-6)
    for s in SEGMENTS_6:
        assert s.length_m > 0.0
        assert s.k_axial_n_per_m > 0.0
        assert s.k_torsion_nm_per_rad >= 0.0
        assert s.contact_area_m2 > 0.0
        assert s.tissue_area_m2 > 0.0
        assert s.tissue_modulus_pa > 0.0


def test_gate_default_none():
    # 门默认关闭（镜像 RIGID_STANDING_CHAIN 范式）：不影响默认路径。
    assert PLANAR_CHAIN is None


# --------------------------------------------------------------------------
# 不变量 ①：旋转→0 收敛回纯压缩
# --------------------------------------------------------------------------
def test_zero_lateral_reduces_to_vertical_chain():
    r = simulate_planar_chain_fall(
        height_m=5.0, mass_kg=80.0, posture="feet-first-stiff",
        on_pad=False, config=PlanarChainConfig(mass_kg=80.0),
    )
    # 无横向扰动 ⇒ 横向坐标与关节转角**精确**为零（纯压缩链）。
    assert np.all(r.node_x == 0.0)
    assert np.all(r.joint_angle_rad == 0.0)
    # 竖向轨迹非平凡（确实发生了塌陷）。
    assert r.node_z[:, 0].max() > 0.0


def test_bending_force_is_gradient_of_quadratic_potential():
    """`_kink_forces` 必须是 V=½kδ² 的负梯度（抗弯弹簧一致性）。"""
    a = np.array([0.0, 0.0])
    b = np.array([0.03, -0.2])
    c = np.array([-0.02, -0.5])
    k = 350.0

    def v(a, b, c):
        p = b - a
        q = c - b
        s = p[0] * q[1] - p[1] * q[0]
        cc = p[0] * q[0] + p[1] * q[1]
        return 0.5 * k * np.arctan2(s, cc) ** 2

    fa, fb, fc = _kink_forces(a, b, c, k)
    eps = 1e-7

    def num_grad(which):
        g = np.zeros(2)
        for d in range(2):
            pp = [a.copy(), b.copy(), c.copy()]
            pp[which][d] += eps
            pm = [a.copy(), b.copy(), c.copy()]
            pm[which][d] -= eps
            g[d] = -(v(*pp) - v(*pm)) / (2.0 * eps)
        return g

    assert np.allclose(fa, num_grad(0), atol=1e-4)
    assert np.allclose(fb, num_grad(1), atol=1e-4)
    assert np.allclose(fc, num_grad(2), atol=1e-4)
    # 内力配平（牛顿第三定律）。
    assert np.allclose(fa + fb + fc, 0.0, atol=1e-9)


# --------------------------------------------------------------------------
# 能量守恒（容差文档化）
# --------------------------------------------------------------------------
def test_energy_conserved_conservative_config():
    cfg = _conservative_cfg()
    r = simulate_planar_chain_fall(
        height_m=5.0, mass_kg=80.0, posture="feet-first-stiff",
        on_pad=False, config=cfg,
    )
    e = _total_energy(r, cfg)
    drift = abs(e.max() - e[0]) / abs(e[0])
    drift = max(drift, abs(e.min() - e[0]) / abs(e[0]))
    assert drift < ENERGY_TOL, f"能量漂移 {drift:.4%} 超过容差 {ENERGY_TOL:.2%}"


# --------------------------------------------------------------------------
# 复现性 / 确定性
# --------------------------------------------------------------------------
def test_two_identical_runs_bit_identical():
    cfg = PlanarChainConfig(mass_kg=80.0)
    r1 = simulate_planar_chain_fall(
        height_m=5.0, mass_kg=80.0, posture="feet-first-stiff",
        on_pad=False, config=cfg,
    )
    r2 = simulate_planar_chain_fall(
        height_m=5.0, mass_kg=80.0, posture="feet-first-stiff",
        on_pad=False, config=cfg,
    )
    assert np.array_equal(r1.t_s, r2.t_s)
    assert np.array_equal(r1.node_z, r2.node_z)
    assert np.array_equal(r1.joint_force_n, r2.joint_force_n)
    assert np.array_equal(r1.contact_force_n, r2.contact_force_n)


# --------------------------------------------------------------------------
# 结果字段
# --------------------------------------------------------------------------
def test_result_fields_present_and_shaped():
    r = simulate_planar_chain_fall(
        height_m=5.0, mass_kg=80.0, posture="feet-first-stiff",
        on_pad=False, config=PlanarChainConfig(mass_kg=80.0),
    )
    t = r.t_s.size
    assert r.node_x.shape == (t, 6)
    assert r.node_z.shape == (t, 6)
    assert r.node_vx.shape == (t, 6)
    assert r.node_vz.shape == (t, 6)
    assert r.joint_force_n.shape == (t, 5)
    assert r.joint_compression_m.shape == (t, 5)
    assert r.joint_angle_rad.shape == (t, 5)
    assert r.contact_force_n.shape == (t, 6)
    assert r.contact_area_m2.shape == (t, 6)
    assert "impact_window_s" in r.meta
    t0, t1 = r.meta["impact_window_s"]
    assert 0.0 <= t0 < t1
    assert r.impact_speed_ms > 0.0
    assert r.peak_force_n > 0.0


def test_contact_force_zero_when_above_ground():
    r = simulate_planar_chain_fall(
        height_m=5.0, mass_kg=80.0, posture="feet-first-stiff",
        on_pad=False, config=PlanarChainConfig(mass_kg=80.0),
    )
    # 非接触节点在窗口前/后接触力应为 0（上面 5 个节点全程离地）。
    assert np.all(r.contact_force_n[:, 1:] == 0.0)


# --------------------------------------------------------------------------
# 不变量 ④：局部组织峰随高度单调提前 + 无塌陷中回弹
# --------------------------------------------------------------------------
def test_tibia_local_peak_monotone_with_height():
    cfg = PlanarChainConfig(mass_kg=80.0)
    tps = []
    for h in HEIGHTS:
        r = simulate_planar_chain_fall(
            height_m=h, mass_kg=80.0, posture="feet-first-stiff",
            on_pad=False, config=cfg,
        )
        tp, val = segment_local_peak(r, 1)   # 1 = 小腿-胫骨段
        tps.append(tp)
        assert val > 0.0
    # 单调不增（更高 ⇒ 更早）。
    for i in range(1, len(tps)):
        assert tps[i] <= tps[i - 1] + 1e-9


def test_force_monotone_through_collapse():
    """塌陷（压缩）相接触力单调上升 —— 修掉 bodies.py 的足部中途回弹。"""
    r = simulate_planar_chain_fall(
        height_m=5.0, mass_kg=80.0, posture="feet-first-stiff",
        on_pad=False, config=PlanarChainConfig(mass_kg=80.0),
    )
    f = r.contact_force_n.sum(axis=1)
    ipk = int(np.argmax(f))
    tol = 1e-6 * max(1.0, f.max())
    assert np.all(np.diff(f[: ipk + 1]) >= -tol)


def test_segment_local_peak_accepts_name_index():
    r = simulate_planar_chain_fall(
        height_m=5.0, mass_kg=80.0, posture="feet-first-stiff",
        on_pad=False, config=PlanarChainConfig(mass_kg=80.0),
    )
    tp_by_int, v_i = segment_local_peak(r, 1)
    tp_by_seg, v_s = segment_local_peak(r, SEGMENTS_6[1])
    assert tp_by_int == tp_by_seg
    assert v_i == v_s
