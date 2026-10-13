"""G4 · 压缩相压溃本构 —— 本构单元（§2.3 五条）+ 不变量 + 复现性测试。

对应冻结设计 ``docs/methods/G4_压缩相本构设计.md`` §5/§7：
1. 本构五条：平台平坦 / 压实斜率 / C⁰ 连续 / 滞回面积精确值 / 永久变形；
2. 能量守恒（保守配置 + 塑性耗散时程闭合）；
3. 不变量：旋转→0 收敛、塑性应变单调不减、接触力非负、复现性；
4. 结果字段 + meta（``xi_r`` / ``crush_active``）。

本文件**只** import 新增的 ``climbing.chain2d_crush``；不改任何既有文件。
**不测 V8a/V8b 的 PASS**：冻结模型如实判定走 ``results/chain_crush_probe.json``
（跑后 PASS/FAIL 由报告/探针如实落账，测试不替判定）。
"""

from __future__ import annotations

import numpy as np
import pytest

from climbing.chain2d import (
    SEGMENTS_6,
    segment_local_peak,
)
from climbing.chain2d_crush import (
    CrushChainConfig,
    CrushParams,
    CrushTissue,
    PLANAR_CRUSH,
    _crush_sigma_and_rate,
    crush_dissipated_energy_j,
    rebound_metrics,
    simulate_crush_chain_fall,
)
from climbing.pad import POSTURES

#: 能量守恒容差（相对）。保守配置（无关节/接触阻尼、无摩擦、无轴向阻尼）
#: 下 Radau 机械能漂移 + 塑性耗散时程闭合误差；实测 ≈0.0x%，留 2% 余量。
ENERGY_TOL = 0.02

P = CrushParams()
T = P.thickness_m
EY, XD, EPL, ED = P.xi_yield, P.xi_dense, P.sigma_plateau_pa, P.E_dense_pa


# --------------------------------------------------------------------------
# 本构单元（§2.3 五条）
# --------------------------------------------------------------------------
def test_constitutive_plateau_is_flat():
    """平台期（ξ_y<ξ≤ξ_d）应力≈σ_pl：斜率 < 0.5%·E。"""
    xs = np.linspace(EY * 1.05, XD * 0.95, 60)
    sigs = np.array([P.flow_pa(x) for x in xs])
    assert np.max(sigs) - np.min(sigs) < 0.005 * P.E_pa
    assert np.allclose(sigs, EPL, rtol=1e-6)


def test_constitutive_dense_slope_is_ed():
    """压实段斜率 = E_d（冻结 20 MPa）。"""
    x1, x2 = XD + 0.05, XD + 0.25
    s1, s2 = P.flow_pa(x1), P.flow_pa(x2)
    assert s1 == pytest.approx(EPL + ED * (x1 - XD), rel=1e-9)
    assert s2 == pytest.approx(EPL + ED * (x2 - XD), rel=1e-9)


def test_constitutive_c0_continuous():
    """加载骨架 C⁰：在 ξ_y 与 ξ_d 处左/右极限连续。"""
    for xc in (EY, XD):
        l = P.flow_pa(xc - 1e-10)
        r = P.flow_pa(xc + 1e-10)
        assert abs(l - r) < 1e-6 * P.E_pa


def test_constitutive_hysteresis_area_exact():
    """完整加载-卸载循环滞回面积 = σ_pl·(ξ_peak−ξ_y)（精确标定 18750 Pa·m）。

    路径：加载 0→ξ_peak=0.20（弹性+平台），卸载沿 E 回 ξ_r=ξ_peak−ξ_y。
    解析：加载功 = ½E·ξ_y² + σ_pl·(ξ_peak−ξ_y)；卸载释放 = ½σ_pl²/E =
    ½E·ξ_y² ⇒ 净耗散 = σ_pl·(ξ_peak−ξ_y)。
    """
    xpeak = 0.20
    xs = np.linspace(0.0, xpeak, 4000)
    sig_load = np.array([P.flow_pa(x) for x in xs])
    xr = xpeak - EY
    xs_u = np.linspace(xr, xpeak, 4000)         # 升序：卸载释放功为正
    sig_unload = P.E_pa * (xs_u - xr)
    w_in = np.trapezoid(sig_load, xs)
    w_out = np.trapezoid(sig_unload, xs_u)
    net = w_in - w_out
    assert net == pytest.approx(EPL * (xpeak - EY), rel=1e-4)
    assert net > 0.0


def test_constitutive_permanent_strain_and_rates():
    """σ 随 ξ_r 卸载沿 E 回残应变；加载 dr=vz/t、卸载 dr=0。

    注意：z 向下为正 ⇒ **vz>0 = 继续压缩（加载）**，vz<0 = 回弹（卸载）。
    """
    xpeak = 0.20
    xr = xpeak - EY
    # 加载（vz>0）在峰顶：σ=σ_pl、dr=vz/t
    sig, dr = _crush_sigma_and_rate(T * xpeak, 0.0, +1.0, P)
    assert sig == pytest.approx(EPL, rel=1e-3)
    assert dr == pytest.approx(1.0 / T, rel=1e-6)
    # 卸载（vz<0）从峰顶：σ 沿 E 回到 ξ_r 处为 0；dr=0
    sig_u, dr_u = _crush_sigma_and_rate(T * xpeak, xr, -1.0, P)
    assert dr_u == 0.0
    assert sig_u == pytest.approx(EPL, rel=1e-3)
    # 完全卸载到 ξ_r 之后（comp < ξ_r·t）：slack（σ=0）
    sig_r, dr_r = _crush_sigma_and_rate(T * xr * 0.5, xr, -1.0, P)
    assert sig_r == pytest.approx(0.0, abs=1e-9)
    assert dr_r == 0.0


def test_crush_tissue_force_interface():
    """CrushTissue.force_n 无状态退化为弹性骨架（诊断接口）。"""
    ct = CrushTissue(params=P)
    f = ct.force_n(compression_m=T * EY, area_m2=0.05, xi_r=0.0, vz=-1.0)
    assert f == pytest.approx(P.E_pa * EY * 0.05, rel=1e-3)


# --------------------------------------------------------------------------
# 能量守恒（保守配置 + 塑性耗散时程闭合）
# --------------------------------------------------------------------------
def _total_energy_crush(r, cfg: CrushChainConfig) -> np.ndarray:
    """KE + 重力势 + 关节弹性势（含止挡）+ 接触弹性势 + 累计塑性耗散。"""
    ctx = _build_like_ctx(cfg)
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
        0.5 * ctx.k_comp_stop[None, :] * np.maximum(comp - ctx.comp_limit[None, :], 0.0) ** 2,
        axis=1,
    )
    pe_j += np.sum(
        0.5 * ctx.k_ext_stop[None, :] * np.maximum(ext - ctx.ext_limit[None, :], 0.0) ** 2,
        axis=1,
    )
    # 接触：可恢复弹性势能 ½E·(ξ−ξ_r)²·A·t（ξ_r 状态驱动）
    z = np.maximum(r.node_z, 0.0)
    xir = r.meta["xi_r"]
    xi = z / T
    A = r.contact_area_m2
    pe_c = 0.5 * P.E_pa * np.sum(np.maximum(xi - xir, 0.0) ** 2 * A * T, axis=1)
    # 累计塑性耗散时程：E_diss[k+1] = E_diss[k] + Σ ½(F[k]+F[k+1])·Δcomp_r
    # （梯形平均力：压溃常发生在单采样步内，F[k]·Δ 会系统性低估耗散）
    d_comp_r = np.diff(xir, axis=0) * T
    F = r.contact_force_n
    e_diss = np.zeros(r.t_s.size)
    for k in range(r.t_s.size - 1):
        dw = float(np.sum(
            0.5 * (F[k] + F[k + 1]) * np.maximum(d_comp_r[k], 0.0)
        ))
        e_diss[k + 1] = e_diss[k] + dw
    return ke + pe_g + pe_j + pe_c + e_diss


def _build_like_ctx(cfg: CrushChainConfig):
    """只取关节势能所需量的上下文（镜像 _build_crush_ctx 的几何部分）。"""
    from climbing.chain2d_crush import _build_crush_ctx
    post = POSTURES["feet-first-stiff"]
    return _build_crush_ctx(post, False, cfg)


def _conservative_crush_cfg(**kw) -> CrushChainConfig:
    base = dict(zeta_joint=0.0, zeta_contact=0.0, c_axial_n_s_per_m=0.0, mu_foot=0.0)
    base.update(kw)
    return CrushChainConfig(mass_kg=80.0, **base)


def test_energy_conserved_conservative_config_crush():
    cfg = _conservative_crush_cfg()
    r = simulate_crush_chain_fall(
        height_m=5.0, mass_kg=80.0, posture="feet-first-stiff",
        on_pad=False, config=cfg,
    )
    e = _total_energy_crush(r, cfg)
    drift = abs(e.max() - e[0]) / abs(e[0])
    drift = max(drift, abs(e.min() - e[0]) / abs(e[0]))
    assert drift < ENERGY_TOL, f"能量漂移 {drift:.4%} 超过容差 {ENERGY_TOL:.2%}"


# --------------------------------------------------------------------------
# 不变量 / 结构
# --------------------------------------------------------------------------
def test_zero_lateral_reduces_to_vertical_chain_crush():
    r = simulate_crush_chain_fall(
        height_m=5.0, mass_kg=80.0, posture="feet-first-stiff",
        on_pad=False, config=CrushChainConfig(mass_kg=80.0),
    )
    assert np.all(r.node_x == 0.0)
    assert np.all(r.joint_angle_rad == 0.0)
    assert r.node_z[:, 0].max() > 0.0


def test_plastic_strain_monotone_non_decreasing():
    """塑性应变 ξ_r 单调不减（压溃不回缩）。"""
    r = simulate_crush_chain_fall(
        height_m=5.0, mass_kg=80.0, posture="feet-first-stiff",
        on_pad=False, config=CrushChainConfig(mass_kg=80.0),
    )
    xir = r.meta["xi_r"]
    assert np.all(np.diff(xir, axis=0) >= -1e-12)


def test_contact_force_nonnegative():
    r = simulate_crush_chain_fall(
        height_m=5.0, mass_kg=80.0, posture="feet-first-stiff",
        on_pad=False, config=CrushChainConfig(mass_kg=80.0),
    )
    assert np.all(r.contact_force_n >= -1e-9)


def test_two_identical_runs_bit_identical_crush():
    cfg = CrushChainConfig(mass_kg=80.0)
    r1 = simulate_crush_chain_fall(
        height_m=5.0, mass_kg=80.0, posture="feet-first-stiff",
        on_pad=False, config=cfg,
    )
    r2 = simulate_crush_chain_fall(
        height_m=5.0, mass_kg=80.0, posture="feet-first-stiff",
        on_pad=False, config=cfg,
    )
    assert np.array_equal(r1.t_s, r2.t_s)
    assert np.array_equal(r1.node_z, r2.node_z)
    assert np.array_equal(r1.joint_force_n, r2.joint_force_n)
    assert np.array_equal(r1.contact_force_n, r2.contact_force_n)
    assert np.array_equal(r1.meta["xi_r"], r2.meta["xi_r"])


def test_result_fields_and_meta_crush():
    r = simulate_crush_chain_fall(
        height_m=5.0, mass_kg=80.0, posture="feet-first-stiff",
        on_pad=False, config=CrushChainConfig(mass_kg=80.0),
    )
    t = r.t_s.size
    assert r.node_z.shape == (t, 6)
    assert r.joint_force_n.shape == (t, 5)
    assert r.contact_force_n.shape == (t, 6)
    assert r.meta["crush_active"] is True
    assert r.meta["xi_r"].shape == (t, 6)
    t0, t1 = r.meta["impact_window_s"]
    assert 0.0 <= t0 < t1
    assert r.impact_speed_ms > 0.0
    assert r.peak_force_n > 0.0


def test_segment_local_peak_and_rebound_metrics_crush():
    r = simulate_crush_chain_fall(
        height_m=5.0, mass_kg=80.0, posture="feet-first-stiff",
        on_pad=False, config=CrushChainConfig(mass_kg=80.0),
    )
    tp, val = segment_local_peak(r, 1)
    assert tp > 0.0
    assert val > 0.0
    rb = rebound_metrics(r)
    assert rb["head_rise_m"] >= 0.0
    assert rb["foot_lift_mm"] >= 0.0
    assert crush_dissipated_energy_j(r) >= 0.0


def test_module_gate_and_config_mirror():
    """模块门（镜像 G1）+ 配置继承（纯增量字段）。"""
    assert PLANAR_CRUSH is None
    cfg = CrushChainConfig(mass_kg=80.0)
    assert cfg.crush.xi_yield == pytest.approx(EPL / P.E_pa)
    assert len(SEGMENTS_6) == 6
