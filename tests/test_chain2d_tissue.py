"""G5 · v2 柔顺层 —— 本构单元 + 不变量 + 复现性。

对应 ``docs/methods/G5_跟垫柔顺层纠错_预注册.md`` §3/§4。
本文件**只** import 新增的 ``climbing.chain2d_tissue``；不改任何既有文件。
**不判 V8b / 骨应力**（预注册声明其不可由本件判定）。
"""

from __future__ import annotations

import numpy as np
import pytest

from climbing.chain2d import SEGMENTS_6
from climbing.chain2d_tissue import (
    PLANAR_TISSUE,
    TissueChainConfig,
    TissueParams,
    TissueV2,
    _build_tissue_ctx,
    _tissue_sigma_and_rate,
    simulate_tissue_chain_fall,
    tissue_dissipated_energy_j,
    tissue_peak_deformation_mm,
    tibia_peak_stress_mpa,
)
from climbing.pad import POSTURES

P = TissueParams()
T = P.thickness_m
E1, XIR, ED = P.E1_pa, P.xi_ref, P.E_dense_pa

#: in-vivo 锚点（PMC8925218；Gefen 2007 佐证）——**参数冻结值，测试只断言不被改**。
ANCHOR_E1_PA = 0.185e6
ANCHOR_T_M = 0.0154
ANCHOR_XI_PEAK = 0.70


# --------------------------------------------------------------------------
# 锚点冻结（防参数漂移）
# --------------------------------------------------------------------------
def test_anchors_frozen():
    assert P.E1_pa == pytest.approx(ANCHOR_E1_PA)
    assert P.thickness_m == pytest.approx(ANCHOR_T_M)
    assert P.xi_ref == pytest.approx(ANCHOR_XI_PEAK)
    # 致密化锁点 = 10.8 mm，正对照 Bennett&Ker in-vivo 8.5–11.3 mm
    assert P.xi_lock_m * 1e3 == pytest.approx(10.78, abs=0.01)


# --------------------------------------------------------------------------
# 本构单元
# --------------------------------------------------------------------------
def test_linear_segment_matches_e1():
    """线性段 σ = E₁·u（in-vivo 锚点直接量）。

    上界取 **ξ=0.55** 而非贴近 ξ_ref：实现用 tanh 过渡带 ``eps_xi=0.05``
    （E_d/E₁=108× 陡过渡需宽带才不塌 Radau 步长），故 ξ_ref 附近 ±0.05 是混合带，
    非纯线性（见 ``test_blend_band_is_regularization_only``）。
    """
    for xi in (0.1, 0.3, 0.5, 0.55):
        sig, _ = _tissue_sigma_and_rate(xi * T, 0.0, 1.0, P)
        assert sig == pytest.approx(E1 * xi, rel=1e-3)


def test_blend_band_is_regularization_only():
    """ξ_ref 附近的 tanh 混合带是**数值正则化**，不改变骨架在带外的取值。

    带内（ξ=0.69，比 ξ_ref 低 0.01 = 0.2 带宽）σ 略高于纯线性；
    带外（ξ ≤ 0.55）严格等于 E₁·ξ（见上）。
    """
    sig_in, _ = _tissue_sigma_and_rate(0.69 * T, 0.0, 1.0, P)
    assert sig_in > E1 * 0.69            # 混合带上翘（致密化贡献）
    assert sig_in < E1 * 0.69 * 1.05    # 但只是 5% 内的正则化偏差


def test_no_plateau_beyond_xi_ref_densifies():
    """越过 ξ_ref 后斜率切到 E_d（**不是**平台）。

    注意：实现用 tanh 过渡带 ``eps_xi=0.05`` 平滑（``E_d/E₁=108×`` 陡过渡会让
    Radau 步长坍缩），故须取**远越过**过渡带的点来量斜率。
    """
    sig_a, _ = _tissue_sigma_and_rate(1.00 * T, 0.0, 1.0, P)
    sig_b, _ = _tissue_sigma_and_rate(1.20 * T, 0.0, 1.0, P)
    slope = (sig_b - sig_a) / 0.20
    assert slope > 10.0 * E1           # ≫ E₁ ⇒ 已致密化，不是平台
    assert sig_a > 5.0 * E1 * XIR      # σ 远超线性段末端 ⇒ 确已致密化


def test_unload_along_e1_to_permanent_strain():
    """卸载沿 E₁ 回永久残余应变（不回弹到 0）。"""
    xir = 0.3
    sig, dr = _tissue_sigma_and_rate(0.30 * T, xir, -1.0, P)   # 恰在永久应变上
    assert sig == pytest.approx(0.0, abs=1e-6)
    assert dr == 0.0                                             # 卸载不增永久
    sig2, _ = _tissue_sigma_and_rate(0.40 * T, xir, -1.0, P)
    assert sig2 == pytest.approx(E1 * (0.40 - 0.30), rel=1e-3)


def test_permanent_strain_saturates_at_xi_ref():
    """ξ_r 饱和规则：**首版让 u≡ξ_ref 恒定 ⇒ 130 kPa 无限塌陷**（首跑 183 mm 的根因）。"""
    # 持续深压：dr 应随 xi_r → ξ_ref 而归零（cap 平滑饱和）
    drs = []
    for xir in (0.0, 0.2, 0.4, 0.6, 0.69):
        _, dr = _tissue_sigma_and_rate(0.20 * T, xir, 1.0, P)
        drs.append(dr)
    assert drs[0] > drs[2] > drs[4]          # 单调衰减
    assert drs[-1] < 0.2 * drs[0]            # 接近饱和
    # 关键：即使持续加载，可恢复弹性应变 u 仍须能越过 ξ_ref（否则退化回无限平台）
    sig_deep, _ = _tissue_sigma_and_rate((0.7316 + 0.925) * T, 0.7316, 1.0, P)
    assert sig_deep > E1 * XIR               # σ 超过线性段末端 ⇒ 致密化确实起作用


def test_xi_r_only_grows_while_linear_loading():
    """深压（u ≫ ξ_ref）时 ξ_r 冻结；卸载时冻结；线性区加载才增长。"""
    _, dr_linear = _tissue_sigma_and_rate(0.20 * T, 0.1, 1.0, P)
    _, dr_dense = _tissue_sigma_and_rate(1.50 * T, 0.5, 1.0, P)
    _, dr_unload = _tissue_sigma_and_rate(0.20 * T, 0.1, -1.0, P)
    assert dr_linear > 0.0
    assert dr_dense < 1e-2 * dr_linear       # 已进致密化 ⇒ 永久压实近似冻结
    assert dr_unload == 0.0


def test_slack_when_no_contact():
    assert _tissue_sigma_and_rate(-0.01, 0.0, -1.0, P) == (0.0, 0.0)


def test_tissue_v2_force_interface():
    ct = TissueV2(params=P)
    f = ct.force_n(compression_m=0.30 * T, area_m2=0.05, xi_r=0.0, vz=1.0)
    assert f == pytest.approx(E1 * 0.30 * 0.05, rel=1e-2)


# --------------------------------------------------------------------------
# 结构 / 门控
# --------------------------------------------------------------------------
def test_module_gate_and_config_mirror():
    assert PLANAR_TISSUE is None
    cfg = TissueChainConfig(mass_kg=80.0)
    assert cfg.tissue.E1_pa == pytest.approx(ANCHOR_E1_PA)
    assert len(SEGMENTS_6) == 6


def test_build_ctx_mirrors_chain():
    """v2 上下文的链几何/关节参数与 G4 crush 上下文逐值一致（只换接触本构）。"""
    from climbing.chain2d_crush import CrushChainConfig, _build_crush_ctx
    post = POSTURES["feet-first-stiff"]
    a = _build_tissue_ctx(post, False, TissueChainConfig(mass_kg=80.0))
    b = _build_crush_ctx(post, False, CrushChainConfig(mass_kg=80.0))
    assert np.allclose(a.k_ax, b.k_ax)
    assert np.allclose(a.rest, b.rest)
    assert np.allclose(a.area_base, b.area_base)
    assert np.allclose(a.m, b.m)
    assert a.n_nodes == b.n_nodes


# --------------------------------------------------------------------------
# ODE 集成
# --------------------------------------------------------------------------
def _sim(h=3.0, **kw):
    return simulate_tissue_chain_fall(
        h, 80.0, "feet-first-stiff", False,
        config=TissueChainConfig(mass_kg=80.0, **kw),
    )


def test_result_fields_and_meta():
    r = _sim()
    t = r.t_s.size
    assert r.node_z.shape == (t, 6)
    assert r.joint_force_n.shape == (t, 5)
    assert r.contact_force_n.shape == (t, 6)
    assert r.meta["tissue_active"] is True
    assert r.meta["xi_r"].shape == (t, 6)
    t0, t1 = r.meta["impact_window_s"]
    assert 0.0 <= t0 < t1


def test_zero_lateral_reduces_to_vertical_chain():
    r = _sim()
    assert np.all(r.node_x == 0.0)
    assert np.all(r.joint_angle_rad == 0.0)
    assert r.node_z[:, 0].max() > 0.0


def test_permanent_strain_monotone_non_decreasing():
    """ξ_r 单调不减（压实不回缩）。

    容差取 1e-6：默认容差 ``rtol=1e-4``（致密化接触刚性所致，见
    ``TissueChainConfig``）下，求解器可在个别步给出 1e-12 以上的微小负步。
    本构侧 ``dr = wv·wlin·cap·max(vz,0)/t`` 四因子皆 ≥0，故负步纯属数值噪声。
    """
    xir = _sim().meta["xi_r"]
    assert np.all(np.diff(xir, axis=0) >= -1e-6)


def test_contact_force_nonnegative():
    assert np.all(_sim().contact_force_n >= -1e-9)


def test_two_identical_runs_bit_identical():
    r1, r2 = _sim(), _sim()
    assert np.array_equal(r1.t_s, r2.t_s)
    assert np.array_equal(r1.node_z, r2.node_z)
    assert np.array_equal(r1.contact_force_n, r2.contact_force_n)
    assert np.array_equal(r1.meta["xi_r"], r2.meta["xi_r"])


def test_energy_conserved_conservative_config():
    """保守配置（无关节/接触阻尼、无摩擦、无轴向阻尼）机械能闭合 < 2%。

    守恒量 = ``KE + PE_g + PE_j + W(t)``，其中 ``W`` 是垫层**路径功时间序列**
    （``tissue_pad_work_t_j``）。

    **2026-10-15 口径更正（替换式，原判被实测证伪）**：原写「**致密化段
    （ξ > ξ_ref）的耗散尚未进入后处理账目** … 补齐需在 ``tissue_dissipated_energy_j``
    增加致密化段项」，并据此取了 h=0.5 m 以「让致密化基本不介入」。

    **实测（本文件 + `temp/g5_*.py`）**：
      1. **h=0.5 m 的峰值变形已 11.40 mm > ``t·ξ_ref`` = 10.8 mm** ⇒
         致密化在**该点上就已介入**，「取 h=0.5 使致密化不介入」的前提不成立。
      2. **按原判补上致密化滞回项后，漂移一字不变**（h=1.5 m：2.060% → 2.059%）
         ⇒ **原判被证伪**。
      3. **真因**：``tissue_dissipated_energy_j`` 返回**标量**，加进逐时刻的
         ``e(t)`` 时**各时刻同增同减、恒被抵消** ⇒ 旧测试**实际只测了 ``KE+PE``
         的漂移**，而那个漂移值 = 未计入的耗散占比。

    **替换后的闭合质量**（rtol=1e-8）：h=0.5 / 1.5 / 2.0 m ⇒
    **0.1709% / 0.2140% / 0.5062%**（旧口径 1.1331% / 2.0642% / 2.4336%）。
    残余为采样网格上的**梯形积分误差**（非物理漏项）。
    """
    from climbing.chain2d_crush import _build_crush_ctx
    from climbing.chain2d_tissue import tissue_pad_work_t_j
    cfg = TissueChainConfig(
        mass_kg=80.0, zeta_joint=0.0, zeta_contact=0.0,
        c_axial_n_s_per_m=0.0, mu_foot=0.0,
        rtol=1e-8, atol=1e-11,
    )
    ctx = _build_tissue_ctx(POSTURES["feet-first-stiff"], False, cfg)
    r = simulate_tissue_chain_fall(0.5, 80.0, "feet-first-stiff", False, config=cfg)
    m = r.meta["node_mass_kg"]
    ke = 0.5 * np.sum(m[None, :] * (r.node_vx**2 + r.node_vz**2), axis=1)
    pe_g = -np.sum(m[None, :] * cfg.g * r.node_z, axis=1)
    dx = r.node_x[:, 1:] - r.node_x[:, :-1]
    dz = r.node_z[:, 1:] - r.node_z[:, :-1]
    ell = np.maximum(np.hypot(dx, dz), 1e-9)
    comp = ctx.rest[None, :] - ell
    ext = ell - ctx.rest[None, :]
    pe_j = 0.5 * np.sum(ctx.k_ax[None, :] * (ell - ctx.rest[None, :]) ** 2, axis=1)
    pe_j += np.sum(
        0.5 * ctx.k_comp_stop[None, :] * np.maximum(comp - ctx.comp_limit[None, :], 0.0) ** 2,
        axis=1)
    pe_j += np.sum(
        0.5 * ctx.k_ext_stop[None, :] * np.maximum(ext - ctx.ext_limit[None, :], 0.0) ** 2,
        axis=1)
    e = ke + pe_g + pe_j + tissue_pad_work_t_j(r)
    drift = max(abs(e.max() - e[0]), abs(e.min() - e[0])) / abs(e[0])
    assert drift < 0.005, f"能量漂移 {drift:.4%}"     # 实测 0.1709%，留 2.9× 余量


def test_pad_work_series_and_dissipation_series():
    """``tissue_pad_work_t_j`` 按定义可独立重算；耗散序列 = 路径功 − 可恢复弹性能。

    ⚠️ **不做单调性断言**（那是臆测、不是事实）：实测耗散序列在**卸载过渡段**
    有最多 −0.044 J 的负增量（`wv`/`wb` 平滑混合期 σ 略高于 ``E₁u``，
    卸载时 ``dW < d(½E₁u²)``）⇒ 允许小幅回落，只验恒等式与端点。
    """
    from climbing.chain2d_tissue import (
        tissue_dissipated_energy_t_j, tissue_pad_work_t_j,
    )
    r = _sim()
    w = tissue_pad_work_t_j(r)
    d = tissue_dissipated_energy_t_j(r)
    assert w.shape == (len(r.t_s),) and d.shape == (len(r.t_s),)
    assert w[0] == 0.0 and abs(d[0]) < 1e-20

    # ---- 独立重算（防重构走样）----
    t = r.meta["tissue_params"].thickness_m
    xi = np.maximum(r.node_z, 0.0) / t
    inc = 0.5 * (r.contact_force_n[:-1] + r.contact_force_n[1:]) * np.diff(xi, axis=0) * t
    w_ref = np.concatenate([[0.0], np.cumsum(inc.sum(axis=1))])
    assert np.allclose(w, w_ref, rtol=0, atol=1e-12)

    # ---- 恒等式：d = w − ½E₁u²A·t ≤ w（可恢复项非负）----
    u = np.maximum(xi - r.meta["xi_r"], 0.0)
    rec = 0.5 * E1 * np.sum(u**2 * r.contact_area_m2 * T, axis=1)
    assert np.allclose(d, w - rec, rtol=0, atol=1e-12)
    assert np.all(d <= w + 1e-12)
    assert np.all(rec >= 0.0)

    # ---- 与旧标量函数的关系（历史口径，**不是**总量的收支项）----
    assert tissue_dissipated_energy_j(r) >= 0.0


def test_report_helpers():
    r = _sim()
    d = tissue_peak_deformation_mm(r)
    assert d > 0.0
    sig, t_ms = tibia_peak_stress_mpa(r)
    assert sig > 0.0 and t_ms >= 0.0
    assert tissue_dissipated_energy_j(r) >= 0.0
