"""G5 · 刚地面软组织柔顺层 v2（跟垫 in-vivo 锚定）—— 6 段链镜像模块。

对照对象 = ``chain2d.TissueLayer``（E=2 MPa，注释称「跟垫量级」，实测偏硬 **10.8×**）。
本模块把该柔顺层换成 **in-vivo 锚定**的 v2 本构，并**只影响刚地面路径**
（``on_pad=True`` 仍走 ``chain2d.CrashPad``，垫侧结论不受影响）。

冻结设计：``docs/methods/G5_跟垫柔顺层纠错_预注册.md``。跑后**如实**报，禁止回解凑靶。

本构（线性 → 致密化 + 永久压实，**无平台**）
------------------------------------------------
in-vivo 实测跟垫在 ε≈0.70 / σ≈145 kPa 前**呈线性**，此后急剧致密
（PMC8925218；Gefen 2007 独立佐证 E=175 kPa）⇒ **没有 Gibson-Ashby 平台**
（G4 的「弹性→平台→压实」平台形态与该测量不一致，本件弃用）。

    u  = ξ - ξ_r                     （可恢复弹性应变，ξ=comp/t）
    σ  = E₁·u                         , u ≤ ξ_ref        （线性，in-vivo）
    σ  = E₁·ξ_ref + E_d·(u - ξ_ref)   , u > ξ_ref        （致密化）
    ξ_r 仅在致密化区且加载时增长 ⇒ 卸载沿 E₁ 回永久残余应变（消除整身回弹）

E_d=20 MPa 沿用 G4 冻结值（**沿用值，非新锚点**）。两处过渡用 tanh 带平滑
（数值正则化，非参数调整）。

**粘弹性已从本件移除**（2026-10-13，替换预注册 §3 的 Maxwell 项，理由见
``G5_跟垫柔顺层纠错_预注册.md`` §3 更正块）：冻结形式 ``σ=E₁·ξ_e`` 建模的是
**应力松弛**，与冲击下**率硬化**方向相反；且 ``τ=η/E₁`` 强制 ``E₀≈E₁`` 使该项退化。
现有锚点（静态模量 185 kPa + 黏性模量 43.9 kPa·s）**识别不出冲击速率的松弛谱** ⇒
按「回到文献重查锚点，不调参」纪律删除该项，登记为**待补锚点**。
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
from scipy.integrate import solve_ivp

from .chain2d import (
    SEGMENTS_6,
    PlanarChainConfig,
    PlanarChainResult,
    _kink_forces,
    segment_local_peak,
)
from .pad import CrashPad, Posture, POSTURES
from .pad2d import sign_reg

__all__ = [
    "TissueParams",
    "TissueV2",
    "TissueChainConfig",
    "PLANAR_TISSUE",
    "_tissue_sigma_and_rate",
    "simulate_tissue_chain_fall",
    "tissue_dissipated_energy_j",
    "tissue_pad_work_t_j",
    "tissue_dissipated_energy_t_j",
    "tissue_peak_deformation_mm",
]


# ==========================================================================
# 参数（in-vivo 锚点，跑前冻结）
# ==========================================================================
@dataclass(frozen=True)
class TissueParams:
    """跟垫/软组织柔顺层 v2 参数（**全部 in-vivo 锚定**，非调参）。

    Attributes
    ----------
    E1_pa
        线性段模量。**0.185 MPa**（in-vivo gait 加载中位，PMC8925218；
        Gefen 2007 in-vivo 175 kPa 独立佐证）。对照值 E=2 MPa **偏硬 10.8×**。
    thickness_m
        垫厚。**15.4 mm**（PMC8925218 in-vivo）。
    xi_ref
        线性段终点 / 致密化起点。**0.70**（PMC8925218 in-vivo 峰值应变）。
        ⇒ 压缩在 ``t·ξ_ref = 10.8 mm`` 处进致密化并锁死（正对照，见 §验收 1）。
    E_dense_pa
        致密化段模量。**沿用 G4 冻结值 20 MPa**（细胞支架压实量级；非新锚点）。
    """

    E1_pa: float = 0.185e6
    thickness_m: float = 0.0154
    xi_ref: float = 0.70
    E_dense_pa: float = 20.0e6

    @property
    def xi_lock_m(self) -> float:
        """致密化起始压缩量 ``t·ξ_ref``（m）—— 峰值变形的饱和值。"""
        return self.thickness_m * self.xi_ref

    def flow_pa(self, u: float) -> float:
        """加载骨架（可恢复弹性应变 → 应力），C⁰ 连续、在 ``u=ξ_ref`` 处斜率跃变。"""
        if u <= self.xi_ref:
            return self.E1_pa * u
        return self.E1_pa * self.xi_ref + self.E_dense_pa * (u - self.xi_ref)


@dataclass(frozen=True)
class TissueV2:
    """v2 柔顺层的无状态力接口（诊断 / 单元测试用；历史量需显式传入）。"""

    params: TissueParams = field(default_factory=TissueParams)

    def force_n(
        self,
        compression_m: float,
        area_m2: float,
        xi_r: float = 0.0,
        vz: float = 0.0,
    ) -> float:
        sig, _ = _tissue_sigma_and_rate(compression_m, xi_r, vz, self.params)
        return float(sig * area_m2)


def _tissue_sigma_and_rate(
    compression_m: float,
    xi_r: float,
    vz: float,
    p: TissueParams,
    eps_xi: float = 0.05,
    eps_vel: float = 0.05,
    eps_cap: float = 0.06,
) -> tuple[float, float]:
    """v2 本构：``(σ, dξ_r/dt)``。

    ``compression_m ≤ 0`` ⇒ ``(0, 0)``（未接触，slack）。
    加载沿骨架、卸载沿 ``E₁``（回弹到永久残余应变 ⇒ 不回弹）。

    **门控用总应变 ξ，不是弹性应变 u**（2026-10-13 规则修正 #4）：
    ``ξ_ref=0.70`` 是 **in-vivo 总应变峰值**（PMC8925218 的 peak strain），
    故致密化须由**总应变**越过 ξ_ref 触发。首版误用 ``u = ξ - ξ_r`` 门控，
    于是 ξ_r→0.70 与 u→0.73 可叠加，总应变跑到 **1.43**（变形 22 mm ≈ 目标带 2×），
    且该错误与输入能量无关（0.05–1.0 m 全部饱和在 22 mm）。
    """
    if compression_m <= 0.0:
        return 0.0, 0.0
    xi = compression_m / p.thickness_m
    xi_ref = p.xi_ref
    u = max(xi - xi_r, 0.0)

    sig_lin = p.E1_pa * u
    # 致密化：σ = E₁·(ξ_ref − ξ_r) + E_d·(ξ − ξ_ref)，在 ξ = ξ_ref 处与线性段衔接
    sig_dense = p.E1_pa * max(xi_ref - xi_r, 0.0) + p.E_dense_pa * max(xi - xi_ref, 0.0)
    wb = 0.5 * (1.0 + np.tanh((xi - xi_ref) / eps_xi))       # 1 = 致密化区（总应变门控）
    wv = 0.5 * (1.0 + np.tanh(vz / eps_vel))                  # 1 = 加载
    sig_load = wb * sig_dense + (1.0 - wb) * sig_lin
    sig = wv * sig_load + (1.0 - wv) * sig_lin

    # ξ_r：仅在**线性区**（总应变 ξ < ξ_ref）且加载时增长，速率 ∝ 剩余容量 → 饱和于 ξ_ref
    wlin = 1.0 - wb
    cap = 0.5 * (1.0 + np.tanh((xi_ref * (1.0 - eps_cap) - xi_r)
                              / (eps_cap * xi_ref)))
    dr = wv * wlin * cap * max(vz, 0.0) / p.thickness_m
    return float(sig), float(dr)


# ==========================================================================
# 配置 + 上下文（镜像 chain2d_crush）
# ==========================================================================
@dataclass(frozen=True)
class TissueChainConfig(PlanarChainConfig):
    """平面链配置 + v2 柔顺层参数（镜像 ``CrushChainConfig``，纯增量）。

    **容差覆盖**：v2 在致密化段的接触刚度（``E_d·A/t ≈ 58 MN/m``，ω ≈ 5150 rad/s）
    比 G1/G4 刚性得多，沿用 ``PlanarChainConfig`` 的 ``rtol=1e-8/atol=1e-11``
    会让 Radau 落到浮点间距（"Required step size is less than spacing between
    numbers"），**在 3 m / 10 m 直接失败**。故本模块自带较宽容差默认。
    **已验收敛**：h=3 m 在 ``rtol=1e-5`` 与 ``1e-4`` 下峰值变形同为 **23.22 mm**
    （逐位一致），松容差不降低精度、只是允许更大步长。
    ��值覆盖是**数值设置**，不是物理参数调整。
    """

    tissue: TissueParams = field(default_factory=TissueParams)
    rtol: float = 1e-4
    atol: float = 1e-6


PLANAR_TISSUE: dict | None = None


@dataclass(frozen=True)
class _TissueCtx:
    n_nodes: int
    m: np.ndarray
    rest: np.ndarray
    k_ax: np.ndarray
    k_tors: np.ndarray
    comp_limit: np.ndarray
    ext_limit: np.ndarray
    k_comp_stop: np.ndarray
    k_ext_stop: np.ndarray
    c_ax: np.ndarray
    area_base: np.ndarray
    posture: Posture
    mu: float
    m_total: float
    g: float
    zeta: float
    tissue: TissueParams
    on_pad: bool


def _forces_tissue(x, vx, z, vz, xir, ctx: _TissueCtx):
    """单状态 (n,) → ``(ax, az, joint_force, contact_force, dxir)``。

    关节/抗弯/摩擦部分与 ``chain2d._forces`` 逐行一致；接触部分走 v2 本构
    （``on_pad=True`` 时用 ``CrashPad``、永久压实不激活）。
    """
    n = ctx.n_nodes
    ax = np.zeros(n)
    az = np.zeros(n) + ctx.g
    joint_force = np.zeros(n - 1)

    for j in range(n - 1):
        d0 = x[j + 1] - x[j]
        d1 = z[j + 1] - z[j]
        ell = float(np.hypot(d0, d1))
        if ell < 1e-9:
            ell = 1e-9
        ux, uz = d0 / ell, d1 / ell
        f = ctx.k_ax[j] * (ell - ctx.rest[j])
        comp = ctx.rest[j] - ell
        if comp > ctx.comp_limit[j]:
            f -= ctx.k_comp_stop[j] * (comp - ctx.comp_limit[j])
        ext = ell - ctx.rest[j]
        if ext > ctx.ext_limit[j]:
            f += ctx.k_ext_stop[j] * (ext - ctx.ext_limit[j])
        drel = (vx[j + 1] - vx[j]) * ux + (vz[j + 1] - vz[j]) * uz
        f += ctx.c_ax[j] * drel
        joint_force[j] = f
        ax[j] += f * ux / ctx.m[j]
        az[j] += f * uz / ctx.m[j]
        ax[j + 1] -= f * ux / ctx.m[j + 1]
        az[j + 1] -= f * uz / ctx.m[j + 1]

    for i in range(1, n - 1):
        kt = ctx.k_tors[i - 1]
        if kt == 0.0:
            continue
        a = np.array([x[i - 1], z[i - 1]])
        b = np.array([x[i], z[i]])
        c = np.array([x[i + 1], z[i + 1]])
        fa, fb, fc = _kink_forces(a, b, c, kt)
        ax[i - 1] += fa[0] / ctx.m[i - 1]; az[i - 1] += fa[1] / ctx.m[i - 1]
        ax[i] += fb[0] / ctx.m[i];         az[i] += fb[1] / ctx.m[i]
        ax[i + 1] += fc[0] / ctx.m[i + 1]; az[i + 1] += fc[1] / ctx.m[i + 1]

    contact_force = np.zeros(n)
    dxir = np.zeros(n)
    for i in range(n):
        comp = z[i]
        if comp > 0.0:
            if i == 0:
                kappa = comp / ctx.posture.collapse_ref_m
                area = ctx.posture.area_at(kappa)
            else:
                area = ctx.area_base[i]
            if ctx.on_pad:
                fn = float(CrashPad().force_n(comp, area))
            else:
                sig, dr = _tissue_sigma_and_rate(comp, xir[i], vz[i], ctx.tissue)
                fn = float(sig * area)
                dxir[i] = dr
            if fn > 0.0 and ctx.zeta > 0.0:
                k_sec = fn / comp
                c = 2.0 * ctx.zeta * float(np.sqrt(k_sec * ctx.m_total))
                if vz[i] < 0.0:
                    fn = max(fn + c * vz[i], 0.0)
            contact_force[i] = fn
            az[i] -= fn / ctx.m[i]
            if ctx.mu != 0.0 and fn > 0.0:
                fx = -ctx.mu * fn * float(sign_reg(vx[i], eps=1e-3))
                ax[i] += fx / ctx.m[i]
    return ax, az, joint_force, contact_force, dxir


def _rhs_tissue(t, y, ctx: _TissueCtx):
    """状态导数。布局 ``[x,vx,z,vz,ξ_r]`` 每节点 5 分量（与 crush 同）。"""
    x = y[0::5]; vx = y[1::5]; z = y[2::5]; vz = y[3::5]; xir = y[4::5]
    ax, az, _, _, dxir = _forces_tissue(x, vx, z, vz, xir, ctx)
    out = np.empty(5 * ctx.n_nodes)
    out[0::5] = vx
    out[1::5] = ax
    out[2::5] = vz
    out[3::5] = az
    out[4::5] = dxir
    return out


def _build_tissue_ctx(
    posture: Posture, on_pad: bool, cfg: TissueChainConfig,
    segments: tuple = SEGMENTS_6,
) -> _TissueCtx:
    n = len(segments)
    m = np.array([cfg.mass_kg * s.mass_frac for s in segments])
    m *= cfg.mass_kg / m.sum()
    lengths = np.array([s.length_m for s in segments])
    rest = 0.5 * (lengths[:-1] + lengths[1:])
    k_ax_seg = np.array([s.k_axial_n_per_m for s in segments])
    k_ax = 1.0 / (1.0 / k_ax_seg[:-1] + 1.0 / k_ax_seg[1:])
    k_tors = np.array([s.k_torsion_nm_per_rad for s in segments[:-1]])
    comp_limit = np.array([s.travel_compress_m for s in segments[:-1]])
    ext_limit = np.array([s.travel_extend_m for s in segments[:-1]])
    k_comp_stop = cfg.k_comp_stop_mult * k_ax
    k_ext_stop = cfg.k_ext_stop_mult * k_ax
    m_above = np.array([m[j + 1:].sum() for j in range(n - 1)])
    c_ax = 2.0 * cfg.zeta_joint * np.sqrt(k_ax * m_above)
    if cfg.c_axial_n_s_per_m > 0.0:
        c_ax = c_ax + cfg.c_axial_n_s_per_m
    area_base = np.array([s.contact_area_m2 for s in segments])
    return _TissueCtx(
        n_nodes=n, m=m, rest=rest, k_ax=k_ax, k_tors=k_tors,
        comp_limit=comp_limit, ext_limit=ext_limit,
        k_comp_stop=k_comp_stop, k_ext_stop=k_ext_stop, c_ax=c_ax,
        area_base=area_base, posture=posture, mu=cfg.mu_foot,
        m_total=float(m.sum()), g=cfg.g, zeta=cfg.zeta_contact,
        tissue=cfg.tissue, on_pad=on_pad,
    )


# ==========================================================================
# 模拟入口（镜像 simulate_crush_chain_fall）
# ==========================================================================
def simulate_tissue_chain_fall(
    height_m: float = 3.0,
    mass_kg: float = 80.0,
    posture: str | Posture = "feet-first-stiff",
    on_pad: bool = False,
    config: TissueChainConfig | None = None,
    segments: tuple = SEGMENTS_6,
) -> PlanarChainResult:
    """模拟一次坠落（6 段平面软关节链 + v2 软组织柔顺层）。

    ``segments`` 是 **G6 旁路入口**（``SEGMENTS_6_SOFT``）；默认 ``SEGMENTS_6``
    ⇒ **既有行为逐位不变**。段定义只影响链的刚度/几何/观测量分母，
    接触本构与其余口径完全不变。

    返回 ``PlanarChainResult``（与 G1 同容器），``meta`` 增补：
    ``tissue_active`` / ``tissue_params`` / ``xi_r``（永久压实时程 (T,6)）。
    """
    post = POSTURES[posture] if isinstance(posture, str) else posture
    cfg = config or TissueChainConfig(mass_kg=mass_kg)
    if PLANAR_TISSUE:
        cfg = replace(cfg, **PLANAR_TISSUE)
    cfg = replace(cfg, mass_kg=float(mass_kg))

    ctx = _build_tissue_ctx(post, on_pad, cfg, segments)
    n = ctx.n_nodes

    z0 = np.empty(n)
    z0[0] = -float(height_m)
    for i in range(n - 1):
        z0[i + 1] = z0[i] - ctx.rest[i]
    x0 = np.zeros(n)

    t_contact = float(np.sqrt(2.0 * height_m / cfg.g))
    fall = 0.5 * cfg.g * t_contact * t_contact
    z0 = z0 + fall
    vz0 = np.full(n, cfg.g * t_contact)

    y0 = np.empty(5 * n)
    y0[0::5] = x0
    y0[1::5] = 0.0
    y0[2::5] = z0
    y0[3::5] = vz0
    y0[4::5] = 0.0

    def settled(_t, y):
        zz = y[2::5]; vv = y[3::5]
        return 1.0 if (zz[0] > 0.0 and np.all(np.abs(vv) < 0.05)) else -1.0

    settled.terminal = True
    settled.direction = 1.0

    t_end = cfg.t_max if cfg.t_max is not None else t_contact + cfg.t_window_s

    sol = solve_ivp(
        lambda t, y: _rhs_tissue(t, y, ctx),
        (t_contact, t_end), y0,
        method="Radau",
        events=[settled],
        max_step=cfg.max_ode_step,
        rtol=cfg.rtol, atol=cfg.atol,
        dense_output=True,
    )
    if not sol.success:
        raise RuntimeError(f"链积分失败: {sol.message}")

    t_out = np.arange(t_contact, sol.t[-1] + cfg.dt_sample, cfg.dt_sample)
    t_out = t_out[t_out <= sol.t[-1] + 1e-12]
    if t_out.size < 4 or t_out[-1] < sol.t[-1] - 1e-9:
        t_out = np.append(t_out, sol.t[-1])

    Y = sol.sol(t_out)
    node_x = Y[0::5].T; node_vx = Y[1::5].T
    node_z = Y[2::5].T; node_vz = Y[3::5].T
    xi_r = Y[4::5].T

    dx = node_x[:, 1:] - node_x[:, :-1]
    dz = node_z[:, 1:] - node_z[:, :-1]
    ell = np.maximum(np.hypot(dx, dz), 1e-9)
    joint_comp = ctx.rest[None, :] - ell
    angle = np.zeros((t_out.size, n - 1))
    if n > 2:
        p = np.stack([dx[:, :-1], dz[:, :-1]], axis=-1)
        q = np.stack([dx[:, 1:], dz[:, 1:]], axis=-1)
        crs = p[:, :, 0] * q[:, :, 1] - p[:, :, 1] * q[:, :, 0]
        dot = np.sum(p * q, axis=-1)
        angle[:, 1:] = np.arctan2(crs, dot)

    joint_force = np.zeros((t_out.size, n - 1))
    contact_force = np.zeros((t_out.size, n))
    for k in range(t_out.size):
        _, _, jf, cf, _ = _forces_tissue(
            node_x[k], node_vx[k], node_z[k], node_vz[k], xi_r[k], ctx
        )
        joint_force[k] = jf
        contact_force[k] = cf

    contact_area = np.zeros((t_out.size, n))
    c0 = np.maximum(node_z[:, 0], 0.0)
    contact_area[:, 0] = np.where(
        c0 > 0, np.array([post.area_at(k / post.collapse_ref_m) for k in c0]), 0.0
    )
    for i in range(1, n):
        contact_area[:, i] = np.where(node_z[:, i] > 0.0, ctx.area_base[i], 0.0)

    f_total = contact_force.sum(axis=1)
    any_c = f_total > 0.0
    if any_c.any():
        i0 = int(np.argmax(any_c))
        gap = np.where(~any_c[i0:])[0]
        i1 = (i0 + int(gap[0])) if gap.size else t_out.size
    else:
        if height_m > 0.0:
            raise RuntimeError(f"积分结束时仍未发生接触：h={height_m:.1f} m。")
        i0, i1 = 0, 0
    i1 = min(max(i1, i0 + 2), t_out.size)
    sl = slice(i0, i1)

    peak_force = float(np.max(f_total[sl])) if f_total[sl].size else 0.0
    max_leg = float(np.max(np.maximum(node_z[sl], 0.0))) if i1 > i0 else 0.0
    a_head = np.gradient(node_vz[:, n - 1], t_out) / cfg.g
    a_trunk = np.gradient(node_vz[:, 4 if n > 4 else n - 1], t_out) / cfg.g
    first = int(np.argmax(np.any(contact_force[sl] > 0, axis=0))) if i1 > i0 else 0

    return PlanarChainResult(
        posture=post.name,
        pad_label=("Crash pad" if on_pad else "硬地面（柔顺层 v2）"),
        height_m=float(height_m), mass_kg=float(mass_kg),
        t_s=t_out,
        node_x=node_x, node_z=node_z, node_vx=node_vx, node_vz=node_vz,
        joint_force_n=joint_force, joint_compression_m=joint_comp,
        joint_angle_rad=angle,
        contact_force_n=contact_force, contact_area_m2=contact_area,
        impact_speed_ms=float(np.sqrt(2.0 * cfg.g * height_m)),
        impact_energy_j=float(mass_kg) * cfg.g * height_m,
        peak_force_n=peak_force,
        max_compression_m=max_leg,
        peak_accel_head_g=float(np.max(np.abs(a_head[sl]))),
        peak_accel_trunk_g=float(np.max(np.abs(a_trunk[sl]))),
        bottomed_out=False,
        first_contact_node=first,
        meta={
            "posture_obj": post,
            "on_pad": on_pad,
            "segments": segments,
            "g": cfg.g,
            "node_mass_kg": ctx.m,
            "impact_window_s": (float(t_out[i0]), float(t_out[i1 - 1])),
            "tissue_active": (not on_pad),
            "tissue_params": ctx.tissue if not on_pad else None,
            "xi_r": xi_r,
        },
    )


# ==========================================================================
# 报告用工具
# ==========================================================================
def tissue_dissipated_energy_j(result: PlanarChainResult) -> float:
    """永久压实路径功 (J)：``Σ F·Δ(ξ_r·t)``（滞回面积）。

    梯形平均力（压实常发生在单采样步内，步首力会系统性低估）。

    ⚠️ **它不是能量账目的收支项** —— 本函数返回**标量**（总量），把它加进
    逐时刻的 ``e(t)`` 时它在**每个时刻都是同一个常数**、自动抵消，
    对能量闭合**零贡献**。能量账目必须用**时间序列** ``tissue_pad_work_t_j``
    （见其 docstring 与 `tests/test_chain2d_tissue.py` 的闭合测试）。
    """
    xi_r = result.meta.get("xi_r")
    if xi_r is None or not result.meta.get("tissue_active"):
        return 0.0
    t = result.meta["tissue_params"].thickness_m
    d_comp_r = np.diff(xi_r, axis=0) * t
    F = result.contact_force_n
    e = 0.0
    for k in range(d_comp_r.shape[0]):
        e += float(np.sum(0.5 * (F[k] + F[k + 1]) * np.maximum(d_comp_r[k], 0.0)))
    return e


def tissue_pad_work_t_j(result: PlanarChainResult) -> np.ndarray:
    """垫层**路径功**的时间序列 (J) —— 能量账目**唯一合法的收支项**。

    ``W(t) = Σ_{k<t} Σ_n ½(F_k,n + F_{k+1},n)·Δ(ξ_k,n·t)``（**有符号**：卸载时回吐能量）。

    **为什么必须用序列**（2026-10-15 定案）：``tissue_dissipated_energy_j`` 返回标量，
    加进 ``e(t)`` 后各时刻同增同减 ⇒ 恒被抵消。旧闭合测试因此**实际只测了
    ``KE+PE`` 的漂移**，漂移值 = 未计入的耗散占比（h=1.5 m 时 2.06%）。

    **守恒量**（保守配置下）：``KE + PE_g + PE_j + W(t)``。
    实测闭合（2026-10-15，rtol=1e-8）：

    =========  ============  ==================
    h (m)      旧测试量误差   用本函数后误差
    =========  ============  ==================
    0.5        1.1331%       **0.1709%**
    1.5        2.0642%       **0.2140%**
    2.0        2.4336%       **0.5062%**
    =========  ============  ==================

    残余 0.17–0.51% 为**采样网格上的梯形积分误差**（非物理漏项；
    在 h=1.5 m 其偏差单向为正、单调累积，符合积分误差特征）。

    ⚠️ **被证伪的旧假设（如实留档）**：原判为「致密化段（ξ>ξ_ref）耗散未进账目」
    ⇒ 补了致密化滞回项后**漂移一字不变**（2.060% → 2.059%），
    因该项同样以标量进入 ``e(t)``。**真因是标量/序列之别，不是漏项。**
    """
    xi_r = result.meta.get("xi_r")
    if xi_r is None or not result.meta.get("tissue_active"):
        return np.zeros(len(result.t_s))
    t = result.meta["tissue_params"].thickness_m
    xi = np.maximum(result.node_z, 0.0) / t
    w_inc = 0.5 * (result.contact_force_n[:-1]
                   + result.contact_force_n[1:]) * np.diff(xi, axis=0) * t
    return np.concatenate([[0.0], np.cumsum(w_inc.sum(axis=1))])


def tissue_dissipated_energy_t_j(result: PlanarChainResult) -> np.ndarray:
    """耗散能量的**时间序列** (J) = 路径功 − 可恢复弹性能。

    ``= tissue_pad_work_t_j(result) − ½E₁·u²·A·t``（``u = max(ξ − ξ_r, 0)``）。

    恒等式：``KE + PE_g + PE_j + PE_c + 本函数 = KE + PE_g + PE_j + W(t)``
    （``PE_c`` 即 ``½E₁u²At``，两项相消）⇒ 用它替换 ``PE_c`` 之外的标量耗散，
    闭合质量与 ``tissue_pad_work_t_j`` 相同。
    """
    p = result.meta["tissue_params"]
    t = p.thickness_m
    xi = np.maximum(result.node_z, 0.0) / t
    u = np.maximum(xi - result.meta["xi_r"], 0.0)
    rec = 0.5 * p.E1_pa * np.sum(u**2 * result.contact_area_m2 * t, axis=1)
    return tissue_pad_work_t_j(result) - rec


def tissue_peak_deformation_mm(result: PlanarChainResult) -> float:
    """足节点峰值压缩 (mm) —— **主判据 1** 的量（Bennett&Ker in-vivo 8.5–11.3 mm）。"""
    return float(np.max(np.maximum(result.node_z[:, 0], 0.0)) * 1e3)


def tibia_peak_stress_mpa(result: PlanarChainResult) -> tuple[float, float]:
    """胫骨段峰值轴向应力 (MPa) 与其时刻（相对冲击窗，ms）—— **判据 2**。

    皮质骨抗压强度 130–180 MPa；> 180 MPa 须显式标注超骨强度。
    """
    from .chain2d import SEGMENTS_6 as _S
    s = _S[1]
    load = np.abs(result.segment_axial_load_n(1))
    sig = load / s.tissue_area_m2
    t0, t1 = result.meta["impact_window_s"]
    m = (result.t_s >= t0) & (result.t_s <= t1)
    tt = (result.t_s - t0)[m]
    ss = sig[m]
    k = int(np.argmax(ss))
    return float(ss[k] / 1e6), float(tt[k] * 1e3)


def tibia_local_peak(result: PlanarChainResult) -> tuple[float, float]:
    """胫骨段首显著局部峰（与 G1/G4 同口径，**仅记录不判定**）。"""
    return segment_local_peak(result, 1)
