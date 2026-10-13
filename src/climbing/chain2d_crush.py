"""G4 · 压缩相渐进压溃本构（G1 6 段平面链的接触层升级）。

为什么需要本模块
----------------
``chain2d.py`` 的接触 = ``_SeriesContact``（``hard_surface(k_n=1e7)`` 串联
``TissueLayer(E=2MPa, t=18mm)``），**纯弹性、保守** ⇒ 冲击窗口内整身弹性
回弹仍在（G1 报告 §5），且 [T22] 4.6 ms 的「接触压溃主导」机制未物化进
材料。本模块把接触层升级为**压缩相渐进压溃本构**（Gibson-Ashby 式
弹性 → 平台 → 压实 + **滞回塑性卸载**），消除整身回弹，并把 V8b 胫骨局部
峰时从 5.8 ms 压进 ≤ 5.0 ms（预注册目标；跑后如实 PASS/FAIL，禁止调参凑靶）。

设计（冻结，见 ``docs/methods/G4_压缩相本构设计.md``）
-----------------------------------------------------
* **本构**：单轴压溃律 ``σ(ξ)``（``ξ = comp/t``）：
  - 加载：弹性 ``σ=E·ξ``（ξ≤ξ_y）→ 平台 ``σ=σ_pl``（ξ_y<ξ≤ξ_d）→
    压实 ``σ=σ_pl+E_d·(ξ−ξ_d)``（ξ>ξ_d）；
  - 卸载/重新加载：沿弹性模量 ``E`` 回到残余应变 ``ξ_r``（永久压溃不回弹）
    ⇒ 滞回面积 = 压溃吸能。
* **叠加**：接触 = **单一压溃本构** ``F=A·σ_crush(comp)``（**整体替换**
  ``_SeriesContact``；``hard_surface``/``TissueLayer`` 并入弹性段）。
* **状态**：ODE 状态 ``[x,vx,z,vz,ξ_r]`` 每节点 5 分量（``ξ_r`` = 压溃塑性
  应变，进 ODE 由 ``_rhs_crush`` 积分、``_forces_crush`` 消费 ⇒ 前后处理
  逐位一致，同 G1 纪律）。``chain2d.py`` / ``pad.py`` 一律不动（字节锁）。
* **参数**：全部跑前冻结（§3 来源，非调参）：E=2.0e6 Pa、t=0.018 m、
  σ_pl=150e3 Pa、ξ_d=0.50、E_d=20e6 Pa；派生 ``ξ_y=σ_pl/E=0.075``。

⚠️ 诚实边界：``ξ_r`` 的演化按冻结公式 ``ξ_r = ξ_peak − σ_pl/E``（平台/压实
统一口径）；σ_pl 下限（``max(ξ_y, ...)``）只在 ξ_peak≈ξ_y 处有差别、对
平台期物理无影响，实现取连续形式（首次屈服 ξ_r 从 0 起）。压实段卸载从
平台应力处回弹（σ 在峰值处向 σ_pl 突降 = 设计意图的耗散 snap-back）。
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
from scipy.integrate import solve_ivp

from . import G
from .chain2d import (
    SEGMENTS_6,
    PlanarChainConfig,
    PlanarChainResult,
    _kink_forces,
)
from .pad import CrashPad, Posture, POSTURES
from .pad2d import sign_reg

__all__ = [
    "CrushParams",
    "CrushTissue",
    "CrushChainConfig",
    "PLANAR_CRUSH",
    "simulate_crush_chain_fall",
    "crush_dissipated_energy_j",
    "rebound_metrics",
]


# ==========================================================================
# 压溃本构（冻结参数，§3）
# ==========================================================================
@dataclass(frozen=True)
class CrushParams:
    """压缩相压溃本构参数（跑前冻结，非调参）。

    Attributes
    ----------
    E_pa
        弹性模量（= G1 冻结 ``TissueLayer.modulus_pa``，跟垫量级）。
    thickness_m
        软组织层厚度（= G1 冻结 ``TissueLayer.thickness_m``）。
    sigma_plateau_pa
        平台应力（跟垫 in-vivo 峰值压应力中位 ~145-156 kPa）。
    xi_dense
        压实应变（组织压实起点）。
    E_dense_pa
        压实段模量（≈10×E，组织压实段量级）。
    """

    E_pa: float = 2.0e6
    thickness_m: float = 0.018
    sigma_plateau_pa: float = 150e3
    xi_dense: float = 0.50
    E_dense_pa: float = 20e6

    @property
    def xi_yield(self) -> float:
        """屈服应变（派生）：ξ_y = σ_pl / E。"""
        return self.sigma_plateau_pa / self.E_pa

    def flow_pa(self, xi: float) -> float:
        """加载骨架曲线（弹性 → 平台 → 压实），C⁰ 连续。"""
        if xi <= self.xi_yield:
            return self.E_pa * xi
        if xi <= self.xi_dense:
            return self.sigma_plateau_pa
        return self.sigma_plateau_pa + self.E_dense_pa * (xi - self.xi_dense)


def _crush_sigma_and_rate(
    compression_m: float, xi_r: float, vz: float, p: CrushParams,
    eps_xi: float = 0.003, eps_vel: float = 0.05,
) -> tuple[float, float]:
    """压溃应力 + 塑性应变率 ``(σ, dξ_r/dt)``。

    分段（§2.1 冻结公式，含两处数值正则化平滑）：
    * **弹性段**（``ξ < ξ_r+ξ_y``）：``σ=E·(ξ−ξ_r)``（低于 ξ_r 夹 0 =
      接触 slack），``dξ_r/dt=0``。
    * **流动区**（``ξ > ξ_r+ξ_y``）：加载 ``σ=flow(ξ)``、``dξ_r/dt=vz/t``
      （ξ_r 天然满足冻结式 ``ξ_r = ξ_peak − σ_pl/E``）；卸载沿弹性模量 E
      回残应变（压实峰顶 snap-back 耗散）。

    数值正则化（**非本构参数**，仅防 Radau 硬不连续；量级 ≪ 冲击量级）：
    1. ``eps_xi``（应变带 0.003，≈0.05 mm 穿透）平滑 弹性段/流动区 边界
       （消除重载越过屈服点时的 σ 突跳）；
    2. ``eps_vel``（速度带 0.05 m/s ≪ 冲击速度 7–20 m/s）平滑 加载/卸载
       切换（消除压实 snap-back 突跳）。
    平台期两分支在峰值处本就连续（均 = σ_pl），滞回面积（吸能）近似不变；
    压实 snap-back 的过渡带被软化（耗散略减，报告如实说明）。
    """
    if compression_m <= 0.0:
        return 0.0, 0.0
    xi = compression_m / p.thickness_m
    boundary = xi_r + p.xi_yield
    wb = 0.5 * (1.0 + np.tanh((xi - boundary) / eps_xi))     # 1=流动区
    wv = 0.5 * (1.0 + np.tanh(vz / eps_vel))                 # 1=加载
    sig_el = p.E_pa * max(xi - xi_r, 0.0)
    sig_flow = p.flow_pa(xi)
    sig = wv * (wb * sig_flow + (1.0 - wb) * sig_el) + (1.0 - wv) * sig_el
    dr = wv * wb * max(vz / p.thickness_m, 0.0)
    return float(sig), float(dr)


@dataclass(frozen=True)
class CrushTissue:
    """压缩相压溃本构（§2.1；整体替换接触层，§2.2/§6 决策 3）。

    接口对齐 ``pad.force_n``，但压溃需要塑性状态 ``ξ_r`` 与速度方向
    ``vz``；纯力接口 ``force_n(comp, area)`` 在无状态下退化为**弹性骨架
    **（E·ξ，无滞回），仅供诊断/兼容。
    """

    params: CrushParams = field(default_factory=CrushParams)

    def force_n(
        self,
        compression_m: float,
        area_m2: float,
        xi_r: float = 0.0,
        vz: float = 0.0,
    ) -> float:
        sig, _ = _crush_sigma_and_rate(compression_m, xi_r, vz, self.params)
        return float(sig * area_m2)

    def sigma_and_rate(
        self, compression_m: float, xi_r: float, vz: float
    ) -> tuple[float, float]:
        return _crush_sigma_and_rate(compression_m, xi_r, vz, self.params)


# ==========================================================================
# 配置（冻结） + 模块门
# ==========================================================================
@dataclass(frozen=True)
class CrushChainConfig(PlanarChainConfig):
    """压溃链配置 = G1 平面链配置 + 压溃本构参数（镜像，纯增量）。

    新增字段：``crush`` = 压溃本构参数（冻结，§3）。刚地面接触一律走
    压溃本构（整体替换 ``_SeriesContact``）；G1 弹性对照用既有
    ``chain2d.py``（本模块不提供回退路径，避免死分支）。
    """

    crush: CrushParams = field(default_factory=CrushParams)


#: 模块门（镜像 ``chain2d.PLANAR_CHAIN``）。默认 ``None`` = 不启用任何覆盖
#: ⇒ 链走 ``simulate_crush_chain_fall`` 的显式 ``config`` / 默认配置。
PLANAR_CRUSH: dict | None = None


# ==========================================================================
# 求解上下文（= _Ctx 镜像 + 压溃）
# ==========================================================================
@dataclass(frozen=True)
class _CrushCtx:
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
    crush: CrushTissue
    on_pad: bool


def _forces_crush(x, vx, z, vz, xir, ctx: _CrushCtx):
    """单状态 (n,) 输入 → ``(ax, az, joint_force, contact_force, dxir)``。

    关节/抗弯/摩擦部分与 ``chain2d._forces`` 逐行一致（字节锁镜像）；
    接触部分用压溃本构（``on_pad=True`` 时用 ``CrashPad``、塑性不激活）。
    ``_rhs_crush`` 与后处理共用本函数 ⇒ 前后处理逐位一致。
    """
    n = ctx.n_nodes
    ax = np.zeros(n)
    az = np.zeros(n) + ctx.g
    joint_force = np.zeros(n - 1)

    # ---- 关节：轴向弹簧 + 行程上下限 + 轴向阻尼（= chain2d._forces） -----
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

    # ---- 抗弯（转角）弹簧（内节点）（= chain2d._forces） -----------------
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

    # ---- 逐节点接触：压溃本构（+ 滞回阻尼 + 摩擦） -----------------------
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
                fn = float(ctx.crush.force_n(comp, area))
            else:
                sig, dr = _crush_sigma_and_rate(
                    comp, xir[i], vz[i], ctx.crush.params
                )
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


def _rhs_crush(t, y, ctx: _CrushCtx):
    """状态导数。状态交错布局 ``[x,vx,z,vz,ξ_r]`` 每节点 5 分量。"""
    x = y[0::5]
    vx = y[1::5]
    z = y[2::5]
    vz = y[3::5]
    xir = y[4::5]
    ax, az, _, _, dxir = _forces_crush(x, vx, z, vz, xir, ctx)
    out = np.empty(5 * ctx.n_nodes)
    out[0::5] = vx
    out[1::5] = ax
    out[2::5] = vz
    out[3::5] = az
    out[4::5] = dxir
    return out


def _build_crush_ctx(
    posture: Posture, on_pad: bool, cfg: CrushChainConfig
) -> _CrushCtx:
    """镜像 ``chain2d._build_ctx``；``on_pad`` 用 ``CrashPad`` 否则压溃本构。"""
    n = len(SEGMENTS_6)
    m = np.array([cfg.mass_kg * s.mass_frac for s in SEGMENTS_6])
    m *= cfg.mass_kg / m.sum()
    lengths = np.array([s.length_m for s in SEGMENTS_6])
    rest = 0.5 * (lengths[:-1] + lengths[1:])
    k_ax_seg = np.array([s.k_axial_n_per_m for s in SEGMENTS_6])
    k_ax = 1.0 / (1.0 / k_ax_seg[:-1] + 1.0 / k_ax_seg[1:])
    k_tors = np.array([s.k_torsion_nm_per_rad for s in SEGMENTS_6[:-1]])
    comp_limit = np.array([s.travel_compress_m for s in SEGMENTS_6[:-1]])
    ext_limit = np.array([s.travel_extend_m for s in SEGMENTS_6[:-1]])
    k_comp_stop = cfg.k_comp_stop_mult * k_ax
    k_ext_stop = cfg.k_ext_stop_mult * k_ax
    m_above = np.array([m[j + 1:].sum() for j in range(n - 1)])
    c_ax = 2.0 * cfg.zeta_joint * np.sqrt(k_ax * m_above)
    if cfg.c_axial_n_s_per_m > 0.0:
        c_ax = c_ax + cfg.c_axial_n_s_per_m
    area_base = np.array([s.contact_area_m2 for s in SEGMENTS_6])

    if on_pad:
        pad: object = CrashPad()
        ctx_pad: object = pad
    else:
        ctx_pad = CrushTissue(params=cfg.crush)
    return _CrushCtx(
        n_nodes=n, m=m, rest=rest, k_ax=k_ax, k_tors=k_tors,
        comp_limit=comp_limit, ext_limit=ext_limit,
        k_comp_stop=k_comp_stop, k_ext_stop=k_ext_stop, c_ax=c_ax,
        area_base=area_base, posture=posture, mu=cfg.mu_foot,
        m_total=float(m.sum()), g=cfg.g, zeta=cfg.zeta_contact,
        crush=ctx_pad, on_pad=on_pad,
    )


# ==========================================================================
# 模拟入口（镜像 simulate_planar_chain_fall，状态 5·n）
# ==========================================================================
def simulate_crush_chain_fall(
    height_m: float = 3.0,
    mass_kg: float = 80.0,
    posture: str | Posture = "feet-first-stiff",
    on_pad: bool = False,
    config: CrushChainConfig | None = None,
) -> PlanarChainResult:
    """模拟一次抱石坠落（6 段平面软关节链 + 压缩相压溃接触）。

    Parameters
    ----------
    height_m
        释放点距地面高度 (m)。足底接触时下落量 = ``height_m``。
    mass_kg
        总质量 (kg)。
    posture
        ``pad.POSTURES`` 的键或 ``Posture``。
    on_pad
        True = 软垫（``CrashPad``，压溃塑性不激活）；False = 刚性地面
        （**压溃本构**，整体替换 ``hard_surface``+``TissueLayer``）。
    config
        ``CrushChainConfig``；None 用默认（``PLANAR_CRUSH`` 门若置非空则覆盖）。

    返回 ``PlanarChainResult``（与 G1 同容器，下游口径不变），``meta`` 增补：
    ``crush_active`` / ``crush_params`` / ``xi_r``（压溃塑性应变时程 (T,6)）。
    """
    post = POSTURES[posture] if isinstance(posture, str) else posture
    cfg = config or CrushChainConfig(mass_kg=mass_kg)
    if PLANAR_CRUSH:
        cfg = replace(cfg, **PLANAR_CRUSH)
    cfg = replace(cfg, mass_kg=float(mass_kg))

    ctx = _build_crush_ctx(post, on_pad, cfg)
    n = ctx.n_nodes

    # ---- 初始几何：竖直链，足底节点在 z = -height ----------------------
    z0 = np.empty(n)
    z0[0] = -float(height_m)
    for i in range(n - 1):
        z0[i + 1] = z0[i] - ctx.rest[i]
    x0 = np.zeros(n)

    # ---- 弹道起点：自由落体解析到**首次接触**那一刻，再启刚体积分 --------
    t_contact = float(np.sqrt(2.0 * height_m / cfg.g))
    fall = 0.5 * cfg.g * t_contact * t_contact
    z0 = z0 + fall
    vz0 = np.full(n, cfg.g * t_contact)

    y0 = np.empty(5 * n)
    y0[0::5] = x0
    y0[1::5] = 0.0
    y0[2::5] = z0
    y0[3::5] = vz0
    y0[4::5] = 0.0                # ξ_r 初始 0（未压溃）

    # ---- 终止：全节点已在接触（节点 0 已穿透）且近静止 ------------------
    def settled(_t, y):
        zz = y[2::5]
        vv = y[3::5]
        touched = zz[0] > 0.0
        return 1.0 if (touched and np.all(np.abs(vv) < 0.05)) else -1.0

    settled.terminal = True
    settled.direction = 1.0

    if cfg.t_max is not None:
        t_end = cfg.t_max
    else:
        t_end = t_contact + cfg.t_window_s

    sol = solve_ivp(
        lambda t, y: _rhs_crush(t, y, ctx),
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
    node_x = Y[0::5].T
    node_vx = Y[1::5].T
    node_z = Y[2::5].T
    node_vz = Y[3::5].T
    xi_r = Y[4::5].T

    # ---- 关节/接触量：用**与 _rhs_crush 同一个 _forces_crush** 重算 -------
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
        _, _, jf, cf, _ = _forces_crush(
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

    # ---- 冲击窗口：首次接触 → 首次分离 ---------------------------------
    any_c = f_total > 0.0
    if any_c.any():
        i0 = int(np.argmax(any_c))
        gap = np.where(~any_c[i0:])[0]
        i1 = (i0 + int(gap[0])) if gap.size else t_out.size
    else:
        if height_m > 0.0:
            raise RuntimeError(
                f"积分结束时仍未发生接触：h={height_m:.1f} m（t={sol.t[-1]:.3f}s）。"
            )
        i0, i1 = 0, 0
    i1 = min(max(i1, i0 + 2), t_out.size)
    sl = slice(i0, i1)

    peak_force = float(np.max(f_total[sl])) if f_total[sl].size else 0.0
    max_leg = float(np.max(np.maximum(node_z[sl], 0.0))) if i1 > i0 else 0.0
    a_head = np.gradient(node_vz[:, n - 1], t_out) / cfg.g
    a_trunk = np.gradient(node_vz[:, 4 if n > 4 else n - 1], t_out) / cfg.g
    pad_max_comp = (
        ctx.crush.max_compression_m if hasattr(ctx.crush, "max_compression_m")
        else 0.20
    )
    bottomed = bool(max_leg > pad_max_comp * 0.999)

    first = int(np.argmax(np.any(contact_force[sl] > 0, axis=0))) if i1 > i0 else 0

    return PlanarChainResult(
        posture=post.name,
        pad_label=("Crash pad" if on_pad else "硬地面（压溃本构）"),
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
        bottomed_out=bottomed,
        first_contact_node=first,
        meta={
            "posture_obj": post,
            "on_pad": on_pad,
            "segments": SEGMENTS_6,
            "g": cfg.g,
            "node_mass_kg": ctx.m,
            "impact_window_s": (float(t_out[i0]), float(t_out[i1 - 1])),
            "crush_active": (not on_pad),
            "crush_params": ctx.crush.params if not on_pad else None,
            "xi_r": xi_r,
        },
    )


# ==========================================================================
# 能量账 / 回弹量化（报告与测试共用）
# ==========================================================================
def crush_dissipated_energy_j(result: PlanarChainResult) -> float:
    """压溃滞回耗能 (J)：``Σ F·Δ(ξ_r·t)``（塑性路径功 = 滞回面积）。

    取 ``meta["xi_r"]`` 与接触面积逐样本累加；单位：F(N)·Δcomp_r(m)。
    仅对 ``crush_active`` 工况有意义（on_pad 时 ξ_r 恒 0 ⇒ 返回 0）。
    """
    xi_r = result.meta.get("xi_r")
    if xi_r is None:
        return 0.0
    t = result.t_s
    A = result.contact_area_m2
    F = result.contact_force_n
    comp_r = xi_r * (
        result.meta["crush_params"].thickness_m
        if result.meta.get("crush_params") is not None else 0.018
    )
    d = np.diff(comp_r, axis=0)
    return float(np.sum(F[:-1] * np.maximum(d, 0.0)))


def rebound_metrics(result: PlanarChainResult) -> dict:
    """冲击窗口内整身回弹量化（G1 §5 诚实短板对比口径）。

    * ``head_rise_m``：头节点相对**首次接触时刻高度**的最大**上升**量
      （z 向下为正，上升 = z 减小）；>0 即存在回弹。
    * ``foot_lift_mm``：足节点最大离地（穿透变负）量 (mm)。
    """
    t0, t1 = result.meta["impact_window_s"]
    m = (result.t_s >= t0) & (result.t_s <= t1)
    z_head = result.node_z[m, -1]
    z_foot = result.node_z[m, 0]
    z0_head = z_head[0]
    head_rise = max(0.0, float(z0_head - z_head.min()))
    foot_lift = max(0.0, float(-z_foot.min())) * 1e3
    return {"head_rise_m": head_rise, "foot_lift_mm": foot_lift}
