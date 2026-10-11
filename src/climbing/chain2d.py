"""G1 · 6 段平面软关节多刚体落地链（P3「二切」）。

为什么需要这个模块
-----------------
1D 双质点链（``pad.py``）的峰时下限被人体等效刚度 ``k_n=1e7`` 钉死在
~14-21 ms；而 [T22] Table 3 报的是**局部组织（胫骨）应变峰时**（20 m
时 4.6 ms）。二者**不是同一个观测量**：1D 的 ``t_torso`` 是链力/整体减速
峰，文献是局部组织应变峰。本模块把人体换成 **6 段平面软关节链**，并
报一个**逐段局部组织峰**（``ε_j = (F_j/A_j)/E_j``），对齐 [T22] 口径。

设计（冻结，见 ``docs/methods/G1_多段链设计.md``）
------------------------------------------------
* **拓扑**：6 段（足 / 小腿-胫骨 / 大腿 / 骨盆-躯干 / 胸廓 / 头颈）/ 5 关节
  （踝 / 膝 / 髋 / 腰 / 颈）；平面软关节串链。
* **自由度**：基座侧向 ``x_f`` + 每关节（轴向压缩 ``u_j``、平面转角 ``φ_j``）。
  本实现用**平面点质量 + 软关节弹簧**（= bodies.py 的结构模式的平面推广）：
  6 个质心节点各带 ``(x_i, z_i)``；相邻节点由带**轴向弹簧（+行程上下限）**
  的「软关节」相连；内节点再加**抗弯（转角）弹簧**。⇒ 12 个平面位置自由度
  = 冻结设计的「11」（基座侧向 + 5×(轴向,转角)）+ 基座**竖直**接触坐标
  （接触柔顺需要一个竖直自由度，冻结文未单列；见报告 §偏差记录）。
* **接触**：逐节点柔顺接触，**逐段穿透激活**（哪一段沉到地面以下哪一段
  受力）+ 面积增长 ``A(κ)``（复用 ``Posture.contact_schedule``）；摩擦/符号
  用 ``sign_reg`` ⇒ C¹、确定性、**种子无关**。刚地面 = 既有
  ``hard_surface(k_n=1e7)`` **串联**人体软组织柔顺层；垫 = ``CrashPad``。
* **观测量**：``segment_local_peak(result, seg)`` 取 ``ε_段`` 在
  ``meta["impact_window_s"]`` 内的**首个显著局部极大**（复用
  ``head_branch._first_significant_peak``）。
* **积分**：scipy **Radau**（``bodies.py`` 已验证、隐式、抗刚）。**不碰**
  RK45 默认路径（字节锁）。

参数来源（跑前冻结，**非调参**）
------------------------------
* 段质量比：de Leva (1996) 男性调整 Zatsiorsky-Seluyanov 数据
  （head 6.94%、trunk 43.46%、thigh 14.16%、shank 4.33%、foot 1.37%、
  上臂 2.71%、前臂 1.62%、手 0.61%）。双腿对称落地 ⇒ 腿部动作合并到
  单一矢状链（见 ``SEGMENTS_6`` 注释）。
* 段长：Drillis & Contini (1966) 节段长比（Winter 2009 转引）× H=1.73 m。
* 关节轴向刚度：远端骨主导（``EA/L``）、近端肌/脊柱主导（量级估计）。
* 组织模量/截面（观测量标度用；峰**时**与之无关）：皮质骨 E≈17 GPa、
  跟垫/软组织 E≈1-2 MPa 等。
* 接触层：``hard_surface`` 的 ``k_n=1e7``（``pad.py`` 文档化）；跟垫厚度/
  模量取文献量级。

⚠️ 诚实边界：本模块是**未标定**的结构性扩展件。观测量与参数跑前冻结，
跑后**如实**报 PASS/FAIL（见 ``results/G1_CHAIN2D_REPORT.md``），不调参凑靶。
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
from scipy.integrate import solve_ivp

from . import G
from .head_branch import _first_significant_peak
from .pad import CrashPad, Posture, POSTURES, hard_surface
from .pad2d import sign_reg

__all__ = [
    "Segment",
    "SEGMENTS_6",
    "PlanarChainConfig",
    "PLANAR_CHAIN",
    "PlanarChainResult",
    "simulate_planar_chain_fall",
    "segment_local_peak",
]


# ==========================================================================
# 段定义（人体测量学，跑前冻结）
# ==========================================================================
@dataclass(frozen=True)
class Segment:
    """链上一段（刚体）。节点 = 该段的质心。

    Attributes
    ----------
    name
        段名（足 / 小腿-胫骨 / 大腿 / 骨盆-躯干 / 胸廓 / 头颈）。
    length_m
        段长（m）。相邻节点（质心）间距取 ``(L_i + L_{i+1})/2``（见求解）。
    mass_frac
        占总质量比（de Leva 1996，双腿合并口径）。
    k_axial_n_per_m
        段自身**轴向**刚度 (N/m)。关节 j 的轴向刚度 = 相邻两段的
        串联和 ``k_j·k_{j+1}/(k_j+k_{j+1})``（两段串联）。
    k_torsion_nm_per_rad
        段**上端**关节的抗弯（转角）刚度 (Nm/rad)。关节 j（连接节点 j
        与 j+1）取 ``SEGMENTS_6[j].k_torsion_nm_per_rad``。
    travel_compress_m
        该段关节的**最大压缩**行程 (m)（超出后加硬止挡）。
    travel_extend_m
        该段关节的**最大拉伸**行程 (m)（超出后加**软**止挡 —— 明确不用
        bodies.py 的过硬抗拉止挡，那是弹跳根因）。
    contact_area_m2
        该段作为**接触节点**时的解剖接触面积 (m²)；节点 0（足）另按
        ``Posture.contact_schedule`` 随塌陷增长。
    tissue_area_m2
        观测量 ``ε=F/(A·E)`` 用的承载截面 (m²)。
    tissue_modulus_pa
        观测量用的组织模量 (Pa)。
    """

    name: str
    length_m: float
    mass_frac: float
    k_axial_n_per_m: float
    k_torsion_nm_per_rad: float
    travel_compress_m: float
    travel_extend_m: float
    contact_area_m2: float
    tissue_area_m2: float
    tissue_modulus_pa: float


#: 6 段链（H = 1.73 m 成年男性；双腿对称落地的矢状合并口径）。
#:
#: 质量比来源 de Leva (1996) 男性：head 6.94%、trunk 43.46%、thigh 14.16%、
#: shank 4.33%、foot 1.37%、上臂 2.71%、前臂 1.62%、手 0.61%。**双腿合并**
#: （双脚同时落地，单矢状链代表全身）：脚 2×1.37、小腿 2×4.33、大腿 2×14.16。
#: 躯干 43.46% 按 ExRx/Zatsiorsky 子段比（胸廓 20.1 / 腹 13.06 / 骨盆 13.66）
#: 拆成骨盆-躯干（腹+骨盆）与胸廓；两臂（3×2）并入胸廓（肩带）。
#: 合计 = 100.0%。段长来源 Drillis & Contini (1966)（Winter 2009 转引），
#: 关节高度比 H：踝 0.039 / 膝 0.285 / 髋 0.530 / 腰 0.650 / 肩 0.818 / 顶 1.0。
SEGMENTS_6: tuple[Segment, ...] = (
    Segment(
        name="足", length_m=0.067, mass_frac=0.0274,
        k_axial_n_per_m=2.0e6, k_torsion_nm_per_rad=120.0,
        travel_compress_m=0.020, travel_extend_m=0.010,
        contact_area_m2=0.045, tissue_area_m2=1.50e-3, tissue_modulus_pa=1.0e6,
    ),
    Segment(
        name="小腿-胫骨", length_m=0.426, mass_frac=0.0866,
        k_axial_n_per_m=1.3e7, k_torsion_nm_per_rad=350.0,
        travel_compress_m=0.012, travel_extend_m=0.008,
        contact_area_m2=0.030, tissue_area_m2=3.20e-4, tissue_modulus_pa=1.7e10,
    ),
    Segment(
        name="大腿", length_m=0.424, mass_frac=0.2832,
        k_axial_n_per_m=8.0e6, k_torsion_nm_per_rad=300.0,
        travel_compress_m=0.012, travel_extend_m=0.008,
        contact_area_m2=0.070, tissue_area_m2=5.00e-4, tissue_modulus_pa=1.7e10,
    ),
    Segment(
        name="骨盆-躯干", length_m=0.208, mass_frac=0.2480,
        k_axial_n_per_m=3.0e6, k_torsion_nm_per_rad=200.0,
        travel_compress_m=0.030, travel_extend_m=0.010,
        contact_area_m2=0.160, tissue_area_m2=1.00e-3, tissue_modulus_pa=1.0e7,
    ),
    Segment(
        name="胸廓", length_m=0.291, mass_frac=0.2854,
        k_axial_n_per_m=2.5e6, k_torsion_nm_per_rad=50.0,
        travel_compress_m=0.030, travel_extend_m=0.010,
        contact_area_m2=0.200, tissue_area_m2=1.20e-3, tissue_modulus_pa=1.0e7,
    ),
    Segment(
        name="头颈", length_m=0.315, mass_frac=0.0694,
        k_axial_n_per_m=2.0e6, k_torsion_nm_per_rad=30.0,
        travel_compress_m=0.020, travel_extend_m=0.008,
        contact_area_m2=0.055, tissue_area_m2=2.00e-4, tissue_modulus_pa=5.0e6,
    ),
)


# ==========================================================================
# 人体软组织柔顺层 + 串联接触（刚地面专用）
# ==========================================================================
@dataclass(frozen=True)
class TissueLayer:
    """人体软组织柔顺层（刚地面冲击时真正变形的那一层）。

    线性本构 ``σ = E·(x/t)``（跟垫量级：E≈1-2 MPa、t≈18 mm）。与
    ``hard_surface`` **串联**（骨/地面 → 软组织 → 足）。纯弹性、C¹。
    """

    thickness_m: float = 0.018
    modulus_pa: float = 2.0e6

    def k_pa_per_m(self) -> float:
        return self.modulus_pa / self.thickness_m


@dataclass(frozen=True)
class _SeriesContact:
    """两个纯弹性接触层**串联**的等效接触（力 = k_eff·comp）。

    串联刚度 ``1/k_eff = 1/(k1·A) + 1/(k2·A)``（A 为接触面积），
    故 ``force = comp·A·k1·k2/(k1+k2)``（k 单位为 Pa/m）。
    """

    k1_pa_per_m: float
    k2_pa_per_m: float

    def force_n(self, compression_m: float, area_m2: float) -> float:
        c = max(compression_m, 0.0)
        k_eff = self.k1_pa_per_m * self.k2_pa_per_m / (
            self.k1_pa_per_m + self.k2_pa_per_m
        )
        return k_eff * c * area_m2


def _contact_force(compression_m: float, area_m2: float, pad) -> float:
    """节点接触反力 (N)。``pad`` 需有 ``force_n(comp, area)``（``CrashPad``
    或 ``_SeriesContact``）。压缩量为负（未接触）时返回 0。
    """
    if compression_m <= 0.0 or area_m2 <= 0.0:
        return 0.0
    return float(pad.force_n(compression_m, area_m2))


# ==========================================================================
# 配置（冻结）
# ==========================================================================
@dataclass(frozen=True)
class PlanarChainConfig:
    """平面链配置。所有权重字段默认值 = 冻结来源，**非调参**。

    Attributes
    ----------
    mass_kg
        总质量 (kg)。
    g
        重力加速度（默认 ``climbing.G``）。
    mu_foot
        地-足摩擦系数（Coulomb，经 ``sign_reg`` 正则化 ⇒ C¹）。默认 0.4
        （``pad2d.Posture2D.mu_foot`` 量级）。
    tissue
        刚地面时串联的人体软组织柔顺层（垫工况不用）。
    rigid_k_n_per_m
        ``hard_surface`` 的 ``k_n``（``pad.py`` 冻结来源，未调）。
    k_comp_stop_mult
        关节压缩行程用尽后的止挡刚度倍数（× 关节轴向刚度）。压缩底是骨/
        组织，可以硬一点。
    k_ext_stop_mult
        关节拉伸行程用尽后的**软**止挡倍数。**明确不为硬止挡**（bodies.py
        的 5e6 过硬抗拉止挡是弹跳根因，此处软化为同量级）。
    c_axial_n_s_per_m
        关节轴向量纲阻尼 (N·s/m)，代表软组织的粘性损耗（量级估计）。
        能量守恒测试可置 0。
    dt_sample, max_ode_step, rtol, atol, t_max
        采样/积分参数（Radau）。
    """

    mass_kg: float = 80.0
    g: float = G
    mu_foot: float = 0.4
    tissue: TissueLayer = field(default_factory=TissueLayer)
    rigid_k_n_per_m: float = 1.0e7
    k_comp_stop_mult: float = 10.0
    k_ext_stop_mult: float = 1.0
    c_axial_n_s_per_m: float = 0.0
    zeta_joint: float = 0.3
    zeta_contact: float = 1.0
    dt_sample: float = 2e-4
    max_ode_step: float = 2e-5
    rtol: float = 1e-8
    atol: float = 1e-11
    t_max: float | None = None
    t_window_s: float = 0.25


#: 模块门（镜像 ``pad.RIGID_STANDING_CHAIN``）。默认 ``None`` = 不启用任何
#: 覆盖 ⇒ 链走 ``simulate_planar_chain_fall`` 的显式 ``config`` / 默认配置。
#: 供未来按工况门控启用（"无调用者默认不启用"，比改既有入口更干净）。
PLANAR_CHAIN: dict | None = None


# ==========================================================================
# 结果容器
# ==========================================================================
@dataclass
class PlanarChainResult:
    """平面链落地结果。字段镜像下游（探针/报告）所需。

    ``meta["impact_window_s"]`` = (首次接触秒, 首次分离秒)（无分离则到末点），
    与 ``pad.py`` / ``head_branch`` 同口径。
    """

    posture: str
    pad_label: str
    height_m: float
    mass_kg: float

    t_s: np.ndarray
    node_x: np.ndarray          # (T, 6)
    node_z: np.ndarray          # (T, 6)  向下为正
    node_vx: np.ndarray         # (T, 6)
    node_vz: np.ndarray         # (T, 6)
    joint_force_n: np.ndarray   # (T, 5)  轴向（拉伸为正）
    joint_compression_m: np.ndarray  # (T, 5) 正 = 压缩
    joint_angle_rad: np.ndarray  # (T, 5) 相邻杆夹角（kink）
    contact_force_n: np.ndarray  # (T, 6) 各节点法向接触力
    contact_area_m2: np.ndarray  # (T, 6)

    impact_speed_ms: float
    impact_energy_j: float
    peak_force_n: float
    max_compression_m: float
    peak_accel_head_g: float
    peak_accel_trunk_g: float
    bottomed_out: bool
    first_contact_node: int
    meta: dict = field(default_factory=dict)

    @property
    def peak_force_kn(self) -> float:
        return self.peak_force_n / 1e3

    @property
    def total_leg_compression_m(self) -> float:
        """踝+膝+髋（关节 0-2）峰值串联压缩。"""
        if self.joint_compression_m.size == 0:
            return 0.0
        t0, t1 = self.meta["impact_window_s"]
        m = (self.t_s >= t0) & (self.t_s <= t1)
        return float(np.max(self.joint_compression_m[m, :3].sum(axis=1)))

    def segment_axial_load_n(self, seg_index: int) -> np.ndarray:
        """段 ``seg_index`` 的轴向载荷时程 (T,)。

        节点 i = 段 i 的质心；**段的载荷 = 其下方关节（关节 i-1）的轴向力**
        （载荷自下传入）。段 0（足）无下方关节 ⇒ 用**接触力**近似其载荷。
        取绝对值（应力 = |F|/A）。
        """
        if seg_index <= 0:
            return np.abs(self.contact_force_n.sum(axis=1))
        return np.abs(self.joint_force_n[:, seg_index - 1])

    def node_accel_g(self, node: int) -> np.ndarray:
        return np.gradient(self.node_vz, self.t_s, axis=0)[:, node] / self.meta["g"]


# ==========================================================================
# 抗弯（kink）弹簧力 —— 内节点处相邻两杆夹角的恢复力矩
# ==========================================================================
def _kink_forces(a, b, c, k: float):
    """内节点处抗弯弹簧对三节点 (a,b,c) 的力（各 2-向量，坐标 (x,z)）。

    势能 ``V = ½·k·δ²``，``δ = atan2(p×q, p·q)`` = 杆 ``p=b-a`` 与 ``q=c-b``
    的相对转角（rest = 直线，δ=0）。力 = ``-∂V/∂节点位置``：
    ``F_c = -k·δ·∂δ/∂q``，``F_a = +k·δ·∂δ/∂p``，``F_b = -(F_a+F_c)``；
    ``∂δ/∂q = (-q_z, q_x)/|q|²``，``∂δ/∂p = (p_z, -p_x)/|p|²``。
    """
    p = b - a
    q = c - b
    pp = float(p[0] * p[0] + p[1] * p[1])
    qq = float(q[0] * q[0] + q[1] * q[1])
    if pp < 1e-18 or qq < 1e-18:
        z2 = np.zeros(2)
        return z2, z2, z2
    s = float(p[0] * q[1] - p[1] * q[0])
    cc = float(p[0] * q[0] + p[1] * q[1])
    delta = float(np.arctan2(s, cc))
    dq = np.array([-q[1], q[0]]) / qq           # ∂δ/∂q
    dp = np.array([p[1], -p[0]]) / pp           # ∂δ/∂p
    fc = -k * delta * dq
    fa = k * delta * dp
    return fa, -(fa + fc), fc


# ==========================================================================
# 求解
# ==========================================================================
@dataclass(frozen=True)
class _Ctx:
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
    pad: object
    zeta: float


def _forces(x, vx, z, vz, ctx: _Ctx):
    """单状态 (n,) 输入 → ``(ax, az, joint_force, contact_force)``。

    ``joint_force`` 带符号（拉伸为正），``contact_force`` 为各节点法向接触力。
    关节含轴向阻尼（``zeta_joint`` 的 Rayleigh 型）、接触含滞回阻尼
    （``zeta_contact``）。**``_rhs`` 与后处理共用本函数** ⇒ 前后处理逐位一致。
    """
    n = ctx.n_nodes
    ax = np.zeros(n)
    az = np.zeros(n) + ctx.g
    joint_force = np.zeros(n - 1)

    # ---- 关节：轴向弹簧 + 行程上下限 + 轴向阻尼 -------------------------
    for j in range(n - 1):
        d0 = x[j + 1] - x[j]
        d1 = z[j + 1] - z[j]
        ell = float(np.hypot(d0, d1))
        if ell < 1e-9:
            ell = 1e-9
        ux, uz = d0 / ell, d1 / ell
        f = ctx.k_ax[j] * (ell - ctx.rest[j])          # >0 拉伸
        comp = ctx.rest[j] - ell                        # >0 压缩
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

    # ---- 抗弯（转角）弹簧（内节点） ------------------------------------
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

    # ---- 逐节点接触（+ 滞回阻尼 + 摩擦） -------------------------------
    contact_force = np.zeros(n)
    for i in range(n):
        comp = z[i]                                     # 地面在 z=0
        if comp > 0.0:
            if i == 0:
                kappa = comp / ctx.posture.collapse_ref_m
                area = ctx.posture.area_at(kappa)
            else:
                area = ctx.area_base[i]
            fn = _contact_force(comp, area, ctx.pad)
            if fn > 0.0 and ctx.zeta > 0.0:
                # 接触滞回/回弹耗散（真实落地**非弹性**的来源）：阻尼器只在
                # **回弹相**（vz<0）作用 —— 加载相不加（避免撞击瞬间 c·v 尖峰），
                # 卸载相强耗散把弹性回弹能量吃掉 ⇒ 力过塌陷后单调、无回弹。
                # 有效质量 = 全链质量（接触承受整体减速）。
                k_sec = fn / comp
                c = 2.0 * ctx.zeta * float(np.sqrt(k_sec * ctx.m_total))
                if vz[i] < 0.0:
                    fn = max(fn + c * vz[i], 0.0)
            contact_force[i] = fn
            az[i] -= fn / ctx.m[i]
            if ctx.mu != 0.0 and fn > 0.0:
                fx = -ctx.mu * fn * float(sign_reg(vx[i], eps=1e-3))
                ax[i] += fx / ctx.m[i]
    return ax, az, joint_force, contact_force


def _rhs(t, y, ctx: _Ctx):
    """状态导数。状态交错布局 ``[x0,vx0,z0,vz0, x1,vx1,z1,vz1, ...]``，
    返回同布局（``[vx, ax, vz, az]``）。``z`` 向下为正。
    """
    x = y[0::4]
    vx = y[1::4]
    z = y[2::4]
    vz = y[3::4]
    ax, az, _, _ = _forces(x, vx, z, vz, ctx)
    out = np.empty(4 * ctx.n_nodes)
    out[0::4] = vx
    out[1::4] = ax
    out[2::4] = vz
    out[3::4] = az
    return out


def _build_ctx(posture: Posture, on_pad: bool, cfg: PlanarChainConfig) -> _Ctx:
    n = len(SEGMENTS_6)
    m = np.array([cfg.mass_kg * s.mass_frac for s in SEGMENTS_6])
    m *= cfg.mass_kg / m.sum()                          # 归一化到总质量
    lengths = np.array([s.length_m for s in SEGMENTS_6])
    rest = 0.5 * (lengths[:-1] + lengths[1:])           # 相邻质心间距
    k_ax_seg = np.array([s.k_axial_n_per_m for s in SEGMENTS_6])
    k_ax = 1.0 / (1.0 / k_ax_seg[:-1] + 1.0 / k_ax_seg[1:])   # 两段串联
    k_tors = np.array([s.k_torsion_nm_per_rad for s in SEGMENTS_6[:-1]])
    comp_limit = np.array([s.travel_compress_m for s in SEGMENTS_6[:-1]])
    ext_limit = np.array([s.travel_extend_m for s in SEGMENTS_6[:-1]])
    k_comp_stop = cfg.k_comp_stop_mult * k_ax
    k_ext_stop = cfg.k_ext_stop_mult * k_ax
    # 关节轴向阻尼（Rayleigh 型：c = 2ζ·√(k·m_above)），代表肌肉**离心收缩**
    # 与软组织的粘性耗能 —— 这是真实落地**非弹性**（不回弹）的物理来源。
    m_above = np.array([m[j + 1:].sum() for j in range(n - 1)])
    c_ax = 2.0 * cfg.zeta_joint * np.sqrt(k_ax * m_above)
    if cfg.c_axial_n_s_per_m > 0.0:
        c_ax = c_ax + cfg.c_axial_n_s_per_m
    area_base = np.array([s.contact_area_m2 for s in SEGMENTS_6])

    if on_pad:
        pad: object = CrashPad()
    else:
        pad = _SeriesContact(
            k1_pa_per_m=cfg.rigid_k_n_per_m,
            k2_pa_per_m=cfg.tissue.k_pa_per_m(),
        )

    return _Ctx(
        n_nodes=n, m=m, rest=rest, k_ax=k_ax, k_tors=k_tors,
        comp_limit=comp_limit, ext_limit=ext_limit,
        k_comp_stop=k_comp_stop, k_ext_stop=k_ext_stop, c_ax=c_ax,
        area_base=area_base, posture=posture, mu=cfg.mu_foot,
        m_total=float(m.sum()), g=cfg.g, pad=pad,
        zeta=cfg.zeta_contact,
    )


def simulate_planar_chain_fall(
    height_m: float = 3.0,
    mass_kg: float = 80.0,
    posture: str | Posture = "feet-first-stiff",
    on_pad: bool = False,
    config: PlanarChainConfig | None = None,
) -> PlanarChainResult:
    """模拟一次抱石坠落（6 段平面软关节链）。

    Parameters
    ----------
    height_m
        释放点距地面高度 (m)。足底接触时下落量 = ``height_m``。
    mass_kg
        总质量 (kg)。
    posture
        ``pad.POSTURES`` 的键或 ``Posture``（接触面积调度 ``A(κ)`` 来源）。
    on_pad
        True = 软垫（``CrashPad``）；False = 刚性地面（``hard_surface`` +
        人体软组织柔顺层，串联）。
    config
        ``PlanarChainConfig``；None 用默认（``PLANAR_CHAIN`` 门若置非空则覆盖）。
    """
    post = POSTURES[posture] if isinstance(posture, str) else posture
    cfg = config or PlanarChainConfig(mass_kg=mass_kg)
    # 模块门（默认 None）：非空时覆盖配置字段（镜像 pad.RIGID_STANDING_CHAIN）。
    if PLANAR_CHAIN:
        cfg = replace(cfg, **PLANAR_CHAIN)
    cfg = replace(cfg, mass_kg=float(mass_kg))

    ctx = _build_ctx(post, on_pad, cfg)
    n = ctx.n_nodes

    # ---- 初始几何：竖直链，足底节点在 z = -height ----------------------
    z0 = np.empty(n)
    z0[0] = -float(height_m)
    for i in range(n - 1):
        z0[i + 1] = z0[i] - ctx.rest[i]
    x0 = np.zeros(n)

    # ---- 弹道起点：自由落体解析到**首次接触**那一刻，再启刚体积分 --------
    # 自由落体段无内力（初始无预应力），可解析：t_c=√(2h/g)、v=g·t_c。
    # 直接以接触瞬间为初值起步，避免 Radau 空跑数百毫秒的无谓刚性步
    # （h=20 时 t_c≈2.02 s，否则仅自由落体就 1e5+ 步）。物理等价。
    t_contact = float(np.sqrt(2.0 * height_m / cfg.g))
    fall = 0.5 * cfg.g * t_contact * t_contact
    z0 = z0 + fall
    vz0 = np.full(n, cfg.g * t_contact)

    y0 = np.empty(4 * n)
    y0[0::4] = x0
    y0[1::4] = 0.0
    y0[2::4] = z0
    y0[3::4] = vz0

    # ---- 终止：全节点已在接触（节点 0 已穿透）且近静止 ------------------
    def settled(_t, y):
        zz = y[2::4]
        vv = y[3::4]
        touched = zz[0] > 0.0
        return 1.0 if (touched and np.all(np.abs(vv) < 0.05)) else -1.0

    settled.terminal = True
    settled.direction = 1.0

    if cfg.t_max is not None:
        t_end = cfg.t_max
    else:
        t_end = t_contact + cfg.t_window_s

    sol = solve_ivp(
        lambda t, y: _rhs(t, y, ctx),
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
    node_x = Y[0::4].T
    node_vx = Y[1::4].T
    node_z = Y[2::4].T
    node_vz = Y[3::4].T

    # ---- 关节/接触量：用**与 _rhs 同一个 _forces** 逐样本重算 ⇒ 前后一致 --
    dx = node_x[:, 1:] - node_x[:, :-1]
    dz = node_z[:, 1:] - node_z[:, :-1]
    ell = np.maximum(np.hypot(dx, dz), 1e-9)
    joint_comp = ctx.rest[None, :] - ell
    # kink 角（相邻杆夹角）
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
        _, _, jf, cf = _forces(
            node_x[k], node_vx[k], node_z[k], node_vz[k], ctx
        )
        joint_force[k] = jf
        contact_force[k] = cf

    # 接触面积（供报告/诊断）
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
    pad_max_comp = ctx.pad.max_compression_m if hasattr(ctx.pad, "max_compression_m") else 0.20
    bottomed = bool(max_leg > pad_max_comp * 0.999)

    first = int(np.argmax(np.any(contact_force[sl] > 0, axis=0))) if i1 > i0 else 0

    return PlanarChainResult(
        posture=post.name,
        pad_label=("Crash pad" if on_pad else "硬地面（离垫）"),
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
        },
    )


def segment_local_peak(result: PlanarChainResult, seg) -> tuple[float, float]:
    """段 ``seg``（``Segment`` 或段下标）的局部组织峰。

    返回 ``(t_peak_s, value)``：``value = ε_seg = (F_seg/A_seg)/E_seg`` 在
    ``meta["impact_window_s"]`` 内**首个显著局部极大**（复用
    ``head_branch._first_significant_peak``，防无阻尼多峰的 argmax 非确定）。
    """
    if isinstance(seg, Segment):
        idx = SEGMENTS_6.index(seg)
        s = seg
    else:
        idx = int(seg)
        s = SEGMENTS_6[idx]
    load = result.segment_axial_load_n(idx)
    eps = (load / s.tissue_area_m2) / s.tissue_modulus_pa
    t = result.t_s
    t0, t1 = result.meta["impact_window_s"]
    m = (t >= t0) & (t <= t1)
    if not m.any():
        raise ValueError("冲击窗口内无数据点，无法定位段局部峰")
    eps_w = np.abs(eps[m])
    t_w = t[m]
    k = _first_significant_peak(eps_w)
    return float(t_w[k] - t0), float(eps_w[k])
