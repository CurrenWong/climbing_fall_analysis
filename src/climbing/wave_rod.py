"""G1 · 预注册备选件 —— 1D 应力波杆（胫骨轴向波，自带观测量）。

为什么需要这个模块
------------------
冻结设计 ``docs/methods/G1_多段链设计.md`` §7 决策 6 预注册：**若 6 段平面
软关节链在冻结参数下 V8b FAIL（不重调），升级到 1D 应力波杆**（自带观测量
定义）。链的实跑结果见 ``results/G1_CHAIN2D_REPORT.md``：20 m 胫骨局部组织
峰 **5.8 ms > 5.0 ms** ⇒ **V8b FAIL** ⇒ 本备选触发。

本模块把胫骨当成**一根一维轴向弹性杆**，用**第一原理**定义波动可观测量，
以检验 [T22] Table 3 的胫骨峰时（5/10/20 m → 11.9/6.5/4.6 ms）是否可能是
**波动主导**的。**不调参**：材料/几何值全部取自文献（见下），跑后如实判定。

观测量定义（第一原理，跑前冻结）
--------------------------------
1. **轴向波速**：``c = sqrt(E/rho)``。取**皮质骨**文档值
   ``E ≈ 17 GPa``、``rho ≈ 1900 kg/m³``（见 ``WaveRod`` 注释与来源），
   得 ``c ≈ 2.99e3 m/s``。
2. **胫骨几何**：``L ≈ 0.38 m``（成人胫骨长度，Drillis & Contini 1966 /
   Winter 2009 转引；与本仓 ``chain2d`` 段长同源族）。
3. **观测点**：**远端胫骨**（impact 端附近，``x = dx`` 的内部节点）的
   局部组织应变 ``ε(x, t) = ∂u/∂x``（压缩为负）。
4. **驱动**：**impact 端（``x = 0``）由接触/GRF 驱动**（Neumann 牵引边界）。
   刚壁撞击的首个 GRF 幅值取 Saint-Venant 结果 ``σ0 = ρ·c·v0``
   （``v0 = sqrt(2 g h)`` 为撞击速度）⇒ 驱动应变 ``ε0 = σ0/E = v0/c``。
   幅值只定标应变**大小**，**不影响峰时**。
5. **近端（``x = L``）**：**固定端**（刚体/大质量近似）。判据：胫骨杆质量
   ``ρ A L ≈ 0.23 kg``，其上人体质量 ``≳ 70 kg``，阻抗比 ``≫ 1``（本模块
   在 ``meta`` 记录该比值）⇒ 在亚毫秒波时标上近端近似不动。
6. **峰时（预注册预测公式）**：撞击在 ``x=0`` 出发的压缩波在 ``x=L`` 固定端
   **同号反射**，回到 ``x=0`` 处使远端应变**叠加/加倍**，往返用时
   ``t_peak = 2L/c``。

**跑前预测值**：``2L/c = 2×0.38/2991.2 = 2.54e-4 s ≈ 0.254 ms``（**亚毫秒**）。
设计 §4.3 预期：1D 波杆会**亚毫秒过冲**，[T22] 的 4.6 ms 应是**接触压溃
主导**而非波动主导。本模块唯一职责是把该预测**干净地算出来并如实报告**，
**不凑 4.6 ms / 5.0 ms**。

诚实边界
--------
* 这是一根**均匀、无耗散、固定近端**的一维弹性杆：不含骨小梁、骨膜、肌肉、
  关节；近端用固定端近似（阻抗比 ≫ 1，见 ``meta``）。
* 观测量（远端应变 ``|ε|`` 的全局极大）**跑前冻结**；峰时对高度**不敏感**
  （``ε0`` 只放大幅值、不改 ``2L/c``）—— 这本身就是本件的结论之一。
* 纯增量：只新增文件，不改任何既有模块/测试/results。

依赖：**numpy**（仅此）。积分用**定步显式蛙跳**（CFL 安全、**确定性、可复现**，
无自适应步噪声）。
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np

from . import G

__all__ = [
    "WaveRod",
    "WAVE_ROD",
    "WaveRodResult",
    "simulate_tibia_wave",
    "wave_observable",
]


# ==========================================================================
# 配置（材料 / 几何，跑前冻结，来源文档化，非调参）
# ==========================================================================
@dataclass(frozen=True)
class WaveRod:
    """胫骨 1D 轴向弹性杆配置。默认值 = **文档化文献值**，**非调参**。

    Attributes
    ----------
    length_m
        胫骨长度 (m)。成人胫骨 ≈ 0.38 m。来源族：Drillis & Contini (1966)
        节段长（Winter 2009 转引），与本仓 ``chain2d.SEGMENTS_6`` 小腿段
        (0.426 m 含足踝) 同源；此处取**胫骨骨干**典型值 0.38 m。
    modulus_pa
        皮质骨**轴向杨氏模量** (Pa)。文献区间 ~15–20 GPa；取 17 GPa
        （骨力学教材典型值，如 Cowin, *Bone Mechanics*；亦为
        ``chain2d`` 小腿段 ``tissue_modulus_pa`` 所用值）。
    density_kg_m3
        皮质骨密度 (kg/m³)。文献 ~1800–2100；取 1900。
    area_m2
        胫骨承载截面 (m²)。~3.2e-4 m²（≈ 直径 2 cm 圆截面，
        与 ``chain2d`` 小腿段 ``tissue_area_m2`` 同值）。
    n_cells
        空间网格数（``dx = L/n_cells``）。
    cfl
        蛙跳 Courant 数 ``c·dt/dx``，取 0.5（< 1 ⇒ 稳定）。
    round_trips
        观测窗 = ``round_trips × (2L/c)``（默认 3 个往返）。
    mass_above_fraction
        "胫骨以上"人体质量占全身比（仅用于**记录近端阻抗比**判据，不影响
        峰时）。默认 0.90。
    """

    length_m: float = 0.38
    modulus_pa: float = 17.0e9
    density_kg_m3: float = 1900.0
    area_m2: float = 3.20e-4
    n_cells: int = 160
    cfl: float = 0.5
    round_trips: float = 3.0
    mass_above_fraction: float = 0.90

    @property
    def wave_speed_ms(self) -> float:
        """轴向波速 ``c = sqrt(E/rho)`` (m/s)。"""
        return float(np.sqrt(self.modulus_pa / self.density_kg_m3))

    @property
    def transit_s(self) -> float:
        """单程过杆时间 ``L/c`` (s)。"""
        return self.length_m / self.wave_speed_ms

    @property
    def round_trip_s(self) -> float:
        """往返时间 ``2L/c`` (s) = 预注册的峰时预测。"""
        return 2.0 * self.length_m / self.wave_speed_ms


#: 模块门（镜像 ``pad.RIGID_STANDING_CHAIN`` / ``chain2d.PLANAR_CHAIN``）。
#: 默认 ``None`` = 不启用任何覆盖 ⇒ 用 ``simulate_tibia_wave`` 的显式
#: ``config`` / 默认 ``WaveRod()``。非空时按字段名覆盖（``dataclasses.replace``）。
WAVE_ROD: dict | None = None


# ==========================================================================
# 结果容器
# ==========================================================================
@dataclass
class WaveRodResult:
    """1D 应力波杆求解结果。

    Attributes
    ----------
    height_m, mass_kg
        工况（释放高度 / 总质量）。
    t_s
        采样时刻 (s)，从撞击起点 0。
    x_m
        节点坐标 (m)，``0`` = 远端（impact），``L`` = 近端（固定）。
    strain
        应变场 ``(T, n_cells+1)``，压缩为负。
    distal_strain
        远端观测点（``x = dx``）的应变时程 ``(T,)``。
    peak_time_s
        观测窗内 ``|distal_strain|`` 的**全局极大时刻**（= 波观测量的峰时）。
    peak_strain
        上述时刻的**带符号**应变值。
    round_trip_s, transit_s
        ``2L/c``、``L/c``（解析参考）。
    analytic_peak_time_s
        预注册预测峰时 ``= 2L/c``。
    impact_speed_ms
        撞击速度 ``v0 = sqrt(2 g h)``（m/s）。
    impact_strain
        撞击驱动应变 ``ε0 = v0/c``（远端边界幅值）。
    window_s
        观测窗（s），``= round_trips × 2L/c``。
    meta
        诊断/溯源字典（波速、阻抗比、网格等）。
    """

    height_m: float
    mass_kg: float
    t_s: np.ndarray
    x_m: np.ndarray
    strain: np.ndarray
    distal_strain: np.ndarray
    peak_time_s: float
    peak_strain: float
    round_trip_s: float
    transit_s: float
    analytic_peak_time_s: float
    impact_speed_ms: float
    impact_strain: float
    window_s: float
    meta: dict = field(default_factory=dict)

    @property
    def wave_speed_ms(self) -> float:
        return float(self.meta["wave_speed_ms"])

    @property
    def peak_time_ms(self) -> float:
        return self.peak_time_s * 1e3

    @property
    def submillisecond(self) -> bool:
        """峰时是否 < 1 ms（预注册的"亚毫秒过冲"判据）。"""
        return self.peak_time_s < 1.0e-3


# ==========================================================================
# 求解（定步显式蛙跳，CFL 安全 ⇒ 确定性、可复现）
# ==========================================================================
def simulate_tibia_wave(
    height_m: float = 20.0,
    mass_kg: float = 80.0,
    config: WaveRod | None = None,
) -> WaveRodResult:
    """模拟一次落地在**胫骨**中激发的 1D 轴向应力波。

    方程：``u_tt = c²·u_xx``（``u`` 轴向位移），边界
    ``x=0`` 远端：``EA·u_x = -F(t)``（Neumann 牵引，F = 接触/GRF，
    刚壁撞击幅值 ``ρc·v0·A``）；``x=L`` 近端：``u=0``（固定端，阻抗比 ≫ 1）。
    初值：静止、无预应力，``t=0`` 起施加 GRF。

    Parameters
    ----------
    height_m
        释放高度 (m)；撞击速度 ``v0 = sqrt(2 g h)``（幅值只定标应变大小）。
    mass_kg
        总质量 (kg)（记录用；亚毫秒波时标上不影响峰时，见 ``meta``）。
    config
        ``WaveRod``；None 用默认（``WAVE_ROD`` 门若置非空则覆盖）。

    Returns
    -------
    WaveRodResult
    """
    cfg = config if config is not None else WaveRod()
    # 模块门（默认 None）：非空时按字段名覆盖（镜像 pad/chain2d 门）。
    if WAVE_ROD:
        cfg = replace(cfg, **WAVE_ROD)

    L = float(cfg.length_m)
    E = float(cfg.modulus_pa)
    rho = float(cfg.density_kg_m3)
    A = float(cfg.area_m2)
    n = int(cfg.n_cells)
    c = cfg.wave_speed_ms
    dx = L / n
    dt = float(cfg.cfl) * dx / c
    cfl = c * dt / dx
    c2 = c * c

    # ---- 驱动 / 撞击量 --------------------------------------------------
    v0 = float(np.sqrt(2.0 * G * float(height_m)))
    eps0 = v0 / c                     # = σ0/E，σ0=ρ c v0（Saint-Venant）
    F0 = rho * c * v0 * A             # 刚壁撞击首个 GRF 幅值 (N)

    window = float(cfg.round_trips) * cfg.round_trip_s
    n_steps = max(int(round(window / dt)), 4)

    # ---- 蛙跳：u_prev, u, u_next ---------------------------------------
    u_prev = np.zeros(n + 1)
    u = np.zeros(n + 1)
    t_s = np.empty(n_steps + 1)
    distal = np.empty(n_steps + 1)
    strain = np.empty((n_steps + 1, n + 1))

    def _accel(ui: np.ndarray) -> np.ndarray:
        a = np.zeros(n + 1)
        a[1:n] = c2 * (ui[2:n + 1] - 2.0 * ui[1:n] + ui[0:n - 1]) / dx ** 2
        # 远端 Neumann：EA·u_x(0,t) = -F0 ⇒ ghost u_{-1} = u_1 + 2 dx F0/(EA)
        a[0] = c2 * (2.0 * ui[1] - 2.0 * ui[0] + 2.0 * dx * F0 / (E * A)) / dx ** 2
        return a

    def _strain_field(ui: np.ndarray) -> np.ndarray:
        e = np.zeros(n + 1)
        e[1:n] = (ui[2:n + 1] - ui[0:n - 1]) / (2.0 * dx)
        # 远端节点应变 = 边界值 (u_1 - u_{-1})/(2dx) = -eps0（驱动施加）
        e[0] = -eps0
        # 近端固定节点应变用单侧差分
        e[n] = (ui[n] - ui[n - 1]) / dx
        return e

    u_next = u + 0.5 * dt * dt * _accel(u)     # 从静止起步的第一步
    t_s[0] = 0.0
    strain[0] = _strain_field(u_next)
    distal[0] = strain[0, 1]

    for k in range(1, n_steps + 1):
        un = 2.0 * u - u_prev + dt * dt * _accel(u)
        un[n] = 0.0                            # 近端固定
        u_prev, u = u, un
        t_s[k] = k * dt
        strain[k] = _strain_field(u)
        distal[k] = strain[k, 1]

    # ---- 观测量：远端 |ε| 全局极大（= 波观测量的峰时） -------------------
    peak_time_s, peak_strain = wave_observable_from_series(t_s, distal)

    # 近端阻抗比判据（记录用）：m_above / (ρ A L)
    rod_mass = rho * A * L
    m_above = float(cfg.mass_above_fraction) * float(mass_kg)

    return WaveRodResult(
        height_m=float(height_m),
        mass_kg=float(mass_kg),
        t_s=t_s,
        x_m=np.linspace(0.0, L, n + 1),
        strain=strain,
        distal_strain=distal,
        peak_time_s=peak_time_s,
        peak_strain=peak_strain,
        round_trip_s=cfg.round_trip_s,
        transit_s=cfg.transit_s,
        analytic_peak_time_s=cfg.round_trip_s,
        impact_speed_ms=v0,
        impact_strain=eps0,
        window_s=window,
        meta={
            "wave_speed_ms": c,
            "modulus_pa": E,
            "density_kg_m3": rho,
            "area_m2": A,
            "length_m": L,
            "n_cells": n,
            "dx_m": dx,
            "dt_s": dt,
            "cfl": cfl,
            "n_steps": n_steps,
            "grf_amplitude_n": F0,
            "rod_mass_kg": rod_mass,
            "mass_above_kg": m_above,
            "impedance_ratio_m_above_over_rod": m_above / rod_mass,
            "proximal_end": "fixed (rigid) — justified by impedance ratio >> 1",
            "distal_probe": "node 1 (x=dx), local tissue strain eps=du/dx",
            "observable": "argmax |distal strain| within window [0, round_trips*(2L/c)]",
            "analytic_peak_time_s": cfg.round_trip_s,
        },
    )


# ==========================================================================
# 波观测量（文档化）
# ==========================================================================
def wave_observable_from_series(
    t_s: np.ndarray, distal_strain: np.ndarray
) -> tuple[float, float]:
    """从远端应变时程取**波观测量**：``(t_peak_s, peak_strain)``。

    ``t_peak_s`` = ``|distal_strain|`` 的**全局极大时刻**（``np.argmax``，
    平顶时取**最早**出现）。按预注册第一原理，该峰对应压缩波在近端固定端
    同号反射、返回远端后应变叠加的**往返时刻 ``2L/c``**。

    这里**不用** ``head_branch._first_significant_peak``：可复现地，行波前沿
    本身即可越过其 0.5×max 相对阈值，会把前沿的数值前沿当作"首峰"；
    全局 ``argmax`` 才是物理上的应变叠加峰，且**确定、种子无关**。
    """
    mag = np.abs(np.asarray(distal_strain, dtype=float))
    k = int(np.argmax(mag))
    return float(t_s[k]), float(distal_strain[k])


def wave_observable(result: WaveRodResult) -> tuple[float, float]:
    """``WaveRodResult`` → 远端局部组织应变的波观测量 ``(t_peak_s, ε_peak)``。

    峰时 = 远端 ``|ε|`` 全局极大时刻（往返叠加）；``ε_peak`` 为其带符号值。
    """
    return wave_observable_from_series(result.t_s, result.distal_strain)
