"""Phase B-3 —— 抱石坠落 2D 模型（z + x + θ），1D 模型严格保留作为子结构。

为什么是 2D 而不是 1D
----------------------
1D 模型的"塌陷进度 κ → 接触面积 A(κ)"机制决定了**峰值力**与**压缩量**，
但漏掉了 [B25] 表格里 36% 的**纵向旋转**和 9% 的**前后向旋转** ——
也就是"踝旋后/旋前"和"躯干前倾/后仰"两类机制。

P3 的 2D 扩展目标：

* 让模型在**默认（所有 2D 参数为 0）** 条件下与 1D 模型逐位一致
  （``np.array_equal``，**不是** ``pytest.approx``）。
* 让 [T22] 提到的"内翻 vs 外翻不对称"和"峰值时刻随高度提前"
  两条不变量成为**可测试的代码**。
* 物理上，把身体当作**绕脚接触点摆动的单刚体**：
  - ``x_f``：脚的横向位移，``vx_f`` = 横向速度；
  - ``theta``：身体相对铅垂的倾角，``omega`` = 角速度；
  - ``L_arm = leg_ext0 = 0.85 m``：脚接触点到身体质心的等效臂长。

设计选择
--------
**纵向子系统不重写**。1D 模型用 ``pad.simulate_boulder_fall`` 解出
``f_pad(t)``，2D 方程的 RHS 把它作为**已知驱动力**线性插值调用。
这是因为本 2D 扩展在结构上**不引入新的纵向自由度** —— ``f_pad`` 的
计算只依赖 ``z_f`` 与 ``kappa``，与 ``x_f``、``theta`` 无关，所以
纵向方程可以独立求解，2D 转动 / 横向方程受其驱动。

**摩擦不写硬 sign**。``f_pad_x = -mu_foot * f_pad * sign_reg(vx_f)``
用 ``sign_reg(v) = v/sqrt(v²+ε²)`` 正则化，避免 RK45 在 ``vx_f=0``
附近的硬不连续处退化为"步长被压到机器精度"。

**默认 → 1D**。所有 2D 字段默认 0.0（``theta0, omega0, lateral_offset,
foot_x_lever, mu_foot, I_body, c_rot``），叠加后 2D RHS 全零，纵向
子系统原样解出 ⇒ 与 1D 完全 bit-identical。
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp

from . import G
from .pad import (
    BoulderFallResult,
    CrashPad,
    Posture,
    simulate_boulder_fall,
)


__all__ = [
    "Posture2D",
    "BoulderFallResult2D",
    "simulate_boulder_fall_2d",
    "sign_reg",
    "LEG_EXT0_M",
]


# 1D pad.py:480 中的 leg_ext0（站立时躯干 CoM 在双足之上的初始高度差）。
# 作本模块默认的 ``L_arm`` —— 脚接触点到身体质心的等效臂长。
# 必须在两模块之间手动同步：任何对 pad.py leg_ext0 的修改都需要同时
# 更新这里的常量（详见 __init__.py 中的约定）。
LEG_EXT0_M: float = 0.85


def sign_reg(v: np.ndarray | float, eps: float = 1e-3) -> np.ndarray | float:
    """正则化 sign —— Coulomb 摩擦的硬 ``sign()`` 是 RK45 的著名刚性陷阱。

    返回 ``v / sqrt(v² + ε²)``，性质：

    * ``|v| >> ε`` 时 → ≈ ``sign(v)``
    * ``|v| << ε`` 时 → ≈ ``v / ε``，平滑过零，绝对值最大为 1
    * 完全可微 → RK45 可大步通过零点

    这把"硬不连续的粘-滑摩擦"换成"软化的等效阻尼"，后者在 ``|v|≈0``
    附近给出连续的横向力，且 ``|f_pad_x| ≤ mu·f_pad``（与 Coulomb
    上限一致）。
    """
    v_arr = np.asarray(v, dtype=float)
    return v_arr / np.sqrt(v_arr * v_arr + eps * eps)


@dataclass(frozen=True)
class Posture2D:
    """2D 姿势扩展 —— 纯数据载体，不影响现有 ``pad.Posture`` 的 1D 行为。

    所有字段默认 0.0 ⇒ 默认与 1D 完全 bit-identical。

    Attributes
    ----------
    theta0_rad
        初始身体倾角（rad）。正负约定与 ``foot_x_lever_m`` 的方向配合
        决定是内翻还是外翻；测试里只需关注**符号敏感**的不对称，
        不假设物理方向。
    omega0_rad_s
        初始角速度。
    lateral_offset_m
        脚的初始横向偏移 ``x_f(0)``。
    foot_x_lever_m
        接触力作用点到脚几何中心的横向偏移（m）。非零时会在
        ``f_pad`` 与力臂 ``L_arm*sin(theta)`` 之外多出一段常数力臂，
        引入**符号敏感**的转动激励 —— 这是内翻 vs 外翻不对称的来源。
    mu_foot
        脚-垫摩擦系数。无量纲。
    I_body_kgm2
        身体绕脚的转动惯量（kg·m²）。**为 0 时跳过转动方程更新**，
        保持 ``theta`` 与 ``omega`` 在初值上不动，避免 0/0。
    c_rot
        转动阻尼（kg·m²/s）。默认 0。
    """

    theta0_rad: float = 0.0
    omega0_rad_s: float = 0.0
    lateral_offset_m: float = 0.0
    foot_x_lever_m: float = 0.0
    mu_foot: float = 0.0
    I_body_kgm2: float = 0.0
    c_rot: float = 0.0


@dataclass
class BoulderFallResult2D(BoulderFallResult):
    """2D 版 ``BoulderFallResult`` —— 继承所有 1D 字段，附加：

    * ``x_foot_m`` —— 脚横向位移时程 (m)
    * ``theta_rad`` —— 身体倾角时程 (rad)
    * ``omega_rad_s`` —— 角速度时程 (rad/s)
    * ``peak_theta_rad`` —— ``max(|theta|)`` 标量
    * ``peak_omega_rad_s`` —— ``max(|omega|)`` 标量

    1D 字段（``t_s``, ``z_foot_m``, ``pad_force_n`` 等）通过共享
    ``BoulderFallResult`` 的存储直接复用，**默认参数下与 1D 轨迹
    ``np.array_equal``**（bit-identicality 验证见
    ``tests/test_pad2d.py::TestInvariant1BitIdentical``）。
    """

    x_foot_m: np.ndarray = field(default_factory=lambda: np.zeros(0))
    theta_rad: np.ndarray = field(default_factory=lambda: np.zeros(0))
    omega_rad_s: np.ndarray = field(default_factory=lambda: np.zeros(0))
    peak_theta_rad: float = 0.0
    peak_omega_rad_s: float = 0.0


def simulate_boulder_fall_2d(
    *,
    height_m: float = 3.0,
    mass_kg: float = 80.0,
    posture: str | Posture = "controlled-drop",
    pad: CrashPad | None = None,
    on_pad: bool = True,
    g: float = G,
    dt_sample: float = 2e-4,
    t_max: float | None = None,
    max_ode_step: float = 1e-4,
    rtol: float = 1e-8,
    atol: float = 1e-11,
    posture2d: Posture2D | None = None,
) -> BoulderFallResult2D:
    """模拟一次抱石坠落的 2D 平面动力学。

    Parameters
    ----------
    height_m, mass_kg, posture, pad, on_pad, g, dt_sample, t_max, max_ode_step, rtol, atol
        与 ``pad.simulate_boulder_fall`` 同名同义。1D 子系统直接调用
        ``pad.simulate_boulder_fall`` 求解，**不允许重写纵向方程**。
    posture2d
        2D 姿势扩展字段集合，默认 ``Posture2D()``（全部 0.0） ⇒
        结果与 1D **bit-identical**（``np.array_equal``）。

    Returns
    -------
    BoulderFallResult2D
        继承 ``BoulderFallResult`` 的全部字段（包括 ``meta`` 里的
        ``posture_obj`` / ``m_up`` / ``m_low`` / ``impact_window_s``
        / ``head_is_contact`` / ``f_flex_peak_n`` 等
        ``injury.py`` / ``bone.py`` 需要的键），附加 2D 状态量。

    Notes
    -----
    **耦合结构**：纵向子系统的 ``f_pad(t)`` 作为已知驱动力进入 2D RHS，
    横向 / 转动子系统之间无耦合（``vx_f_dot`` 不依赖 ``theta``，
    ``omega_dot`` 中 ``theta`` 项虽存在但与 ``vx_f`` 无关），所以 2D
    4 维状态可以一起解。默认 ``Posture2D()`` 下 2D RHS 全 0，2D 状态
    保持初值 0 ⇒ 1D 轨迹原样返回。
    """
    p2d = posture2d or Posture2D()

    # ---- Step 1: 解 1D 子系统（调用 pad.simulate_boulder_fall，bit-identical）----
    r1d = simulate_boulder_fall(
        height_m=height_m, mass_kg=mass_kg,
        posture=posture, pad=pad, on_pad=on_pad,
        g=g, dt_sample=dt_sample, t_max=t_max,
        max_ode_step=max_ode_step, rtol=rtol, atol=atol,
    )

    # ---- Step 2: 解 2D 子系统（受 f_pad(t) 驱动）----
    t_1d = r1d.t_s
    f_pad_1d = r1d.pad_force_n

    L_arm = LEG_EXT0_M
    mu = p2d.mu_foot
    foot_x_lever = p2d.foot_x_lever_m
    I_body = p2d.I_body_kgm2
    c_rot = p2d.c_rot
    m_total = float(mass_kg)

    # 检测"全零 2D"短路：所有 2D 参数都为零 ⇒ RHS 全零 ⇒ 状态保持初值。
    # 防止 RK45 在零域里产生 1e-17 量级的浮点噪声污染 bit-identical 测试。
    trivially_zero_2d = (
        mu == 0.0
        and I_body == 0.0
        and c_rot == 0.0
        and foot_x_lever == 0.0
        and p2d.theta0_rad == 0.0
        and p2d.omega0_rad_s == 0.0
        and p2d.lateral_offset_m == 0.0
    )

    if trivially_zero_2d:
        # 短路：所有 2D 状态保持初值（零），返回与 r1d 同样形状的全零数组。
        x_foot = np.zeros_like(t_1d)
        vx_foot = np.zeros_like(t_1d)
        theta = np.zeros_like(t_1d)
        omega = np.zeros_like(t_1d)
    else:
        def rhs_2d(_t, y):
            xf, vxf, th, om = y
            # 从 1D 轨迹线性插值取 f_pad(t)
            f_pad = float(np.interp(_t, t_1d, f_pad_1d))
            # ---- 横向 ----
            # sign_reg 避免 RK45 在 vx=0 处被不连续卡死（粘-滑刚性陷阱）
            f_pad_x = -mu * f_pad * sign_reg(vxf, eps=1e-3)
            vxf_dot = f_pad_x / m_total
            # ---- 转动 ----
            # I_body=0 ⇒ 跳过（默认情形下被 trivially_zero_2d 短路，
            # 这里作为对"I_body 单独自零但其它参数非零"的兜底），避免 0/0。
            if I_body == 0.0 and c_rot == 0.0:
                om_dot = 0.0
            else:
                d_lever = foot_x_lever + L_arm * float(np.sin(th))
                om_dot = (
                    f_pad * d_lever
                    - m_total * g * L_arm * float(np.sin(th))
                    - c_rot * om
                ) / I_body
            return [vxf, vxf_dot, om, om_dot]

        y0_2d = [
            p2d.lateral_offset_m,
            0.0,
            p2d.theta0_rad,
            p2d.omega0_rad_s,
        ]
        sol_2d = solve_ivp(
            rhs_2d,
            (0.0, float(t_1d[-1])),
            y0_2d,
            method="RK45",
            t_eval=t_1d,            # 与 1D 同一采样栅格 ⇒ z-轨迹 bit-identical
            max_step=max_ode_step,
            rtol=rtol,
            atol=atol,
        )
        if not sol_2d.success:
            raise RuntimeError(f"2D 积分失败: {sol_2d.message}")
        x_foot = sol_2d.y[0]
        vx_foot = sol_2d.y[1]
        theta = sol_2d.y[2]
        omega = sol_2d.y[3]

    # ---- Step 3: 构造结果对象 ----
    peak_theta = float(np.max(np.abs(theta))) if theta.size else 0.0
    peak_omega = float(np.max(np.abs(omega))) if omega.size else 0.0

    # 浅拷贝 1D meta（包含 posture_obj, m_up, m_low, impact_window_s 等所有
    # injury.py / bone.py 需要的字段），不修改 r1d.meta 原对象。
    meta_2d = dict(r1d.meta)

    return BoulderFallResult2D(
        posture=r1d.posture,
        pad_label=r1d.pad_label,
        height_m=r1d.height_m,
        mass_kg=r1d.mass_kg,
        t_s=r1d.t_s,
        z_foot_m=r1d.z_foot_m,
        z_torso_m=r1d.z_torso_m,
        pad_compression_m=r1d.pad_compression_m,
        contact_area_m2=r1d.contact_area_m2,
        pad_force_n=r1d.pad_force_n,
        accel_leg_g=r1d.accel_leg_g,
        accel_torso_g=r1d.accel_torso_g,
        impact_speed_ms=r1d.impact_speed_ms,
        impact_energy_j=r1d.impact_energy_j,
        peak_force_n=r1d.peak_force_n,
        peak_accel_torso_g=r1d.peak_accel_torso_g,
        peak_accel_leg_g=r1d.peak_accel_leg_g,
        max_compression_m=r1d.max_compression_m,
        bottomed_out=r1d.bottomed_out,
        energy_into_pad_j=r1d.energy_into_pad_j,
        energy_into_flex_j=r1d.energy_into_flex_j,
        energy_residual_j=r1d.energy_residual_j,
        hic=r1d.hic,
        meta=meta_2d,
        x_foot_m=x_foot,
        theta_rad=theta,
        omega_rad_s=omega,
        peak_theta_rad=peak_theta,
        peak_omega_rad_s=peak_omega,
    )
