"""Phase A-1 —— 绳索坠落 1D 模型。

模型结构（自上而下三层）::

    固定锚点/保护器  ──[ 保护器：锁死 / 滑动 两态 ]──  自由绳段 L(t), 伸长 e(t)
                                                              │
                                                    非线性绳本构 T(ε)
                                                              │
                                                             人体 m

符号约定
--------
* 竖直向下为正方向。
* ``x(t)``  坠落者自释放点起的向下位移（m）
* ``v(t)``  向下速度（m/s）
* ``L(t)``  保护器与坠落者之间的**自由绳长**（m），保护器滑动时增长
* ``e(t)``  该自由绳段的**弹性伸长**（m），几何约束 ``L0 + x = L + e``

本构模型
--------
动态绳不是线性弹簧。真实动态绳的力-伸长曲线是 *硬化* 的：
低应变时较软（编织层松散），高应变时急剧变硬（纱线拉直）。
用幂律拟合::

    T(ε) = T_ref · (ε / ε_ref)^n ,  ε = e / L

参数标定在 UIAA 101 坠落试验条件（80 kg, FF = 1.77）下完成，
见 ``scripts/calibrate_rope.py``。

已知简化
--------
1. 单轴垂直接触，无姿态、无岩壁碰撞。
2. 保护器滑动态下取 ``ε_hold`` 常数（等价于滑动态绳近似弹性）。
3. 绳索阻尼用单一粘性系数近似（真实绳是粘弹性的：加载/卸载路径不同）。
4. 不区分绳股、编织方向（各向异性）。
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
from scipy.integrate import solve_ivp

from . import G

__all__ = [
    "RopeModel",
    "RopeFallResult",
    "simulate_rope_fall",
    "default_rope",
]


# --------------------------------------------------------------------------
# 绳索本构
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class RopeModel:
    """动态绳的力-伸长本构（幂律硬化）。

    Attributes
    ----------
    t_ref_n
        参考张力 (N)。在 ``eps_ref`` 应变下的张力。
    eps_ref
        参考应变（无量纲，= 伸长率）。动态绳典型 0.25-0.35。
    exponent
        硬化指数。1.0 = 线性弹簧；>1 硬化。
    c_damp
        等效粘性阻尼 (N·s/m)。真实绳的滞后损耗 + 保护器摩擦的粗略等效。

    标定状态（2026-10-02）
    ---------------------
    UIAA 101 条件（80 kg, FF=1.77, L0=2 m）下本默认值给出
    **峰值 8.17 kN / 9.4 g / 伸长率 28.8%**，对比公开实测的
    8-9 kN / 30-40%。由 ``scripts/calibrate_rope.py`` 用加权 RSS
    目标函数（峰值 0.7 / 伸长 0.3）标定。

    ⚠️ **峰值力和伸长率没法同时对上，这是幂律本构的结构性局限**，
    不是标定没调好。标定网格（t_ref 4-12 kN × exponent 1.5-3.0 ×
    c_damp 0-350）里，峰值力随 ``t_ref`` 单调上升、伸长率随 ``t_ref``
    单调下降；要伸长率到 35% 就得把峰值压到 7 kN 附近，反之亦然。

    根因：真实动态绳有一条明显的**趾部**（编织层 + 护套先被拉紧，
    载荷-伸长曲线开头非常软，然后纱线芯才逐渐咬合），之后才硬化。
    幂律从零点就硬化，画不出"先软后硬"这个拐点，只能用指数去凑，
    于是标定出的 exponent=2.3 高于文献的动态绳量级（1.3-1.7）。

    **取舍**：这里优先保峰值力（UIAA 12 kN 判据直接作用在它上面，
    低估冲击力是安全分析里更危险的偏置方向），牺牲伸长率
    （28.8% vs 实测 30-40%）。

    修好它的正解是换成**趾部 + 硬化**的分段本构（见 Phase A-2），
    而不是继续调这四个数。旧的默认值 (8000/1.40/250) 给出的 5.58 kN
    偏软 30%，而且因为够不到任何保护器的 ``t_hold``，保护器对比表里
    ATC / Reverso / 8字结 四行的输出**完全相同**。
    """

    t_ref_n: float = 9_000.0
    eps_ref: float = 0.30
    exponent: float = 2.40
    c_damp: float = 100.0

    def tension(self, eps: float | np.ndarray) -> float | np.ndarray:
        """给定应变返回张力 (N)。``eps <= 0`` 时绳松弛，张力为 0。"""
        eps_arr = np.asarray(eps, dtype=float)
        out = np.where(
            eps_arr > 0.0,
            self.t_ref_n * np.power(np.maximum(eps_arr, 0.0) / self.eps_ref, self.exponent),
            0.0,
        )
        return float(out) if out.ndim == 0 else out

    def eps_at_tension(self, t: float) -> float:
        """给定张力返回应变（``tension`` 的反函数，单调）。

        张力必须 > 0。
        """
        if t <= 0.0:
            return 0.0
        return self.eps_ref * (t / self.t_ref_n) ** (1.0 / self.exponent)

    @property
    def stiffness(self) -> float:
        """初始刚度 (N/m)，用于对比线性弹簧近似。"""
        return self.t_ref_n / self.eps_ref / self.exponent

    def energy_per_length(self, eps: float) -> float:
        """单位绳长储存的弹性能 (J/m)，∫T de / L。"""
        if eps <= 0.0:
            return 0.0
        # ∫_0^eps T_ref (x/eps_ref)^n dx = T_ref eps_ref / (n+1) * (eps/eps_ref)^(n+1)
        return self.t_ref_n * self.eps_ref / (self.exponent + 1.0) * (eps / self.eps_ref) ** (self.exponent + 1.0)


def default_rope() -> RopeModel:
    """UIAA 标准条件（80 kg, FF=1.77, 单段 2 m）下峰值约 8-9 kN 的标定动态绳。"""
    return RopeModel()


# --------------------------------------------------------------------------
# 结果容器
# --------------------------------------------------------------------------
@dataclass
class RopeFallResult:
    """一次绳索坠落仿真的全部输出。"""

    # 输入回显
    mass_kg: float
    rope_length_m: float          # 保护器到坠落者的绳长 L0
    drop_height_m: float         # 释放点距保护器的高度 = FF * L0
    free_fall_m: float            # 绳绷紧前的自由落差 = max(0, FF-1) * L0
    fall_factor: float
    device: str
    t_hold_n: float

    # 时间序列（已按 dt 重采样）
    t_s: np.ndarray
    x_m: np.ndarray              # 向下位移
    v_ms: np.ndarray             # 向下速度
    accel_g: np.ndarray          # 加速度，正 = 向上（悬挂感）
    tension_n: np.ndarray        # 坠落者侧绳张力
    elongation_m: np.ndarray     # 自由绳段伸长
    free_length_m: np.ndarray    # 保护器到坠落者的自由绳长
    slip_m: np.ndarray           # 保护器累计放绳量

    # 标量指标
    peak_tension_n: float
    peak_accel_g: float
    max_drop_m: float            # 相对释放点的最大向下位移
    total_rope_slide_m: float    # 绳索相对岩壁/锚点的总滑动量 ≈ 落差 + 放绳
    max_elongation: float
    max_elongation_pct: float    # 最大伸长率 %（相对 L0）
    max_payout_m: float          # 保护器放出的绳量
    final_suspension_m: float    # 悬挂后相对释放点的高度差
    device_sat_time_s: float     # 保护器处于滑动态的累计时长
    hic: float = 0.0

    meta: dict = field(default_factory=dict)

    @property
    def peak_tension_kn(self) -> float:
        return self.peak_tension_n / 1e3

    @property
    def is_uiaa_pass(self) -> bool:
        """UIAA 101：单次坠落峰值力 < 12 kN。"""
        return self.peak_tension_n < 12_000.0


# --------------------------------------------------------------------------
# 求解器
# --------------------------------------------------------------------------
def simulate_rope_fall(
    *,
    mass_kg: float = 80.0,
    rope_length_m: float = 2.0,
    fall_factor: float = 1.77,
    t_hold_n: float = 1e9,
    device: str = "ATC-braked",
    rope: RopeModel | None = None,
    g: float = G,
    dt_sample: float = 0.0005,
    t_impact_max: float = 3.0,
    t_rebound: float = 0.6,
    max_ode_step: float = 2e-4,
    rtol: float = 1e-8,
    atol: float = 1e-10,
    t_peak_window_s: float = 0.15,
) -> RopeFallResult:
    """模拟一次绳索坠落。

    Parameters
    ----------
    mass_kg
        坠落者质量（含装备），kg。
    rope_length_m
        保护器（最后一个保护点）到坠落者的绳长 L0，m。
    fall_factor
        FF = 落差 / L0。0.3（低位）~ 1.77（UIAA 标准）~ 2.0（理论上限）。
    t_hold_n
        保护器可承受的**坠落者侧**最大张力 (N)，超过即放绳滑动。
        1e9 表示"等效锁死"（带刹停的 ATC 在常规坠落中不滑动）。
    device
        设备名，仅用于标注。
    rope
        绳索本构，默认用 ``default_rope()``。
    dt_sample
        结果重采样步长 (s)。默认 0.5 ms，远细于 50 ms 级的冲击过程。
    t_peak_window_s
        峰值统计窗口 (s)：从绳绷紧（首次 ``v=0``）前后各取一段。
        峰值力 / 峰值加速度只在这个窗口内取 max。

        **为什么必须限定窗口**：绳绷紧后把人拉停，随即回弹、绳再松弛、
        再绷紧……本模型的绳本构是纯弹性的（加载/卸载同一条曲线，没有
        滞回），阻尼项又只在拉伸时耗能，于是这个"自由振荡"几乎不衰减，
        3 s 内会反复重绷。**第二次重绷的瞬时冲击速度高于第一次**
        （第一次的重力势能被阻尼吃掉一部分），阻尼力 `c·ė` 又正比于
        伸长速率，于是重绷峰值可以超过首次冲击峰值。

        实测：标定后的较硬绳本构在 FF=1.5 下全程 max 给出 **17.3 kN**，
        而同一根绳 FF=1.77（UIAA 标准条件）只有 **10.8 kN** —— 峰值对
        FF 非单调，物理上不可能。限定窗口后 FF 单调性恢复。

    Returns
    -------
    RopeFallResult
    """
    rope = rope or default_rope()
    l0 = float(rope_length_m)
    h_total = fall_factor * l0        # 释放点距保护器的总高度
    # 几何：坠落者距保护器 h_total + x；保护器到坠落者的绳只有 l0。
    # 所以有 (FF-1)·L0 的绳是"松弛待用"的 —— 这段是**自由落体**。
    # FF <= 1 时绳在释放瞬间就绷紧，坠落完全由保护器摩擦承担。
    h_free = max(0.0, h_total - l0)
    eps_hold = rope.eps_at_tension(t_hold_n) if t_hold_n < 1e8 else 0.0

    # ---- 状态 y = [x, v, L_free] -----------------------------------------
    # 伸长 e = (x - h_free) + (L_free - L0)   （松弛时截断为 0）
    def elongation(y):
        return max(0.0, y[0] - h_free + (y[2] - l0))

    def tension_of(e, l_free):
        """给定伸长和自由绳长，返回 (绳张力, 保护器是否滑动)。"""
        if e <= 0.0 or l_free <= 0.0:
            return 0.0, False
        t_rope = rope.tension(e / l_free)
        return t_rope, t_rope > t_hold_n

    def rhs(_t, y):
        x, v, l_free = y
        e = elongation(y)
        t_rope, slipping = tension_of(e, l_free)

        xdot = v
        if slipping:
            # 滑动态：坠落者受到的是摩擦限幅张力 T_hold，与放绳速率无关
            t_on_body = t_hold_n
            # 放绳速率 L̇ = v / (1 + ε_hold)，见模块 docstring 推导
            ldot = v / (1.0 + eps_hold)
        else:
            t_on_body = t_rope
            ldot = 0.0

        # 阻尼：绳索阻尼器抵抗的是**伸长速率**，对坠落者而言张力方向朝上。
        #   ė > 0（拉伸中）→ 阻尼力朝上，叠加到张力（耗能）
        #   ė < 0（回缩中）→ 阻尼力朝下，抵消部分张力
        # 符号写反会变成"负阻尼"（加载时反而卸载），会向系统注入能量。
        edot = v - ldot
        t_damp = rope.c_damp * edot if e > 0.0 else 0.0
        t_total = max(0.0, t_on_body + t_damp)

        vdot = g - t_total / mass_kg
        return [xdot, vdot, ldot]

    # ---- 事件：第一次 v=0（冲击极点）--------------------------------------
    def v_zero(_t, y):
        return y[1]

    v_zero.terminal = False
    v_zero.direction = -1.0

    # ---- 单次积分（事件非终止，跑满 t_impact_max 后自然悬挂）--------------
    sol = solve_ivp(
        rhs,
        (0.0, t_impact_max),
        [0.0, 0.0, l0],
        events=[v_zero],
        max_step=max_ode_step,
        rtol=rtol,
        atol=atol,
        dense_output=True,
    )
    if not sol.success:
        raise RuntimeError(f"积分失败: {sol.message}")

    t_peak = float(sol.t_events[0][0]) if len(sol.t_events[0]) else float(sol.t[-1])

    # ---- 重采样到均匀时间栅格 --------------------------------------------
    dt_sample = min(dt_sample, max(1e-4, (t_impact_max - t_peak) / 2))
    t_out = np.arange(0.0, sol.t[-1], dt_sample)
    t_out = np.append(t_out, sol.t[-1])
    from scipy.interpolate import CubicSpline

    cs = CubicSpline(sol.t, sol.y, axis=1)
    x_out, v_out, l_out = cs(t_out)
    e_out = np.maximum(0.0, x_out - h_free + (l_out - l0))

    tension_out = np.where(e_out > 0.0, rope.tension(e_out / np.maximum(l_out, 1e-9)), 0.0)
    # 保护器限幅后的实际受力
    slip_now = tension_out > t_hold_n
    tension_body = np.where(slip_now, t_hold_n, tension_out)
    # 悬挂者感受到的向上减速加速度（g）
    a_out = (tension_body / mass_kg - g) / g

    # ---- 标量指标 --------------------------------------------------------
    # 峰值只在"首次绷紧"窗口内取 max，排除回弹重绷的伪峰（见 t_peak_window_s）
    win = (t_out >= t_peak - t_peak_window_s) & (t_out <= t_peak + t_peak_window_s)
    if not win.any():
        win = t_out <= t_peak + t_peak_window_s
    peak_tension = float(np.max(tension_body[win]))
    peak_accel = float(np.max(a_out[win]))
    max_drop = float(np.max(x_out))
    payout = float(np.max(l_out) - l0)
    slide = float(np.max(x_out)) + payout

    # 滑动态时长（保护器处于放绳状态的时间）—— 只统计首次冲击窗口，
    # 否则回弹期间的保护器反复滑/锁会被算进去
    sat = float(np.sum(slip_now[win]) * (t_out[1] - t_out[0])) if t_out.size > 1 else 0.0

    return RopeFallResult(
        mass_kg=mass_kg,
        rope_length_m=l0,
        drop_height_m=h_total,
        free_fall_m=h_free,
        fall_factor=fall_factor,
        device=device,
        t_hold_n=t_hold_n,
        t_s=t_out,
        x_m=x_out,
        v_ms=v_out,
        accel_g=a_out,
        tension_n=tension_body,
        elongation_m=e_out,
        free_length_m=l_out,
        slip_m=np.maximum(0.0, l_out - l0),
        peak_tension_n=peak_tension,
        peak_accel_g=peak_accel,
        max_drop_m=max_drop,
        total_rope_slide_m=slide,
        max_elongation=float(np.max(e_out)),
        max_elongation_pct=float(np.max(e_out) / l0 * 100.0),
        max_payout_m=payout,
        final_suspension_m=float(x_out[-1] - l0),
        device_sat_time_s=sat,
        hic=0.0,  # 由 metrics.hic() 填充
        meta={"t_peak_s": t_peak, "eps_hold": eps_hold, "rope": rope,
              "peak_window_s": (float(t_peak - t_peak_window_s), float(t_peak + t_peak_window_s))},
    )
