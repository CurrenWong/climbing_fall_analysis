"""Phase A-1 —— 保护器（belay device）模型。

方案里的设计选择：*把保护器简化为「滑动摩擦 + 锁死」两态*。
本模块把这个两态具体化为一个参数：**坠落者侧保持张力** ``t_hold``。

物理依据
--------
保护器（和绳结）通过绳与挂点的**摩擦**来限制绳索的相对滑动。
静止时能保持的最大张力近似由欧拉绞盘公式给出::

    T_hold = T_bight · exp(μ · θ)

其中 ``T_bight`` 是挂片处"绳兜"的张力，``θ`` 是绳与挂片的等效包角，
``μ`` 是绳-金属/绳-绳摩擦系数。制动式保护器（ATC）额外叠加**带刹停的
保护者手力**，这一项通常远大于纯几何摩擦，因而在常规坠落中不滑动。

数据性质
--------
下表的 ``t_hold_n`` 与 ``mu`` / ``theta`` 都是**量级估计**，来源是
运动医学/绳索动力学的公开文献与厂商标称值区间，不是实测常数。
用途是驱动参数扫描与相对比较，不作为绝对预测。
模型可用 ``BelayDevice(...)`` 覆盖任意自定义取值。
"""

from __future__ import annotations

import math
from dataclasses import dataclass

__all__ = ["BelayDevice", "DEVICES", "get_device", "belayer_side_force"]


@dataclass(frozen=True)
class BelayDevice:
    """保护器参数集。"""

    name: str
    t_hold_n: float              # 坠落者侧保持张力 (N)
    mu: float = 0.20             # 绳-挂片摩擦系数（用于算保护者侧力）
    theta_rad: float = 1.20      # 保护器处的等效包角 (rad)
    slips_under_brake: bool = False   # 带刹停时是否仍可能滑动
    note: str = ""

    @property
    def e_friction(self) -> float:
        """欧拉绞盘放大系数 exp(μ·θ)。"""
        return math.exp(self.mu * self.theta_rad)

    def t_hold_kn(self) -> float:
        return self.t_hold_n / 1e3


# --------------------------------------------------------------------------
# 保护器库
# --------------------------------------------------------------------------
#: 保护器参数库。t_hold_n 为坠落者侧保持张力估计 (N)。
DEVICES: dict[str, BelayDevice] = {
    # 带刹停的传统保护：几何摩擦 + 保护者手力，常规坠落中几乎不滑动
    "ATC-braked": BelayDevice(
        "ATC + 带刹停", t_hold_n=15_000.0, mu=0.22, theta_rad=1.30,
        slips_under_brake=True,
        note="保护者手力主导；极限坠落或保护者被拉起时会放绳。",
    ),
    # 不刹停：绳直接从保护器里滑过，保护者只能靠握持
    "ATC-nobrake": BelayDevice(
        "ATC 无刹停", t_hold_n=6_500.0, mu=0.22, theta_rad=1.30,
        note="保护者未刹停；绳在挂片处大量滑动，绳磨损最快。",
    ),
    # 自锁辅助：靠凸轮 + 摩擦在低张力下放绳，把冲击压在较低水平
    "Grigri": BelayDevice(
        "Grigri", t_hold_n=4_000.0, mu=0.25, theta_rad=0.55,
        note="低张力放绳，坠落者侧峰值力最低；但放绳量大、绳在器内摩擦大。",
    ),
    "Reverso": BelayDevice(
        "Reverso", t_hold_n=7_000.0, mu=0.23, theta_rad=1.10,
        note="中等；锁止模式接近 ATC-braked，放绳模式介于 ATC 与 Grigri。",
    ),
    # 8 字结备份：被动时靠结的高摩擦缓降，主动锁止后几乎不滑动
    "Figure8-passive": BelayDevice(
        "8 字结（被动）", t_hold_n=9_000.0, mu=0.45, theta_rad=0.90,
        note="结的高摩擦提供缓慢降绳；未主动锁止。",
    ),
    "Figure8-locked": BelayDevice(
        "8 字结（锁止）", t_hold_n=1e9, mu=0.45, theta_rad=0.90,
        slips_under_brake=False,
        note="已锁止，等效刚性；此时 8 字结是纯备份，主保护仍由主系统承担。",
    ),
}


def get_device(name: str) -> BelayDevice:
    """按名称取保护器。找不到时抛 ``KeyError`` 并列出可用名称。"""
    if name not in DEVICES:
        raise KeyError(f"未知保护器 {name!r}；可选: {sorted(DEVICES)}")
    return DEVICES[name]


# --------------------------------------------------------------------------
# 保护者侧（锚点侧）张力
# --------------------------------------------------------------------------
def belayer_side_force(
    climber_tension_n: float | list[float] | tuple[float, ...],
    mu: float = 0.20,
    theta_rad: float = 0.60,
) -> float | list[float]:
    """把坠落者侧张力换算成保护者（锚点）侧张力。

    沿绳向上每经过一处摩擦（保护器、每一个快挂/主锁），张力按欧拉绞盘
    公式放大。**这正是"保护者被猛地拽起"的来源，也是多段攀岩中
    "最后一挂承担绝大部分力"的原因。**

    Parameters
    ----------
    climber_tension_n
        坠落者侧张力 (N)。
    mu
        绳-金属摩擦系数。干燥尼龙绳-铝合金 0.15-0.25。
    theta_rad
        系统总等效包角 (rad)。运动式攀登每个快挂约贡献 0.4-0.7 rad。
        - 顶绳/单段：0.2-0.4
        - 运动式两挂：0.8-1.4
        - 多段/传统排：2.0+

    Returns
    -------
    保护者侧张力 (N)，与输入同类型。
    """
    k = math.exp(mu * theta_rad)
    if isinstance(climber_tension_n, (list, tuple)):
        return [t * k for t in climber_tension_n]
    return climber_tension_n * k
