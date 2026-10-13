"""G6 · 载荷段软组织串联柔顺 —— 旁路段定义（**真根因修复**）。

冻结设计：``docs/methods/G6_载荷段软组织串联柔顺_预注册.md``。

缺陷
----
``chain2d.SEGMENTS_6`` 的 ``k_axial_n_per_m`` 取**纯骨柱值**（如小腿段
``1.3e7 ≡ EA/L = 17e9×3.2e-4/0.426``），**零软组织串联柔顺**，也不含膝屈曲等
制动分散机制。实测后果：σ_胫骨 **631.5 / 1013.9 MPa = 皮质骨抗压强度
（130–180 MPa）的 3.5–5.6×**——**模型不自洽**。

修复
----
逐段加**软组织串联元件**（而非整体缩放 ``k_axial``，以保留「远端骨主导、
近端肌主导」的解剖分层）::

    k_soft,i = E_st · A_st,i / L_i        （A_st,i = ALPHA_ST · A_cort,i）
    k_eff,i  = 1 / (1/k_bone,i + 1/k_soft,i)

**由公式推导，不手抄数值**（预注册表曾出现 `足` 行算术错 10 倍，本函数杜绝该类错误）。

纪律
----
* ``SEGMENTS_6`` 是跑前冻结的解剖参数 ⇒ **本模块一字不改它**，只提供旁路元组
  ``SEGMENTS_6_SOFT``。默认路径仍走 ``SEGMENTS_6`` ⇒ 既有结果**全部可复现**。
* 锚点是 §2 的 **A1 整链等效轴向刚度 30–170 kN/m** 与 **A2 步行胫骨骨应变
  237–1250 µε**，两者均**独立可测、与 V8b 无关**。
* 判据失败时**先回文献/解剖重查 ``ALPHA_ST`` / ``E_ST_PA``**，**不得反向调参凑区间**。
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from .chain2d import SEGMENTS_6, Segment

__all__ = [
    "E_ST_PA",
    "ALPHA_ST",
    "SEGMENTS_6_SOFT",
    "segments_soft",
    "chain_effective_stiffness_n_per_m",
    "arrest_stiffness_n_per_m",
]

#: 软组织模量 (Pa)。来源：PMC7856106（胫腓关节/生长板等软组织取 10 MPa）。
#: **冻结值**（预注册 §3）。
E_ST_PA: float = 10.0e6

#: 软组织横截面积 / 皮质承载面积 之比。**解剖估计值，非文献实测**（预注册 §3 已声明）。
#: **冻结值**；判据失败时先回解剖重估，不得反向凑区间。
ALPHA_ST: float = 3.0


def _k_series(k_bone: float, k_soft: float) -> float:
    """两弹簧串联等效刚度。"""
    return 1.0 / (1.0 / k_bone + 1.0 / k_soft)


def segments_soft(
    segments: tuple[Segment, ...] = SEGMENTS_6,
    E_st_pa: float = E_ST_PA,
    alpha: float = ALPHA_ST,
) -> tuple[Segment, ...]:
    """逐段加软组织串联柔顺后返回新段元组（不改原 ``Segment``）。

    ``tissue_area_m2`` / ``tissue_modulus_pa`` **保持不变**——它们是**观测量的分母**
    （骨应变 ``ε=(F/A_cort)/E_bone``），改了就不是同一个观测量了。
    只改 ``k_axial_n_per_m``。
    """
    out = []
    for s in segments:
        a_st = alpha * s.tissue_area_m2
        k_soft = E_st_pa * a_st / s.length_m
        out.append(replace(s, k_axial_n_per_m=_k_series(s.k_axial_n_per_m, k_soft)))
    return tuple(out)


#: 旁路段定义（G6 默认入口）。既有路径仍用 ``SEGMENTS_6``。
SEGMENTS_6_SOFT: tuple[Segment, ...] = segments_soft()


def chain_effective_stiffness_n_per_m(
    segments: tuple[Segment, ...] = SEGMENTS_6_SOFT,
) -> float:
    """整链等效轴向刚度 (N/m) = 各关节刚度的**串联**。

    与 ``chain2d._build_ctx`` 的 ``k_ax`` 同一算法（相邻段刚度串联后再串联各关节），
    故与模拟器内部实际用的轴向刚度一致。
    """
    k_seg = np.array([s.k_axial_n_per_m for s in segments])
    k_joint = 1.0 / (1.0 / k_seg[:-1] + 1.0 / k_seg[1:])
    return float(1.0 / np.sum(1.0 / k_joint))


def arrest_stiffness_n_per_m(
    mass_kg: float, t_arrest_s: float
) -> float:
    """由制动半周期反推的等效刚度 (N/m)：``k = m·(π/(2·T_arrest))²``。

    **判据 1 的观测量定义**（预注册 §4，跑前冻结）：``T_arrest`` = 首次接触时刻
    到**节点 0 速度首次过零**（动量反转）的时长。半正弦制动假设下
    ``ω = π/(2T)`` ⇒ ``k = mω²``。
    """
    if t_arrest_s <= 0.0:
        return float("nan")
    omega = np.pi / (2.0 * t_arrest_s)
    return float(mass_kg * omega * omega)
