"""G4 降阶损伤本构：双线性弹塑 + σ(h) 后处理映射（纯 numpy，无 FE）。

预注册冻结设计：``docs/methods/G4_损伤演化设计.md`` §3（架构）/§4（本构）/§5（验收）。
本模块**只做后处理**：把冻结的严格线弹性 σ(h) 矩阵经双线性弹塑律映射为新的 σ(h)
矩阵，并给出 rank-2 的数据驱动 null（平行分析 / 断棍）。它**不调用** FEBio /
OpenSim、**不引入新依赖**（仅 numpy）、**不触碰任何冻结文件**。

双线性弹塑（降阶）::

    σ(ε) = E·ε                        , ε ≤ σ_y/E
    σ(ε) = σ_y + E_t·(ε − σ_y/E)      , ε > σ_y/E

等价于对线性弹性应力 ``σ_lin = E·ε`` 的逐点映射（后处理核心）::

    σ = σ_lin                         , σ_lin ≤ σ_y
    σ = σ_y + (E_t/E)·(σ_lin − σ_y)   , σ_lin > σ_y

一个 ``E_t`` 即跨越三种形状：**饱和**（E_t→0）/ **阶跃**（快速软化）/ **渐进**（缓降）。
逐 PID 的 σ_y、E_t 直接取自 THUMS 材料卡（``bone_feb.read_material_damage``，
卡布局 ``[mid, rho, E, nu, yield, tan_mod, …]``），**非调参**（设计 §6 D2/D3）。

诚实边界（设计 §4，必须随结果记录）
------------------------------------
* 皮质骨 ``E_t/E`` 实测 ≈ 0.82–0.90（如胫骨 81000700 = 14800/18000 = 0.822）
  ⇒ **只给膝拐（单调），非平台**。
* 许多壳/松质域卡 ``E_t = 0`` ⇒ **完美塑性平台**（饱和端）。
* ⇒ 长骨形状分化**可能温和**，rank-2 脱离 0.024% 但**不保证** clears 数据驱动 null。
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "BoneDamageLaw",
    "damage_map",
    "sigma_h_from_law",
    "classify_step_progressive",
    "broken_stick",
    "parallel_analysis_null",
]


@dataclass(frozen=True, slots=True)
class BoneDamageLaw:
    """双线性弹塑（降阶）本构：逐 PID 卡参数 ``E``、``σ_y``、``E_t``。"""

    E_mpa: float
    yield_stress_mpa: float
    tangent_modulus_mpa: float

    @property
    def yield_strain(self) -> float:
        """σ_y/E（屈服应变）。"""
        return float(self.yield_stress_mpa) / float(self.E_mpa)

    @property
    def tangent_ratio(self) -> float:
        """``E_t/E`` ∈ [0, ∞)（后屈服斜率 / 弹性斜率）。饱和端 = 0。"""
        return float(self.tangent_modulus_mpa) / float(self.E_mpa)

    def sigma(self, eps) -> np.ndarray:
        """σ(ε)：双线性弹塑（ε 可为标量 / 任意形状数组）。"""
        eps = np.asarray(eps, dtype=float)
        eps_y = self.yield_strain
        return np.where(
            eps <= eps_y,
            float(self.E_mpa) * eps,
            float(self.yield_stress_mpa)
            + float(self.tangent_modulus_mpa) * (eps - eps_y),
        )

    def sigma_from_linear_stress(self, sigma_lin) -> np.ndarray:
        """把严格线弹性应力 ``σ_lin = E·ε`` 映射为双线性 σ（后处理核心）。"""
        return damage_map(sigma_lin, self.yield_stress_mpa,
                          self.tangent_modulus_mpa, self.E_mpa)


def damage_map(sigma_lin, yield_stress_mpa, tangent_modulus_mpa, E_mpa) -> np.ndarray:
    """逐点双线性映射（numpy 广播）。

    ``σ_lin ≤ σ_y ⇒ σ = σ_lin``；否则 ``σ = σ_y + (E_t/E)·(σ_lin − σ_y)``。

    参数可按部位广播：例如把 ``(n_parts, 1)`` 的 σ_y/E_t/E 与 ``(n_parts, n_h)``
    的 σ_lin 矩阵相乘。**不改**入参。
    """
    x = np.asarray(sigma_lin, dtype=float)
    sy = np.asarray(yield_stress_mpa, dtype=float)
    ratio = np.asarray(tangent_modulus_mpa, dtype=float) / np.asarray(E_mpa, dtype=float)
    # σ_y 可为 +inf（"永不屈服" 极限）⇒ 塑性支产生 inf−inf；该支不会被选中，静默之。
    with np.errstate(invalid="ignore"):
        return np.where(x <= sy, x, sy + ratio * (x - sy))


def sigma_h_from_law(load_h, law: BoneDamageLaw, *, sigma_ref_mpa: float,
                     load_ref_n: float) -> np.ndarray:
    """由载荷 ``load(h)`` 与参考应力构造双线性 σ(h)。

    ``σ_lin(h) = σ_ref · load(h)/load_ref``（冻结的严格线弹性），再经 `law` 映射。
    """
    lin = float(sigma_ref_mpa) * np.asarray(load_h, dtype=float) / float(load_ref_n)
    return law.sigma_from_linear_stress(lin)


def classify_step_progressive(sigma_lin_h, law: BoneDamageLaw, *,
                              heights=None, plateau_ratio: float = 0.05,
                              gradual_ratio: float = 0.5) -> dict:
    """按**膝点锐度**把 σ(h) 曲线分类（设计 §5 V7 口径）。

    膝点锐度 = ``1 − E_t/E``（双线性律下即后屈服斜率相对弹性斜率的跌幅）。
    分类（阈值显式、可审计，非拟合）：

    * 区间内未达到 σ_y ⇒ ``elastic_no_knee``（纯弹性，无演化）；
    * ``E_t/E ≤ plateau_ratio`` ⇒ ``saturation``（E_t→0，完美塑性平台）；
    * ``plateau_ratio < E_t/E ≤ gradual_ratio`` ⇒ ``step``（快速软化）；
    * ``E_t/E > gradual_ratio`` ⇒ ``gradual``（缓降；皮质骨膝拐、单调）。
    """
    x = np.asarray(sigma_lin_h, dtype=float)
    gate = x > float(law.yield_stress_mpa)
    yielded = bool(np.any(gate))
    knee_h = None
    if yielded and heights is not None:
        knee_h = float(np.asarray(heights, dtype=float)[int(np.argmax(gate))])
    r = law.tangent_ratio
    if not yielded:
        cls = "elastic_no_knee"
    elif r <= plateau_ratio:
        cls = "saturation"
    elif r <= gradual_ratio:
        cls = "step"
    else:
        cls = "gradual"
    return {
        "class": cls,
        "tangent_ratio_Et_over_E": float(r),
        "knee_sharpness": float(1.0 - r),
        "yielded_in_range": yielded,
        "knee_h_m": knee_h,
        "max_sigma_lin_mpa": float(x.max()),
        "plateau_at_sigma_y_mpa": float(law.yield_stress_mpa),
    }


def broken_stick(p: int) -> np.ndarray:
    """断棍 null 的逐分量期望解释方差比（p 个变量）。第 k 分量 = (1/p)Σ_{j=k}^{p}1/j。"""
    p = int(p)
    return np.array([(1.0 / p) * float(np.sum(1.0 / np.arange(k, p + 1)))
                     for k in range(1, p + 1)], dtype=float)


def parallel_analysis_null(scales, *, n_samples: int, n_trials: int = 2000,
                           seed: int = 0, percentile: float = 95.0):
    """平行分析（Horn）null。

    生成与观测**同维、同方差、独立同分布**的正态数据（``scales`` = 各变量观测标准差），
    每 trial 取协方差特征值降序。返回 ``(pct, all)``：

    * ``pct``  : shape ``(p,)``，逐分量的 ``percentile`` 分位（默认 95）；
    * ``all``  : shape ``(n_trials, p)``，全部 trial 特征值（供 ratio 复算）。

    用固定 ``seed`` ⇒ 逐位可复现。
    """
    scales = np.asarray(scales, dtype=float)
    p = int(scales.size)
    rng = np.random.default_rng(int(seed))
    ev = np.empty((int(n_trials), p), dtype=float)
    for i in range(int(n_trials)):
        z = rng.standard_normal((int(n_samples), p)) * scales
        ev[i] = np.linalg.eigvalsh(np.cov(z, rowvar=False))[::-1]
    return np.percentile(ev, float(percentile), axis=0), ev
