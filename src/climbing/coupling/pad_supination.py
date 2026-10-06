"""pad_supination.py -- R4: 垫子刚度 k ↔ 踝旋后(β) 耦合的两难模型。

为什么需要这一层
----------------
``pad_stiffness``（S1）只回答"软垫把力 **降** 给骨"；``ankle_supination``
（R1）只回答"已知 β → 侧向风险"。把两条独立路线串起来才能看清
**[Beurienne2025] Discussion** 描述的两难：

    "Increasing pad rigidity could prevent excessive ankle supination, thereby
    reducing force on ligament structures and lowering injury risk. However,
    this solution may be ineffective or even worse for high impact energy
    scenarios."

即：垫子越软 → F_peak 越低（**好**）但脚沉得越深 → 踝越内翻（**坏**）；
垫子越硬 → F_peak 越高（**坏**）但踝不内翻（**好**）。
本文档把这两条 **耦合** 起来，给出 (h, k) 平面上的腓骨远端 **总风险**，
并寻找 k* = argmin_total risk(h, k)。

设计原则（全部 opt-in）
------------------------
* **不动任何既有模块**——``pad_stiffness.r_fe`` / ``ankle_supination`` 的签名
  与返回值一律不变。
* **所有量都标注 measured / modeled / assumed**：
  - R_FE(k)：**measured**（FE 跖面 BC 探针，``pad_stiffness.json``）。
  - 腓骨 A / I / c：measured（``bc_robust_metric_route1.json`` +
    ``risk_1d_bending.json``）；主口径 = 等效圆。
  - σ_c = 70 MPa：assumed 用作弯曲极限（[Y25] ``fibula_ends`` **压缩**强度）。
  - 力臂 d = 30 mm：assumed（[Y25] S5 R1 已用）。
  - β(k) 几何：**modeled**（见下文 §β(k) 模型）。
* **β(k) 不当真理**：本模块只在几何上把"脚沉 → 踝内翻"做成可计算的量；
  β(k) 的绝对值、α(k) 系数都用灵敏度扫描给出区间。
* **总风险定义**（保守线性叠加）：

      risk_total = σ_aesisa / σ_c + σ_lat / σ_c
                 = (F_aesisa · R_FE(k)) / (A·σ_c)
                 + (M·c/I + F_lat/A) / σ_c

β(k) 模型（**modeIled**）
------------------------
[Heck2024] §5.2.3 / [Beurienne2025] Discussion：脚沉入垫 → 踝旋后 →
达垫子模量后载荷上升。把这句几何化成：

    δ(k) = F_peak(k) / k                              # mm，Winkler 凹陷深度
    β(k) = clip(β_0 + arctan(δ(k) / L_roll), 0, β_max)   # deg

    - δ(k)：把垫子看成 Winkler 弹簧（与 [S1] FE 探针同源），峰值力下脚的
      凹陷深度。
    - L_roll：从跟骨铰点到足外侧缘的横向半宽度（foot half-width）。
    - arctan(δ/L_roll)：足外缘相对跟骨铰点的几何倾角。
    - β_0：基线旋后角（**默认 0°**——几何模型仅建模"垫子引起的"额外旋转，
      不叠加自然不对称）。
    - β_max：几何饱和角（**默认 45°**——脚能翻到的上限）。

参数（中央 / 灵敏区间）：

| 参数 | 中央 | 低 | 高 | 来源 / 备注 |
|------|-----:|----:|----:|---|
| L_roll | **40 mm** | 25 mm | 60 mm | 足部横向半宽量级：跟骨 ~35 mm、前足 ~50 mm |
| β_0 | **0°** | 0° | 5° | 默认 0° = "无垫子 = 无几何旋后"基线 |
| β_max | 45° | 30° | 60° | 几何上限；扭伤 / 脱位前 |

灵敏度扫描见 :func:`sweep_roll_lever_arm` 与脚本 §3.2。

    总风险 = axial + lateral（**保守线性叠加**）
---------------------------------------------
    axial_risk  = axial_risk_rigid · R_FE(k)         # R1 纯轴向风险 × S1 折减比
                  （或备用：mult · F_peak(k) / (A·σ_c)）
    F_peak(k)   = F_aesisa_per_foot · R_FE(k)        # S1 转移函数（实测）
    F_lat(k)    = F_peak(k) · sin(β(k))              # β 重新指向横向
    M_inv(k)    = F_lat(k) · d                       # 距下轴力臂 assumed 30 mm
    σ_bend      = M_inv · c / I
    τ_shear     = F_lat / A
    σ_lat       = σ_bend + τ_shear                   # 主口径（线性叠加）
    lateral_risk = σ_lat / σ_c
    total_risk  = axial_risk + lateral_risk

⚠️ **诚实边界（同 §S5 R1）**：上界式 1D 筛选，绝对值不可信；只看相对与
"是否越阈"。临床 #1 踝伤是 **韧带扭伤（71%）**——本骨模型只能对应
**踝/外踝骨折（27%）**。

耦合退化
--------
当 β_0 = 0 且 k→∞ 时：
    - R_FE → 1.0、F_peak → F_aesisa、δ → 0、β → 0
    - lateral_risk → 0、total_risk → axial_risk_only（**与今天逐位一致**）。
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .ankle_supination import (
    DEFAULT_SUBTALAR_LEVER_ARM_MM,
    FIBULA_ENDS_SIGMA_C_MPA,
    fibula_lateral_risk,
    frontal_lateral_force_n,
    inversion_moment_nmm,
)

__all__ = [
    "BETA_MAX_DEG",
    "BETA_BASELINE_DEG_DEFAULT",
    "BETA_ROLL_LEVER_ARM_MM_DEFAULT",
    "CoupledPadResult",
    "pad_reduction_ratio",
    "load_fe_transfer_table",
    "indentation_depth_mm",
    "supination_angle_deg",
    "coupled_pad_risk",
    "coupled_pad_risk_at_grid",
    "sweep_roll_lever_arm",
    "find_optimal_stiffness_n_per_mm",
    "required_axial_reduction_for_soft_optimum",
]


# ---------------------------------------------------------------------------
# β(k) 几何模型参数（默认值；调用方可显式覆盖）
# ---------------------------------------------------------------------------
#: β 几何饱和角（deg）。脚可旋转的最大角度；超过此值已被扭伤 / 脱位拦截。
BETA_MAX_DEG: float = 45.0

#: 默认基线旋后角（deg）。中央取 **0**——只把"垫子引起"的旋后建模进来；
#: 不与自然落地的细微不对称叠加。可由 :func:`coupled_pad_risk` 覆盖。
BETA_BASELINE_DEG_DEFAULT: float = 0.0

#: 默认 roll 几何杠杆（mm）。足部外侧缘到跟骨铰点的横向半宽。
#: 中央 40 mm 在跟骨 ~35 mm 与前足 ~50 mm 之间。
BETA_ROLL_LEVER_ARM_MM_DEFAULT: float = 40.0

#: 灵敏度扫描默认杠杆档位（mm）。低 = 跟骨处（更陡），高 = 前足处（更缓）。
_LEVER_ARM_SENSITIVITY_MM: tuple[float, ...] = (25.0, 40.0, 60.0)


# ---------------------------------------------------------------------------
# S1 FE 实测折减比（measured）—— 只读复用 pad_stiffness.json
# ---------------------------------------------------------------------------
def load_fe_transfer_table(pad_stiffness_json: str | "Path") -> dict:
    """从 ``results/opensim_fe/pad_stiffness.json`` 读出 FE 转移函数表。

    返回 ``fe_transfer_table`` 字典（k_N_per_mm / gauge_max_mpa / ratio_vs_fixed /
    r_min / k_at_r_min_N_per_mm 等）。
    """
    from pathlib import Path

    p = Path(pad_stiffness_json)
    import json as _json_mod

    data = _json_mod.loads(p.read_text(encoding="utf-8"))
    return data["fe_transfer_table"]


def pad_reduction_ratio(
    k_n_per_mm: float | None,
    fe_table: dict,
) -> float:
    """``R_FE(k)``：垫子相对刚性地面的 **折减比**（**measured**，来自 S1 FE 探针）。

    在 ``log10(k)`` 上对 ``fe_table['k_N_per_mm']`` / ``ratio_vs_fixed`` 做线性插值；
    ``k >= max(k)`` 视为刚性极限 ``R = 1.0``；``k <= min(k)`` 平推软端最小 ``R``；
    ``k is None`` ⇒ ``R = 1.0``。

    与 :func:`scripts.opensim_fe.pad_stiffness.r_fe` 同语义；为避免脚本 ↔ src 耦合，
    在这里重写一次（输入只依赖 JSON 字典，不依赖脚本路径）。
    """
    if k_n_per_mm is None:
        return 1.0
    k = float(k_n_per_mm)
    ks = np.asarray(fe_table["k_N_per_mm"], dtype=float)
    rs = np.asarray(fe_table["ratio_vs_fixed"], dtype=float)
    if ks.size == 0 or rs.size == 0:
        raise ValueError("fe_table 必须含非空的 k_N_per_mm / ratio_vs_fixed。")
    if k >= ks[-1]:
        return 1.0
    if k <= ks[0]:
        return float(rs[0])
    return float(np.interp(math.log10(k), np.log10(ks), rs))


# ---------------------------------------------------------------------------
# β(k) 几何模型（modeIled）
# ---------------------------------------------------------------------------
def indentation_depth_mm(F_peak_n: float, k_n_per_mm: float) -> float:
    """Winkler 凹陷深度 ``δ = F_peak / k`` (mm)。

    当 ``k ≤ 0`` 或 ``F_peak ≤ 0`` 时返回 ``0.0``（无凹陷 → 无几何旋后）。
    """
    f = float(F_peak_n)
    k = float(k_n_per_mm)
    if k <= 0.0 or f <= 0.0:
        return 0.0
    return f / k


def supination_angle_deg(
    F_peak_n: float,
    k_n_per_mm: float,
    *,
    lever_arm_mm: float = BETA_ROLL_LEVER_ARM_MM_DEFAULT,
    beta_baseline_deg: float = BETA_BASELINE_DEG_DEFAULT,
    beta_max_deg: float = BETA_MAX_DEG,
) -> float:
    """**ModeIled** 旋后角 ``β(k)`` (deg)。

    几何：``β = clip(β_0 + arctan(δ / L_roll), 0, β_max)``，
    其中 ``δ = F_peak / k`` 为 Winkler 凹陷深度（mm）。

    边界：
      * ``k → ∞``（刚性）：``δ → 0``、``β → β_0``。
      * ``k → 0+``（极软）：``δ → ∞``、``β → β_max``（饱和）。
      * ``F_peak ≤ 0``：``δ = 0``、``β = β_0``（无峰值力 → 无旋后）。
    """
    delta = indentation_depth_mm(F_peak_n, k_n_per_mm)
    if delta <= 0.0:
        return float(beta_baseline_deg)
    if lever_arm_mm <= 0.0:
        raise ValueError(f"lever_arm_mm 必须为正，得到 {lever_arm_mm!r}")
    beta_pure = math.degrees(math.atan(delta / float(lever_arm_mm)))
    out = float(beta_baseline_deg) + beta_pure
    return float(min(max(out, 0.0), float(beta_max_deg)))


# ---------------------------------------------------------------------------
# 耦合：axial + lateral → total
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class CoupledPadResult:
    """一个 (h, k) 格的耦合解。

    量纲 mm-N-MPa-s；``*_mpa`` 与 ``*_risk`` 均为 scalar。
    """

    height_m: float
    k_n_per_mm: float
    R_fe: float
    F_peak_n: float              # 折减后的单足峰值（per-foot）
    indentation_mm: float
    beta_deg: float              # modeIled
    F_lat_n: float               # 横向力（N，|F_y·sinβ|）
    M_inv_nmm: float             # 内翻力矩（N·mm）
    sigma_bend_mpa: float
    sigma_shear_mpa: float
    sigma_lateral_mpa: float     # 主口径 = σ_bend + τ
    sigma_vonmises_mpa: float    # 交叉校核
    axial_risk: float            # 折减后 F_aesisa·R / (A·σ_c)
    lateral_risk: float          # σ_lat / σ_c
    total_risk: float            # axial + lateral（保守）
    fracture: bool               # total_risk >= 1

    # --- 几何 / σ_c / 力臂 / 模型来源 ------------------------------------
    A_section_mm2: float
    I_mm4: float
    c_mm: float
    sigma_c_mpa: float
    lever_arm_mm: float
    beta_baseline_deg: float
    lever_arm_geom_mm: float
    beta_max_deg: float

    def as_dict(self) -> dict:
        return {
            "height_m": self.height_m,
            "k_n_per_mm": self.k_n_per_mm,
            "R_fe": self.R_fe,
            "F_peak_n": self.F_peak_n,
            "indentation_mm": self.indentation_mm,
            "beta_deg": self.beta_deg,
            "F_lat_n": self.F_lat_n,
            "M_inv_nmm": self.M_inv_nmm,
            "M_inv_nm": self.M_inv_nm,
            "sigma_bend_mpa": self.sigma_bend_mpa,
            "sigma_shear_mpa": self.sigma_shear_mpa,
            "sigma_lateral_mpa": self.sigma_lateral_mpa,
            "sigma_vonmises_mpa": self.sigma_vonmises_mpa,
            "axial_risk": self.axial_risk,
            "lateral_risk": self.lateral_risk,
            "total_risk": self.total_risk,
            "fracture": self.fracture,
            "A_section_mm2": self.A_section_mm2,
            "I_mm4": self.I_mm4,
            "c_mm": self.c_mm,
            "sigma_c_mpa": self.sigma_c_mpa,
            "lever_arm_mm": self.lever_arm_mm,
            "beta_baseline_deg": self.beta_baseline_deg,
            "lever_arm_geom_mm": self.lever_arm_geom_mm,
            "beta_max_deg": self.beta_max_deg,
        }

    @property
    def M_inv_nm(self) -> float:
        return self.M_inv_nmm / 1000.0


def coupled_pad_risk(
    height_m: float,
    k_n_per_mm: float | None,
    *,
    F_peak_rigid_per_foot_n: float,
    fe_table: dict,
    A_section_mm2: float,
    I_mm4: float,
    c_mm: float,
    sigma_c_mpa: float = FIBULA_ENDS_SIGMA_C_MPA,
    axial_force_multiplier: float = 0.0,
    axial_risk_rigid: float | None = None,
    lever_arm_mm: float = DEFAULT_SUBTALAR_LEVER_ARM_MM,
    lever_arm_geom_mm: float = BETA_ROLL_LEVER_ARM_MM_DEFAULT,
    beta_baseline_deg: float = BETA_BASELINE_DEG_DEFAULT,
    beta_max_deg: float = BETA_MAX_DEG,
) -> CoupledPadResult:
    """计算一个 (h, k) 格的耦合解。

    Parameters
    ----------
    height_m
        跌落高度（m）。仅记录 / 用于报告。
    k_n_per_mm
        垫子 / 地面 FE Winkler 刚度（``None`` ⇒ 刚性极限 ``R = 1``）。
    F_peak_rigid_per_foot_n
        刚性地面下、**单足**峰值地面反力（N）。通常取
        ``max(ground_reaction.f_right_n)``。**整个模型的"载荷侧"基线**。
    fe_table
        :func:`load_fe_transfer_table` 返回的字典（**measured**，S1）。
    A_section_mm2, I_mm4, c_mm
        腓骨远端截面（等效圆 / 真实中段）；见 :mod:`ankle_supination`。
    sigma_c_mpa
        弯曲极限 (MPa)。默认 70 = [Y25] ``fibula_ends`` **压缩**强度，
        用作弯曲极限是近似（拉伸侧仅 30）。
    axial_force_multiplier
        **备用**轴向口径：``axial_risk = mult · F_peak(k) / (A·σ_c)``。
        默认 0 = 禁用。仅当 ``axial_risk_rigid is None`` 时生效。
    axial_risk_rigid
        **推荐**轴向口径：R1 在刚性地面下的纯轴向腓骨风险
        ``axial_risk_ends``（``nonvertical_s5.json`` 已算好的 OpenSim 口径）。
        若给出，则 ``axial_risk(k) = axial_risk_rigid · R_FE(k)``
        （与 S1 的"垫子按同一折减比折减整条链"假设一致）。
        ``k→∞`` 时 ``axial_risk = axial_risk_rigid``（**逐位复现 R1**）。
        优先级高于 ``axial_force_multiplier``。
    lever_arm_mm
        CoP→距下轴力臂 (mm)，assumed 30。``F_lat`` 经此力臂折算到关节轴。
    lever_arm_geom_mm
        β 几何杠杆 (mm)，脚外侧缘到跟骨铰点。中央 40，灵敏度 {25, 60}。
    beta_baseline_deg
        基线旋后角（deg）。默认 0 = 只建模"垫子引起"的旋后。
    beta_max_deg
        几何饱和角 (deg)。默认 45。

    Returns
    -------
    CoupledPadResult
        包含 ``F_peak`` / ``δ`` / ``β`` / ``F_lat`` / ``M_inv`` / ``σ_lat`` /
        ``axial_risk`` / ``lateral_risk`` / ``total_risk`` / ``fracture``。
    """
    # ---- S1 measured 折减比 → 折减后 F_peak -------------------------------
    R = pad_reduction_ratio(k_n_per_mm, fe_table)
    F_peak = float(F_peak_rigid_per_foot_n) * float(R)

    # ---- β(k) modeIled ----------------------------------------------------
    beta = supination_angle_deg(
        F_peak_n=F_peak,
        k_n_per_mm=k_n_per_mm if k_n_per_mm is not None else float("inf"),
        lever_arm_mm=lever_arm_geom_mm,
        beta_baseline_deg=beta_baseline_deg,
        beta_max_deg=beta_max_deg,
    )
    delta = indentation_depth_mm(F_peak, k_n_per_mm if k_n_per_mm is not None else float("inf"))

    # ---- R1 复用 supination_load / fibula_lateral_risk --------------------
    fvec = (0.0, F_peak, 0.0)            # 地面水平：β 把竖直载荷重新指向横向
    F_lat = abs(frontal_lateral_force_n(fvec, beta))
    M_inv = inversion_moment_nmm(F_lat, lever_arm_mm)
    lat = fibula_lateral_risk(
        F_lat, M_inv,
        a_section_mm2=float(A_section_mm2),
        i_mm4=float(I_mm4),
        c_mm=float(c_mm),
        sigma_c_mpa=float(sigma_c_mpa),
    )

    # ---- axial risk（按 S1 折减比一致缩放）--------------------------------
    A = float(A_section_mm2)
    sc = float(sigma_c_mpa)
    if axial_risk_rigid is not None:
        # 推荐口径：R1 的刚性纯轴向风险 × S1 折减比（k→∞ 时逐位复现 R1）
        axial = float(axial_risk_rigid) * float(R)
    elif axial_force_multiplier > 0.0 and A > 0.0 and sc > 0.0:
        axial = float(axial_force_multiplier) * F_peak / (A * sc)
    else:
        axial = 0.0

    total = axial + lat["risk"]
    return CoupledPadResult(
        height_m=float(height_m),
        k_n_per_mm=(None if k_n_per_mm is None else float(k_n_per_mm)),
        R_fe=float(R),
        F_peak_n=float(F_peak),
        indentation_mm=float(delta),
        beta_deg=float(beta),
        F_lat_n=float(F_lat),
        M_inv_nmm=float(M_inv),
        sigma_bend_mpa=float(lat["sigma_bend_mpa"]),
        sigma_shear_mpa=float(lat["sigma_shear_mpa"]),
        sigma_lateral_mpa=float(lat["sigma_lateral_mpa"]),
        sigma_vonmises_mpa=float(lat["sigma_vonmises_mpa"]),
        axial_risk=float(axial),
        lateral_risk=float(lat["risk"]),
        total_risk=float(total),
        fracture=bool(total >= 1.0),
        A_section_mm2=A,
        I_mm4=float(I_mm4),
        c_mm=float(c_mm),
        sigma_c_mpa=sc,
        lever_arm_mm=float(lever_arm_mm),
        beta_baseline_deg=float(beta_baseline_deg),
        lever_arm_geom_mm=float(lever_arm_geom_mm),
        beta_max_deg=float(beta_max_deg),
    )


# ---------------------------------------------------------------------------
# 扫描：(h, k) 网格 + 杠杆灵敏度
# ---------------------------------------------------------------------------
def coupled_pad_risk_at_grid(
    heights_m: Sequence[float],
    k_grid_n_per_mm: Sequence[float | None],
    *,
    F_peak_rigid_per_foot_by_h: dict[float, float],
    fe_table: dict,
    A_section_mm2: float,
    I_mm4: float,
    c_mm: float,
    sigma_c_mpa: float = FIBULA_ENDS_SIGMA_C_MPA,
    axial_force_multiplier: float = 0.0,
    axial_risk_rigid_by_h: dict[float, float] | None = None,
    lever_arm_mm: float = DEFAULT_SUBTALAR_LEVER_ARM_MM,
    lever_arm_geom_mm: float = BETA_ROLL_LEVER_ARM_MM_DEFAULT,
    beta_baseline_deg: float = BETA_BASELINE_DEG_DEFAULT,
    beta_max_deg: float = BETA_MAX_DEG,
) -> list[CoupledPadResult]:
    """在 ``(h, k)`` 矩形网格上调用 :func:`coupled_pad_risk`，返回平铺结果列表。

    ``axial_risk_rigid_by_h``（可选）：``{h: R1 纯轴向 risk}``，逐高度传入
    :func:`coupled_pad_risk` 的 ``axial_risk_rigid``。
    """
    out: list[CoupledPadResult] = []
    for h in heights_m:
        F_rigid = float(F_peak_rigid_per_foot_by_h[float(h)])
        ax_rigid = None
        if axial_risk_rigid_by_h is not None:
            ax_rigid = float(axial_risk_rigid_by_h[float(h)])
        for k in k_grid_n_per_mm:
            out.append(
                coupled_pad_risk(
                    height_m=float(h),
                    k_n_per_mm=k,
                    F_peak_rigid_per_foot_n=F_rigid,
                    fe_table=fe_table,
                    A_section_mm2=A_section_mm2,
                    I_mm4=I_mm4,
                    c_mm=c_mm,
                    sigma_c_mpa=sigma_c_mpa,
                    axial_force_multiplier=axial_force_multiplier,
                    axial_risk_rigid=ax_rigid,
                    lever_arm_mm=lever_arm_mm,
                    lever_arm_geom_mm=lever_arm_geom_mm,
                    beta_baseline_deg=beta_baseline_deg,
                    beta_max_deg=beta_max_deg,
                )
            )
    return out


def sweep_roll_lever_arm(
    heights_m: Sequence[float],
    k_grid_n_per_mm: Sequence[float | None],
    *,
    F_peak_rigid_per_foot_by_h: dict[float, float],
    fe_table: dict,
    A_section_mm2: float,
    I_mm4: float,
    c_mm: float,
    sigma_c_mpa: float = FIBULA_ENDS_SIGMA_C_MPA,
    axial_force_multiplier: float = 0.0,
    axial_risk_rigid_by_h: dict[float, float] | None = None,
    lever_arm_mm: float = DEFAULT_SUBTALAR_LEVER_ARM_MM,
    lever_arm_geom_grid_mm: Sequence[float] = _LEVER_ARM_SENSITIVITY_MM,
    beta_baseline_deg: float = BETA_BASELINE_DEG_DEFAULT,
    beta_max_deg: float = BETA_MAX_DEG,
) -> dict[float, list[CoupledPadResult]]:
    """对 ``lever_arm_geom`` 做灵敏度扫描（默认 ``(25, 40, 60)`` mm）。

    返回 ``{lever_arm_geom_mm: [CoupledPadResult ...]}``，键 = ``float`` 杠杆值。
    """
    out: dict[float, list[CoupledPadResult]] = {}
    for L in lever_arm_geom_grid_mm:
        out[float(L)] = coupled_pad_risk_at_grid(
            heights_m,
            k_grid_n_per_mm,
            F_peak_rigid_per_foot_by_h=F_peak_rigid_per_foot_by_h,
            fe_table=fe_table,
            A_section_mm2=A_section_mm2,
            I_mm4=I_mm4,
            c_mm=c_mm,
            sigma_c_mpa=sigma_c_mpa,
            axial_force_multiplier=axial_force_multiplier,
            axial_risk_rigid_by_h=axial_risk_rigid_by_h,
            lever_arm_mm=lever_arm_mm,
            lever_arm_geom_mm=float(L),
            beta_baseline_deg=beta_baseline_deg,
            beta_max_deg=beta_max_deg,
        )
    return out


# ---------------------------------------------------------------------------
# k* 寻找（argmin total_risk）—— 处理"单调→边界极值"的情形
# ---------------------------------------------------------------------------
def find_optimal_stiffness_n_per_mm(
    cells: Sequence[CoupledPadResult],
    *,
    prefer: str = "interior_or_boundary",
) -> dict:
    """在一组同 ``height_m`` 的 :class:`CoupledPadResult` 中寻找 k* = argmin total_risk。

    Parameters
    ----------
    cells
        同一高度的耦合结果列表（按 ``k`` 排序；``None`` 表示刚性）。
    prefer
        * ``"interior"`` —— 只在内部格（``k`` 两端有限）找最小；若等值则取较软。
        * ``"interior_or_boundary"``（默认）—— 在 **所有** 格上找最小；
          内部存在显式极小 → 返回该 k；否则取边界（``k → ∞`` / ``k → 0``）中
          的较小者。

    Returns
    -------
    dict
        ``{"k_star_n_per_mm", "total_risk_star", "interior_min", "boundary_min",
           "all_totals": [...]}``。``interior_min`` / ``boundary_min`` 各自给出
        ``k``、``total_risk``、``fracture`` 与解释文本。
    """
    if not cells:
        raise ValueError("cells 不能为空")

    # ---- 排序（None = 刚性末端最大） -------------------------------------
    def _k_key(c: CoupledPadResult) -> float:
        return 1e30 if c.k_n_per_mm is None else float(c.k_n_per_mm)

    sorted_cells = sorted(cells, key=_k_key)

    # ---- 内部极小（两端都是有限 k） --------------------------------------
    interior: list[CoupledPadResult] = []
    for i in range(1, len(sorted_cells) - 1):
        prev_t = sorted_cells[i - 1].total_risk
        cur_t = sorted_cells[i].total_risk
        nxt_t = sorted_cells[i + 1].total_risk
        if cur_t <= prev_t and cur_t <= nxt_t:
            interior.append(sorted_cells[i])

    interior_min = None
    if interior:
        best = min(interior, key=lambda c: c.total_risk)
        interior_min = {
            "k_n_per_mm": best.k_n_per_mm,
            "total_risk": best.total_risk,
            "fracture": best.fracture,
            "explanation": (
                "内部极小：total_risk 在该 k 处对两侧同时不增；"
                f"k={best.k_n_per_mm} N/mm, total={best.total_risk:.3f}"
            ),
        }

    # ---- 边界极小 ---------------------------------------------------------
    soft_end = sorted_cells[0]
    rigid_end = sorted_cells[-1]
    if soft_end.total_risk <= rigid_end.total_risk:
        boundary_min = {
            "k_n_per_mm": soft_end.k_n_per_mm,
            "total_risk": soft_end.total_risk,
            "fracture": soft_end.fracture,
            "explanation": (
                "软端边界（k→0）极小：total_risk 从软端到硬端单调递增。"
                f"k={soft_end.k_n_per_mm} N/mm, total={soft_end.total_risk:.3f}"
            ),
        }
    else:
        boundary_min = {
            "k_n_per_mm": rigid_end.k_n_per_mm,
            "total_risk": rigid_end.total_risk,
            "fracture": rigid_end.fracture,
            "explanation": (
                "硬端边界（k→∞）极小：total_risk 从软端到硬端单调递减。"
                f"k={rigid_end.k_n_per_mm} N/mm, total={rigid_end.total_risk:.3f}"
            ),
        }

    # ---- 报告 -------------------------------------------------------------
    if interior_min is None:
        # 无内部极小 → 选边界极小
        k_star = boundary_min["k_n_per_mm"]
        total_star = boundary_min["total_risk"]
        verdict_kind = "boundary_min"
        verdict_text = (
            "无内部极小：total_risk 在 k 网格上单调；k* 落在边界（"
            f"{'硬端 (rigid)' if boundary_min['k_n_per_mm'] in (None, 1e30) else '软端'}）。"
            f" {boundary_min['explanation']}"
        )
    else:
        if prefer == "interior":
            k_star = interior_min["k_n_per_mm"]
            total_star = interior_min["total_risk"]
            verdict_kind = "interior_min"
            verdict_text = interior_min["explanation"]
        else:
            # interior_or_boundary：若 internal 比 boundary 更小就用 internal；
            # 否则用 boundary（典型情况：单调时 ∈ interior 实际是"局部平段"）。
            if (interior_min["total_risk"] < boundary_min["total_risk"]
                    or abs(interior_min["total_risk"] - boundary_min["total_risk"]) < 1e-9):
                k_star = interior_min["k_n_per_mm"]
                total_star = interior_min["total_risk"]
                verdict_kind = "interior_min"
                verdict_text = interior_min["explanation"]
            else:
                k_star = boundary_min["k_n_per_mm"]
                total_star = boundary_min["total_risk"]
                verdict_kind = "boundary_min"
                verdict_text = boundary_min["explanation"]

    return {
        "k_star_n_per_mm": k_star,
        "total_risk_star": float(total_star),
        "verdict_kind": verdict_kind,
        "verdict_text": verdict_text,
        "interior_min": interior_min,
        "boundary_min": boundary_min,
        "all_totals": [
            {
                "k_n_per_mm": c.k_n_per_mm,
                "total_risk": c.total_risk,
                "lateral_risk": c.lateral_risk,
                "axial_risk": c.axial_risk,
                "fracture": c.fracture,
            }
            for c in sorted_cells
        ],
    }


# ---------------------------------------------------------------------------
# 内部极小存在性的定量门槛（"软端能否胜过硬端"）
# ---------------------------------------------------------------------------
def required_axial_reduction_for_soft_optimum(
    axial_risk_rigid: float,
    lateral_risk_at_soft_k: float,
) -> float:
    """使"软端 k"成为最优，所需的轴向折减比例 ``(1 − R)``。

    推导（刚性端侧向 = 0）：

        total(soft) < total(rigid)
        ⟺  axial_rigid·R + lateral_soft < axial_rigid
        ⟺  1 − R > lateral_soft / axial_rigid

    返回值 ``> 1`` 表示**物理上不可能**（等价于要求负的 R）。用实测 R_FE 的
    最大折减比例 ``1 − R_min`` 与之比较，即可定量判断软端能否胜出。
    """
    a = float(axial_risk_rigid)
    lat = float(lateral_risk_at_soft_k)
    if a <= 0.0:
        return float("inf")
    return lat / a