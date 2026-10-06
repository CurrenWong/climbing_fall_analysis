"""R6 —— 踝韧带扭伤（sprain）判据（**opt-in，默认零**）。

为什么需要这一层
----------------
现有 OpenSim→FEBio 链路 + Phase-S5 ``ankle_supination`` 只把踝载荷视为**骨**
承载（腓骨两端压缩 / 远端弯曲 / 剪切）。但临床**#1 抱石踝伤是韧带扭伤**：

    [Heck2024 / Müller2022] 急诊 430 例 / 447 伤：
    - 踝伤：扭伤 **71.3%** vs 骨折 27.4%。
    [Beurienne2025] 问卷 245 人 / 301 伤：
    - 踝扭伤 #1（28%）、踝骨折 8%。

我们的骨模型只能覆盖踝骨骨折（27%）——扭伤这条主线完全留白。本模块把
**踝韧带扭伤判据**作为 R1 旋后/内翻载荷模型的下游判据补上。

模型（全部标注为 **measured / modeled / assumed**）
--------------------------------------------------
输入：``M_inv`` (N·mm) = 旋后力矩（来自 ``ankle_supination.inversion_moment_nmm``）。

1. **关节级判据（主口径，避开了 r_lig 力臂假设）**

       risk = M_inv / M_fail
       sprain = risk ≥ 1.0

   * ``M_fail`` 来自尸体文献：
     - 无预载 ATFL: **21–34 N·m** (Funk 2002 = 21 ± 5 @ 38 ± 8°;
       Parenteau 1998 = 34.1 ± 14.5 @ 34.3 ± 7.5°)
     - 2 kN 轴向预载: **77 ± 27 N·m** @ 40 ± 12° (Funk 2002)
     - 临床实测反向支持：Bahr 1998 报告 ATFL 步态期峰值力
       **76 ± 23 N** ≈ 极限载荷的一半 → 安全系数 ≈ 2。

2. **预载效应外推（两口径，避免单阈值错觉）**

   文献数据点：
   - (0 N, 21 N·m)        — Funk 2002 无预载均值
   - (2000 N, 77 N·m)     — Funk 2002 2 kN 预载均值

   ``preload_capacity_interp(P)`` 在 [0, 2 kN] 上做线性插值；
   ``P > 2 kN`` 的外推用 ``extrapolation="clamp"`` 钳到 77 N·m
   或 ``extrapolation="linear"`` 继续外推（**UNCERTAIN**）：
   抱石冲击轴向 ≈ 21 kN 是 2 kN 的 10.5×，远超文献外推范围。

3. **力 / 应变判据（r_lig assumed）**

       F_lig = M_inv / r_lig                (N)
       ε_lig ≈ F_lig / (k · L_lig)          (linear, secant stiffness)
       sprain = (F_lig / F_fail) ≥ 1 OR (ε_lig / ε_fail) ≥ 1

   * ``r_lig`` 是韧带力臂（**assumed**，不在 OpenSim 里）；中央 22 mm
     （接近 ATFL 静息长度 18–25 mm 的中点），灵敏度 {20, 25} mm。
   * 极限参数 (ATFL)：F_fail=200 N, ε_fail=14%, k=14 N/mm, L=22 mm,
     CSA=13, Bahr 1998 步态期峰值力 76 ± 23 N。
   * 极限参数 (CFL)：F_fail=340 N, ε_fail=13%, k=12 N/mm, L=28.5 mm,
     CSA=10。

4. **[Tochigi2006] 参考带**（**不是失效阈值**）
   - 正常步态 ATFL 应变 ≤ 6.2 %、CFL ≤ 2.1 %。
   - "taut" 阈值 = 1.0 % 应变。
   - 用作"非生理"参照线与"开始受力"下限，**不当作失效判据**。

单位 mm–N–MPa–s。
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

__all__ = [
    # ATFL 参数
    "ATFL_ULTIMATE_LOAD_N",
    "ATFL_FAILURE_STRAIN",
    "ATFL_FAILURE_STRAIN_LOW",
    "ATFL_FAILURE_STRAIN_HIGH",
    "ATFL_STIFFNESS_N_PER_MM",
    "ATFL_RESTING_LENGTH_MM",
    "ATFL_RESTING_LENGTH_LOW_MM",
    "ATFL_RESTING_LENGTH_HIGH_MM",
    "ATFL_CSA_MM2",
    "ATFL_E_MPA",
    "ATFL_PEAK_FORCE_STANCE_N",
    "ATFL_PEAK_FORCE_STANCE_STD_N",
    # CFL 参数
    "CFL_ULTIMATE_LOAD_N",
    "CFL_FAILURE_STRAIN",
    "CFL_FAILURE_STRAIN_LOW",
    "CFL_FAILURE_STRAIN_HIGH",
    "CFL_STIFFNESS_N_PER_MM",
    "CFL_RESTING_LENGTH_MM",
    "CFL_RESTING_LENGTH_LOW_MM",
    "CFL_RESTING_LENGTH_HIGH_MM",
    "CFL_CSA_MM2",
    "CFL_E_MPA",
    # 关节级 M_fail（无预载 / 2 kN 预载）
    "ATFL_INVERSION_FAIL_MOMENT_NO_PRELOAD_NM",
    "ATFL_INVERSION_FAIL_MOMENT_NO_PRELOAD_HIGH_NM",
    "ATFL_INVERSION_FAIL_MOMENT_2KN_PRELOAD_NM",
    "ATFL_INVERSION_FAIL_ANGLE_DEG",
    "ATFL_INVERSION_FAIL_ANGLE_HIGH_DEG",
    "ATFL_INVERSION_FAIL_ANGLE_2KN_DEG",
    "CFL_INVERSION_FAIL_ANGLE_DEG",
    # Tochigi 参考带
    "TOCHIGI_ATFL_NORMAL_STRAIN_MAX",
    "TOCHIGI_CFL_NORMAL_STRAIN_MAX",
    "TOCHIGI_TAUT_STRAIN",
    # 默认 r_lig
    "DEFAULT_LIGAMENT_LEVER_ARM_MM",
    "R_LIG_SENSITIVITY_MM",
    # API
    "LigamentCriterion",
    "preload_capacity_interp",
    "atfl_sprain_risk",
    "cfl_sprain_risk",
    "ligament_force_n",
    "ligament_strain",
    "sprain_risk",
    "critical_inversion_angle_deg",
    "critical_inversion_moment_nmm",
    "ankle_ligament_criterion",
    "ankle_ligament_criterion_at_grid",
]


# ---------------------------------------------------------------------------
# ATFL 极限 / 力学参数（**measured**，来源见 [L1]）
# ---------------------------------------------------------------------------
#: ATFL 极限载荷 (N) —— 均值 (St.Pierre 1983; Chapman 2019; Siegler 1988)
ATFL_ULTIMATE_LOAD_N: float = 200.0
#: ATFL 失效应变（线性，secant）—— 均值 14 % (11–15 %，Siegler 1988)
ATFL_FAILURE_STRAIN: float = 0.14
ATFL_FAILURE_STRAIN_LOW: float = 0.11
ATFL_FAILURE_STRAIN_HIGH: float = 0.15
#: ATFL 割线刚度 (N/mm)
ATFL_STIFFNESS_N_PER_MM: float = 14.0
#: ATFL 静息长度 (mm) —— 区间中点
ATFL_RESTING_LENGTH_MM: float = 22.0
ATFL_RESTING_LENGTH_LOW_MM: float = 18.0
ATFL_RESTING_LENGTH_HIGH_MM: float = 25.0
#: ATFL 截面积 (mm²)
ATFL_CSA_MM2: float = 13.0
#: ATFL 弹性模量 (MPa)
ATFL_E_MPA: float = 260.0
#: 步态期 ATFL 峰值力 (N) —— Bahr 1998 步态测量（safety factor ≈ 2）
ATFL_PEAK_FORCE_STANCE_N: float = 76.0
ATFL_PEAK_FORCE_STANCE_STD_N: float = 23.0

# ---------------------------------------------------------------------------
# CFL 极限 / 力学参数
# ---------------------------------------------------------------------------
#: CFL 极限载荷 (N)
CFL_ULTIMATE_LOAD_N: float = 340.0
#: CFL 失效应变
CFL_FAILURE_STRAIN: float = 0.13
CFL_FAILURE_STRAIN_LOW: float = 0.09
CFL_FAILURE_STRAIN_HIGH: float = 0.13
#: CFL 割线刚度 (N/mm)
CFL_STIFFNESS_N_PER_MM: float = 12.0
#: CFL 静息长度 (mm)
CFL_RESTING_LENGTH_MM: float = 28.5
CFL_RESTING_LENGTH_LOW_MM: float = 25.0
CFL_RESTING_LENGTH_HIGH_MM: float = 32.0
#: CFL 截面积 (mm²)
CFL_CSA_MM2: float = 10.0
#: CFL 弹性模量 (MPa)
CFL_E_MPA: float = 510.0

# ---------------------------------------------------------------------------
# 关节级内翻失效力矩（**measured**，Funk 2002 + Parenteau 1998）
# ---------------------------------------------------------------------------
#: 无预载 ATFL 内翻失效力矩 (N·m) —— Funk 2002 均值 (21 ± 5)
ATFL_INVERSION_FAIL_MOMENT_NO_PRELOAD_NM: float = 21.0
#: 无预载 ATFL 上限参考（Parenteau 1998 均值）—— 给"区间"上沿
ATFL_INVERSION_FAIL_MOMENT_NO_PRELOAD_HIGH_NM: float = 34.1
#: 2 kN 轴向预载 ATFL 内翻失效力矩 (N·m) —— Funk 2002 均值 (77 ± 27)
ATFL_INVERSION_FAIL_MOMENT_2KN_PRELOAD_NM: float = 77.0
#: 无预载 ATFL 内翻失效角 (deg) —— Funk 2002 均值
ATFL_INVERSION_FAIL_ANGLE_DEG: float = 38.0
ATFL_INVERSION_FAIL_ANGLE_HIGH_DEG: float = 34.3  # Parenteau 1998 均值
#: 2 kN 预载 ATFL 内翻失效角 (deg) —— Funk 2002 均值
ATFL_INVERSION_FAIL_ANGLE_2KN_DEG: float = 40.0
#: CFL 在组合载荷下失效角 (deg) —— 晚于 ATFL
CFL_INVERSION_FAIL_ANGLE_DEG: float = 27.0

# ---------------------------------------------------------------------------
# Tochigi 2006 正常步态参考带（**measured**，不是失效阈值）
# ---------------------------------------------------------------------------
#: Tochigi 2006 Table 1：ATFL 步态应变上限（dynamic loading）= 6.2 %
TOCHIGI_ATFL_NORMAL_STRAIN_MAX: float = 0.062
#: Tochigi 2006 Table 1：CFL 步态应变上限 = 2.1 %
TOCHIGI_CFL_NORMAL_STRAIN_MAX: float = 0.021
#: Tochigi 2006 "taut" 阈值 = 1.0 %
TOCHIGI_TAUT_STRAIN: float = 0.010

# ---------------------------------------------------------------------------
# r_lig（韧带力臂，**assumed**）—— 默认 / 灵敏度
# ---------------------------------------------------------------------------
#: 默认 r_lig (mm) —— ATFL 静息长度中点附近
DEFAULT_LIGAMENT_LEVER_ARM_MM: float = 22.0
#: r_lig 灵敏度档位 (mm)
R_LIG_SENSITIVITY_MM: tuple[float, ...] = (20.0, 25.0)


# ---------------------------------------------------------------------------
# 关节级预载插值
# ---------------------------------------------------------------------------
def preload_capacity_interp(
    axial_preload_n: float,
    *,
    m_fail_no_preload_nm: float = ATFL_INVERSION_FAIL_MOMENT_NO_PRELOAD_NM,
    m_fail_2kn_nm: float = ATFL_INVERSION_FAIL_MOMENT_2KN_PRELOAD_NM,
    extrapolation: Literal["clamp", "linear"] = "clamp",
) -> float:
    """轴向预载 P → 内翻承载力 ``M_fail(P)`` (N·m)。

    锚点（[Funk 2002]）：
        P = 0       → 21 N·m        （无预载）
        P = 2000 N  → 77 N·m        （2 kN 预载）

    在 ``[0, 2000] N`` 上做线性插值。``P > 2000`` 时按 ``extrapolation`` 处理：

    * ``"clamp"`` (默认)  → 钳到 ``m_fail_2kn_nm``。**保守**。
    * ``"linear"``        → 沿同一斜率外推。**外推不确定**（抱石冲击 ≈ 21 kN
      是 2 kN 的 10.5×，外推远超出文献范围）。

    Returns
    -------
    float
        M_fail (N·m)。
    """
    p = float(axial_preload_n)
    p_lo, p_hi = 0.0, 2000.0
    m_lo, m_hi = float(m_fail_no_preload_nm), float(m_fail_2kn_nm)
    if p_hi == p_lo:
        return m_hi
    if p <= p_lo:
        return m_lo
    if p >= p_hi:
        if extrapolation == "clamp":
            return m_hi
        if extrapolation == "linear":
            slope = (m_hi - m_lo) / (p_hi - p_lo)
            return m_hi + slope * (p - p_hi)
        raise ValueError(
            f"extrapolation 必须是 'clamp'/'linear'，得到 {extrapolation!r}"
        )
    # 线性插值
    return m_lo + (m_hi - m_lo) * (p - p_lo) / (p_hi - p_lo)


# ---------------------------------------------------------------------------
# 关节级判据（ATFL / CFL）—— risk = M_inv / M_fail
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class LigamentCriterion:
    """一次韧带扭伤判据解析。"""

    ligament: str                 # "ATFL" / "CFL"
    M_inv_nmm: float              # 输入（N·mm）
    M_fail_nmm: float             # 内翻承载力（N·mm）
    preload_model: str            # "conservative" / "elevated"
    axial_preload_n: float        # 用于插值的轴向预载（N）
    risk: float                   # M_inv / M_fail
    sprain: bool                  # risk ≥ 1

    @property
    def M_inv_nm(self) -> float:
        return self.M_inv_nmm / 1000.0

    @property
    def M_fail_nm(self) -> float:
        return self.M_fail_nmm / 1000.0


def _select_m_fail(
    ligament: str,
    *,
    preload_model: str,
    axial_preload_n: float,
    m_fail_nm: float | None,
) -> tuple[float, float, str]:
    """根据韧带 + 预载模型选择 ``M_fail`` (N·m)。

    返回 ``(M_fail_Nm, axial_preload_N_used, model_label)``。
    """
    if m_fail_nm is not None:
        # 调用者显式给定 → 跳过模型
        return float(m_fail_nm), float(axial_preload_n), "explicit"

    # 默认映射：CFL 用与 ATFL 同一口径（保守：沿用 ATFL M_fail，
    # 现实中 CFL 极限 ≈ 1.6× ATFL 但失效晚于 ATFL——文献口径不一；
    # 用同一 M_fail 让"CFL 越阈 = 至少 ATFL 也已越阈"成立）。
    if preload_model == "conservative":
        # 无预载口径：用 ATFL_INVERSION_FAIL_MOMENT_NO_PRELOAD_NM
        m = ATFL_INVERSION_FAIL_MOMENT_NO_PRELOAD_NM
        ax = 0.0
        return float(m), float(ax), "conservative_no_preload"

    if preload_model == "elevated":
        # 2 kN 预载口径：用 2 kN 实测值
        m = ATFL_INVERSION_FAIL_MOMENT_2KN_PRELOAD_NM
        ax = 2000.0
        return float(m), float(ax), "elevated_2kn_preload"

    raise ValueError(
        f"preload_model 必须是 'conservative' / 'elevated'，得到 {preload_model!r}"
    )


def atfl_sprain_risk(
    inversion_moment_nmm: float,
    *,
    preload_model: Literal["conservative", "elevated"] = "conservative",
    axial_preload_n: float = 0.0,
    m_fail_nm: float | None = None,
) -> LigamentCriterion:
    """ATFL 关节级扭伤判据：``risk = M_inv / M_fail``。

    Parameters
    ----------
    inversion_moment_nmm
        内翻力矩 ``M_inv`` (N·mm)，来自 ``ankle_supination.inversion_moment_nmm``。
    preload_model
        ``"conservative"`` (默认) → 无预载 ``M_fail = 21 N·m``。
        ``"elevated"``            → 2 kN 预载 ``M_fail = 77 N·m``。
    axial_preload_n
        当 ``preload_model`` 在内部再细分时使用的轴向预载 (N)。
        当 ``preload_model="conservative"`` 时此值忽略（钳 0）。
    m_fail_nm
        调用者显式给定的 M_fail (N·m)。非 None 时**覆盖**模型选择。

    Notes
    -----
    默认 ``preload_model="conservative"`` ⇒ 在 β=0（``M_inv=0``）时
    risk = 0 ⇒ **逐位复现**今天的行为（无任何越阈）。
    """
    m_inv = abs(float(inversion_moment_nmm))
    m_fail, ax_used, label = _select_m_fail(
        "ATFL",
        preload_model=str(preload_model),
        axial_preload_n=float(axial_preload_n),
        m_fail_nm=m_fail_nm,
    )
    risk = m_inv / (m_fail * 1000.0)
    return LigamentCriterion(
        ligament="ATFL",
        M_inv_nmm=m_inv,
        M_fail_nmm=m_fail * 1000.0,
        preload_model=str(preload_model),
        axial_preload_n=float(ax_used),
        risk=float(risk),
        sprain=bool(risk >= 1.0),
    )


def cfl_sprain_risk(
    inversion_moment_nmm: float,
    *,
    preload_model: Literal["conservative", "elevated"] = "conservative",
    axial_preload_n: float = 0.0,
    m_fail_nm: float | None = None,
) -> LigamentCriterion:
    """CFL 关节级扭伤判据（与 ATFL 同口径）。

    文献上 CFL 极限力 ≈ 1.6× ATFL 但**失效晚于 ATFL**（组合载荷 ~27°，
    ATFL 先失效）。这里**沿用同一 ``M_fail`` 口径**作为保守下限——任何
    ATFL 越阈必然伴 CFL 越阈；任何 CFL 越阈必然伴 ATFL 也越阈。

    Parameters 与 :func:`atfl_sprain_risk` 同。
    """
    m_inv = abs(float(inversion_moment_nmm))
    m_fail, ax_used, label = _select_m_fail(
        "CFL",
        preload_model=str(preload_model),
        axial_preload_n=float(axial_preload_n),
        m_fail_nm=m_fail_nm,
    )
    risk = m_inv / (m_fail * 1000.0)
    return LigamentCriterion(
        ligament="CFL",
        M_inv_nmm=m_inv,
        M_fail_nmm=m_fail * 1000.0,
        preload_model=str(preload_model),
        axial_preload_n=float(ax_used),
        risk=float(risk),
        sprain=bool(risk >= 1.0),
    )


# ---------------------------------------------------------------------------
# 力 / 应变判据（r_lig assumed）
# ---------------------------------------------------------------------------
def ligament_force_n(
    inversion_moment_nmm: float,
    r_lig_mm: float,
) -> float:
    """韧带内力 ``F_lig = M_inv / r_lig`` (N)。

    ``r_lig`` 是韧带到内翻轴的力臂，**assumed**（无 CoP）。
    """
    if r_lig_mm <= 0.0:
        raise ValueError(f"r_lig_mm 必须为正，得到 {r_lig_mm!r}")
    return abs(float(inversion_moment_nmm)) / float(r_lig_mm)


def ligament_strain(
    force_n: float,
    *,
    stiffness_n_per_mm: float,
    resting_length_mm: float,
) -> float:
    """线弹性应变 ``ε = F / (k · L)``（secant stiffness 线性假设）。

    当 ``k`` 或 ``L`` ≤ 0 时返回 0（防零除）。
    """
    f = abs(float(force_n))
    k = float(stiffness_n_per_mm)
    L = float(resting_length_mm)
    if k <= 0.0 or L <= 0.0:
        return 0.0
    return f / (k * L)


def sprain_risk(
    inversion_moment_nmm: float,
    *,
    ligament: Literal["ATFL", "CFL"] = "ATFL",
    r_lig_mm: float = DEFAULT_LIGAMENT_LEVER_ARM_MM,
    stiffness_n_per_mm: float | None = None,
    resting_length_mm: float | None = None,
    ultimate_load_n: float | None = None,
    failure_strain: float | None = None,
) -> dict:
    """力 / 应变双判据：``F_lig / F_fail`` 与 ``ε / ε_fail``。

    返回 ``F_lig``、``ε``、``risk_force``、``risk_strain``、``sprain_force``、
    ``sprain_strain``、``sprain``（任一越阈）。

    默认参数为 ATFL；显式切换 ligament 自动套用 CFL 文献值（也可显式覆盖）。
    """
    lig = str(ligament).upper()
    if lig == "ATFL":
        k = stiffness_n_per_mm if stiffness_n_per_mm is not None else ATFL_STIFFNESS_N_PER_MM
        L = resting_length_mm if resting_length_mm is not None else ATFL_RESTING_LENGTH_MM
        Ff = ultimate_load_n if ultimate_load_n is not None else ATFL_ULTIMATE_LOAD_N
        ef = failure_strain if failure_strain is not None else ATFL_FAILURE_STRAIN
    elif lig == "CFL":
        k = stiffness_n_per_mm if stiffness_n_per_mm is not None else CFL_STIFFNESS_N_PER_MM
        L = resting_length_mm if resting_length_mm is not None else CFL_RESTING_LENGTH_MM
        Ff = ultimate_load_n if ultimate_load_n is not None else CFL_ULTIMATE_LOAD_N
        ef = failure_strain if failure_strain is not None else CFL_FAILURE_STRAIN
    else:
        raise ValueError(f"ligament 必须是 'ATFL'/'CFL'，得到 {ligament!r}")

    F_lig = ligament_force_n(inversion_moment_nmm, r_lig_mm)
    eps = ligament_strain(F_lig, stiffness_n_per_mm=k, resting_length_mm=L)
    r_force = F_lig / Ff
    r_strain = eps / ef
    sprain = (r_force >= 1.0) or (r_strain >= 1.0)
    # 透明度诊断：两条子判据的"触发力"（文献参数在 toe 区 vs 失效区不自洽，
    # 线性弹簧下 F@ε_fail = ε_fail·k·L 常远小于 F_fail）。
    F_at_failure_strain = ef * k * L
    M_thr_force = Ff * r_lig_mm
    M_thr_strain = F_at_failure_strain * r_lig_mm
    return {
        "ligament": lig,
        "r_lig_mm": float(r_lig_mm),
        "stiffness_n_per_mm": float(k),
        "resting_length_mm": float(L),
        "ultimate_load_n": float(Ff),
        "failure_strain": float(ef),
        "F_lig_n": float(F_lig),
        "strain": float(eps),
        "risk_force": float(r_force),
        "risk_strain": float(r_strain),
        "sprain_force": bool(r_force >= 1.0),
        "sprain_strain": bool(r_strain >= 1.0),
        "sprain": bool(sprain),
        # --- 透明度 / 自洽性诊断 -------------------------------------------
        "force_trigger_n": float(F_at_failure_strain),
        "force_threshold_n": float(Ff),
        "M_inv_threshold_force_nmm": float(M_thr_force),
        "M_inv_threshold_strain_nmm": float(M_thr_strain),
        "params_consistent": bool(abs(F_at_failure_strain - Ff) <= 0.05 * Ff),
        "criterion_note": (
            "力判据触发力 = F_fail；应变判据触发力 = ε_fail·k·L。"
            "文献的 secant 刚度（toe 区）与极限载荷（失效区）在纯线弹性下不自洽："
            "ε_fail·k·L 通常 ≪ F_fail。应变判据因此更保守；"
            "这也是主口径选关节级 M_inv/M_fail 的原因。"
        ),
    }


# ---------------------------------------------------------------------------
# 临界角 / 力矩
# ---------------------------------------------------------------------------
def critical_inversion_moment_nmm(
    *,
    ligament: Literal["ATFL", "CFL"] = "ATFL",
    preload_model: Literal["conservative", "elevated"] = "conservative",
    axial_preload_n: float = 0.0,
    m_fail_nm: float | None = None,
) -> float:
    """使该韧带首次扭伤（risk=1）的最小 ``M_inv`` (N·mm)。"""
    m_fail, _, _ = _select_m_fail(
        str(ligament),
        preload_model=str(preload_model),
        axial_preload_n=float(axial_preload_n),
        m_fail_nm=m_fail_nm,
    )
    return float(m_fail * 1000.0)


def critical_inversion_angle_deg(
    force_vec_n,
    *,
    ligament: Literal["ATFL", "CFL"] = "ATFL",
    preload_model: Literal["conservative", "elevated"] = "conservative",
    axial_preload_n: float = 0.0,
    lever_arm_mm: float = DEFAULT_LIGAMENT_LEVER_ARM_MM,  # CoP→距下轴力臂
    m_fail_nm: float | None = None,
) -> float | None:
    """使该韧带首次扭伤（risk=1）的最小踝旋后角 (度)。

    在 ``[0, 60]°`` 范围二分（保守口径下阈值 < 2°，故下限用 0° 起点）。
    若 60° 仍未越阈返回 ``None``。``F_y`` 为 0 时返回 ``None``。

    与 :func:`ankle_supination.critical_inversion_angle_deg` 的差别：
    这里是 **韧带扭伤** 阈值（M_inv/M_fail_ATFL），后者是 **腓骨远端
    弯曲** 阈值（σ_lat/σ_c）。
    """
    f = np.asarray(force_vec_n, dtype=float)
    if f.shape != (3,) or not np.all(np.isfinite(f)):
        raise ValueError(f"force_vec_n 必须是 3 个有限数，得到 {force_vec_n!r}")
    from .ankle_supination import frontal_lateral_force_n, inversion_moment_nmm

    m_thr = critical_inversion_moment_nmm(
        ligament=str(ligament),
        preload_model=str(preload_model),
        axial_preload_n=float(axial_preload_n),
        m_fail_nm=m_fail_nm,
    )

    def _m_inv(beta_deg: float) -> float:
        fl = frontal_lateral_force_n(f, beta_deg)
        return inversion_moment_nmm(fl, lever_arm_mm)

    if abs(_m_inv(0.0)) >= m_thr:
        return 0.0
    if abs(_m_inv(60.0)) < m_thr:
        return None
    lo, hi = 0.0, 60.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if _m_inv(mid) >= m_thr:
            hi = mid
        else:
            lo = mid
    return float(hi)


# ---------------------------------------------------------------------------
# 单格联合判据（M_inv → ATFL + CFL 双扭伤 + r_lig 灵敏度）
# ---------------------------------------------------------------------------
def ankle_ligament_criterion(
    m_inv_nmm: float,
    *,
    r_lig_mm: float = DEFAULT_LIGAMENT_LEVER_ARM_MM,
    preload_model: Literal["conservative", "elevated"] = "conservative",
    axial_preload_n: float = 0.0,
) -> dict:
    """在 ``M_inv`` 上同时跑 ATFL / CFL 关节级 + 力/应变判据。

    返回 ``m_inv_nmm``、``atfl_joint``、``cfl_joint``、``atfl_force_strain``、
    ``cfl_force_strain``，及 ``r_lig_mm`` 灵敏度口径下的 ATFL 力/应变重跑。

    Notes
    -----
    * 不调用 OpenSim；纯模块函数；与 ``ankle_supination`` 解耦。
    * 默认 ``preload_model="conservative"`` ⇒ 在 β=0（M_inv=0）时所有 risk
      均为 0 ⇒ **逐位复现今天行为**。
    """
    m_inv = abs(float(m_inv_nmm))

    atfl = atfl_sprain_risk(
        m_inv,
        preload_model=preload_model,
        axial_preload_n=axial_preload_n,
    )
    cfl = cfl_sprain_risk(
        m_inv,
        preload_model=preload_model,
        axial_preload_n=axial_preload_n,
    )
    atfl_fs = sprain_risk(
        m_inv, ligament="ATFL", r_lig_mm=r_lig_mm,
    )
    cfl_fs = sprain_risk(
        m_inv, ligament="CFL", r_lig_mm=r_lig_mm,
    )

    # r_lig 灵敏度（默认 + 20 + 25 mm）
    sensitivity = {}
    for r in R_LIG_SENSITIVITY_MM:
        if abs(r - r_lig_mm) < 1e-12:
            continue  # 默认已算
        sensitivity[r] = sprain_risk(
            m_inv, ligament="ATFL", r_lig_mm=r,
        )

    return {
        "m_inv_nmm": m_inv,
        "m_inv_nm": m_inv / 1000.0,
        "r_lig_mm": float(r_lig_mm),
        "preload_model": str(preload_model),
        "axial_preload_n": float(axial_preload_n),
        "atfl_joint": atfl,
        "cfl_joint": cfl,
        "atfl_force_strain": atfl_fs,
        "cfl_force_strain": cfl_fs,
        "r_lig_sensitivity": sensitivity,
    }


def ankle_ligament_criterion_at_grid(
    heights_m: tuple[float, ...],
    angles_deg: tuple[float, ...],
    measured_heights_m: list[dict],
    *,
    r_lig_mm: float = DEFAULT_LIGAMENT_LEVER_ARM_MM,
    lever_arm_mm: float = 30.0,
) -> list[dict]:
    """在 ``(h, β)`` 网格上跑韧带判据；输入 ``measured_heights_m`` 来自 R1 脚本。

    ``measured_heights_m`` 是 ``results/opensim_fe/nonvertical_s5.json`` 中
    ``measured_heights`` 的列表；这里只取 ``on_pad=False`` 的硬格。
    """
    from .ankle_supination import frontal_lateral_force_n, inversion_moment_nmm

    hard = [m for m in measured_heights_m if not m.get("on_pad", True)]
    by_h: dict[float, dict] = {float(m["height_m"]): m for m in hard}
    rows: list[dict] = []
    for h in heights_m:
        m = by_h.get(float(h))
        if m is None:
            continue
        f_vert = float(m["grf_peak_right_n"])
        fvec = (0.0, f_vert, 0.0)
        for beta in angles_deg:
            fl = abs(frontal_lateral_force_n(fvec, float(beta)))
            m_inv = inversion_moment_nmm(fl, lever_arm_mm)
            for pm in ("conservative", "elevated"):
                row = {
                    "height_m": float(h),
                    "supination_deg": float(beta),
                    "preload_model": pm,
                    "F_vert_n": f_vert,
                    "F_lat_n": float(fl),
                    "M_inv_nmm": float(m_inv),
                    "M_inv_nm": float(m_inv) / 1000.0,
                    "r_lig_mm": float(r_lig_mm),
                    "lever_arm_mm": float(lever_arm_mm),
                }
                # 关节级
                row["atfl_joint"] = atfl_sprain_risk(
                    m_inv, preload_model=pm,
                )
                row["cfl_joint"] = cfl_sprain_risk(
                    m_inv, preload_model=pm,
                )
                # 力 / 应变
                row["atfl_force_strain"] = sprain_risk(
                    m_inv, ligament="ATFL", r_lig_mm=r_lig_mm,
                )
                row["cfl_force_strain"] = sprain_risk(
                    m_inv, ligament="CFL", r_lig_mm=r_lig_mm,
                )
                rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# 自检：python -m climbing.coupling.ankle_ligament
# ---------------------------------------------------------------------------
def _main() -> int:
    print("=== ankle_ligament self-check ===")
    # 默认：β=0 → risk=0（逐位复现旧行为）
    z = atfl_sprain_risk(0.0)
    print(f"ATFL M_inv=0 → risk={z.risk:.3e} sprain={z.sprain}")
    z = cfl_sprain_risk(0.0)
    print(f"CFL  M_inv=0 → risk={z.risk:.3e} sprain={z.sprain}")

    # 临界力矩
    print(f"ATFL M_fail (conservative) = {critical_inversion_moment_nmm(preload_model='conservative')/1000.0:.1f} N·m")
    print(f"ATFL M_fail (elevated)     = {critical_inversion_moment_nmm(preload_model='elevated')/1000.0:.1f} N·m")

    # 插值
    for p in (0.0, 1000.0, 2000.0, 5000.0, 21000.0):
        print(f"  preload {p:6.0f} N → M_fail = {preload_capacity_interp(p):.1f} N·m")

    # 临界角
    fvec = (0.0, 21000.0, 0.0)
    for pm in ("conservative", "elevated"):
        b = critical_inversion_angle_deg(fvec, preload_model=pm)
        print(f"  critical_inversion_angle ({pm}) = {b}")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())