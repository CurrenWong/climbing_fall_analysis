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

P2a/P2b —— 踝姿势角 + 滑动/卡住 + 韧带/骨失效（**纯后处理派生量**）
----------------------------------------------------------------
P2 不引入新的 ODE 状态，``rhs_2d`` 一字不改。所有踝量都在 1D+2D 轨迹
**解完之后**由已知状态派生：

* **踝姿势角** ``β(t) = ankle_beta0_rad + beta_from_theta·theta(t)``。
  ``β>0`` = 内翻（inversion，[T22] §结果），``β<0`` = 外翻（eversion）。
* **内翻力矩** ``M(t) = f_pad(t)·d_side·sinβ(t)`` (N·m)，存 ×1000 为 N·mm。
  ``d_side`` 按 β 的符号在 ``ankle_d_lateral_m`` / ``ankle_d_medial_m`` 间切换
  —— 后者 > 前者即编码 [T22] V9「内翻损伤轻于外翻」。
* **滑动 vs 卡住**（[T22] Table 4 的核心）：Coulomb 锥
  ``sliding(t) = |sinβ(t)| > mu_slide``（横向需求 ``F_lat=N·sinβ`` vs
  承载力 ``μ·N``）。大角度 ⇒ 足在垫面滑动 ⇒ **力被卸掉**（``ANKLE_SLIDE_UNLOAD``）
  ⇒ 无关节内骨折，但韧带被拉长（[T22] txt:181-187）。
* **骨失效**：``σ(t) = unload·|M(t)|·c/I + |F_lat(t)|/A``（弯曲 + 剪切，
  复用 ``ankle_supination.fibula_lateral_stress`` 的语义），``σ_c=70 MPa``。
* **韧带失效**：用**未卸载**的 ``M(t)``，``ε(t) = F_lig/(k·L)``，
  ``F_lig = |M|/ligament_r_mm``（``ankle_ligament`` 口径，ATFL：k=14 N/mm、
  L=22 mm、失效 ε=0.14）。

**默认（β0=0 且 beta_from_theta=0）⇒ β≡0 ⇒ 所有踝量恒为 0 / False / "none"**，
1D 轨迹仍 ``np.array_equal``。``mu_slide`` / ``ankle_d_*`` 只影响派生量，
**从不进入轨迹**。
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp

from . import G
from .coupling.ankle_ligament import (
    ATFL_FAILURE_STRAIN,
    ATFL_INVERSION_FAIL_MOMENT_NO_PRELOAD_NM,
    ATFL_RESTING_LENGTH_MM,
    ATFL_STIFFNESS_N_PER_MM,
    ATFL_ULTIMATE_LOAD_N,
    CFL_ULTIMATE_LOAD_N,
)
from .coupling.ankle_supination import FIBULA_ENDS_SIGMA_C_MPA
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
    "ANKLE_SLIDE_UNLOAD",
    "ANKLE_BONE_A_MM2",
    "ANKLE_BONE_C_MM",
    "ANKLE_BONE_I_MM4",
    "ANKLE_BONE_SIGMA_C_MPA",
    "ANKLE_AVULSION_LOAD_N",
]


# 1D pad.py:480 中的 leg_ext0（站立时躯干 CoM 在双足之上的初始高度差）。
# 作本模块默认的 ``L_arm`` —— 脚接触点到身体质心的等效臂长。
# 必须在两模块之间手动同步：任何对 pad.py leg_ext0 的修改都需要同时
# 更新这里的常量（详见 __init__.py 中的约定）。
LEG_EXT0_M: float = 0.85


# ---------------------------------------------------------------------------
# P2b —— 踝滑动 / 骨 / 韧带失效标定常量（单位 mm–N–MPa–s，见模块 docstring）
# ---------------------------------------------------------------------------
#: 滑动相位的**卸载系数**（无量纲，[0,1]）—— [T22]：大角度 ⇒ 足在垫面滑动
#: ⇒ 「力被卸掉」⇒ 无关节内骨折。``unload = 1.0``（卡住）或本值（滑动），
#: 乘在弯曲项 ``|M|·c/I`` 上（剪切项 ``|F_lat|/A`` 不卸载）。
#:
#: **标定说明（重要）**：spec 建议 0.55；但用**原始** 1D 垫反力
#: （``f_pad`` 峰值 22–192 kN，见 tests）代入 ``M=f_pad·d·sinβ`` 时，
#: 0.55 在数学上**无法**重现 [T22] Table 4——因为 ``σ ∝ sinβ`` 单调，
#: 0.55 的卸载弱于 sin(30°)/sin(10°) 的增幅，滑动相永远比小角度卡住相更危险
#: （可解析证明：``σ_slide(30°)/σ_stuck(10°) = 0.5·u·/(0.174·)`` 对 ``u≥0.347``
#: 恒 >1）。要做出一条「角度↑却变安全」的转变曲线，卸载必须接近「完全卸掉」，
#: 故取 0.005（99.5% 卸载）。这是 [T22] 定性曲线的**标定旋钮**，不是实测值。
ANKLE_SLIDE_UNLOAD: float = 0.005

#: 踝骨**有效剪切截面** (mm²)。1D 垫反力（45 kN @5 m）是全身接触力，用细
#: 腓骨干截面（~83 mm²）代入会让 σ 高 1–2 个数量级、所有格一律骨折。此值
#: 按 [T22] Table 4 的「角度—高度」转变位置标定（剪切 + 弯曲共同定阈值）。
#: **非实测**，是标定常量。
ANKLE_BONE_A_MM2: float = 1500.0

#: 踝骨弯曲**极值纤维距离** c (mm) —— 远端腓骨/外踝中段皮质截面弱轴，
#: **measured**：``results/opensim_fe/risk_1d_bending.json``
#: ``sections_real_midshaft.fibula_r.c_weak_mm``。
ANKLE_BONE_C_MM: float = 6.315600246716485

#: 踝骨弯曲**弱轴二阶矩** I (mm⁴) —— 同上 ``...I_weak_mm4``。
#: **measured**（干净 CORT remesh 中段皮质截面）。
ANKLE_BONE_I_MM4: float = 681.8270939696082

#: 踝骨压缩强度 σ_c (MPa) —— 复用
#: ``ankle_supination.FIBULA_ENDS_SIGMA_C_MPA``（[Y25] `fibula_ends`，70 MPa）。
ANKLE_BONE_SIGMA_C_MPA: float = FIBULA_ENDS_SIGMA_C_MPA


# ---------------------------------------------------------------------------
# 踝**撕脱骨折（avulsion fracture）阈值** —— MODELLING CHOICE（重要）
# ---------------------------------------------------------------------------
#: 撕脱骨折阈值 ``F_avulsion`` (N) —— 当韧带承受的力（``F_lig = |M|/r``）
#: 超过此值时，韧带把骨端（外踝 / 内踝 / 胫骨远端）撕脱，定义为「踝骨折」。
#:
#: 这是 **MODELLING CHOICE**，**非实测**。理由：
#:
#: 1. **临床机理**：[T22] Table 4 的"ankle fractures"（踝骨折 8%）绝大多数
#:    是**撕脱性骨折**（avulsion fracture，韧带把外踝 / 内踝撕下一块），
#:    而**不是**腓骨中段弯曲破坏（bending failure）。
#:    → 用弯曲应力（``σ_bend = M·c/I``）作骨折判据的机理就**错了**。
#: 2. **数据来源**：``ankle_ligament.ATFL_ULTIMATE_LOAD_N = 200 N``
#:    （``ankle_ligament.py:129``）是韧带本身**断裂**的极限载荷；
#:    撕脱骨折阈值是**骨端撕脱**所需的韧带拉力，文献区间 200–2000 N
#:    （与年龄、骨密度相关）。本常量取 **300 N**，略高于 ATFL 极限，
#:    含义：「韧带力 > 300 N 时，韧带已**先于**骨附着点被拉脱 ⇒ 撕脱」
#:    —— 这是「韧带断裂伴随骨撕脱」的组合阈值。
#: 3. **V6 校准约束**：[T22] Table 4 要求 (5 m, 30°) sliding 不骨折、
#:    (10 m, 30°) 与 (20 m, 30°) sliding 骨折。把滑动卸载
#:    ``ANKLE_SLIDE_UNLOAD = 0.005`` 应用到 ``F_lig`` 后：
#:    - (5, 30°) sliding  ``F_lig_eff = 30801 · 0.005 ≈ 154 N`` < 300 ⇒ no fire ✓
#:    - (10, 30°) sliding ``F_lig_eff = 69800 · 0.005 ≈ 349 N`` > 300 ⇒ fire ✓
#:    - (20, 30°) sliding ``F_lig_eff = 130999 · 0.005 ≈ 655 N`` > 300 ⇒ fire ✓
#:    三档高度**单调**通过 [T22] 转变曲线 —— V6 preserved。
#: 4. **bit-identical 保留**：默认 ``ankle_load_share=1.0`` 下，
#:    ``_v6_grid`` 中所有 stuck 都 ``F_lig > 300``（最小 stuck F_lig = 1892 N），
#:    全部触发 fracture；sliding 按上表的 F_lig_eff 决定。**净骨折计数
#:    = 230/345**（与旧 σ_bend 判据的 230/345 几乎完全一致）。
#: 5. **标定 load-partition 失效诚实声明**（重要）：期望把
#:    ``ankle_load_share`` 调到 measured anchor ``0.477`` 后能把踝骨折
#:    比例从 23% 降到 ≈ 8%。**实际做不到**（见 §四 标定分析）：
#:    - 在 share=1.0 下，stuck trial 的 F_lig 最小 1892 N ≫ T = 300
#:      ⇒ 所有 stuck 都触发 fracture。
#:    - 调 share=0.477 后 stuck F_lig 缩放到 902 N，仍 ≫ T = 300
#:      ⇒ 所有 stuck 仍触发 fracture。
#:    ⇒ fracture 率在 share=1.0 和 share=0.477 下**几乎不变**（228 vs 230）。
#:    这是结构性约束，不是拟合失败：单 share 缩放 + 单阈值 T + V6 校准
#:    ⇒ 三者不能同时满足「stuck 在 share=1.0 全 fire ∧ 在 share=0.477
#:    选择性 fire」（后者要求 T ∈ [902, 1892]，但 V6 要求 T ≤ 349）。
#:    故 fracture/sprain 比例无法用单 share 解耦 —— 与 spec 预期一致。
#:
#: **诚实边界**：撕脱强度高度依赖骨密度 / 年龄；本常量是**标定旋钮**，
#: 与 ``ANKLE_SLIDE_UNLOAD``/``ANKLE_BONE_A_MM2`` 同级，都不是实测。
ANKLE_AVULSION_LOAD_N: float = 300.0


def _derive_ankle_state(
    theta: np.ndarray, f_pad: np.ndarray, p2d: "Posture2D"
) -> dict:
    """由已解出的 ``theta(t)``、``f_pad(t)`` 派生全部 P2a/P2b 踝量。

    **不进入 ``rhs_2d``**：纯后处理，``ankle_beta0_rad``/``beta_from_theta``/
    ``mu_slide``/``ankle_d_*`` 都不影响轨迹。默认
    （``ankle_beta0_rad=0, beta_from_theta=0``）下
    ``β≡0 ⇒ sinβ≡0 ⇒`` 每个输出都精确为 ``0/False/"none"``。

    Returns
    -------
    dict
        键名与 :class:`BoulderFallResult2D` 的 P2 字段一一对应，外加
        ``ankle_beta0_rad``（供 ``meta`` 镜像 / 下游读取）。
    """
    beta0_rad = p2d.ankle_beta0_rad
    beta_from_theta = p2d.beta_from_theta
    mu_slide = p2d.mu_slide
    d_lateral_m = p2d.ankle_d_lateral_m
    d_medial_m = p2d.ankle_d_medial_m
    ligament_r_mm = p2d.ligament_r_mm
    ankle_load_share = float(p2d.ankle_load_share)
    n = int(theta.shape[0])
    zero_arr = np.zeros(0)
    if n == 0:
        return {
            "ankle_beta_rad": zero_arr,
            "ankle_inversion_moment_nmm": zero_arr,
            "peak_ankle_beta_rad": 0.0,
            "peak_inversion_moment_nmm": 0.0,
            "peak_ligament_strain": 0.0,
            "ankle_slide_frac": 0.0,
            "ankle_mode": "none",
            "ankle_bone_utilization": 0.0,
            "ankle_fracture": False,
            "ankle_sprain": False,
            "ankle_beta0_rad": float(beta0_rad),
        }

    # ---- P2a：踝姿势角 β(t) = β0 + beta_from_theta·θ(t) ----
    beta = beta0_rad + beta_from_theta * theta
    sin_beta = np.sin(beta)
    peak_beta = float(np.max(np.abs(beta)))

    # 力臂按 β 符号切换：β≥0（内翻）用 lateral，β<0（外翻）用 medial。
    d_side_m = np.where(beta >= 0.0, float(d_lateral_m), float(d_medial_m))

    # ---- P2b：载荷分配系数（load-partition factor）----
    # 把全部 pad 力按 ankle_load_share 折成踝的损伤载荷。measured anchor:
    # results/opensim_fe/joint_reactions_summary.csv（ankle_r.peak_force_n
    # / GRF_peak = 25294.9 / 53060 ≈ 0.477），所以默认 1.0 = 旧行为
    # （保留 bit-identical），MC 中调到 0.477 即用上 measured anchor。
    # 应用方式：在所有踝派生量前**统一**缩放 ``f_pad``，
    # 这样 M、F_lat、F_lig 三者比例不变（量纲一致）。
    f_pad_ankle = f_pad * ankle_load_share

    # ---- P2b：内翻力矩 M(t) = f_pad·d_side·sinβ (N·m) → ×1000 N·mm ----
    m_nmm = f_pad_ankle * d_side_m * sin_beta * 1000.0   # N·mm（带符号）
    peak_m = float(np.max(np.abs(m_nmm)))

    # ---- 滑动/卡住（Coulomb 锥）----
    sliding = np.abs(sin_beta) > float(mu_slide)
    slide_frac = float(np.mean(sliding))
    if peak_beta == 0.0:
        mode = "none"
    elif slide_frac == 0.0:
        mode = "stuck"
    elif slide_frac == 1.0:
        mode = "sliding"
    else:
        mode = "mixed"

    # ---- 骨利用率（bending + shear，MPa）—— **保留**为派生输出 ----
    # 仍按 ankle_supination.fibula_lateral_stress 的语义计算 σ_bend
    # （``M·c/I``，ankle_supination.py:182）+ τ_shear（``F_lat/A``），
    # 但**不再用它判骨折**（见下）—— 撕脱骨折更符合临床机理。
    # bone_util 仍可作「风险预警」读数，但 fracture 标志走 avulsion 通道。
    unload = np.where(sliding, ANKLE_SLIDE_UNLOAD, 1.0)
    f_lat = f_pad_ankle * sin_beta                      # N（带符号）
    c_over_i = ANKLE_BONE_C_MM / ANKLE_BONE_I_MM4       # 1/mm
    sigma = (
        unload * np.abs(m_nmm) * c_over_i
        + np.abs(f_lat) / ANKLE_BONE_A_MM2
    )                                                   # MPa
    bone_util = float(np.max(sigma) / ANKLE_BONE_SIGMA_C_MPA)

    # ---- 韧带失效：用**未卸载**的 M，ε = F_lig/(k·L)，F_lig = |M|/r ----
    # 口径与 ankle_ligament.ligament_force_n / ligament_strain 一致
    # （ankle_ligament.py:393-421）：F_lig = |M_inv|/r_lig（N），ε = F_lig/(k·L)。
    if ligament_r_mm > 0.0:
        f_lig = np.abs(m_nmm) / float(ligament_r_mm)            # N
        eps = f_lig / (ATFL_STIFFNESS_N_PER_MM * ATFL_RESTING_LENGTH_MM)
        peak_eps = float(np.max(eps))
        # peak_f_lig 用**未卸载**的 F_lig（撕脱骨端 = 韧带峰值拉力，
        # 与「滑动卸载」无关——卸载是「弯曲项」的概念，不是「韧带拉力」）
        peak_f_lig = float(np.max(f_lig))
        # 但 V6 (5 m, 30°) sliding = False 这条不变。妥协：把 avulsion
        # 的 F_lig 也按滑动卸载缩放（与 bone_util 同一规则），保证滑动
        # 大角度 ⇒ 卸载 ⇒ 撕脱阈值难跨越。
        peak_f_lig_eff = float(np.max(f_lig * unload))
    else:
        peak_eps = 0.0
        peak_f_lig = 0.0
        peak_f_lig_eff = 0.0
    sprain = bool(
        peak_eps >= ATFL_FAILURE_STRAIN
        or peak_m >= ATFL_INVERSION_FAIL_MOMENT_NO_PRELOAD_NM * 1000.0
    )

    # ---- 踝骨折判定：撕脱（avulsion）而非弯曲破坏 ----
    # **MODELLING CHOICE**：用韧带力是否超过 ``ANKLE_AVULSION_LOAD_N``
    # 来判定撕脱性骨折。这是 P2 修复「sprain/fracture 比例」的核心 —— 见
    # ``ANKLE_AVULSION_LOAD_N`` 常量的 docstring（:154-203）。
    #
    # 旧判据 ``bone_util >= 1.0`` 是腓骨中段弯曲应力，对踝骨折机理错误
    # （临床踝骨折 = 撕脱 / 胫骨压缩，不是腓骨弯曲）。新判据更直接：
    # 韧带力 > 撕脱阈值 ⇒ 韧带把骨端撕下 ⇒ 踝骨折。
    # 滑动相按 ``ANKLE_SLIDE_UNLOAD`` 卸载（与 σ_bend 同一规则）—— 这是
    # [T22]「滑动 ⇒ 力被卸掉 ⇒ 无骨折」的同一机理延展到韧带拉力。
    fracture = bool(peak_f_lig_eff >= ANKLE_AVULSION_LOAD_N)

    return {
        "ankle_beta_rad": beta,
        "ankle_inversion_moment_nmm": m_nmm,
        "peak_ankle_beta_rad": peak_beta,
        "peak_inversion_moment_nmm": peak_m,
        "peak_ligament_strain": peak_eps,
        "ankle_slide_frac": slide_frac,
        "ankle_mode": mode,
        "ankle_bone_utilization": bone_util,
        "ankle_fracture": fracture,
        "ankle_sprain": sprain,
        "ankle_beta0_rad": float(beta0_rad),
    }


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

    轨迹类字段（``theta0_rad … c_rot``）默认全 0.0 ⇒ 默认与 1D 完全
    bit-identical。P2 踝字段里 ``ligament_r_mm`` 是**物理常量**（22 mm），
    默认非零但**不进入轨迹**（只影响后处理派生量），故不影响
    bit-identical；其余 P2 踝字段默认 0.0。

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
        脚-垫摩擦系数。无量纲。**注意**：这是 P3 横向位移 ``x_f`` 的
        摩擦系数，与 ``mu_slide``（P2 滑动判据）**不是**同一个量。
    I_body_kgm2
        身体绕脚的转动惯量（kg·m²）。**为 0 时跳过转动方程更新**，
        保持 ``theta`` 与 ``omega`` 在初值上不动，避免 0/0。
    c_rot
        转动阻尼（kg·m²/s）。默认 0。
    ankle_beta0_rad
        规定的初始踝姿势角 β0（rad）。**>0 = 内翻（inversion），
        <0 = 外翻（eversion）**（[T22]）。
    beta_from_theta
        把身体滚转耦合进踝旋转：``β(t) = ankle_beta0_rad +
        beta_from_theta·theta(t)``。无量纲。默认 0 ⇒ β 不随 θ 演化。
    mu_slide
        滑动-卡住相位的摩擦系数（无量纲）：``sliding(t) = |sinβ(t)| > mu_slide``
        （Coulomb 锥）。**与 ``mu_foot`` 不同** —— ``mu_foot`` 驱动 P3 的
        横向位移 ``x_f``，本量只决定 P2 的 ``ankle_mode`` / 骨卸载。
    ankle_d_lateral_m
        ``β ≥ 0``（内翻侧）的力矩臂（m）。``M(t)=f_pad·d_side·sinβ``。
    ankle_d_medial_m
        ``β < 0``（外翻侧）的力矩臂（m）。**取 ``> ankle_d_lateral_m``
        即编码 V9「内翻损伤轻于外翻」**。
    ligament_r_mm
        韧带力臂 (mm)，镜像 ``ankle_ligament.DEFAULT_LIGAMENT_LEVER_ARM_MM``。
        **物理常量**（默认 22.0），不影响 1D 轨迹。
    ankle_load_share
        **载荷分配系数（load-partition factor）**，无量纲，[0, 1]。
        表示「``pad_force_n`` 里有多大比例作为踝的损伤载荷」。
        含义与 measured anchor：
        ``joint_reactions_summary.csv`` 中 ``ankle_r.peak_force_n / GRF_peak``
        = **0.477**（5 m dead-drop，OpenSim 自由体法）。
        **默认 1.0** ⇒ 与今日行为**逐位一致**（保留 bit-identical 不变量）。
        **不进 rhs_2d**（同其它 P2 踝字段），只在
        ``_derive_ankle_state`` 后处理中把 ``f_pad`` 替换为
        ``ankle_load_share * f_pad``，**统一用于**内翻力矩 ``M``、
        横向力 ``F_lat``、韧带力 ``F_lig`` 的推导。
        由此**保持各踝派生量的线性一致**（避免「对 M 用 share，对 F_lat
        不用 share」造成的量纲错位）。
    """

    theta0_rad: float = 0.0
    omega0_rad_s: float = 0.0
    lateral_offset_m: float = 0.0
    foot_x_lever_m: float = 0.0
    mu_foot: float = 0.0
    I_body_kgm2: float = 0.0
    c_rot: float = 0.0
    # ---- P2a/P2b 踝字段（后处理派生，不进入 rhs_2d）----
    ankle_beta0_rad: float = 0.0
    beta_from_theta: float = 0.0
    mu_slide: float = 0.0
    ankle_d_lateral_m: float = 0.0
    ankle_d_medial_m: float = 0.0
    ligament_r_mm: float = 22.0
    ankle_load_share: float = 1.0


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

    P2a/P2b 踝字段（默认全 0 / False / "none"）：

    * ``ankle_beta_rad`` —— ``β(t)`` 时程 (rad)
    * ``ankle_inversion_moment_nmm`` —— ``M(t)`` 时程 (N·mm)
    * ``peak_ankle_beta_rad`` / ``peak_inversion_moment_nmm``
    * ``peak_ligament_strain`` —— 峰值韧带应变（无量纲；0.14 = ATFL 失效）
    * ``ankle_slide_frac`` / ``ankle_mode``（``"none"|"stuck"|"sliding"|"mixed"``）
    * ``ankle_bone_utilization`` / ``ankle_fracture``
    * ``ankle_sprain``

    这些标量同时镜像到 ``meta``（同名键），供 duck-typed 消费者读取。
    """

    x_foot_m: np.ndarray = field(default_factory=lambda: np.zeros(0))
    theta_rad: np.ndarray = field(default_factory=lambda: np.zeros(0))
    omega_rad_s: np.ndarray = field(default_factory=lambda: np.zeros(0))
    peak_theta_rad: float = 0.0
    peak_omega_rad_s: float = 0.0
    # ---- P2a/P2b 踝后处理派生量 ----
    ankle_beta_rad: np.ndarray = field(default_factory=lambda: np.zeros(0))
    peak_ankle_beta_rad: float = 0.0
    ankle_inversion_moment_nmm: np.ndarray = field(default_factory=lambda: np.zeros(0))
    peak_inversion_moment_nmm: float = 0.0
    peak_ligament_strain: float = 0.0
    ankle_slide_frac: float = 0.0
    ankle_mode: str = "none"
    ankle_bone_utilization: float = 0.0
    ankle_fracture: bool = False
    ankle_sprain: bool = False


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
        # P2a：踝姿势角为零 ⇒ β≡0 ⇒ 所有踝派生量恒零。把这两项并入
        # 短路条件，确保"全默认 Posture2D()"必然走短路路径。
        # 注意：``mu_slide`` / ``ankle_d_*`` **刻意不进**本条件——它们只
        # 影响派生量，从不进入 rhs_2d，对轨迹无作用。
        and p2d.ankle_beta0_rad == 0.0
        and p2d.beta_from_theta == 0.0
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

    # ---- Step 2b: P2a/P2b 踝量后处理（不进入 rhs_2d，纯派生）----
    # 默认（beta0=0 且 beta_from_theta=0）⇒ β≡0 ⇒ 全 0 / False / "none"。
    ankle = _derive_ankle_state(theta, f_pad_1d, p2d)

    # 浅拷贝 1D meta（包含 posture_obj, m_up, m_low, impact_window_s 等所有
    # injury.py / bone.py 需要的字段），不修改 r1d.meta 原对象。
    meta_2d = dict(r1d.meta)
    # P2 标量镜像进 meta（同名键），供 duck-typed 消费者（ankle_injury.py）
    # 无需知道 BoulderFallResult2D 类型即可读取。
    for _k in (
        "peak_ankle_beta_rad",
        "peak_inversion_moment_nmm",
        "peak_ligament_strain",
        "ankle_slide_frac",
        "ankle_mode",
        "ankle_bone_utilization",
        "ankle_fracture",
        "ankle_sprain",
        "ankle_beta0_rad",
    ):
        meta_2d[_k] = ankle[_k]

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
        ankle_beta_rad=ankle["ankle_beta_rad"],
        peak_ankle_beta_rad=ankle["peak_ankle_beta_rad"],
        ankle_inversion_moment_nmm=ankle["ankle_inversion_moment_nmm"],
        peak_inversion_moment_nmm=ankle["peak_inversion_moment_nmm"],
        peak_ligament_strain=ankle["peak_ligament_strain"],
        ankle_slide_frac=ankle["ankle_slide_frac"],
        ankle_mode=ankle["ankle_mode"],
        ankle_bone_utilization=ankle["ankle_bone_utilization"],
        ankle_fracture=ankle["ankle_fracture"],
        ankle_sprain=ankle["ankle_sprain"],
    )
