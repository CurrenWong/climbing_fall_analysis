"""P4 · 头/颈惯性分支（路径 B：自上而下的应力传递）。

为什么需要这个模块
------------------
1D/2D 模型把躯干当作质点（``BoulderFallResult.z_torso_m`` 是躯干 CoM 的轨迹），
意味着"头"和"躯干"是一回事，**结构性漏掉了** [Y25] §2.3 的关键发现：

* **颅骨 / 颈椎在 t < 0.0125 s 就出现应力** —— 这**不可能**是从足部
  应力波传导上来（足到脑的距离 ~1.5 m，应力波传导时间 ~1 ms 倒是
  够快，但这条路径在刚体链里根本不存在）。
* 真正的机制是 **惯性**：足骤停 → 躯干骤停 → **头因惯性继续向下** →
  颈椎被压缩。

这就是 [Y25] 区分的两条路径：

* **路径 A（自下而上）**：足 → 胫骨/腓骨 → 股骨 → 骨盆
* **路径 B（自上而下）**：颅骨 + 颈椎因惯性继续下压 → 胸椎 → 腰椎

两路径在脊柱汇集，导致脊柱在冲击**后期**才出现应力集中
（[docs/方案更新_v2.md:548-553]）。

**脊柱后期的汇集峰 out of scope。** 本模块**只建路径 B**（头部惯性
分支）。两路径在脊柱叠加形成的"后期应力集中"需要带脊柱节的**多刚体链**
（``bodies.py`` 有结构但验收判据未通过、未标定）才能表达；把它硬塞进
这个 1-DOF 派生量只会得到一个没有物理依据的假峰。因此本模块的职责边界是：
把路径 B 的**早期**特征（头部应力峰）算准，路径 A→脊柱的后期汇集
留给 Phase B-3 的多刚体模型。

建模方案
--------
**严格的下游派生量**，**不进入** pad / pad2d 的 ODE 状态（这是 1D/2D
bit-identicality 的前提）。读 ``BoulderFallResult`` 的 ``z_torso_m(t)``
与 ``accel_torso_g(t)``，在结果**自己的**时间网格 ``t_s`` 上解一个
**单自由度基础激励问题**：

    相对坐标 delta(t) = z_head(t) - z_torso(t)（向下为正）
    颈部内力 F_neck = k_neck*delta + c_neck*delta_dot（标量，带符号）
    头部运动方程（向下为正）：
        m_head·z_head_ddot = m_head·g − k_neck·delta − c_neck·delta_dot
    其中 z_torso(t) 是**已知**基础激励，由 ``result.z_torso_m`` /
    ``result.accel_torso_g`` 提供。

**头部不反馈到躯干** —— 这是本模块与全耦合求解（多刚体）的关键区别。
模型把颈部当作"只在头端有质量、躯干端是刚性基础"的一根弹簧，与
[Y25] 的"**脊柱-头段有足够刚度，否则头部响应会被延迟**"对应：
``k_neck`` 越大 → 头越接近刚性跟随躯干（``delta → 0``）；
``k_neck`` 越小 → 头越接近自由质点（``delta`` 延迟最严重）。

V3 / V8 验证判据
----------------
[docs/方案更新_v2.md:861] V3：头部应力峰值时刻 ≤ 0.0125 s（与足部冲击
几乎同时），而不是延迟到载荷传导之后。

[docs/方案更新_v2.md:866] V8：峰值时刻随高度**提前**，20 m 时 ≤ 5 ms。
（更高 ⇒ 冲击速度更大 ⇒ 减速脉冲更窄 ⇒ 头更早达到峰值应力。）

单位约定（mm–N–MPa–s）
--------------------
所有 ``*_n_per_mm`` / ``*_n_s_per_mm`` 字段用 N 与 mm 组合；进入 ODE 时
换算到 SI（m、N、s）。``head_section_mm2`` 是把 ``F_neck`` 折算成
应力 σ = F/A 的承载截面（默认 1.0 mm² 把 F [N] 直接读作 MPa，
参见 ``HeadBranchParams.head_section_mm2`` 的注释）。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator

from . import G

__all__ = [
    "HeadBranchParams",
    "HeadVerdict",
    "simulate_head_branch",
    "head_band",
    "V3_THRESHOLD_S",
    "V8_THRESHOLD_AT_20M_S",
]


# ==========================================================================
# 验证阈值（直接 cite 自 docs/方案更新_v2.md）
# ==========================================================================
#: V3 判定阈值（s）。颅骨 / 颈椎应力峰值距冲击起点的延迟上界。
#: cite: docs/方案更新_v2.md:861 —— "V3 | [Y25] §2.3 | 头部应力峰值时刻
#: | ≤ 0.0125 s | ⏳ 需 P4"
V3_THRESHOLD_S: float = 0.0125

#: V8 判定阈值（s）。20 m 落差时头部峰值距冲击起点应 ≤ 5 ms。
#: cite: docs/方案更新_v2.md:866 —— "V8 | [T22] Table 3 | 峰值时刻随高度提前
#: | 20 m 时 ≤ 5 ms | ⏳ 需 P3/P4"
V8_THRESHOLD_AT_20M_S: float = 0.005


# ==========================================================================
# 参数容器
# ==========================================================================
@dataclass(frozen=True)
class HeadBranchParams:
    """头/颈惯性分支的标定参数。

    Attributes
    ----------
    m_head_kg
        头部质量 (kg)。**ASSUMPTION**：默认 5.0 kg。
        文献基础：de Leva (1996) "Adjustments to Zatsiorsky-Seluyanov's
        segment inertia parameters" —— 头部 CoM 质量约占全身 8.1%，
        80 kg 人对应 ~6.5 kg；本模型保守取 5.0 kg（占全身 6.25%），
        与 biomechanics 教材中"头占全身 6-8%"的下沿一致。
    k_neck_n_per_mm
        颈椎轴向刚度 (N/mm)。**ASSUMPTION**：默认 200 N/mm。
        文献基础：颈椎韧带柱（C2-C7）在亚失效轴向载荷下的准静态刚度
        量级 50-150 N/mm（Nightingale et al., 1996, J. Biomechanics
        Eng.）；动态冲击下等效刚度上浮 2-4× 至 ~200 N/mm。
        量级验证：80 kg 全身对应的颈部椎体极限压缩 ~3-4 kN，
        极限位移 1.5-2.5 mm ⇒ 极限刚度 1.5-2.5 kN/mm；本默认取
        极限刚度的 ~10%，对应**亚失效**区间，与 [Y25] §2.3 描述的
        "亚临界颅/颈应力" 一致。
        ``k_neck`` 太小 → 头响应延迟超出 V3 的 12.5 ms 上界。
    c_neck_n_s_per_mm
        颈椎粘性阻尼系数 (N·s/mm)。**ASSUMPTION**：默认 0（无阻尼，
        最坏情况保守估计 —— 实测难，文献缺）。
        临界阻尼比 ``ζ = c / (2·sqrt(k·m))``：k=200 N/mm, m=5 kg
        ⇒ ``c_critical ≈ 2.0 N·s/mm``；本默认 0 ⇒ ζ = 0（无阻尼，
        振荡峰最大）。
    head_section_mm2
        头部 / 颈部承载截面 (mm²)。**ASSUMPTION**：默认 1.0 mm²
        把 ``σ = F/A`` 折算成一个**无量纲数值上接近 F [kN] 的应力度量**
        —— 数值上 ``F [N] / 1.0 mm² = F [N/mm²] = F [MPa]``。
        物理意义：实际颈椎 C1-C2 关节面 ~50-150 mm²，
        椎体截面 ~150-300 mm²；要换算成真实应力 (MPa)，
        应把 ``head_section_mm2`` 设为 ~200 mm²。
        分档阈值 ``HEAD_STRESS_BANDS_MPA`` 在默认截面下解读为
        "F [N] 对应的 MPa 数值"；切换截面时分档会自动等比例缩放。
    g
        重力加速度 (m/s²)。默认 ``climbing.G`` = 9.80665，
        允许低重力场景在测试中覆盖。
    """

    m_head_kg: float = 5.0
    k_neck_n_per_mm: float = 200.0
    c_neck_n_s_per_mm: float = 0.0
    head_section_mm2: float = 1.0
    g: float = G


# ==========================================================================
# 结果容器
# ==========================================================================
@dataclass(frozen=True)
class HeadVerdict:
    """头/颈惯性分支的派生判读。

    Attributes
    ----------
    peak_stress_mpa
        头部应力度量峰值 (MPa)，``σ_peak = |F_neck_peak| / head_section``。
    peak_time_s
        峰值距**冲击起点** ``t0 = meta['impact_window_s'][0]`` 的延迟 (s)。
        V3 判据：``peak_time_s <= 0.0125 s``。
    peak_force_n
        颈部内力绝对值的峰值 (N)。
    delta_peak_m
        峰值时刻头-躯干相对位移 (m)。``>0`` 表示头"赶上"躯干
        （颈部压缩）；``<0`` 表示头被颈部拉离躯干（颈部拉伸）。
    within_y25_window
        ``True`` iff ``peak_time_s <= V3_THRESHOLD_S``。
    band
        四档分档：``{"低", "中", "高", "极高"}``，与 ``injury.VERDICTS``
        同源。
    note
        物理机制说明（中文），含分档依据与 V3 通过状态。
    """

    peak_stress_mpa: float
    peak_time_s: float
    peak_force_n: float
    delta_peak_m: float
    within_y25_window: bool
    band: str
    note: str


# ==========================================================================
# 头/颈应力分档（mm–N–MPa 单位下）
# ==========================================================================
#: 头/颈应力度量分档 (MPa)。默认 ``head_section_mm2 = 1.0 mm²`` 时，
#: ``σ [MPa]`` 数值上**等于** ``F [N] / 1 mm² = F [MPa]``，因此
#: 本表 1 MPa 对应 1 kN 力，6 MPa 对应 6 kN，与
#: ``injury.BONE_BANDS`` 的骨骼轴向分档（6/10/16 kN）口径一致。
#:
#: **标定依据（ASSUMPTION）**：
#: * "低" 上界 1 MPa ⇒ F < 1 kN —— 远低于任何颈椎骨折阈值；
#: * "中" 上界 3 MPa ⇒ F < 3 kN —— 头部惯性响应的常规峰值；
#: * "高" 上界 6 MPa ⇒ F < 6 kN —— ``injury.BONE_BANDS`` 的"低→中"分界，
#:   此处作为"高"分档是因本应力度量用的是 mm² 单位，比 BONE_BANDS 的
#:   kN 单位更敏感；
#: * "极高" 6 MPa+ ⇒ F ≥ 6 kN —— 进入"显著颈椎/颅骨压缩"区间。
HEAD_STRESS_BANDS_MPA: tuple[tuple[float, float, str], ...] = (
    (0.0, 1.0, "低"),
    (1.0, 3.0, "中"),
    (3.0, 6.0, "高"),
    (6.0, float("inf"), "极高"),
)


def head_band(peak_stress_mpa: float) -> str:
    """把 ``σ_peak`` 分到四档：``{"低","中","高","极高"}``。

    与 ``injury.VERDICTS = ["低","中","高","极高"]`` 同源；阈值表见
    ``HEAD_STRESS_BANDS_MPA``。
    """
    for lo, hi, label in HEAD_STRESS_BANDS_MPA:
        if lo <= peak_stress_mpa < hi:
            return label
    return HEAD_STRESS_BANDS_MPA[-1][2]


def _first_significant_peak(values: np.ndarray, rel_threshold: float = 0.5) -> int:
    """返回 ``values`` 中**首个显著局部极大**的下标。

    为什么不能用全局 ``argmax``
    ---------------------------
    默认 ``c_neck = 0``（无阻尼）时，基础激励结束后头-颈系统进入
    **等幅自由振荡**：``|F_neck|`` 的每个峰都几乎一样高（实测 100% /
    100% / 100%）。全局 ``argmax`` 会被浮点噪声在几个等幅峰之间随机
    挑一个，让 ``peak_time_s`` **不确定**（同一输入两次可得到不同结果）。

    物理上"头部应力峰值时刻"指的是**首次**应力集中（[Y25] §2.3：
    颅骨/颈椎在 t < 0.0125 s 就出现应力），所以这里取**首个**超过
    全局最大值 ``rel_threshold`` 倍的局部极大。有阻尼时首个峰本就是
    最大峰，结果不变；无局部极大（单调上升 / 数组过短）时退回 argmax。
    """
    n = int(values.size)
    if n < 3:
        return int(np.argmax(values)) if n else 0
    vmax = float(np.max(values))
    if not np.isfinite(vmax) or vmax <= 0.0:
        return 0
    for i in range(1, n - 1):
        if (values[i] >= values[i - 1] and values[i] > values[i + 1]
                and values[i] >= rel_threshold * vmax):
            return i
    return int(np.argmax(values))


# ==========================================================================
# 主函数
# ==========================================================================
def simulate_head_branch(
    result: Any,
    params: HeadBranchParams | None = None,
) -> HeadVerdict:
    """在 ``BoulderFallResult`` 的时间网格上解头/颈 1-DOF 基础激励。

    Parameters
    ----------
    result
        ``climbing.pad.BoulderFallResult``（或任意 duck-typed 对象，
        需有 ``t_s``, ``z_torso_m``, ``meta['impact_window_s']``，
        以及可选的 ``accel_torso_g``）。不修改 ``result`` 的任何字段。
    params
        ``HeadBranchParams``。默认参数即可用（80 kg 人 / 标准颈椎刚度）。

    Returns
    -------
    HeadVerdict
        头/颈应力峰值及其时刻，以及 V3 / band 判读。

    Raises
    ------
    ValueError
        若 ``result.meta['impact_window_s']`` 缺失，或 ``t_s`` 与
        ``z_torso_m`` 长度不一致 / 长度 < 4（无法解 ODE）。

    Notes
    -----
    **严格的下游派生量** —— ``result`` 的数组不被修改，且本函数**不**
    把头部反作用力回灌到躯干的 ODE 中。这是 1D/2D 模型 bit-identicality
    的前提。
    """
    p = params or HeadBranchParams()

    # ---- 取输入时序 -----------------------------------------------------
    t = np.asarray(getattr(result, "t_s"), dtype=float)
    z_torso = np.asarray(getattr(result, "z_torso_m"), dtype=float)
    if t.shape != z_torso.shape:
        raise ValueError(
            f"result.t_s (len={t.size}) 与 result.z_torso_m "
            f"(len={z_torso.size}) 长度不一致；头颈分支需要两者等长。"
        )
    if t.size < 4:
        raise ValueError(
            f"result.t_s 长度仅 {t.size}（<4），无法解 1-DOF ODE。"
        )

    meta = getattr(result, "meta", None) or {}
    win = meta.get("impact_window_s")
    if win is None:
        # 镜像 src/climbing/bone.py:301-307 的同款守卫 —— 中文消息同 style。
        raise ValueError(
            "result.meta['impact_window_s'] 缺失：头颈分支必须在冲击窗口内"
            "定位峰值时刻，否则会被自由落体段 / 反弹重峰污染"
            "（参考 src/climbing/bone.py:301-307 同款守卫）。"
        )
    t0_impact = float(win[0])
    t1_impact = float(win[1])

    # ---- 躯干加速度（向下为正 = -accel_torso_g * g） --------------------
    # accel_torso_g 的口径来自 src/climbing/pad.py:614：
    #   a_torso = (f_flex/m_up - g) / g
    # "正 = 向上减速"（注释明示）。z 向下为正 ⇒ z_ddot (向下为正) = -a_torso_g*g。
    a_torso_g_raw = getattr(result, "accel_torso_g", None)
    if a_torso_g_raw is None:
        # 兜底：只给 z_torso_m 时，用**二阶**导数的插值求躯干加速度。
        # 注意必须是**二阶**：np.gradient(z_torso, t) 只是速度。
        # 用 CubicSpline 的二阶导（C² 连续）而不是连续两次 np.gradient ——
        # 后者在 5 ms 级的减速脉冲上会产生明显的数值噪声，把峰值力低估
        # 一个量级（实测 135 N vs 7542 N）；样条二阶导稳得多。
        from scipy.interpolate import CubicSpline

        cs = CubicSpline(t, z_torso)
        z_torso_ddot_arr = cs(t, 2)
        a_torso_g_arr = -z_torso_ddot_arr / p.g
    else:
        a_torso_g_arr = np.asarray(a_torso_g_raw, dtype=float)
        z_torso_ddot_arr = -a_torso_g_arr * p.g

    # ---- 单位换算：N/mm → N/m ------------------------------------------
    k_si = p.k_neck_n_per_mm * 1.0e3   # N/mm → N/m
    c_si = p.c_neck_n_s_per_mm * 1.0e3  # N·s/mm → N·s/m

    # ---- 1-DOF ODE：y = [delta, delta_dot] ------------------------------
    # PCHIP 插值：保留单调性，一阶导数连续 —— 比线性插值的尖角更适合
    # RK45 在基础激励曲率大的位置（脉冲峰值附近）保持大步长。
    try:
        z_torso_ddot_interp = PchipInterpolator(
            t, z_torso_ddot_arr, extrapolate=True,
        )
    except Exception:
        # 退化（极小数组 / 全常数）：线性插值兜底
        z_torso_ddot_interp = lambda x: np.interp(x, t, z_torso_ddot_arr)

    m_inv = 1.0 / p.m_head_kg
    k_over_m = k_si * m_inv
    c_over_m = c_si * m_inv

    def rhs(_t: float, y: np.ndarray) -> list[float]:
        delta, delta_dot = y
        # 头运动方程（向下为正）：
        #   m_head * z_head_ddot = m_head*g - k*delta - c*delta_dot
        # 在相对坐标下：
        #   delta_ddot = g - z_torso_ddot(t) - (k/m)*delta - (c/m)*delta_dot
        # 物理直观：自由落体段 z_torso_ddot = +g ⇒ delta_ddot = 0
        # （头与躯干共动，相对位移恒为 0）；冲击段 z_torso_ddot < 0
        # ⇒ delta_ddot > 0（头被"留在空中"），delta 上升，颈部被压缩。
        a_drive = float(z_torso_ddot_interp(_t))
        return [
            delta_dot,
            p.g - a_drive - k_over_m * delta - c_over_m * delta_dot,
        ]

    y0 = [0.0, 0.0]   # 自由落体初始：相对位移 0，相对速度 0

    sol = solve_ivp(
        rhs,
        (float(t[0]), float(t[-1])),
        y0,
        method="RK45",
        rtol=1e-8,
        atol=1e-11,
        # 默认颈椎自振周期 ~31 ms（k=200 N/mm, m=5 kg），必须用足够细的
        # 步长才能解析头部的惯性振荡。2e-4 s = 200 µs，约周期的 1/150。
        max_step=2.0e-4,
        dense_output=True,   # 让我们能在结果自己的 t_s 上精确采样
    )
    if not sol.success:
        raise RuntimeError(f"头颈分支 ODE 积分失败: {sol.message}")

    # ---- 采样回 result 的时间网格 --------------------------------------
    delta_t = sol.sol(t)[0]
    delta_dot_t = sol.sol(t)[1]
    # 颈部内力（向下为正的标量约定）。正值 = 颈部"被向下拉"
    # （头被向下拽，躯干被向上顶），即 head 在 torso 之下、弹簧压缩；
    # 负值 = 颈部"被向上拉"。取 |F_neck| 是物理应力。
    f_neck_t = k_si * delta_t + c_si * delta_dot_t
    f_mag = np.abs(f_neck_t)

    # ---- 在冲击窗口内取峰（避免自由落体段与反弹重峰污染） ---------------
    in_window = (t >= t0_impact) & (t <= t1_impact)
    if not in_window.any():
        raise ValueError(
            f"冲击窗口 [{t0_impact:.4f}, {t1_impact:.4f}] s 内没有有效数据点"
            f"（t_s 范围 [{float(t[0]):.4f}, {float(t[-1]):.4f}] s）。"
        )
    f_in = f_mag[in_window]
    t_in = t[in_window]
    # 首个显著局部极大（见 _first_significant_peak 的说明）：
    # 无阻尼时全局 argmax 落在多个等幅自由振荡峰之间是不确定的，
    # 取首个峰既是物理上的"首次应力集中"，也让 peak_time 可复现。
    idx_local = _first_significant_peak(f_in)
    t_peak = float(t_in[idx_local])
    peak_force_n = float(f_in[idx_local])
    delta_peak_m = float(delta_t[in_window][idx_local])

    # peak_time_s 距冲击起点 t0_impact（自由落体段被扣除）
    peak_time_s = t_peak - t0_impact
    within_y25 = peak_time_s <= V3_THRESHOLD_S

    # ---- 应力度量 (MPa) -------------------------------------------------
    # head_section_mm2: mm² → m²；F [N] / A [m²] = Pa → /1e6 = MPa
    a_section_m2 = p.head_section_mm2 * 1.0e-6
    peak_stress_mpa = peak_force_n / a_section_m2 / 1.0e6

    band = head_band(peak_stress_mpa)

    note = (
        f"delta_peak = {delta_peak_m*1e3:+.2f} mm；"
        f"F_neck_peak = {peak_force_n:.1f} N；"
        f"σ_peak = {peak_stress_mpa:.2f} MPa（A = {p.head_section_mm2:g} mm²）；"
        f"V3 {'通过' if within_y25 else '未通过'}"
        f"（阈值 ≤ {V3_THRESHOLD_S*1e3:.1f} ms，实测 {peak_time_s*1e3:.2f} ms）；"
        f"分档 [{band}]。"
    )

    return HeadVerdict(
        peak_stress_mpa=peak_stress_mpa,
        peak_time_s=peak_time_s,
        peak_force_n=peak_force_n,
        delta_peak_m=delta_peak_m,
        within_y25_window=within_y25,
        band=band,
        note=note,
    )
