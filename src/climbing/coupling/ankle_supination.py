"""S5 —— 踝旋后 / 内翻（冠状面）载荷模型（**opt-in，默认零**）。

为什么需要这一层
----------------
现有 OpenSim→FEBio 链路只把地面反力当**纯轴向竖直力**处理（``subtalar_dir=(0,-1,0)``）。
但两篇临床论文指出，抱石踝伤的**第一机制是冠状面内的旋后 / 内翻**：

    [Heck2024 / Müller2022, §5.2.3] "The soft surface initially allows the foot to
    sink in, accompanied by a **supination movement of the ankle joint**. Once the
    modulus of elasticity of the soft-floor mat is reached, the weight forces on the
    supinated joint increase until ligamentous structures or bones give way."

    [Beurienne2025, Discussion] "Increasing pad rigidity could prevent excessive
    ankle supination, thereby reducing force on ligament structures ..."

即：脚沉入垫 → 踝**旋后（invert）**→ 达垫子模量后载荷上升 → 韧带 / 骨失效。
这条轴对应[Y25]里**腓骨异常敏感**、以及临床 **踝骨折 27.4%** 的现实。
本模块把这条"被写死的轴向"补成**冠状面 1D 载荷**，作为纯轴向口径的 opt-in 扩展。

模型（全部标注为 **modeled / assumed**）
----------------------------------------
坐标系（OpenSim 全局系，y 向上）：

* ``x`` —— 前后（矢状面），≈ 距下轴的投影方向；
* ``y`` —— 竖直（up）；
* ``z`` —— 左右（冠状面），**旋后 / 内翻力矩轴为 x**。

给定峰值 3D 地面反力矢量 ``F⃗=(F_x,F_y,F_z)``、踝旋后角 ``β``、趾下 CoP 到距下轴的
力臂 ``d``：

1. **横向力**（冠状面内、垂直于旋后后的距下轴的分量）：

       ``F_lat = F_y·sinβ + F_z·cosβ``

   * ``β=0`` 时退化为 ``F_z``（纯冠状地面倾斜产生的横向剪切）；
   * ``F_z=0``（地面水平）时退化为 ``F_y·sinβ``（竖直载荷被旋后角重新指向横向）。

2. **内翻力矩**（绕距下轴 / 踝关节 x 轴）：

       ``M_inv = F_lat · d``

   ``d`` = CoP 到距下轴的力臂（mm）。这是把 ``F_lat`` 以力臂 ``d`` 折算到关节轴的
   标准力臂口径（与跑步生物力学里"GRF 对距下轴的力矩臂"同源）。

3. **腓骨远端（外踝）1D 侧向应力** —— 弯曲来自 ``M_inv``，剪切来自 ``F_lat``：

       ``σ_bend = M_inv · c / I``        （梁弯曲，c=极值纤维距离）
       ``τ_shear = F_lat / A_section``   （名义平均剪应力）
       ``σ_lat  = σ_bend + τ_shear``     （保守线性叠加，主口径）
       ``σ_vm   = sqrt(σ_bend² + 3·τ_shear²)``（von Mises，交叉校核）

4. **风险** = ``σ_lat / σ_c``，``σ_c = 70 MPa``（[Y25] `fibula_ends` 压缩强度，
   任务给定阈值）。

几何来源（**state the geometry source**）
-----------------------------------------
* ``A_section``：``results/opensim_fe/bc_robust_metric_route1.json`` 的 ``fibula_r``
  —— **measured**（THUMS AM50 CORT 实体 hex 的 V/L 截面）。
* ``I`` / ``c``：``results/opensim_fe/risk_1d_bending.json`` 的
  ``sections_real_midshaft.fibula_r``（干净 CORT remesh 中段皮质截面的弱轴二阶矩与
  极值纤维距离，**measured**）或 ``sections_equivalent_circle.fibula_r``（由
  ``A_section`` 反推的等效圆，**modeled**）。**骨库中没有独立的"远端腓骨 / 外踝"截面**，
  故用**中段皮质截面**作代理 —— 外踝略粗，用中段截面偏保守（高估风险）。
* ``d``（力臂）：默认 **assumed** ``30 mm``（跟骨 / 后足横向量级 50–70 mm 的半宽，
  实测 bbox 见报告），灵敏度覆盖 {15, 30, 45} mm。

⚠️ 诚实边界：这是**上界式 1D 筛选**——忽略韧带 / 腱的卸载、应力集中、接触摩擦与
能量吸收，且用**材料强度**（非整体骨失效载荷）作阈值。绝对值不可信，**只看相对与阈值
是否跨越**。临床第一大伤是**韧带扭伤（71%）**，本骨模型只能对应**骨折（27%）**。

单位 mm–N–MPa–s。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "FIBULA_ENDS_SIGMA_C_MPA",
    "DEFAULT_SUBTALAR_LEVER_ARM_MM",
    "SupinationLoad",
    "frontal_lateral_force_n",
    "inversion_moment_nmm",
    "supination_load",
    "fibula_lateral_stress",
    "fibula_lateral_risk",
    "critical_inversion_angle_deg",
]

#: [Y25] `fibula_ends` 压缩材料强度（MPa）—— 任务给定的腓骨两端阈值，另见 `bone.py`。
FIBULA_ENDS_SIGMA_C_MPA: float = 70.0

#: 默认 CoP→距下轴力臂 (mm)。**assumed**：后足横向半宽量级（跟骨 bbox 50–70 mm）。
DEFAULT_SUBTALAR_LEVER_ARM_MM: float = 30.0


def _as_vec3(force_vec_n) -> np.ndarray:
    f = np.asarray(force_vec_n, dtype=float)
    if f.shape != (3,) or not np.all(np.isfinite(f)):
        raise ValueError(f"force_vec_n 必须是 3 个有限数，得到 {force_vec_n!r}")
    return f


def frontal_lateral_force_n(force_vec_n, supination_deg: float) -> float:
    """冠状面横向力 ``F_lat = F_y·sinβ + F_z·cosβ`` (N)。

    ``force_vec_n`` 是全局系 3D 力 (N)，``supination_deg`` 是踝旋后 / 内翻角 (度)。
    ``β=0`` 时退化为全局横向分量 ``F_z``；纯轴向默认（``F_z=0``）在 ``β=0`` 时给 **0**。
    返回**带符号**分量（内翻 / 外翻方向由符号区分）；判风险时取绝对值。
    """
    f = _as_vec3(force_vec_n)
    b = np.deg2rad(float(supination_deg))
    return float(f[1] * np.sin(b) + f[2] * np.cos(b))


def inversion_moment_nmm(f_lat_n: float, lever_arm_mm: float) -> float:
    """内翻力矩 ``M_inv = |F_lat| · |d|`` (N·mm)，绕距下 / 踝 x 轴。"""
    return abs(float(f_lat_n)) * abs(float(lever_arm_mm))


@dataclass(frozen=True)
class SupinationLoad:
    """一次旋后载荷解析（峰值帧）。力 (N)，力矩 (N·mm)，力臂 (mm)，角 (度)。"""

    supination_deg: float
    lever_arm_mm: float
    force_vec_n: tuple[float, float, float]
    f_lat_n: float
    m_inv_nmm: float
    m_inv_nm: float

    @property
    def is_zero(self) -> bool:
        """默认口径判定：无横向力 → 无内翻力矩（等价于旧轴向行为）。"""
        return abs(self.f_lat_n) < 1e-12 and abs(self.m_inv_nmm) < 1e-12


def supination_load(
    force_vec_n,
    supination_deg: float = 0.0,
    lever_arm_mm: float = DEFAULT_SUBTALAR_LEVER_ARM_MM,
) -> SupinationLoad:
    """由 3D GRF 矢量 + 旋后角 + 力臂 → :class:`SupinationLoad`。

    默认 ``supination_deg=0.0`` 且若 ``F_z=0`` ⇒ ``F_lat=M_inv=0``（**旧行为**）。
    """
    f = _as_vec3(force_vec_n)
    fl = frontal_lateral_force_n(f, supination_deg)
    m = inversion_moment_nmm(fl, lever_arm_mm)
    return SupinationLoad(
        supination_deg=float(supination_deg),
        lever_arm_mm=float(lever_arm_mm),
        force_vec_n=(float(f[0]), float(f[1]), float(f[2])),
        f_lat_n=float(fl),
        m_inv_nmm=float(m),
        m_inv_nm=float(m) / 1000.0,
    )


def fibula_lateral_stress(
    f_lat_n: float,
    m_inv_nmm: float,
    *,
    a_section_mm2: float,
    i_mm4: float,
    c_mm: float,
) -> dict:
    """腓骨远端 1D 侧向应力：弯曲（``M_inv``）+ 剪切（``F_lat``）。

    返回 ``sigma_bend_mpa`` / ``sigma_shear_mpa`` / ``sigma_lateral_mpa``
    （= 弯曲 + 剪切，主口径）/ ``sigma_vonmises_mpa``。
    """
    f = abs(float(f_lat_n))
    m = abs(float(m_inv_nmm))
    a = float(a_section_mm2)
    i = float(i_mm4)
    c = float(c_mm)
    if a <= 0.0 or i <= 0.0 or c < 0.0:
        raise ValueError(f"截面参数非法：A={a}, I={i}, c={c}")
    sig_bend = m * c / i
    tau = f / a
    return {
        "f_lat_n": f,
        "m_inv_nmm": m,
        "sigma_bend_mpa": float(sig_bend),
        "sigma_shear_mpa": float(tau),
        "sigma_lateral_mpa": float(sig_bend + tau),
        "sigma_vonmises_mpa": float(np.sqrt(sig_bend ** 2 + 3.0 * tau ** 2)),
    }


def fibula_lateral_risk(
    f_lat_n: float,
    m_inv_nmm: float,
    *,
    a_section_mm2: float,
    i_mm4: float,
    c_mm: float,
    sigma_c_mpa: float = FIBULA_ENDS_SIGMA_C_MPA,
) -> dict:
    """腓骨远端侧向风险 = ``σ_lat / σ_c``（``σ_c`` 默认 70 MPa = [Y25] fibula_ends）。"""
    sig = fibula_lateral_stress(
        f_lat_n, m_inv_nmm, a_section_mm2=a_section_mm2, i_mm4=i_mm4, c_mm=c_mm
    )
    sc = float(sigma_c_mpa)
    if sc <= 0.0:
        raise ValueError(f"sigma_c_mpa 必须为正，得到 {sigma_c_mpa!r}")
    out = dict(sig)
    out.update(
        sigma_c_mpa=sc,
        risk=float(sig["sigma_lateral_mpa"] / sc),
        risk_vonmises=float(sig["sigma_vonmises_mpa"] / sc),
        fracture=bool(sig["sigma_lateral_mpa"] / sc >= 1.0),
    )
    return out


def critical_inversion_angle_deg(
    force_vec_n,
    *,
    a_section_mm2: float,
    i_mm4: float,
    c_mm: float,
    lever_arm_mm: float = DEFAULT_SUBTALAR_LEVER_ARM_MM,
    sigma_c_mpa: float = FIBULA_ENDS_SIGMA_C_MPA,
) -> float | None:
    """使侧向风险首次达到 1 的最小旋后角 (度)；若 [0,45]° 内不越阈则返回 ``None``。

    单调性：``σ_lat(β)`` 在 ``β∈[0,45]`` 单调增（``F_lat`` 单调增），故可用二分。
    若 ``β=0`` 已越阈（地面横向倾斜很大）返回 ``0.0``。
    """
    f = _as_vec3(force_vec_n)

    def _risk(b_deg: float) -> float:
        fl = frontal_lateral_force_n(f, b_deg)
        m = inversion_moment_nmm(fl, lever_arm_mm)
        return fibula_lateral_risk(
            fl, m, a_section_mm2=a_section_mm2, i_mm4=i_mm4, c_mm=c_mm,
            sigma_c_mpa=sigma_c_mpa,
        )["risk"]

    lo, hi = 0.0, 45.0
    if _risk(lo) >= 1.0:
        return 0.0
    if _risk(hi) < 1.0:
        return None
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if _risk(mid) >= 1.0:
            hi = mid
        else:
            lo = mid
    return float(hi)
