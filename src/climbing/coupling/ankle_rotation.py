"""S6 —— 踝外旋（横截面 / 水平面，external rotation / SER）载荷模型（**opt-in，默认零**）。

为什么需要这一层
----------------
现有 OpenSim→FEBio 链路 + Phase-S5 ``ankle_supination`` 已经覆盖了**冠状面内翻**对
远端腓骨（外踝）的侧向弯曲 / 剪切风险。但 **trimalleolar 三踝骨折的支配机制不是
内翻**，而是 **external rotation / supination-external-rotation（SER）** —— 临床
文献 [Lauge-Hansen 1950; Yde 1980] 显示 SER 占踝骨折的 **57–85%**。

力学（Lauge-Hansen SER + 临床解剖）：
脚以胫骨长轴为轴**向外旋转** → 距骨在踝穴（mortise）内旋转 → 距骨外侧楔形面把
**腓骨远端（外踝）**沿水平面**扭转**（torsion about fibular long axis） +
**前后方向**弯曲 → **斜形 / 螺旋骨折**（oblique / spiral fibular fracture）。

临床/实验文献锚点（**measured**）：
* [Hirsch & Lewis, 1960s]：踝关节**能承受很大竖直压缩**，但**很小的旋转力矩即
  失效**。这是本模型只敢称"上界"的根本原因。
* [Funk et al. 2000]：轴向冲击 + <10° 旋转即可在外踝产生骨折。
* [Haraguchi & Armiger 2020; Ma et al. 2024]：距骨楔形 + 后踝受力的现代生物力学
  综述，支撑"距骨推动外踝"为本模型机制源头。

模型（全部标注为 **modeled / assumed**）
----------------------------------------
坐标系（OpenSim 全局系，y 向上）：``x`` 前后 / ``y`` 竖直 / ``z`` 左右。
外旋角 ``θ`` 是绕 **y 轴** 的足部旋转（**向外为正**）。

给定峰值 3D 地面反力 ``F⃗=(F_x, F_y, F_z)``、踝外旋角 ``θ``、CoP→腓骨长轴的
**水平面杠杆** ``d``（mm）：

1. **扭转力矩**（绕腓骨长轴 / 踝 y 轴，由外旋把竖直载荷横向化产生）：

       ``T_ext = F_y · d · sinθ``

   * 这是把 "GRF 沿胫骨轴线" 通过 CoP 到腓骨轴线的水平偏距 ``d`` 折算为绕
     腓骨长轴的扭矩，标准生物力学口径（与跑步生物力学 "GRF 对踝 y 轴的力矩"
     同源）。
   * ``θ=0`` 时 ``T_ext=0`` → 默认零行为。

2. **外踝 1D 扭转 / 弯曲应力**（远端腓骨外踝段，薄壁梁端部）：
   * 极惯性矩：对等效圆 ``J = 2·I``；``c`` = 极值纤维距离。
   * 扭转剪应力：``τ = T_ext · c / J``
   * 弯曲应力（外旋时距骨楔形面把外踝前后推；为 1D 上界**保守**地略去，纯扭转
     下 von Mises ≈ ``√3·τ``，偏大但属上限）：

       ``σ_torsion = T_ext · c / J``
       ``σ_vm      = sqrt(3) · σ_torsion``

3. **风险** = ``σ_vm / σ_c``，``σ_c = 70 MPa``（[Y25] `fibula_ends` 压缩强度，
   任务给定阈值；外旋主要走扭转 + 弯曲拉伸侧，但统一用同一材料阈值作上界）。

几何来源（**state the geometry source**）
-----------------------------------------
* ``A_section``：``results/opensim_fe/bc_robust_metric_route1.json`` 的 ``fibula_r``
  —— **measured**（THUMS AM50 CORT 实体 hex 的 V/L 截面）。
* ``I`` / ``c``：``results/opensim_fe/risk_1d_bending.json`` 的
  ``sections_equivalent_circle.fibula_r``（由 ``A_section`` 反推的等效圆，
  **modeled**），主口径。等效圆下 ``J = 2·I``，``c = r``。骨库中没有独立的
  "远端腓骨 / 外踝"截面，故用**中段皮质截面**作代理 —— 外踝略粗，用中段截面
  偏保守（高估扭转应力）。
* ``d``（杠杆）：默认 **assumed** ``30 mm``（CoP→胫骨长轴的水平偏距，量级
  30–50 mm；灵敏度覆盖 {15, 30, 45} mm）。

⚠️ 诚实边界：这是**上界式 1D 筛选**——忽略韧带 / 腱的卸载、应力集中、接触摩擦、
能量吸收与前后弯曲的复合作用，且用**材料强度**（非整体骨失效载荷）作阈值。
绝对值不可信，**只看相对与阈值是否跨越**。临床第一大伤是**韧带扭伤（71%）**，
本骨模型只能对应**骨折（27%）**。

单位 mm–N–MPa–s。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "FIBULA_ENDS_SIGMA_C_MPA",
    "DEFAULT_TIBIAL_LEVER_ARM_MM",
    "ExternalRotationLoad",
    "horizontal_torque_nmm",
    "external_rotation_load",
    "fibula_torsional_stress",
    "fibula_external_rotation_risk",
    "critical_external_rotation_angle_deg",
]

#: [Y25] `fibula_ends` 压缩材料强度（MPa）—— 任务给定的腓骨两端阈值；与
#: ``ankle_supination.FIBULA_ENDS_SIGMA_C_MPA`` 同源（70 MPa）。
FIBULA_ENDS_SIGMA_C_MPA: float = 70.0

#: 默认 CoP→胫骨长轴的水平偏距 (mm)。**assumed**：足部 CoP 落在距骨体上，距
#: 腓骨长轴的水平距离 ≈ 30 mm（足部半宽 50–70 mm 的中点偏内）。
DEFAULT_TIBIAL_LEVER_ARM_MM: float = 30.0


def _as_vec3(force_vec_n) -> np.ndarray:
    f = np.asarray(force_vec_n, dtype=float)
    if f.shape != (3,) or not np.all(np.isfinite(f)):
        raise ValueError(f"force_vec_n 必须是 3 个有限数，得到 {force_vec_n!r}")
    return f


def horizontal_torque_nmm(
    force_vec_n,
    external_rotation_deg: float,
    lever_arm_mm: float,
) -> float:
    """外旋扭矩 ``T_ext = |F_y| · |d| · |sinθ|`` (N·mm)，绕腓骨长轴（y 轴）。

    ``force_vec_n`` 是 3D 力 (N)，``external_rotation_deg`` 是踝外旋角 (度)，
    ``lever_arm_mm`` 是 CoP→胫骨长轴的水平偏距 (mm)。返回**非负**绝对值（与
    :func:`ankle_supination.inversion_moment_nmm` 同口径）。
    """
    f = _as_vec3(force_vec_n)
    th = np.deg2rad(float(external_rotation_deg))
    return abs(float(f[1])) * abs(float(lever_arm_mm)) * abs(float(np.sin(th)))


@dataclass(frozen=True)
class ExternalRotationLoad:
    """一次外旋载荷解析（峰值帧）。力 (N)，力矩 (N·mm)，力臂 (mm)，角 (度)。"""

    external_rotation_deg: float
    lever_arm_mm: float
    force_vec_n: tuple[float, float, float]
    t_ext_nmm: float
    t_ext_nm: float

    @property
    def is_zero(self) -> bool:
        """默认口径判定：无外旋角 → 无扭矩（等价于旧轴向行为）。"""
        return abs(self.t_ext_nmm) < 1e-12


def external_rotation_load(
    force_vec_n,
    external_rotation_deg: float = 0.0,
    lever_arm_mm: float = DEFAULT_TIBIAL_LEVER_ARM_MM,
) -> ExternalRotationLoad:
    """由 3D GRF 矢量 + 外旋角 + 力臂 → :class:`ExternalRotationLoad`。

    默认 ``external_rotation_deg=0.0`` ⇒ ``T_ext=0``（**旧行为**）。
    """
    f = _as_vec3(force_vec_n)
    t = horizontal_torque_nmm(f, external_rotation_deg, lever_arm_mm)
    return ExternalRotationLoad(
        external_rotation_deg=float(external_rotation_deg),
        lever_arm_mm=float(lever_arm_mm),
        force_vec_n=(float(f[0]), float(f[1]), float(f[2])),
        t_ext_nmm=float(t),
        t_ext_nm=float(t) / 1000.0,
    )


def fibula_torsional_stress(
    t_ext_nmm: float,
    *,
    i_mm4: float,
    c_mm: float,
) -> dict:
    """外踝 1D 扭转 / 弯曲应力：等效圆 ``J=2·I``，``τ=T·c/J``，von Mises = ``√3·τ``。

    返回 ``t_ext_nmm`` / ``sigma_torsion_mpa`` / ``sigma_shear_mpa`` /
    ``sigma_vonmises_mpa``。``t_ext_nmm≤0`` 时全部归 0。
    """
    t = abs(float(t_ext_nmm))
    i = float(i_mm4)
    c = float(c_mm)
    if i <= 0.0 or c < 0.0:
        raise ValueError(f"截面参数非法：I={i}, c={c}")
    if t <= 0.0:
        return {
            "t_ext_nmm": 0.0,
            "j_mm4": float(2.0 * i),
            "sigma_torsion_mpa": 0.0,
            "sigma_shear_mpa": 0.0,
            "sigma_vonmises_mpa": 0.0,
        }
    j = 2.0 * i  # 等效圆：J = πr⁴/2 = 2·I
    tau = t * c / j
    sigma_vm = float(np.sqrt(3.0) * tau)
    return {
        "t_ext_nmm": float(t),
        "j_mm4": float(j),
        "sigma_torsion_mpa": float(tau),
        "sigma_shear_mpa": float(tau),  # τ = T·c/J 是剪应力，等价命名
        "sigma_vonmises_mpa": float(sigma_vm),
    }


def fibula_external_rotation_risk(
    t_ext_nmm: float,
    *,
    i_mm4: float,
    c_mm: float,
    sigma_c_mpa: float = FIBULA_ENDS_SIGMA_C_MPA,
) -> dict:
    """腓骨远端外旋风险 = ``σ_vm / σ_c``（``σ_c`` 默认 70 MPa = [Y25] fibula_ends）。

    ``t_ext_nmm=0`` ⇒ risk=0 ⇒ **逐位复现**今天默认行为。
    """
    sig = fibula_torsional_stress(t_ext_nmm, i_mm4=i_mm4, c_mm=c_mm)
    sc = float(sigma_c_mpa)
    if sc <= 0.0:
        raise ValueError(f"sigma_c_mpa 必须为正，得到 {sigma_c_mpa!r}")
    out = dict(sig)
    out.update(
        sigma_c_mpa=sc,
        risk=float(sig["sigma_vonmises_mpa"] / sc),
        fracture=bool(sig["sigma_vonmises_mpa"] / sc >= 1.0),
    )
    return out


def critical_external_rotation_angle_deg(
    force_vec_n,
    *,
    i_mm4: float,
    c_mm: float,
    lever_arm_mm: float = DEFAULT_TIBIAL_LEVER_ARM_MM,
    sigma_c_mpa: float = FIBULA_ENDS_SIGMA_C_MPA,
) -> float | None:
    """使外旋风险首次达到 1 的最小外旋角 (度)；若 ``[0, 60]°`` 内不越阈则返回 ``None``。

    单调性：``T_ext(θ)=|F|·|d|·|sinθ|`` 在 ``θ∈[0, π/2]`` 单调增，故可用二分。
    若 ``θ=0`` 已越阈返回 ``0.0``。
    """
    f = _as_vec3(force_vec_n)

    def _risk(theta_deg: float) -> float:
        t = horizontal_torque_nmm(f, theta_deg, lever_arm_mm)
        return fibula_external_rotation_risk(
            t, i_mm4=i_mm4, c_mm=c_mm, sigma_c_mpa=sigma_c_mpa,
        )["risk"]

    lo, hi = 0.0, 60.0
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


# --------------------------------------------------------------------------
# 自检：python -m climbing.coupling.ankle_rotation
# --------------------------------------------------------------------------
def _main() -> int:
    print("=== ankle_rotation (external / SER) self-check ===")
    # 默认：θ=0 → risk=0（逐位复现旧行为）
    z = external_rotation_load((0.0, 21000.0, 0.0))
    print(f"θ=0  T_ext={z.t_ext_nmm:.3e} N·mm  is_zero={z.is_zero}")
    r = fibula_external_rotation_risk(
        z.t_ext_nmm, i_mm4=826.8918876962899, c_mm=5.696258799381752,
    )
    print(f"      σ_vm={r['sigma_vonmises_mpa']:.3e} MPa  risk={r['risk']:.3e}  "
          f"fracture={r['fracture']}")
    for th in (5.0, 10.0, 15.0, 30.0):
        z2 = external_rotation_load((0.0, 21000.0, 0.0), external_rotation_deg=th)
        r2 = fibula_external_rotation_risk(
            z2.t_ext_nmm, i_mm4=826.8918876962899, c_mm=5.696258799381752,
        )
        print(f"θ={th:5.1f}°  T_ext={z2.t_ext_nmm/1000.0:7.2f} N·m  "
              f"σ_vm={r2['sigma_vonmises_mpa']:8.2f} MPa  risk={r2['risk']:6.2f}  "
              f"fracture={r2['fracture']}")
    crit = critical_external_rotation_angle_deg(
        (0.0, 21000.0, 0.0),
        i_mm4=826.8918876962899, c_mm=5.696258799381752,
    )
    print(f"critical_external_rotation_angle = {crit}")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
