"""S1.2 —— OpenSim 正动力学 "dead-drop"：把 GRF 注入 Rajagopal2015 并积分。

方案依据（``docs/OpenSim_FE复现方案.md``）
-----------------------------------------
* §3.1 D1：正动力学 + ExternalForce（**不用**逆动力学）；
* §3.3 D3：肌肉 dead-drop（关激活动力学、激活=0，仅被动）。

为什么用 :class:`PrescribedForce` 而不是 :class:`ExternalForce`
-------------------------------------------------------------
OpenSim 4.6 的 Python 绑定**没有** `ExternalForce.setForceFunction`，
时变外力无法通过它施加。``PrescribedForce`` 支持把力/作用点写成
**时间的函数**（`setForceFunctions` / `setPointFunctions`），正好用于
已知 GRF 时程。力用**全局系**（`setForceIsInGlobalFrame(True)`），
作用点用**体内坐标系**（`setPointIsInGlobalFrame(False)`，取 calcn 原点）。

作用点约定
----------
``calcn_r`` 的体坐标系原点 = 距下关节中心（``subtalar_r`` 的 child 位置
为 `0 0 0`）。所以把 GRF 施加在 calcn 原点 = 施加在距下关节中心，
**对该关节无力矩**。S1 是"直立轴向"场景（[Y25]），忽略足底 CoP 产生的
踝力矩是可接受的；力矩通路留到 S3 的 ``load_transfer``。

已知限制
--------
* 数据源是 ``pad.py`` 的 1D 等效双质点 GRF，套到 22 刚体模型上是
  **载荷边界条件的近似**（方案 §3.1）；S1 只验证管线。
* 默认模型未缩放（Rajagopal 通用 75.34 kg），方案 §八已列为风险；
  需要 [Y25] AM50 (77 kg) 口径时用 ``run_dead_drop(grf, scale_mass_kg=77.0)``
  （见 :func:`scale_model_mass`；默认 ``None`` 时行为不变）。

⚠️ 注意：``Model.getWorkingState()`` / ``StatesTrajectory`` 的实现细节
以 4.6 绑定为准，见 ``_state_from_table_row``。
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import opensim as os

from .opensim_grf import GroundReaction

__all__ = [
    "FallResult",
    "LegPose",
    "PassiveStiffness",
    "run_dead_drop",
    "scale_model_mass",
    "model_total_mass",
    "DEFAULT_MODEL",
    "FOOT_BODIES",
    "PELVIS_TY_SPEED_STATE",
    "PELVIS_TX_SPEED_STATE",
    "PELVIS_TZ_SPEED_STATE",
]

ROOT = Path(__file__).resolve().parents[3]
DEFAULT_MODEL = ROOT / "model" / "opensim" / "FullBodyModel-4.0" / "Rajagopal2015.osim"

#: 左右足（承受 GRF 的刚体）
FOOT_BODIES: tuple[str, str] = ("calcn_r", "calcn_l")

#: 根平移坐标的速度状态（向下 = 负 y）
PELVIS_TY_SPEED_STATE = "/jointset/ground_pelvis/pelvis_ty/speed"
#: 根水平平移速度状态（L2：初始速度加水平分量，默认 0 = 纯竖直）
PELVIS_TX_SPEED_STATE = "/jointset/ground_pelvis/pelvis_tx/speed"
PELVIS_TZ_SPEED_STATE = "/jointset/ground_pelvis/pelvis_tz/speed"

#: 下肢关节坐标。S1 把它们**锁死**，以近似 [Y25] 的"直立、腿绷直"坠落。
#:
#: 为什么必须锁：OpenSim 的关节是**理想无摩擦铰链，零被动刚度**。dead-drop
#: 又关闭了肌肉，于是 26 kN 的足底冲击会把小腿"甩"起来 —— 实测踝关节
#: 冲到 −4.4 rad（范围只有 [−0.7, 0.52]），膝关节 −1.55 rad，整个模型四肢
#: 乱飞。这不是物理，是"零刚度关节 + 无肌肉"的模型缺陷。
#: [Y25] 的直立绷直腿在 ~20 ms 冲击内近似刚性柱体，故 S1 先用**刚性腿**
#: 打通管线；关节顺应性/肌肉（两路径机制）留到 S3。
LEG_JOINTS: tuple[str, ...] = (
    "ankle_angle_r", "ankle_angle_l",
    "knee_angle_r", "knee_angle_l",
    "hip_flexion_r", "hip_adduction_r", "hip_rotation_r",
    "hip_flexion_l", "hip_adduction_l", "hip_rotation_l",
)

#: 仅参与主矢状面"屈膝 / 屈髋 / 踝背屈"的三个下肢关节坐标（S3 姿势接口用）。
#: ``knee_angle_r`` 与 ``knee_angle_l`` 对称；只设一侧时另一侧自动镜像。
LEG_FLEX_JOINTS: tuple[str, ...] = (
    "ankle_angle_r", "ankle_angle_l",
    "knee_angle_r", "knee_angle_l",
    "hip_flexion_r", "hip_flexion_l",
)

_G = 9.80665  # OpenSim 默认重力大小


@dataclass(frozen=True)
class LegPose:
    """S3 —— OPT-IN 目标下肢姿势 (rad, 关节坐标系)。

    Path (a) **POSE-LOCKED POSTURE** 的输入：把 ``LEG_FLEX_JOINTS`` 的默认
    值改为 ``self`` 给出的值，再 ``set_locked(True)``。锁的是**目标姿势**而不是
    直立（默认 0）；力学上等价于"刚性腿 + 弯曲几何"。

    **默认值** = 直立、腿绷直，与 ``LEG_JOINTS`` 的旧锁定路径**逐位一致**。

    Parameters
    ----------
    knee_flexion_rad
        单膝屈曲 (rad, 0 = 直立；≈1.05 = 60°)。双膝对称施加。
    hip_flexion_rad
        单侧髋屈 (rad, 0 = 直立；≈0.52 = 30°)。双髋对称施加。
    ankle_rad
        单侧踝背屈/跖屈 (rad, 0 = 中立；>0 = 背屈)。模型踝关节范围
        ``[-0.698, 0.524]`` rad = ``[-40°, 30°]``。
    """

    knee_flexion_rad: float = 0.0
    hip_flexion_rad: float = 0.0
    ankle_rad: float = 0.0
    #: **可选**躯干前倾 (rad, 作用在 ``lumbar_extension``；负 = 前倾/屈，正 = 后仰)。
    #: 默认 0 = **不动躯干**（保持改前"躯干自由"的动力学，逐位一致）。
    torso_lean_rad: float = 0.0

    @classmethod
    def flexed_landing(
        cls,
        knee_deg: float = 30.0,
        *,
        hip_deg: float = 0.0,
        ankle_deg: float = 0.0,
        torso_lean_deg: float = 0.0,
    ) -> "LegPose":
        """便捷构造：``flexed_landing(knee_deg=30)`` → 膝屈 30°，余 0。

        ``torso_lean_deg`` 非 0 时**额外**把 ``lumbar_extension`` 锁在该角度
        （正 = 后仰 / 负 = 前倾）。
        """
        return cls(
            knee_flexion_rad=math.radians(knee_deg),
            hip_flexion_rad=math.radians(hip_deg),
            ankle_rad=math.radians(ankle_deg),
            torso_lean_rad=math.radians(torso_lean_deg),
        )

    def as_default_value_dict(self) -> dict[str, float]:
        """把姿势展开成 ``{coord_name: value}``，双侧对称。"""
        return {
            "ankle_angle_r": self.ankle_rad,
            "ankle_angle_l": self.ankle_rad,
            "knee_angle_r": self.knee_flexion_rad,
            "knee_angle_l": self.knee_flexion_rad,
            "hip_flexion_r": self.hip_flexion_rad,
            "hip_flexion_l": self.hip_flexion_rad,
        }


@dataclass(frozen=True)
class PassiveStiffness:
    """S3 —— OPT-IN 下肢被动刚度 / 阻尼 (Path b: UNLOCK + 被动刚度)。

    把 ``LEG_FLEX_JOINTS`` 的默认锁定**解除**后，必须补上**被动抗弯刚度**，
    否则 26 kN 的足底冲击会让关节冲到范围外（实测踝 −4.4 rad，超范围）
    —— 模型不是物理飞，而是**零刚度铰链 + 无肌肉**的缺陷。

    这里同时挂两种 OpenSim 力元件：
    * ``SpringGeneralizedForce`` —— **线性**扭转弹簧 + 阻尼，绕 ``rest_length``
      （= 起始姿势）抗弯。模拟"组织在中立位的本构抗力"。
    * ``CoordinateLimitForce`` —— **非线性**边壁刚度（``upper`` / ``lower``
      stiffness × C²-continuous transition），限制坐标不冲出解剖范围。

    **作用范围**：只解锁 :data:`LEG_FLEX_JOINTS`（矢状面的踝 / 膝 / 髋屈曲）；
    髋内收 / 旋转（非矢状面，无生理被动弹簧标定）**仍锁死**，避免侧向失控。

    **所有数值默认 = 0（完全关闭）**。旧 S1/S2 路径不传此参数时，模型行为
    **逐位不变**。

    Units
    -----
    * ``*_stiffness`` → N·m/rad（扭转弹簧 / 边壁刚度）
    * ``*_damping``   → N·m·s/rad（粘性阻尼）
    * ``*_limit_*_rad`` → rad（坐标范围；``None`` 表示用模型 ``range_min/max``）
    * ``transition_rad`` → rad（C²-continuous 过渡区半宽；≈3°）

    Caveats / Data sources
    ----------------------
    **以下数值均为 ASSUMED / MODEL，不是实测数据**：
    * 中立位线性刚度参照散打文献量级（~5–10 N·m/deg in mid-range；
      ~25–50 N·m/deg near end-range；参见 Stein 等 1996, Silder 2008）。
    * 边壁刚度足以在 0.5 ms 内阻止越界（实际值按"无过分形变"标定）。
    * 没有标定到具体受试者 / 没有肌肉反射 / 没有任何主动控制。
    若后续有实测数据，**应**改用真实曲线并去掉 ``ASSUMED`` 标签。
    """

    # 中立位线性扭转弹簧（N·m/rad）：膝 / 髋 / 踝
    knee_stiffness: float = 0.0
    knee_damping: float = 0.0
    hip_stiffness: float = 0.0
    hip_damping: float = 0.0
    ankle_stiffness: float = 0.0
    ankle_damping: float = 0.0

    # 边壁刚度（N·m/rad）+ 阻尼（N·m·s/rad）
    knee_limit_stiffness: float = 0.0
    knee_limit_damping: float = 0.0
    hip_limit_stiffness: float = 0.0
    hip_limit_damping: float = 0.0
    ankle_limit_stiffness: float = 0.0
    ankle_limit_damping: float = 0.0

    # 边壁坐标（rad）。``None`` = 用 ``Coordinate.getRangeMin/Max``（模型自带范围）
    knee_limit_lower_rad: float | None = None
    knee_limit_upper_rad: float | None = None
    hip_limit_lower_rad: float | None = None
    hip_limit_upper_rad: float | None = None
    ankle_limit_lower_rad: float | None = None
    ankle_limit_upper_rad: float | None = None

    # C²-continuous 过渡区半宽（rad）。≈3° = 0.0524 rad；太大积分器慢，
    # 太小接触不光滑 → 力跳变。OpenSim 文档建议"越大越光滑"。
    transition_rad: float = 0.05

    @classmethod
    def nominal_landing(cls) -> "PassiveStiffness":
        """一组**数值稳定**的 "标定起点"（全部 ASSUMED，非实测）。

        * 膝/髋/踝中立位线性刚度 = 500 N·m/rad（≈ 8.7 N·m/deg；近端范围量级）
        * 阻尼 = 500 N·m·s/rad（重阻尼，避免在 20 ms 冲击下激发欠阻尼振铃）
        * 边壁刚度 = 20000 N·m/rad（硬壁，约 0.1 ms 内阻止 1 rad 越界）
        * 边壁阻尼 = 200 N·m·s/rad
        * 过渡区 = 0.01 rad (≈ 0.57°)

        **为什么用这组偏硬的数值**（诚实说明）
        -----------------------------------------
        实测发现：**软弹簧 + 轻阻尼**（如 k=50 N·m/rad, c=5 N·m·s/rad，虽是
        文献量级的"生理中段被动刚度"）在这 20 ms 冲击下**欠阻尼振铃**，数值上把
        髋/腰椎峰值放大 2–5 倍——那是**数值伪影**，不是物理。这组偏硬、重阻尼的
        数值能给出**稳定、无振铃**的解（关节仍只弯 ~2°，见 S3 报告），从而干净地
        回答"被动刚度能不能吸能"这个问题。

        ⚠️ **本组数值不是受试者标定值**；``scripts/opensim_fe/nonvertical_s3_posture.py``
        的 §4.4 敏感性扫描会显式展示"软 → 振铃"与"硬 → 稳定"两个 regime。
        """
        return cls(
            knee_stiffness=500.0, knee_damping=500.0,
            hip_stiffness=500.0, hip_damping=500.0,
            ankle_stiffness=500.0, ankle_damping=500.0,
            knee_limit_stiffness=20000.0, knee_limit_damping=200.0,
            hip_limit_stiffness=20000.0, hip_limit_damping=200.0,
            ankle_limit_stiffness=20000.0, ankle_limit_damping=200.0,
            transition_rad=0.01,
        )

    @classmethod
    def soft_landing(cls) -> "PassiveStiffness":
        """**文献量级的中段被动刚度**（ASSUMED）：k=50 N·m/rad、c=5 N·m·s/rad。

        保留它作为**反例**：S3 报告用它证明"生理中段刚度 → 欠阻尼振铃 → 数值放大"。
        这不是失败模式，而是**该模型在短促冲击下的固有局限**（见报告 §5）。
        """
        return cls(
            knee_stiffness=50.0, knee_damping=5.0,
            hip_stiffness=50.0, hip_damping=5.0,
            ankle_stiffness=50.0, ankle_damping=5.0,
            knee_limit_stiffness=2000.0, knee_limit_damping=50.0,
            hip_limit_stiffness=2000.0, hip_limit_damping=50.0,
            ankle_limit_stiffness=2000.0, ankle_limit_damping=50.0,
            transition_rad=0.05,
        )

    def is_active(self) -> bool:
        """是否真有任何一个刚度 / 阻尼 > 0（决定要不要走 unlock 路径）。"""
        for v in (
            self.knee_stiffness, self.knee_damping,
            self.hip_stiffness, self.hip_damping,
            self.ankle_stiffness, self.ankle_damping,
            self.knee_limit_stiffness, self.knee_limit_damping,
            self.hip_limit_stiffness, self.hip_limit_damping,
            self.ankle_limit_stiffness, self.ankle_limit_damping,
        ):
            if float(v) > 0.0:
                return True
        return False


@dataclass
class FallResult:
    """一次 dead-drop 正动力学的产物。"""

    model: "os.Model"
    states: "os.TimeSeriesTable"
    t_s: np.ndarray
    peak_applied_n: float
    applied_total_n: np.ndarray
    v0_ms: float
    final_time_s: float
    model_path: Path
    #: 积分时模型的实际总质量 (kg)；未缩放时 = 75.337，缩放后 = 目标值。
    model_mass_kg: float = 0.0
    #: 质量缩放因子 f = target_kg / 原始总质量；未缩放 = 1.0。
    mass_scale: float = 1.0

    @property
    def n_states(self) -> int:
        return int(self.states.getNumRows())


def model_total_mass(model: "os.Model") -> float:
    """返回模型全部 body 的质量之和 (kg)（不依赖 state，加载后即可调用）。"""
    bodies = model.getBodySet()
    return float(sum(bodies.get(i).getMass() for i in range(bodies.getSize())))


def scale_model_mass(model: "os.Model", target_kg: float = 77.0) -> float:
    """把模型各 body 的质量/惯量按比例缩放到目标总质量。

    单参数 (uniform) 缩放：每个 body 的质量乘 ``f = target_kg / m_total``，
    其**对质心的惯量**同样乘 ``f``（固定几何下密度均匀缩放 ⇒ 惯量与质量
    成正比），质心位置不变。这正对应 [Y25] 的 AM50 (77 kg, 1.75 m) 与
    Rajagopal2015 通用模型 (75.337 kg) 之间的"同一身材不同总质量"折算。

    Parameters
    ----------
    model
        已 ``os.Model(...)`` 载入但**尚未** ``initSystem()`` 的模型（就地修改）。
    target_kg
        目标总质量 (kg)。必须为正有限值。

    Returns
    -------
    float
        实际施加的缩放因子 ``f``。

    Raises
    ------
    ValueError
        原始总质量非正/非有限，或 ``target_kg`` 非正/非有限。

    Notes
    -----
    **默认不调用**：``run_dead_drop`` 的 ``scale_mass_kg=None`` 时模型逐位不动，
    现有 75.337 kg 结果保持不变。仅在显式传入目标质量时才缩放。
    """
    total = model_total_mass(model)
    if not np.isfinite(total) or total <= 0.0:
        raise ValueError(f"模型原始总质量非法：{total!r} kg，拒绝缩放。")
    if not np.isfinite(target_kg) or target_kg <= 0.0:
        raise ValueError(f"目标质量必须为正有限值，得到 {target_kg!r} kg。")
    f = float(target_kg) / total

    bodies = model.getBodySet()
    for i in range(bodies.getSize()):
        b = bodies.get(i)
        b.setMass(b.getMass() * f)
        inert = b.get_inertia()  # Vec6: (xx, yy, zz, xy, xz, yz)，对质心
        mom = os.Vec3(inert.get(0) * f, inert.get(1) * f, inert.get(2) * f)
        prod = os.Vec3(inert.get(3) * f, inert.get(4) * f, inert.get(5) * f)
        b.setInertia(os.Inertia(mom, prod))
    return f


def _pwl(t: np.ndarray, f: np.ndarray) -> "os.PiecewiseLinearFunction":
    """把 ``(t, f)`` 采样做成 OpenSim 分段线性函数。"""
    p = os.PiecewiseLinearFunction()
    for x, y in zip(np.asarray(t, float), np.asarray(f, float)):
        p.addPoint(float(x), float(y))
    return p


def _time_function(t: np.ndarray, v):
    """把标量或时程做成 OpenSim 时间函数（``None``/标量→Constant，数组→PWL）。"""
    if v is None:
        return os.Constant(0.0)
    arr = np.asarray(v, dtype=float)
    if arr.ndim == 0:
        return os.Constant(float(arr))
    return _pwl(t, arr)


def _add_foot_force(
    model: "os.Model",
    body: str,
    t: np.ndarray,
    f: np.ndarray,
    *,
    f_x=None,
    f_z=None,
    point=None,
) -> None:
    """把一条 GRF 时程施加到 ``body``（力在全局系，向上为正）。

    Parameters
    ----------
    f
        **y 分量**时程（向上为正）—— 保持旧位置参数语义。
    f_x, f_z
        L2 可选水平/前后分量时程。``None``（默认）→ 该轴恒为 0，
        与旧的纯竖直注入**逐位一致**。
    point
        L2 可选作用点 ``(x, y, z)``（**体内坐标系**，默认 calcn 原点 =
        距下中心、对关节无力矩）。每项可为标量或时程；``None``（默认）→
        原点，等价于旧行为。
    """
    pf = os.PrescribedForce()
    pf.setName(f"grf_{body}")
    # 4.x 用 frame socket 解析受力体（setBodyName 会解析失败）
    pf.connectSocket_frame(model.getBodySet().get(body))
    pf.setForceFunctions(_time_function(t, f_x), _pwl(t, f), _time_function(t, f_z))
    if point is None:
        pf.setPointFunctions(os.Constant(0.0), os.Constant(0.0), os.Constant(0.0))
    else:
        px, py, pz = point
        pf.setPointFunctions(
            _time_function(t, px), _time_function(t, py), _time_function(t, pz)
        )
    pf.setForceIsInGlobalFrame(True)
    pf.setPointIsInGlobalFrame(False)
    model.addForce(pf)


def _disable_muscles(model: "os.Model") -> int:
    """dead-drop：关闭激活动力学/腱顺应，激活置 0，返回肌肉数。"""
    nmus = model.getMuscles().getSize()
    for i in range(nmus):
        mus = model.getMuscles().get(i)
        mus.set_ignore_activation_dynamics(True)
        mus.set_ignore_tendon_compliance(True)
    return nmus


def _apply_leg_pose_locked(model: "os.Model", pose: LegPose) -> None:
    """Path (a)：把 ``LEG_FLEX_JOINTS`` 的默认改为 ``pose``，再 ``set_locked(True)``。

    必须先 ``setDefaultValue`` 再 ``set_locked(True)`` —— 锁定的是 *当时*
    的坐标值（实现细节：锁定后再 ``setStateVariableValue`` 在 OpenSim 4.6
    上会被忽略）。

    **锁全部** :data:`LEG_JOINTS`（屈曲关节锁在 ``pose``，其余锁在默认 0）。
    这样 ``posture=LegPose()``（全 0）与改前 ``lock_leg_joints=True`` 路径
    **逐位一致**（`tests/test_nonvertical_s3.py` 的默认回归锁定此契约）。
    """
    cs = model.getCoordinateSet()
    for name, val in pose.as_default_value_dict().items():
        cs.get(name).setDefaultValue(float(val))
        cs.get(name).setDefaultSpeedValue(0.0)
    for name in LEG_JOINTS:
        cs.get(name).set_locked(True)
    # 可选躯干前倾：仅在非 0 时锁 lumbar_extension（0 时保持模型默认"躯干自由"，
    # 从而 posture=LegPose() 与改前逐位一致）。
    if pose.torso_lean_rad != 0.0:
        cs.get("lumbar_extension").setDefaultValue(float(pose.torso_lean_rad))
        cs.get("lumbar_extension").setDefaultSpeedValue(0.0)
        cs.get("lumbar_extension").set_locked(True)


def _apply_passive_stiffness(
    model: "os.Model",
    pose: LegPose,
    passive: PassiveStiffness,
) -> None:
    """Path (b)：**解锁** ``LEG_FLEX_JOINTS`` 并挂上被动刚度 / 阻尼。

    * 先把默认坐标值改为 ``pose``（让 ``SpringGeneralizedForce.rest_length``
      围绕此值抗弯），再 ``set_locked(False)``（解锁）。
    * 然后挂三类力：
      - 线性扭转弹簧（``SpringGeneralizedForce``）：F = −k·(q − q₀) − c·qdot
      - 边壁刚度（``CoordinateLimitForce``）：用 C²-continuous 过渡区
      - 边壁阻尼已包含在 CoordinateLimitForce 的 ``damping`` 项

    **注意** :func:`_disable_muscles` 已关激活动力学；本函数**不重新激活**
    任何肌肉 —— 完全是被动抗力（关节囊 + 韧带 + 肌腱被动弹性）。
    """
    if not passive.is_active():
        raise ValueError(
            "PassiveStiffness 没有任何刚度/阻尼 > 0；要走 unlock+stiffness 路径"
            "至少需要把 *_stiffness 或 *_damping 设非零（例如 PassiveStiffness.nominal_landing()）。"
        )
    cs = model.getCoordinateSet()
    # 先把默认改为姿势（让 rest_length / 边壁下界一致），再解锁
    for name, val in pose.as_default_value_dict().items():
        cs.get(name).setDefaultValue(float(val))
        cs.get(name).setDefaultSpeedValue(0.0)
    # 屈曲关节 = 矢状面自由度 → 解锁（被动刚度作用于此）
    for name in LEG_FLEX_JOINTS:
        cs.get(name).set_locked(False)
    # 其余下肢关节（髋内收 / 旋转）非矢状面 → 仍锁死（避免侧向失控）
    for name in LEG_JOINTS:
        if name not in LEG_FLEX_JOINTS:
            cs.get(name).set_locked(True)

    # 每关节：(扭转弹簧 + 阻尼) + (边壁刚度 + 阻尼)
    spec_table = (
        # (name, stiffness, damping, limit_low, limit_high, limit_k, limit_c)
        ("knee_angle_r", passive.knee_stiffness, passive.knee_damping,
         passive.knee_limit_lower_rad, passive.knee_limit_upper_rad,
         passive.knee_limit_stiffness, passive.knee_limit_damping),
        ("knee_angle_l", passive.knee_stiffness, passive.knee_damping,
         passive.knee_limit_lower_rad, passive.knee_limit_upper_rad,
         passive.knee_limit_stiffness, passive.knee_limit_damping),
        ("hip_flexion_r", passive.hip_stiffness, passive.hip_damping,
         passive.hip_limit_lower_rad, passive.hip_limit_upper_rad,
         passive.hip_limit_stiffness, passive.hip_limit_damping),
        ("hip_flexion_l", passive.hip_stiffness, passive.hip_damping,
         passive.hip_limit_lower_rad, passive.hip_limit_upper_rad,
         passive.hip_limit_stiffness, passive.hip_limit_damping),
        ("ankle_angle_r", passive.ankle_stiffness, passive.ankle_damping,
         passive.ankle_limit_lower_rad, passive.ankle_limit_upper_rad,
         passive.ankle_limit_stiffness, passive.ankle_limit_damping),
        ("ankle_angle_l", passive.ankle_stiffness, passive.ankle_damping,
         passive.ankle_limit_lower_rad, passive.ankle_limit_upper_rad,
         passive.ankle_limit_stiffness, passive.ankle_limit_damping),
    )
    trans = float(passive.transition_rad)
    for name, k_lin, c_lin, lo_lim, hi_lim, k_lim, c_lim in spec_table:
        coord = cs.get(name)
        # ---- 线性扭转弹簧（围绕 rest_length = 当前默认 = pose 值） -------
        if k_lin > 0.0 or c_lin > 0.0:
            sp = os.SpringGeneralizedForce()
            sp.setName(f"passive_{name}")
            sp.set_coordinate(name)
            sp.set_rest_length(float(coord.getDefaultValue()))
            sp.set_stiffness(float(k_lin))
            sp.set_viscosity(float(c_lin))
            model.addForce(sp)
        # ---- 边壁刚度（C²-continuous，限范围） --------------------------
        if k_lim > 0.0 or c_lim > 0.0:
            lo = (float(lo_lim) if lo_lim is not None
                  else float(coord.getRangeMin()))
            hi = (float(hi_lim) if hi_lim is not None
                  else float(coord.getRangeMax()))
            # 防退化：lower ≥ upper ⇒ 退化为单边壁
            clf = os.CoordinateLimitForce()
            clf.setName(f"limit_{name}")
            clf.set_coordinate(name)
            clf.set_lower_limit(float(lo))
            clf.set_lower_stiffness(float(k_lim))
            clf.set_upper_limit(float(hi))
            clf.set_upper_stiffness(float(k_lim))
            clf.set_damping(float(c_lim))
            clf.set_transition(trans)
            model.addForce(clf)


def run_dead_drop(
    grf: GroundReaction,
    *,
    model_path: str | Path = DEFAULT_MODEL,
    final_time_s: float | None = None,
    accuracy: float = 1e-4,
    max_step_s: float = 2e-5,
    lock_leg_joints: bool = True,
    scale_mass_kg: float | None = None,
    vx0_ms: float = 0.0,
    vz0_ms: float = 0.0,
    posture: LegPose | None = None,
    passive: PassiveStiffness | None = None,
) -> FallResult:
    """把 ``grf`` 注入模型并做 Forward Dynamics 积分。

    Parameters
    ----------
    grf
        :func:`ground_reaction` 的产物（每只脚一条时程）。
    final_time_s
        积分终止时间；默认 = GRF 结束 + 3 ms 余量。
    accuracy, max_step_s
        积分器精度 / 最大步长。``max_step_s`` 默认 20 µs，
        匹配冲击载荷的最快特征。
    lock_leg_joints
        锁死下肢关节（见 :data:`LEG_JOINTS` 的说明）。S1 刚性腿近似，默认 True。
    scale_mass_kg
        **可选** body 质量缩放目标 (kg)。``None``（默认）**不缩放**，模型保持
        Rajagopal2015 原始 75.337 kg，现有结果逐位不变。传入数值（如 77.0）
        时先 :func:`scale_model_mass` 把各 body 质量/惯量等比缩放到该总质量，
        用于 [Y25] 的 AM50 (77 kg) 质量标定实验。
    vx0_ms, vz0_ms
        **可选** 骨盆水平初速度分量 (m/s)，全局系 x / z。默认均为 ``0`` =
        纯竖直（现有行为逐位不变）。非零时给出撞击速度的水平分量（L2，方案
        §2 自由度①）。**GRF 方向**由 ``grf.ground_normal`` 决定（非竖直时
        自动按三分量注入；见 :func:`_add_foot_force`）。
    posture
        **OPT-IN** S3 目标下肢姿势（见 :class:`LegPose`）。Path (a)：
        把 ``LEG_FLEX_JOINTS`` 的默认改为 ``posture`` 后 ``set_locked(True)``
        —— 力学上等价于"刚性腿 + 弯曲几何"。``None``（默认）= 直立（与旧
        行为**逐位一致**）。**互斥于** ``passive``。
    passive
        **OPT-IN** S3 被动刚度（见 :class:`PassiveStiffness`）。Path (b)：
        解锁 ``LEG_FLEX_JOINTS`` 并挂上扭转弹簧 + 边壁刚度，让关节在冲击下
        有限地弯。``None``（默认）= 不加刚度（与旧行为**逐位一致**）。
        **互斥于** ``posture``（同时给 → ``ValueError``）；起始姿势固定为
        中立位（``LegPose()``），只显示被动刚度的纯效应。

    Returns
    -------
    FallResult
        ``model_mass_kg`` / ``mass_scale`` 记录积分时模型的实际总质量与缩放因子。

    Notes
    -----
    * **默认行为契约** (``posture=None`` 且 ``passive=None``)：
      ``posture = LegPose()`` ≡ 全 0，``_apply_leg_pose_locked`` 把默认改为 0
      后再锁 → 与改前 ``for c in LEG_JOINTS: set_locked(True)`` **逐位一致**
      （坐标默认就是 0）。``tests/test_nonvertical_s3.py`` 的
      ``test_default_reproduces_cached_axial`` 会硬门槛验证。
    * ``passive`` 设了但 ``is_active() is False`` → 抛 ``ValueError``（防止
      "什么都没做还跑了 30 秒"）。
    * ``posture`` 与 ``passive`` 同时给 → 抛 ``ValueError``（语义冲突）。
    """
    if posture is not None and passive is not None:
        raise ValueError(
            "posture (Path a: lock-at-pose) 与 passive (Path b: unlock+stiffness)"
            "语义冲突；同时给时无法决定走哪条路径。请只传一个。"
        )
    if passive is not None and not passive.is_active():
        raise ValueError(
            "passive 给定了但没有任何刚度/阻尼 > 0；用 PassiveStiffness.nominal_landing()"
            "或显式设 *_stiffness / *_damping > 0。"
        )

    model = os.Model(str(model_path))
    applied_scale = 1.0
    if scale_mass_kg is not None:
        applied_scale = scale_model_mass(model, scale_mass_kg)
    _disable_muscles(model)

    if grf.is_vertical:
        # 默认路径：与原实现逐位一致（x/z 恒为 0、作用于原点）
        _add_foot_force(model, FOOT_BODIES[0], grf.t_s, grf.f_right_n)
        _add_foot_force(model, FOOT_BODIES[1], grf.t_s, grf.f_left_n)
    else:
        # L1/L2 矢量路径：三分量注入（方向 = grf.ground_normal）
        for body, side in zip(FOOT_BODIES, ("r", "l")):
            vec = grf.force_vec_n(side)
            _add_foot_force(
                model, body, grf.t_s, vec[:, 1], f_x=vec[:, 0], f_z=vec[:, 2]
            )

    # ---- S3 姿势 / 被动刚度（OPT-IN）---------------------------------------
    if posture is not None:
        # Path (a)：把默认改为目标姿势，再锁
        _apply_leg_pose_locked(model, posture)
    elif passive is not None:
        # Path (b)：解锁 + 挂被动刚度。起始姿势固定为 LegPose()（= 全 0）——
        # 用户若想看"屈膝 30° 起步的有限弯"，应手动跑两个脚本（一个 LegPose 锁定
        # baseline，一个 passive 路径对比）。这里保持中立位以显示被动刚度的纯效应。
        _apply_passive_stiffness(model, LegPose(), passive)
    elif lock_leg_joints:
        # 旧路径：LEG_JOINTS 全部锁死（默认 = 直立，逐位复现改前）
        for c in LEG_JOINTS:
            model.getCoordinateSet().get(c).set_locked(True)

    state = model.initSystem()

    # 初速度：整体以 v0 向下（y 向上为正，故取负）；水平分量可选（默认 0）
    model.setStateVariableValue(state, PELVIS_TY_SPEED_STATE, -grf.impact_speed_ms)
    if vx0_ms != 0.0:
        model.setStateVariableValue(state, PELVIS_TX_SPEED_STATE, float(vx0_ms))
    if vz0_ms != 0.0:
        model.setStateVariableValue(state, PELVIS_TZ_SPEED_STATE, float(vz0_ms))

    # 激活清零（ignore_activation_dynamics 下 excitation 直接当 activation）
    sv = model.getStateVariableNames()
    for i in range(sv.size()):
        nm = sv.get(i)
        if nm.endswith("/activation"):
            model.setStateVariableValue(state, nm, 0.0)

    tf = float(final_time_s) if final_time_s is not None else float(grf.t_s[-1]) + 3e-3

    man = os.Manager(model)
    man.setIntegratorMethod(os.Manager.IntegratorMethod_RungeKuttaMerson)
    man.setIntegratorAccuracy(accuracy)
    man.setIntegratorMaximumStepSize(max_step_s)
    man.initialize(state)
    man.setIntegratorFinalTime(tf)
    man.integrate(tf)

    table = man.getStatesTable()
    t = np.asarray(table.getIndependentColumn(), dtype=float)
    total = grf.f_left_n + grf.f_right_n

    return FallResult(
        model=model,
        states=table,
        t_s=t,
        peak_applied_n=float(np.max(total)),
        applied_total_n=total,
        v0_ms=grf.impact_speed_ms,
        final_time_s=tf,
        model_path=Path(model_path),
        model_mass_kg=model_total_mass(model),
        mass_scale=applied_scale,
    )
