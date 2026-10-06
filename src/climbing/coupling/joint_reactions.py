"""S1.4 —— 把"距下关节反力"泛化到**全身各关节反力**（载荷端）。

背景
----
``joint_loads.subtalar_reaction`` 用**牛顿自由体**算距下关节 wrench（不用
``os.JointReaction``，它在 4.6 里经 Manager 积分后不回填 storage，跑不通，见
``joint_loads`` docstring）。本模块把同一套自由体数学推广到**任意关节**：
把隔离体从"距下关节的远端"换成"该关节 child body 的**远端刚体子树**"。

自由体方程（与 ``joint_loads`` 逐字一致）
----------------------------------------
对关节 ``J`` 的远端子系统 ``D``（child body + 其全部子代刚体）：

    R + Σ F_ext + Σᵢ mᵢ·g = Σᵢ mᵢ·aᵢ
    ⇒ R = Σᵢ mᵢ·(aᵢ − g) − Σ F_ext                 (力，全局系)

    M = Σᵢ (pᵢ − p_joint) × (mᵢ·aᵢ)                (动量矩率)
        − Σᵢ (pᵢ − p_joint) × (mᵢ·g)                (重力矩)
        − Σ (p_apply − p_joint) × F_ext             (外力矩)

关节中心 ``p_joint`` = **child body 的体坐标系原点**（与距下关节取 calcn
原点完全一致）。所有量在**全局系**下、对 ``p_joint`` 取矩。

关节 → 远端子树（Rajagopal2015 拓扑，实测）
-------------------------------------------
模型 joint 树（``*_offset`` 帧去后缀即 body）：

    ground_pelvis : ground → pelvis
    hip_{s}       : pelvis → femur_{s}
    walker_knee_{s}: femur_{s} → tibia_{s}
    patellofemoral_{s}: femur_{s} → patella_{s}
    ankle_{s}     : tibia_{s} → talus_{s}
    subtalar_{s}  : talus_{s} → calcn_{s}
    mtp_{s}       : calcn_{s} → toes_{s}
    back          : pelvis → torso
    acromial_{s}  : torso → humerus_{s}
    elbow_{s}     : humerus_{s} → ulna_{s}
    radioulnar_{s}: ulna_{s} → radius_{s}
    radius_hand_{s}: radius_{s} → hand_{s}

由此得到的子树（``s ∈ {r,l}``）：

    subtalar_{s} : calcn, toes
    ankle_{s}    : talus, calcn, toes
    knee_{s}     : tibia, talus, calcn, toes
    hip_{s}      : femur, tibia, patella, talus, calcn, toes
    lumbar       : torso, humerus, ulna, radius, hand（左右）

⚠️ 两点与直觉不同，均为**模型 joint 树的必然结果**，T2 会逐项自检：
* ``hip`` 子树**包含 patella**（patellofemoral 关节的 parent 是 femur，
  patella 在髋切面远端）—— 任务列出的 "femur, tibia, talus, calcn, toes"
  是其子集，patella 必须计入才是闭合子树。
* ``knee`` 子树**不含 patella**（patella 挂在 femur 侧，是膝切面近端）。
* 本模型**没有 head/neck 刚体**：``torso`` 是躯干链顶，故 ``lumbar`` 的
  "躯干/上肢/头" 退化为 torso + 双上肢。

不变量 T2（逐关节）
-------------------
静止、无外力时，每个关节的纵向支撑力应 = **该子树远端自重**：
``R_y = m_d · g``（``gravity_baseline_joint``）。残差应 ≪ 冲击峰值。

⚠️ 状态时间修正（本模块相对既有 ``joint_loads._realized_states`` 的唯一差异）
----------------------------------------------------------------------------
``joint_loads._realized_states`` 逐帧写入 q/u 并 ``realizeAcceleration``，但
**从不写 state 的时间**（实测 ``ws.getTime()==0``）。OpenSim 的外力
（``PrescribedForce``）按 ``state.getTime()`` 求值，于是所有帧都用
``GRF(t=0)`` —— 得到与真实轨迹**不一致**的加速度。

实测（``temp/opensim_fe/joint_reactions/_probe_fixtime.py``，@5 m dead-drop）：
把 ``findStationAccelerationInGround``（骨盆 COM 的 a_y）与由
``pelvis_ty/speed`` 列做数值微分得到的 ``d(speed)/dt`` 对比：

    | 口径                    | max|a_y| | corr  |
    |-------------------------|----------|-------|
    | 旧（不设时间）          |    7.7   | 0.13  |
    | 新（每帧 setTime(row_t)）| 2382.4   | 0.9997 |

旧口径的骨盆加速度几乎是"自由落体"（≈ −g），完全没反映 ~2000 m/s² 的
减速；新口径与状态轨迹一致。对**距下**等含足底外力的关节，旧口径因
``R = m_d(a−g) − F_ext`` 中 ``−F_ext`` 项（用插值 GRF，正确）占主导而
"看起来对"；但对**腰椎**这类子树内**无外力**的关节，``R`` 完全由被污染的
加速度决定 → 旧口径**不可信**。

因此本模块新增 ``set_state_time=True``（默认）的**修正口径**：每帧
``ws.setTime(row_t)`` 后再 realize。这**不改动** ``joint_loads.py``；旧口径仍
可通过 ``set_state_time=False`` 复现（``legacy`` 值），用于回归对照。
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import opensim as os

from .opensim_fall import FallResult
from .joint_loads import _G, _as_vec3, _mass_center, _realized_states

__all__ = [
    "JointReactionResult",
    "JointSpec",
    "JOINT_SPECS",
    "DEFAULT_JOINT_ORDER",
    "joint_subtree_bodies",
    "joint_reaction",
    "gravity_baseline_joint",
    "locked_static_model",
    "subtree_mass_kg",
    "subtree_weight_n",
    "t2_residual",
    "verify_subtree_topology",
]

#: 承受 GRF 的足刚体 → ``GroundReaction`` 上的力时程字段
_FOOT_SERIES = {"calcn_r": "f_right_n", "calcn_l": "f_left_n"}


def _realized_states_timed(fall: FallResult):
    """同 ``joint_loads._realized_states``，但**逐帧显式写入 state 时间**。

    为什么必须：``prescribed force``（足底 GRF）按 ``state.getTime()`` 求值；
    不写时间 ⇒ 所有帧都用 ``GRF(0)`` ⇒ 加速度与真实轨迹不符（见模块 docstring
    的实测表）。这里在 ``realizeAcceleration`` 前 ``ws.setTime(row_t)``。
    """
    model = fall.model
    tab = fall.states
    labels = [tab.getColumnLabel(j) for j in range(tab.getNumColumns())]
    cols = [tab.getDependentColumnAtIndex(j).to_numpy() for j in range(tab.getNumColumns())]
    times = np.asarray(tab.getIndependentColumn(), dtype=float)
    ws = model.getWorkingState()
    for i in range(tab.getNumRows()):
        for j, l in enumerate(labels):
            try:
                model.setStateVariableValue(ws, l, float(cols[j][i]))
            except Exception:
                pass  # locked 坐标不可写，保持初值即可
        ws.setTime(float(times[i]))
        model.realizeAcceleration(ws)
        yield ws


@dataclass(frozen=True)
class JointSpec:
    """一个关节的几何/拓扑定义。"""

    key: str
    model_joint: str
    child_body: str
    bodies: tuple[str, ...]
    side: str | None = None

    @property
    def distal_body(self) -> str:
        """关节中心所在的 child body（取体坐标原点）。"""
        return self.child_body


def _leg_specs(side: str) -> dict[str, JointSpec]:
    s = side
    return {
        f"subtalar_{s}": JointSpec(
            key=f"subtalar_{s}", model_joint=f"subtalar_{s}",
            child_body=f"calcn_{s}", bodies=(f"calcn_{s}", f"toes_{s}"), side=s,
        ),
        f"ankle_{s}": JointSpec(
            key=f"ankle_{s}", model_joint=f"ankle_{s}",
            child_body=f"talus_{s}",
            bodies=(f"talus_{s}", f"calcn_{s}", f"toes_{s}"), side=s,
        ),
        f"knee_{s}": JointSpec(
            key=f"knee_{s}", model_joint=f"walker_knee_{s}",
            child_body=f"tibia_{s}",
            bodies=(f"tibia_{s}", f"talus_{s}", f"calcn_{s}", f"toes_{s}"), side=s,
        ),
        f"hip_{s}": JointSpec(
            key=f"hip_{s}", model_joint=f"hip_{s}",
            child_body=f"femur_{s}",
            bodies=(f"femur_{s}", f"tibia_{s}", f"patella_{s}",
                    f"talus_{s}", f"calcn_{s}", f"toes_{s}"), side=s,
        ),
    }


_JL = _leg_specs("l")
_JR = _leg_specs("r")
_LUMBAR = JointSpec(
    key="lumbar", model_joint="back", child_body="torso",
    bodies=("torso",
            "humerus_r", "ulna_r", "radius_r", "hand_r",
            "humerus_l", "ulna_l", "radius_l", "hand_l"),
    side=None,
)

#: 关节键 → 定义（含右/左腿链 + 腰椎链）
JOINT_SPECS: dict[str, JointSpec] = {
    **_JR, **_JL, "lumbar": _LUMBAR,
}

#: 报告的默认输出顺序（近端→远端单调）
DEFAULT_JOINT_ORDER: tuple[str, ...] = (
    "subtalar_r", "ankle_r", "knee_r", "hip_r",
    "subtalar_l", "ankle_l", "knee_l", "hip_l",
    "lumbar",
)


@dataclass
class JointReactionResult:
    """一个关节的反力 wrench 时程（全局系，对关节中心取矩）。"""

    t_s: np.ndarray
    force_n: np.ndarray        # (N,3)
    moment_nm: np.ndarray      # (N,3)
    joint: str
    side: str | None
    bodies: tuple[str, ...] = ()
    subtree_mass_kg: float = 0.0

    @property
    def peak_force_n(self) -> float:
        return float(np.max(np.linalg.norm(self.force_n, axis=1)))

    @property
    def peak_force_time_s(self) -> float:
        return float(self.t_s[int(np.argmax(np.linalg.norm(self.force_n, axis=1)))])

    @property
    def peak_vertical_n(self) -> float:
        """峰值纵向（y）分量的绝对值 —— 传给单骨 FE 的主载荷。"""
        return float(np.max(np.abs(self.force_n[:, 1])))

    @property
    def peak_vertical_time_s(self) -> float:
        return float(self.t_s[int(np.argmax(np.abs(self.force_n[:, 1])))])

    @property
    def peak_moment_nm(self) -> float:
        return float(np.max(np.linalg.norm(self.moment_nm, axis=1)))

    @property
    def peak_moment_time_s(self) -> float:
        return float(self.t_s[int(np.argmax(np.linalg.norm(self.moment_nm, axis=1)))])


# --------------------------------------------------------------------------
# 拓扑解析 / 校验
# --------------------------------------------------------------------------
def _body_of_frame(frame_name: str) -> str:
    return frame_name[:-len("_offset")] if frame_name.endswith("_offset") else frame_name


def _derived_subtree(model: "os.Model", root_body: str) -> list[str]:
    """按模型 joint 树求 ``root_body`` 的全部子代刚体（不含 root）。"""
    js = model.getJointSet()
    parents: dict[str, list[str]] = {}
    for i in range(js.getSize()):
        j = js.get(i)
        p = _body_of_frame(j.getParentFrame().getName())
        c = _body_of_frame(j.getChildFrame().getName())
        parents.setdefault(p, []).append(c)
    seen: dict[str, bool] = {}
    stack = [root_body]
    while stack:
        b = stack.pop()
        for c in parents.get(b, []):
            if c not in seen:
                seen[c] = True
                stack.append(c)
    return [root_body] + list(seen.keys())


def _resolve_spec(joint: str, side: str | None) -> JointSpec:
    """把 ``(joint, side)`` 解析成 :class:`JointSpec`（支持全名或基名+side）。"""
    key = joint
    if key not in JOINT_SPECS:
        if side is not None and f"{joint}_{side}" in JOINT_SPECS:
            key = f"{joint}_{side}"
        else:
            raise ValueError(
                f"未知关节 {joint!r}（side={side!r}）；可选：{sorted(JOINT_SPECS)}"
            )
    spec = JOINT_SPECS[key]
    if spec.side is None:
        if side is not None:
            raise ValueError(f"关节 {joint!r} 无侧别（side 必须为 None），得到 {side!r}")
    else:
        if side is not None and side != spec.side:
            raise ValueError(f"关节 {joint!r} 是 {spec.side} 侧，与 side={side!r} 冲突")
    return spec


def joint_subtree_bodies(joint: str, side: str | None = None) -> tuple[str, ...]:
    """返回该关节远端子树的刚体名（含 child body）。"""
    return _resolve_spec(joint, side).bodies


def verify_subtree_topology(model: "os.Model") -> dict[str, list[str]]:
    """用模型 joint 树重算各关节子树，返回 ``key -> [差异描述]``（空 = 一致）。"""
    diffs: dict[str, list[str]] = {}
    for key, spec in JOINT_SPECS.items():
        derived = set(_derived_subtree(model, spec.distal_body))  # 含 root
        declared = set(spec.bodies)
        missing = sorted(derived - declared)
        extra = sorted(declared - derived)
        msgs: list[str] = []
        if missing:
            msgs.append(f"声明缺 {missing}")
        if extra:
            msgs.append(f"声明多 {extra}")
        if msgs:
            diffs[key] = msgs
    return diffs


# --------------------------------------------------------------------------
# 质量 / T2
# --------------------------------------------------------------------------
def subtree_mass_kg(model: "os.Model", *, joint: str, side: str | None = None) -> float:
    spec = _resolve_spec(joint, side)
    bs = model.getBodySet()
    return float(sum(bs.get(n).getMass() for n in spec.bodies))


def subtree_weight_n(model: "os.Model", *, joint: str, side: str | None = None) -> float:
    return subtree_mass_kg(model, joint=joint, side=side) * _G


def gravity_baseline_joint(
    model: "os.Model", *, joint: str, side: str | None = None
) -> float:
    """T2：静止、无外力时，该关节纵向支撑力应 = 子树远端自重 (N)。"""
    spec = _resolve_spec(joint, side)
    bs = model.getBodySet()
    bodies = [bs.get(n) for n in spec.bodies]
    masses = np.array([b.getMass() for b in bodies], dtype=float)
    m_d = float(masses.sum())
    model.initSystem()
    st = model.initializeState()
    model.realizeAcceleration(st)
    acc = np.zeros(3)
    for m_i, b in zip(masses, bodies):
        acc += m_i * _as_vec3(b.findStationAccelerationInGround(st, _mass_center(b)))
    a_d = acc / m_d
    r = m_d * (a_d - np.array([0.0, -_G, 0.0]))
    return float(r[1])


# --------------------------------------------------------------------------
# 主计算
# --------------------------------------------------------------------------
def joint_reaction(
    fall: FallResult,
    grf,
    *,
    joint: str,
    side: str | None = None,
    applied_series: np.ndarray | None = None,
    set_state_time: bool = True,
) -> JointReactionResult:
    """计算关节 ``joint``（可选 ``side``）的反力 wrench。

    与 :func:`joint_loads.subtalar_reaction` **同一自由体数学**，仅把隔离体
    换成该关节的远端刚体子树，并把子树内所有足底外力（GRF）计入。

    Parameters
    ----------
    fall
        :func:`opensim_fall.run_dead_drop` 结果（含 model + states）。
    grf
        :func:`opensim_grf.ground_reaction` 结果。
    joint
        关节键，如 ``"ankle_r"`` / ``"knee_l"`` / ``"hip_r"`` / ``"lumbar"``
        （也接受基名 ``"ankle"`` + ``side="r"``）。
    side
        侧别 ``"r"`` / ``"l"``；``lumbar`` 必须为 ``None``。默认从 ``joint`` 后缀解析。
    applied_series
        可选：覆盖**子树内唯一足**的竖直力时程（N，向上为正）。若子树含 0 或
        >1 只足则报错（避免静默忽略）。
    set_state_time
        ``True``（默认）= 修正口径：逐帧 ``ws.setTime(row_t)`` 后再 realize，
        使加速度与真实轨迹一致（见模块 docstring）。``False`` = 复现既有
        ``joint_loads._realized_states`` 的旧口径（不写时间），仅用于回归对照。
    """
    spec = _resolve_spec(joint, side)
    resolved_side = spec.side

    model = fall.model
    bs = model.getBodySet()
    bodies = [bs.get(n) for n in spec.bodies]
    masses = np.array([b.getMass() for b in bodies], dtype=float)
    m_d = float(masses.sum())
    mcs = [_mass_center(b) for b in bodies]
    child_body = bs.get(spec.distal_body)

    # 子树内的足底外力（默认取 grf 的对应侧）
    feet = [n for n in spec.bodies if n in _FOOT_SERIES]
    if applied_series is not None:
        if len(feet) != 1:
            raise ValueError(
                f"关节 {spec.key!r} 子树含 {len(feet)} 只足，applied_series 只能覆盖单足"
            )
    ext_bodies: list[tuple[str, str]] = []   # (body_name, series_key)
    for n in feet:
        ext_bodies.append((n, _FOOT_SERIES[n]))

    # 只在**外力窗口内**分析：PrescribedForce 的 PWL 在窗口外会线性外推成
    # 巨大负值（实测 t=24 ms 时约 −88 kN），污染 run_dead_drop 的 3 ms 余量段。
    t_all = np.asarray(fall.t_s, dtype=float)
    n_valid = int(np.searchsorted(t_all, float(np.max(grf.t_s)) + 1e-12))
    t = t_all[:n_valid]
    g_vec = np.array([0.0, -_G, 0.0])

    # 预先插值外力时程到（窗口内的）fall.t_s。
    # 默认 grf 为竖直法向 → force_vec_n 给出 (0, f(t), 0)，与旧行为逐位一致；
    # 非竖直法向时按完整 3D GRF 计入自由体（L1/L2 矢量通路）。
    ext_series: list[tuple["os.Body", np.ndarray]] = []
    for body_name, series_attr in ext_bodies:
        if applied_series is not None:
            y = np.interp(t, grf.t_s, applied_series)
            vec = np.column_stack([np.zeros_like(y), y, np.zeros_like(y)])
        else:
            fv = grf.force_vec_n(body_name[-1])
            vec = np.column_stack(
                [np.interp(t, grf.t_s, fv[:, k]) for k in range(3)]
            )
        ext_series.append((bs.get(body_name), vec))

    forces = np.full((t.size, 3), np.nan)
    moments = np.full((t.size, 3), np.nan)
    realize = _realized_states_timed(fall) if set_state_time else _realized_states(fall)

    for i, st in enumerate(realize):
        if i >= n_valid:
            break
        p_joint = _as_vec3(child_body.findStationLocationInGround(st, os.Vec3(0, 0, 0)))

        sum_ma = np.zeros(3)
        rate = np.zeros(3)        # Σ (p_i−p_j)×m_i a_i
        grav = np.zeros(3)        # Σ (p_i−p_j)×m_i g
        for m_i, b, c in zip(masses, bodies, mcs):
            a_i = _as_vec3(b.findStationAccelerationInGround(st, c))
            p_i = _as_vec3(b.findStationLocationInGround(st, c))
            r_i = p_i - p_joint
            sum_ma += m_i * a_i
            rate += np.cross(r_i, m_i * a_i)
            grav += np.cross(r_i, m_i * g_vec)
        a_d = sum_ma / m_d

        f_ext_total = np.zeros(3)
        ext_moment = np.zeros(3)
        for body, vec in ext_series:
            f = vec[i]
            p_apply = _as_vec3(body.findStationLocationInGround(st, os.Vec3(0, 0, 0)))
            f_ext_total += f
            ext_moment += np.cross(p_apply - p_joint, f)

        forces[i] = m_d * (a_d - g_vec) - f_ext_total
        moments[i] = rate - grav - ext_moment

    return JointReactionResult(
        t_s=t, force_n=forces, moment_nm=moments,
        joint=spec.key, side=resolved_side,
        bodies=spec.bodies, subtree_mass_kg=m_d,
    )


def locked_static_model(model_path: str | None = None) -> "os.Model":
    """构建一个**全部坐标锁死**的静止模型，供 T2 支撑力自检使用。

    为什么全锁：T2 的物理含义是"身体被支撑住、静止"——此时每个关节的纵向
    反力 = 该子树的远端自重。若不锁根（ground_pelvis 的 6 个坐标）与上体
    坐标，模型会**自由落体**，``Σm·(a−g)=0``，得到 R_y≈0（这是 ``joint_loads
    .gravity_baseline`` 的隐式口径）。要复现"支撑=自重"，必须把模型完全约束。
    """
    from .opensim_fall import DEFAULT_MODEL as _DM

    m = os.Model(str(model_path if model_path is not None else _DM))
    cs = m.getCoordinateSet()
    for i in range(cs.getSize()):
        cs.get(i).set_locked(True)
    return m


def t2_residual(model_locked: "os.Model", key: str) -> tuple[float, float, float]:
    """返回 ``(subtree_mass_kg, weight_n, residual_n)``；``residual = R_y − weight``。"""
    w = subtree_weight_n(model_locked, joint=key)
    r = gravity_baseline_joint(model_locked, joint=key)
    return subtree_mass_kg(model_locked, joint=key), w, r - w


if __name__ == "__main__":  # pragma: no cover - 手跑自检
    base_model = locked_static_model()

    diffs = verify_subtree_topology(base_model)
    print("子树拓扑校验：", "✅ 全部一致" if not diffs else f"❌ {diffs}")

    print("\n关节          子树质量(kg)   自重(N)    T2残差(N)")
    for key in DEFAULT_JOINT_ORDER:
        msub, w, resid = t2_residual(base_model, key)
        print(f"  {key:12s} {msub:9.4f} {w:9.2f}  {resid:+.3e}")
