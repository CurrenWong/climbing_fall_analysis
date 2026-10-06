"""S1.3 —— 距下关节反力（wrench）的透明计算。

为什么不用 ``os.JointReaction``
------------------------------
OpenSim 4.6 Python 绑定里 ``JointReaction`` **经 Manager 积分后不回填
``getStorageList()``**（实测 size=0，``printResults`` 也写不出文件），
即分析根本没有被 Manager 执行。与其依赖一个跑不通的黑箱，这里用
**牛顿自由体**直接算 —— 更透明、可测，也天然挂得上不变量 T2。

受力隔离体
----------
距下关节 ``subtalar_r`` 的**远端**是 ``calcn_r``（+ 其远端的 ``toes_r``，
mtp 锁死 → 与 calcn 刚性相连）。对该远端子系统列牛顿第二定律：

    R + F_ext + m_d·g = m_d·a_d
    ⇒ R = m_d·(a_d − g) − F_ext            (力, 全局系)

    M = Σᵢ (pᵢ − p_joint) × (mᵢ·aᵢ)          (动量矩率)
        − Σᵢ (pᵢ − p_joint) × (mᵢ·g)         (重力矩)
        − (p_apply − p_joint) × F_ext         (外力矩)

其中 ``a_d = Σᵢ mᵢ·aᵢ / m_d`` 是远端子系统质心加速度，由每个刚体
``mass_center`` 站的 ``findStationAccelerationInGround`` 得到
（``p_apply = p_joint``，因为 GRF 施加在 calcn 原点 = 距下关节中心）。

不变量 T2：无外力、静止时 ``R`` 应恰为远端自重的**向上支撑力**
（= ``m_d·g``）。``gravity_baseline()`` 就是这个自检。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import opensim as os

from .opensim_fall import FallResult

__all__ = ["JointReactionResult", "subtalar_reaction", "gravity_baseline"]

_G = 9.80665

#: 距下关节远端的刚体（mtp 锁死 → calcn 与 toes 刚性相连）
_DISTAL = {"r": ("calcn_r", "toes_r"), "l": ("calcn_l", "toes_l")}


@dataclass
class JointReactionResult:
    t_s: np.ndarray
    force_n: np.ndarray       # (N,3) 全局系，talus 对 calcn 的作用力
    moment_nm: np.ndarray     # (N,3) 对距下关节中心
    side: str

    @property
    def peak_force_n(self) -> float:
        return float(np.max(np.linalg.norm(self.force_n, axis=1)))

    @property
    def peak_vertical_n(self) -> float:
        """峰值纵向（y）分量的绝对值 —— 传给单骨 FE 的主载荷。"""
        return float(np.max(np.abs(self.force_n[:, 1])))

    @property
    def peak_moment_nm(self) -> float:
        return float(np.max(np.linalg.norm(self.moment_nm, axis=1)))

    @property
    def peak_vertical_time_s(self) -> float:
        return float(self.t_s[int(np.argmax(np.abs(self.force_n[:, 1])))])


def _mass_center(body) -> "os.Vec3":
    for name in ("get_mass_center", "getMassCenter"):
        if hasattr(body, name):
            return getattr(body, name)()
    raise AttributeError(f"{body.getName()}: 找不到 mass_center 访问器")


def _as_vec3(result) -> np.ndarray:
    return np.array([result.get(0), result.get(1), result.get(2)], dtype=float)


def _realized_states(fall: FallResult):
    """把 states 表逐帧写进模型的 working state，并 realize 到 Acceleration。

    不用 ``StatesTrajectory``：它遇到 locked 坐标会**静默产出 NaN 状态**
    （实测 ``findStationLocationInGround`` 全 NaN，而表里的 q 都是有限值）。
    """
    model = fall.model
    tab = fall.states
    labels = [tab.getColumnLabel(j) for j in range(tab.getNumColumns())]
    cols = [tab.getDependentColumnAtIndex(j).to_numpy() for j in range(tab.getNumColumns())]
    ws = model.getWorkingState()
    for i in range(tab.getNumRows()):
        for j, l in enumerate(labels):
            try:
                model.setStateVariableValue(ws, l, float(cols[j][i]))
            except Exception:
                pass  # locked 坐标不可写，保持初值即可
        model.realizeAcceleration(ws)
        yield ws


def subtalar_reaction(
    fall: FallResult,
    grf,
    *,
    side: str = "r",
    applied_series: np.ndarray | None = None,
) -> JointReactionResult:
    """计算距下关节 ``subtalar_{side}`` 的反力 wrench。

    Parameters
    ----------
    fall
        :func:`opensim_fall.run_dead_drop` 的结果（含 model + states）。
    grf
        :func:`ground_reaction` 结果，用于取外力时程。
    side
        ``"r"`` / ``"l"``。
    applied_series
        施加在该足的竖直力时程（N，向上为正）；默认按 ``grf`` 的左右分配。
    """
    if side not in _DISTAL:
        raise ValueError(f"side 必须是 'r'/'l'，得到 {side!r}")
    if applied_series is None:
        applied_series = grf.f_right_n if side == "r" else grf.f_left_n

    model = fall.model
    bodies = [model.getBodySet().get(n) for n in _DISTAL[side]]
    masses = np.array([b.getMass() for b in bodies], dtype=float)
    m_d = float(masses.sum())
    mc = [_mass_center(b) for b in bodies]

    t = np.asarray(fall.t_s, dtype=float)
    ext_y = np.interp(t, grf.t_s, applied_series)

    forces = np.full((t.size, 3), np.nan)
    moments = np.full((t.size, 3), np.nan)
    g_vec = np.array([0.0, -_G, 0.0])

    for i, st in enumerate(_realized_states(fall)):
        p_joint = _as_vec3(bodies[0].findStationLocationInGround(st, os.Vec3(0, 0, 0)))
        acc_terms = np.zeros(3)
        pos_terms = np.zeros(3)      # Σ m_i a_i
        for m_i, b, c in zip(masses, bodies, mc):
            a_i = _as_vec3(b.findStationAccelerationInGround(st, c))
            p_i = _as_vec3(b.findStationLocationInGround(st, c))
            acc_terms += m_i * a_i
            pos_terms += m_i * (p_i - p_joint)
        a_d = acc_terms / m_d

        f_ext = np.array([0.0, ext_y[i], 0.0])
        forces[i] = m_d * (a_d - g_vec) - f_ext

        # 动量矩率 − 重力矩（GRF 作用在关节中心，力偶为 0）
        # Σ (p_i−p_j)×m_i a_i
        rate = np.zeros(3)
        grav = np.zeros(3)
        for m_i, b, c in zip(masses, bodies, mc):
            a_i = _as_vec3(b.findStationAccelerationInGround(st, c))
            p_i = _as_vec3(b.findStationLocationInGround(st, c))
            r = p_i - p_joint
            rate += np.cross(r, m_i * a_i)
            grav += np.cross(r, m_i * g_vec)
        moments[i] = rate - grav

    return JointReactionResult(t_s=t, force_n=forces, moment_nm=moments, side=side)


def gravity_baseline(model: "os.Model", *, side: str = "r") -> float:
    """不变量 T2：静止、无外力时，距下关节纵向支撑力应 = 远端自重 (N)。"""
    bodies = [model.getBodySet().get(n) for n in _DISTAL[side]]
    m_d = float(sum(b.getMass() for b in bodies))
    model.initSystem()
    st = model.initializeState()
    model.realizeAcceleration(st)
    p_joint = bodies[0].findStationLocationInGround(st, os.Vec3(0, 0, 0))
    acc = np.zeros(3)
    for b in bodies:
        a_i = _as_vec3(b.findStationAccelerationInGround(st, _mass_center(b)))
        acc += b.getMass() * a_i
    a_d = acc / m_d
    r = m_d * (a_d - np.array([0.0, -_G, 0.0]))
    return float(r[1])


if __name__ == "__main__":
    from .opensim_grf import ground_reaction
    from .opensim_fall import run_dead_drop, LEG_JOINTS, DEFAULT_MODEL

    base_model = os.Model(str(DEFAULT_MODEL))
    for c in LEG_JOINTS:
        base_model.getCoordinateSet().get(c).set_locked(True)
    base = gravity_baseline(base_model, side="r")
    m_d = sum(base_model.getBodySet().get(n).getMass() for n in _DISTAL["r"])
    print(f"[T2] 自由落体残差 = {base:+.2f} N "
          f"(远端自重 {m_d*_G:.1f} N；应远小于冲击峰值)")

    grf = ground_reaction(height_m=5.0)
    print(f"GRF: peak={grf.peak_total_n/1e3:.1f} kN, 窗口 {grf.window_s[-1]*1e3:.1f} ms, "
          f"T1 err={grf.impulse_rel_err*100:+.1f}%")
    fall = run_dead_drop(grf)
    print(f"FD: {fall.n_states} 帧, 末时 {fall.t_s[-1]*1e3:.1f} ms")
    for side in ("r", "l"):
        jr = subtalar_reaction(fall, grf, side=side)
        print(f"  距下关节 {side}: peak|F|={jr.peak_force_n/1e3:.1f} kN "
              f"(纵向 {jr.peak_vertical_n/1e3:.1f} kN @ {jr.peak_vertical_time_s*1e3:.1f} ms), "
              f"peak|M|={jr.peak_moment_nm:.2f} N·m")
