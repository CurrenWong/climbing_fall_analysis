"""S1.1 —— 由 ``pad.py`` 的 1D 冲击模型生成地面反力时程 GRF(t)。

为什么需要这一层
----------------
``pad.py`` 只有**双质点 + 渐进接触**的 1D 模型，给的是整条通路的地面反力
时程 ``pad_force_n``。OpenSim 的 ``ExternalForce`` 需要的是
**每只脚**的时间序列 ``(t, F)``，且必须落在**单次冲击窗口**
（首触 → 首分离）内 —— 否则回弹重绷会污染峰值（见 ``pad.py`` 的窗口注释）。

本模块把"1D 结果"翻译成"OpenSim 可直接吃的 GRF"：

    输入：h, m, posture, on_pad
      → pad.simulate_boulder_fall(...)
      → 取 meta['impact_window_s'] 窗口
      → 平移到 t0 = 0
      → 左右各半（S1 假设对称）
    输出：GroundReaction(t_s, f_total_n, f_left_n, f_right_n, 冲量诊断)

不变量 T1
---------
``∫ F dt`` 应 ≈ 动量 ``m·√(2gh)``。
严格地说，接触窗口内地面对身体的冲量 = ``m·v0 + m·g·t_contact``
（重力在接触期间也在往下的方向积累动量，也必须被地面抵消），
故比值会略大于 1（对 5 m/20 ms 约 +2%）。这里如实报告相对误差，
阈值取 10%（方案 §七 T1）。

⚠️ 局限：``pad.py`` 的身体是 80 kg 的等效双质点，套到 OpenSim 的
75.337 kg 刚体模型上是**载荷边界条件的近似**，不是同一次动力学。
S1 只用它做"接触物理的唯一诚实来源"（方案 §3.1 D1）。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..pad import POSTURES, simulate_boulder_fall

__all__ = [
    "GroundReaction",
    "ground_reaction",
    "single_foot_split",
    "DEFAULT_POSTURE",
    "IMPULSE_TOL",
]

G_ACC = 9.81

#: S1 默认姿势：直立、双足先落、腿几乎绷直 —— 对应 [Y25] 的"成人直立位双足坠地"。
DEFAULT_POSTURE = "feet-first-stiff"

#: T1 冲量守恒相对误差阈值（方案 §七）。
IMPULSE_TOL = 0.10


@dataclass(frozen=True)
class GroundReaction:
    """一条 OpenSim 可直接施加的地面反力时程。

    Attributes
    ----------
    t_s
        时间 (s)，**以首次接触为 0**，单调递增。
    f_total_n, f_left_n, f_right_n
        总 / 左 / 右地面反力 (N)，向上为正。
    impulse_ns, momentum_ns, impulse_rel_err
        T1 诊断：``∫F dt``、``m·√(2gh)``、相对误差。
    """

    t_s: np.ndarray
    f_total_n: np.ndarray
    f_left_n: np.ndarray
    f_right_n: np.ndarray
    height_m: float
    mass_kg: float
    posture: str
    on_pad: bool
    impact_speed_ms: float
    peak_total_n: float
    impulse_ns: float
    momentum_ns: float
    impulse_rel_err: float
    window_s: tuple[float, float]
    split: tuple[float, float]
    #: 标量 GRF 时程作用于何种方向（地面/垫面法向，全局系、归一化、向上为正）。
    #: 默认 ``(0,1,0)`` = 竖直 —— 与非垂直扩展前的轴向行为**逐位一致**。
    ground_normal: tuple[float, float, float] = (0.0, 1.0, 0.0)
    #: S2：``True`` 表示左右力时程由 ``foot_forces`` 显式给出（>1 自由度，
    #: 允许一足恒 0 = 单脚）；``False``（默认）= 由 ``split`` 比例拆分标量 GRF。
    per_foot_forces: bool = False

    @property
    def duration_s(self) -> float:
        return float(self.t_s[-1] - self.t_s[0])

    @property
    def loaded_sides(self) -> tuple[str, ...]:
        """峰值力 > 0 的足侧（``'l'``/``'r'``）—— S2 单脚/不对称落地的判定。"""
        out = []
        if float(np.max(self.f_left_n)) > 0.0:
            out.append("l")
        if float(np.max(self.f_right_n)) > 0.0:
            out.append("r")
        return tuple(out)

    @property
    def peak_share(self) -> tuple[float, float]:
        """峰值意义下的 (左, 右) 力占比；总和为 0 时返回 ``(0,0)``。"""
        pl = float(np.max(self.f_left_n))
        pr = float(np.max(self.f_right_n))
        tot = pl + pr
        if tot <= 0.0:
            return (0.0, 0.0)
        return (pl / tot, pr / tot)

    @property
    def impulse_ok(self) -> bool:
        return abs(self.impulse_rel_err) <= IMPULSE_TOL

    @property
    def is_vertical(self) -> bool:
        """默认（竖直）方向判定 —— 供上行只走旧路径、保证逐位复现。"""
        return bool(np.allclose(np.asarray(self.ground_normal, float),
                                np.array([0.0, 1.0, 0.0])))

    def force_vec_n(self, side: str) -> np.ndarray:
        """每脚 3D GRF 时程 ``(N,3)``：标量大小 × ``ground_normal``。

        非垂直扩展的 L1 最小口径（方案 §3）：沿用现有标量大小 ``f(t)`` 作为
        沿 ``ground_normal`` 的力，不改接触本构。
        """
        if side == "r":
            scalar = self.f_right_n
        elif side == "l":
            scalar = self.f_left_n
        else:
            raise ValueError(f"side 必须是 'r'/'l'，得到 {side!r}")
        n = np.asarray(self.ground_normal, dtype=float)
        return np.asarray(scalar, dtype=float)[:, None] * n[None, :]

    def resample(self, n: int | None = None) -> "GroundReaction":
        """线性重采样到 ``n`` 点（默认 ≤ 2000，供 OpenSim 的 PWL 函数降采样）。"""
        if n is None:
            n = min(len(self.t_s), 2000)
        t = np.linspace(self.t_s[0], self.t_s[-1], n)
        return self._replace(t=t,
                             ft=np.interp(t, self.t_s, self.f_total_n),
                             fl=np.interp(t, self.t_s, self.f_left_n),
                             fr=np.interp(t, self.t_s, self.f_right_n))

    def _replace(self, t, ft, fl, fr) -> "GroundReaction":
        return GroundReaction(
            t_s=t, f_total_n=ft, f_left_n=fl, f_right_n=fr,
            height_m=self.height_m, mass_kg=self.mass_kg, posture=self.posture,
            on_pad=self.on_pad, impact_speed_ms=self.impact_speed_ms,
            peak_total_n=float(np.max(ft)),
            impulse_ns=float(np.trapezoid(ft, t)),
            momentum_ns=self.momentum_ns,
            impulse_rel_err=float(np.trapezoid(ft, t) / self.momentum_ns - 1.0),
            window_s=(float(t[0]), float(t[-1])), split=self.split,
            ground_normal=self.ground_normal,
            per_foot_forces=self.per_foot_forces,
        )


def _resolve_ground_normal(
    ground_normal: tuple[float, float, float],
    tilt_deg: float | None,
    roll_deg: float | None = None,
) -> tuple[float, float, float]:
    """把 ``ground_normal`` / ``tilt_deg`` / ``roll_deg`` 解析成归一化的全局系法向。

    * ``tilt_deg``（矢状面 / 前后向）：绕 z 轴把竖直法向倾斜 ``tilt_deg`` 度，
      即 ``n̂ = (sinθ, cosθ, 0)``（水平分量落在 +x，模拟前后向斜接触）。
    * ``roll_deg``（**冠状面 / 左右向**，S5 新增）：再绕 x 轴倾斜 ``roll_deg`` 度，
      给出**横向（±z）**分量，是踝旋后 / 内翻的接触侧来源。两者组合时先矢状后冠状：

          ``n̂ ∝ (sinθ, cosθ·cosφ, cosθ·sinφ)``   （θ=tilt, φ=roll）

      ``roll_deg=None``（默认）时 φ≡0，退化为扩展前的行为（**逐位一致**）。
    * 两者都 ``None`` 时直接用 ``ground_normal``。
    """
    if tilt_deg is not None or roll_deg is not None:
        th = np.deg2rad(float(tilt_deg) if tilt_deg is not None else 0.0)
        ph = np.deg2rad(float(roll_deg) if roll_deg is not None else 0.0)
        n = np.array(
            [np.sin(th), np.cos(th) * np.cos(ph), np.cos(th) * np.sin(ph)],
            dtype=float,
        )
    else:
        n = np.asarray(ground_normal, dtype=float)
    if n.shape != (3,) or not np.all(np.isfinite(n)):
        raise ValueError(f"ground_normal 必须是 3 个有限数，得到 {ground_normal!r}")
    norm = float(np.linalg.norm(n))
    if norm <= 0.0:
        raise ValueError("ground_normal 不能是零向量")
    return (float(n[0] / norm), float(n[1] / norm), float(n[2] / norm))


#: 单脚落地的 ``split`` 值（承力足占全部 GRF，另一足恒 0）。
_SINGLE_FOOT_SPLIT = {"l": (1.0, 0.0), "r": (0.0, 1.0)}


def single_foot_split(side: str) -> tuple[float, float]:
    """返回单脚落地的 ``split`` 值：``'l'`` → ``(1,0)``，``'r'`` → ``(0,1)``。

    ``ground_reaction(split=single_foot_split('r'))`` = **只有右脚承全部 GRF**
    （左脚恒 0），即方案 §2 自由度④的"单脚"最小实现。默认 ``(0.5,0.5)``
    仍是双足对称，行为不变。
    """
    if side not in _SINGLE_FOOT_SPLIT:
        raise ValueError(f"side 必须是 'l'/'r'，得到 {side!r}")
    return _SINGLE_FOOT_SPLIT[side]


def _resolve_foot_series(tw: np.ndarray, spec, side: str) -> np.ndarray:
    """把 ``foot_forces`` 单侧规格解析成与 ``tw`` 对齐的非负 GRF 时程 (N)。

    ``spec`` 三种形式（S2 泛化接口）：

    * **标量** → 整窗常数；
    * **一维数组** → 视作在 ``[tw[0], tw[-1]]`` 上均匀采样，线性重采样到 ``tw``；
    * **``(t, f)`` 二元组** → 按 ``t`` 线性插值到 ``tw``（窗外观测钳到端点）。

    仅做输入校验（有限、非负），不改接触物理。
    """
    if isinstance(spec, tuple) and len(spec) == 2:
        t_i = np.asarray(spec[0], dtype=float)
        f_i = np.asarray(spec[1], dtype=float)
        if t_i.ndim != 1 or f_i.ndim != 1 or t_i.size != f_i.size or t_i.size < 2:
            raise ValueError(
                f"foot_forces[{side!r}] 的 (t, f) 必须是等长一维数组（≥2 点）"
            )
        y = np.interp(tw, t_i, f_i)
    else:
        a = np.asarray(spec, dtype=float)
        if a.ndim == 0:
            y = np.full(tw.shape, float(a), dtype=float)
        elif a.ndim == 1:
            if a.size < 2:
                raise ValueError(f"foot_forces[{side!r}] 一维数组至少 2 点")
            t_u = np.linspace(float(tw[0]), float(tw[-1]), a.size)
            y = np.interp(tw, t_u, a)
        else:
            raise ValueError(f"foot_forces[{side!r}] 只接受标量 / 一维数组 / (t,f)")
    if not np.all(np.isfinite(y)):
        raise ValueError(f"foot_forces[{side!r}] 含 NaN/Inf")
    if np.any(y < 0.0):
        raise ValueError(f"foot_forces[{side!r}] 出现负力（地面反力不可为拉）")
    return y


def ground_reaction(
    *,
    height_m: float = 5.0,
    mass_kg: float = 75.337,
    posture: str = DEFAULT_POSTURE,
    on_pad: bool = False,
    dt_sample: float = 5e-5,
    split: tuple[float, float] = (0.5, 0.5),
    foot_forces: tuple | None = None,
    stop_at_zero_momentum: bool = True,
    ramp_s: float = 1e-3,
    ground_normal: tuple[float, float, float] = (0.0, 1.0, 0.0),
    tilt_deg: float | None = None,
    roll_deg: float | None = None,
) -> GroundReaction:
    """生成冲击窗口内的地面反力时程。

    Parameters
    ----------
    height_m
        释放高度 (m)。
    mass_kg
        身体质量 (kg)。**必须传 OpenSim 模型质量**，否则 T1 标定错位。
    posture
        ``pad.POSTURES`` 的键。默认 ``feet-first-stiff``。
    on_pad
        False = 刚性地面（S1 验证阶段）；True = 软垫（S5）。
    dt_sample
        采样步长 (s)。默认 5e-5，保证 ExternalForce 的分段线性足够光滑。
    split
        左右脚**分摊比例** ``(左, 右)``，非负、之和为 1。默认对称 ``(0.5, 0.5)``
        （与扩展前逐位一致）。**单脚**用 ``(1.0, 0.0)``（仅左）/ ``(0.0, 1.0)``
        （仅右）—— 见 :func:`single_foot_split`。S2 起该参数从"固定对称"推广为
        任意每足占比（含一足为 0）。
    foot_forces
        **可选**（默认 ``None``）逐足力时程 ``(左, 右)``（N），**完全覆盖** ``split``。
        每侧可为：标量（整窗常数）、一维数组（在窗口上均匀采样）、或 ``(t, f)``
        二元组（按 ``t`` 插值）。这是方案 §3 L1 的"两条独立时程"完全版：允许
        左右形状 / 峰值时刻不同，也允许一足恒 0（单脚）。``None`` → 走 ``split``
        路径（逐位复现旧行为）。总力 ``f_total_n = f_left + f_right``。
    stop_at_zero_momentum
        把时程截断在**全身动量首次归零**处（= 身体减速到静止），
        而不是弹回分离处。理由见下。默认 True。
    ramp_s
        截断后把残余力线性降到 0 的时长 (s)，避免 ExternalForce 出现阶跃。
    ground_normal
        标量 GRF 时程的作用方向（全局系、向上为正）。默认 ``(0, 1, 0)`` =
        竖直，**与非垂直扩展前的轴向行为逐位一致**。传入非竖直法向即启用
        L1 最小矢量通路：``F⃗(t) = f(t) · n̂``（方案 §3）。
    tilt_deg
        地面/垫面**矢状面**倾角（度）。非 ``None`` 时覆盖 ``ground_normal``，令
        ``n̂ = (sinθ, cosθ, 0)``。默认 ``None`` = 竖直。
    roll_deg
        **冠状面（左右向）**地面/垫面倾角（度），S5 新增的 opt-in 接口。非 ``None``
        时在 ``tilt_deg`` 之后**再绕 x 轴倾斜**，给法向一个横向（±z）分量：
        ``n̂ ∝ (sinθ, cosθ·cosφ, cosθ·sinφ)``（φ=roll）。这是踝旋后 / 内翻的
        **接触侧来源**（配合 ``ankle_supination`` 的显式旋后角使用）。
        **默认 ``None`` = 不加冠状分量 → 与扩展前逐位一致。**

    Notes
    -----
    **为什么默认截断在"动量归零"而不是"力归零（分离）"**：

    ``pad.hard_surface()`` 是**纯弹性的**（`linear_modulus_pa`，无接触阻尼），
    所以 5 m 落地会把身体**以约 7–10 m/s 弹回** —— 实测整窗
    ``∫F dt ≈ 2·m·v0``（速度 +9.9 → −7.1 m/s）。这是标定缺陷（见 README
    「硬地面刚度标定最弱」），不是物理。

    而 [Y25] 的场景是**人摔在刚性地面上并停住**（软组织耗散），不是反弹。
    因此 S1 取"首触 → 全身动量首次归零"这一段：此时
    ``∫(F − m·g)dt = m·v0`` **精确成立**，且 ``∫F dt / (m·v0) ≈ 1 + g·t/v0``
    （5 m 时约 +2.5%），T1 的 10% 阈值才有意义。

    注意：峰值力发生在动量归零**之前**，故截断不影响峰值。
    """
    if posture not in POSTURES:
        raise KeyError(f"未知姿势 {posture!r}；可选：{sorted(POSTURES)}")
    if foot_forces is None:
        if not np.isclose(sum(split), 1.0):
            raise ValueError(f"split 之和必须为 1，得到 {split}")
        if split[0] < 0.0 or split[1] < 0.0:
            raise ValueError(f"split 不可为负，得到 {split}")
    else:
        if not (isinstance(foot_forces, tuple) and len(foot_forces) == 2):
            raise ValueError("foot_forces 必须是 (左, 右) 二元组")

    normal = _resolve_ground_normal(ground_normal, tilt_deg, roll_deg)

    res = simulate_boulder_fall(
        height_m=height_m, mass_kg=mass_kg, posture=posture,
        on_pad=on_pad, dt_sample=dt_sample,
    )

    win = res.meta.get("impact_window_s")
    if win is None:
        raise RuntimeError("pad 结果缺少 meta['impact_window_s']，无法定位冲击窗口。")
    t0, t1 = float(win[0]), float(win[1])

    t = np.asarray(res.t_s, dtype=float)
    f = np.asarray(res.pad_force_n, dtype=float)
    # 取冲击窗口，并平移到 t=0
    sel = (t >= t0 - 1e-12) & (t <= t1 + 1e-12)
    tw = t[sel] - t0
    fw = np.maximum(f[sel], 0.0)

    # ---- 截断到"全身动量首次归零"（= 减速到静止，见 docstring Notes）---------
    if stop_at_zero_momentum:
        ml = res.meta.get("m_low")
        mu = res.meta.get("m_up")
        if ml and mu:
            zf = np.asarray(res.z_foot_m, dtype=float)[sel]
            zt = np.asarray(res.z_torso_m, dtype=float)[sel]
            # pad.py 约定：z 向下为正 → 动量向下为正，冲击中 p 由 +m·v0 降到 0 再变负
            vzf = np.gradient(zf, tw)
            vzt = np.gradient(zt, tw)
            p = ml * vzf + mu * vzt
            crossed = np.nonzero(p <= 0.0)[0]
            if crossed.size and crossed[0] > 1:
                k = int(crossed[0])
                tw, fw = tw[:k], fw[:k]
                # 残余力在 ramp_s 内线性降到 0，避免阶跃
                if ramp_s > 0 and fw[-1] > 0.0:
                    n_ramp = 10
                    tw = np.concatenate([tw, tw[-1] + np.linspace(ramp_s / n_ramp, ramp_s, n_ramp)])
                    fw = np.concatenate([fw, np.linspace(fw[-1], 0.0, n_ramp)])

    if tw.size < 4:
        raise RuntimeError(
            f"冲击窗口只有 {tw.size} 个采样点（t0={t0:.4f}, t1={t1:.4f}）—— "
            "窗口定位可疑，拒绝静默返回。"
        )
    if tw[0] > 0.0:
        tw = np.concatenate([[0.0], tw])
        fw = np.concatenate([[fw[0]], fw])

    v0 = float(res.impact_speed_ms)
    momentum = mass_kg * v0

    # ---- 左右力：默认 = split 比例拆分（逐位复现）；foot_forces = 独立逐足时程 ----
    if foot_forces is None:
        f_total = fw
        f_left = fw * split[0]
        f_right = fw * split[1]
        stored_split = split
        per_foot = False
    else:
        f_left = _resolve_foot_series(tw, foot_forces[0], "l")
        f_right = _resolve_foot_series(tw, foot_forces[1], "r")
        f_total = f_left + f_right
        peak = float(np.max(f_total)) if f_total.size else 0.0
        stored_split = (
            (float(np.max(f_left)) / peak, float(np.max(f_right)) / peak)
            if peak > 0.0 else (0.0, 0.0)
        )
        per_foot = True

    impulse = float(np.trapezoid(f_total, tw))

    return GroundReaction(
        t_s=tw,
        f_total_n=f_total,
        f_left_n=f_left,
        f_right_n=f_right,
        height_m=height_m,
        mass_kg=mass_kg,
        posture=posture,
        on_pad=on_pad,
        impact_speed_ms=v0,
        peak_total_n=float(np.max(f_total)),
        impulse_ns=impulse,
        momentum_ns=momentum,
        impulse_rel_err=impulse / momentum - 1.0,
        window_s=(float(tw[0]), float(tw[-1])),
        split=stored_split,
        ground_normal=normal,
        per_foot_forces=per_foot,
    )


# --------------------------------------------------------------------------
# 自检：python -m climbing.coupling.opensim_grf
# --------------------------------------------------------------------------
def _main() -> int:
    m_model = 75.337
    for h in (1.0, 5.0, 10.0):
        grf = ground_reaction(height_m=h, mass_kg=m_model)
        flag = "✅" if grf.impulse_ok else "❌"
        print(
            f"h={h:5.1f}m  v0={grf.impact_speed_ms:5.2f} m/s  "
            f"t_window={grf.window_s[0]*1e3:5.2f}→{grf.window_s[1]*1e3:5.2f} ms  "
            f"peak={grf.peak_total_n/1e3:6.1f} kN  "
            f"∫Fdt={grf.impulse_ns:6.1f} N·s  m·v0={grf.momentum_ns:6.1f} N·s  "
            f"err={grf.impulse_rel_err*100:+5.1f}% {flag}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
