"""S3.2 —— 关节 wrench + 跟腱力 → FEBio 单骨载荷 BC。

把「多体段产出的关节力」与「跟腱力」翻译成 :func:`febio_model.build_calcaneus_feb`
的入参，集中管理方向/单位/符号，避免散落在各脚本里。单位 mm–N–MPa–s。

为什么跟腱力不能从当前多体模型直接取（S3.2b 的诚实边界）
------------------------------------------------------
S1 的 dead-drop：**下肢关节锁死 + 肌肉关闭**；GRF 是 1D 竖直力、作用于距下关节
中心（``pad.py`` 无 CoP）。因此模型自身**不含踝跖屈力矩**，也就没有跟腱力的来源
（``subtalar_reaction`` 只给距下关节反力）。在 S3 恢复肌肉/关节顺应性之前，
S3.2b 先用**参数化锚定**：

    F_achilles = a · F_max(triceps surae),   a ∈ [0, 1]

``a`` = 三头肌激活水平（a=0 无跟腱、a=1 最大等长）。F_max 取 Rajagopal2015
``soleus_r + gasmed_r + gaslat_r`` 的 ``getMaxIsometricForce()`` 之和 = **10885 N**
（实测，见 ``temp/opensim_fe/probe_achilles.py``）。这是唯一有**生理上界**、
可复现、不依赖具体激活时程的口径 —— 用于检验「补上跟腱能否把首次骨折高度
从 ~26.5 m 拉回 [Y25] 的 7–9 m」这个假设。
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = [
    "FeLoadSpec",
    "achilles_force",
    "transfer",
    "SUBTALAR_DIR",
    "ACHILLES_DIR",
    "F_MAX_TRICEPS_N",
]

#: 距下关节（talus）对跟骨的合力方向：向下。
SUBTALAR_DIR: tuple[float, float, float] = (0.0, -1.0, 0.0)
#: 跟腱对跟骨结节（后上）的力方向：向上、略偏前（≈ 腱线作用方向）。
#: 由 Rajagopal 三头肌几何路径实测：soleus (+0.145,+0.989,+0.022)、gastro 类似。
ACHILLES_DIR: tuple[float, float, float] = (0.15, 0.985, 0.0)
#: 三头肌最大等长力之和（N）——Rajagopal2015 实测。
F_MAX_TRICEPS_N: float = 10885.0


def _unit(v: tuple[float, float, float]) -> tuple[float, float, float]:
    """把 3 向量归一化（零向量原样返回，交由调用方决定兜底）。"""
    n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) ** 0.5
    if n <= 0.0:
        return (0.0, 0.0, 0.0)
    return (v[0] / n, v[1] / n, v[2] / n)


@dataclass(frozen=True)
class FeLoadSpec:
    """一次单骨 FE 的载荷规格（可直接喂给 :func:`febio_model.build_calcaneus_feb`）。"""

    subtalar_n: float
    achilles_n: float
    subtalar_dir: tuple[float, float, float] = SUBTALAR_DIR
    achilles_dir: tuple[float, float, float] = ACHILLES_DIR
    #: L3（非垂直扩展）：完整 3D 距下力 (N,3)。``None``（默认）时由
    #: ``subtalar_n × subtalar_dir`` 复现旧的标量行为。
    subtalar_force: tuple[float, float, float] | None = None
    #: L3：完整 3D 距下力矩 (N·m,3)，对距下中心。``None`` = 无力矩（旧行为）。
    subtalar_moment: tuple[float, float, float] | None = None

    @property
    def force_vector_n(self) -> tuple[float, float, float]:
        """保留完整 3D 力（N）。默认路径 = ``subtalar_n × unit(subtalar_dir)``。"""
        if self.subtalar_force is not None:
            return self.subtalar_force
        d = _unit(self.subtalar_dir)
        return (self.subtalar_n * d[0], self.subtalar_n * d[1], self.subtalar_n * d[2])

    @property
    def moment_vector_nm(self) -> tuple[float, float, float]:
        """保留完整 3D 力矩 (N·m)；默认路径 = 零力偶。"""
        if self.subtalar_moment is not None:
            return self.subtalar_moment
        return (0.0, 0.0, 0.0)

    def build(self, mesh, out_feb, **kwargs):
        """用本规格写出 ``.feb``（默认压力 + 跟腱牵引，静态 1 步）。

        默认路径（``subtalar_force is None``）逐位复用旧标量接口。
        传入完整 3D 力时，退回用其**大小 + 方向**喂当前 FE 适配器
        （力矩 / 偏心映射是 S4，见方案 §4；此处不静默旋转方向）。
        """
        from .febio_model import build_calcaneus_feb

        if self.subtalar_force is None:
            load_n, load_dir = self.subtalar_n, self.subtalar_dir
        else:
            import numpy as np

            f = np.asarray(self.subtalar_force, dtype=float)
            mag = float(np.linalg.norm(f))
            load_n = mag
            load_dir = tuple(f / mag) if mag > 0.0 else self.subtalar_dir

        kwargs.setdefault("use_rigid", False)
        kwargs.setdefault("time_steps", 1)
        return build_calcaneus_feb(
            mesh,
            out_feb,
            load_n=load_n,
            load_dir=load_dir,
            achilles_n=self.achilles_n,
            achilles_dir=self.achilles_dir,
            **kwargs,
        )


def achilles_force(activation: float, *, f_max_n: float = F_MAX_TRICEPS_N) -> float:
    """跟腱力 = 激活 × F_max（N）。``activation`` 必须在 [0, 1]。"""
    a = float(activation)
    if not 0.0 <= a <= 1.0:
        raise ValueError(f"activation 必须在 [0, 1]，得到 {activation!r}")
    return a * float(f_max_n)


def transfer(
    subtalar_n: float,
    *,
    achilles_activation: float = 0.0,
    f_max_n: float = F_MAX_TRICEPS_N,
    subtalar_dir: tuple[float, float, float] = SUBTALAR_DIR,
    achilles_dir: tuple[float, float, float] = ACHILLES_DIR,
    subtalar_force: tuple[float, float, float] | None = None,
    subtalar_moment: tuple[float, float, float] | None = None,
) -> FeLoadSpec:
    """把距下关节峰值力 + 跟腱激活水平 → :class:`FeLoadSpec`。

    Parameters
    ----------
    subtalar_n
        距下关节峰值纵向力（N），来自 ``joint_loads.subtalar_reaction(...).peak_vertical_n``。
    achilles_activation
        三头肌激活 ``a``（0=无跟腱，1=最大等长）。见模块 docstring 的 S3.2b 说明。
    subtalar_force
        **可选** 完整 3D 距下力 (N,3)。默认 ``None`` → 走旧标量口径
        ``subtalar_n × subtalar_dir``，行为逐位不变。传入时 spec 保留该矢量
        （见 :attr:`FeLoadSpec.force_vector_n`），不再 scalarize。
    subtalar_moment
        **可选** 完整 3D 力距 (N·m,3)。默认 ``None`` = 零力偶（旧行为）。
    """
    return FeLoadSpec(
        subtalar_n=float(subtalar_n),
        achilles_n=achilles_force(achilles_activation, f_max_n=f_max_n),
        subtalar_dir=tuple(subtalar_dir),
        achilles_dir=tuple(achilles_dir),
        subtalar_force=None if subtalar_force is None else tuple(subtalar_force),
        subtalar_moment=None if subtalar_moment is None else tuple(subtalar_moment),
    )


def f_max_triceps(*, model_path=None) -> float:
    """从 Rajagopal2015 实测三头肌 F_max 之和（N）——可追溯，不用常数硬编码时调用。"""
    import opensim as os

    from .opensim_fall import DEFAULT_MODEL

    m = os.Model(str(model_path or DEFAULT_MODEL))
    total = 0.0
    for name in ("soleus_r", "gasmed_r", "gaslat_r"):
        total += m.getMuscles().get(name).getMaxIsometricForce()
    return float(total)
