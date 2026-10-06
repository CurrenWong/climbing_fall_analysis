"""Phase-S3 非垂直落地（自由度③：落地姿势）回归测试。

覆盖两条路径（方案 §3 L2）：

* **Path (a)** —— 锁关节到目标位姿（``LegPose``）；刚性腿 + 弯曲几何。
* **Path (b)** —— 解锁 + 被动刚度（``PassiveStiffness``）；关节有限弯曲吸能。

**硬门槛（方案 §7）**：默认 ``posture=None, passive=None`` 必须数值复现扩展前
轴向结果。快测试锁定接口契约（无需正动力学积分）；慢测试（``RUN_OPEN_FE=1``）
直接对比改前缓存 ``results/opensim_fe/s1_h5_opensim.json``，并实测 (a)/(b)。

慢测试放在 ``RUN_OPEN_FE`` 门后，与 ``tests/test_nonvertical_s1.py`` /
``test_nonvertical_s2.py`` 同一约定。
"""

from __future__ import annotations

import inspect
import json
import math
import os
import sys
from pathlib import Path

import numpy as np
import pytest

from climbing.coupling.opensim_fall import (
    LEG_FLEX_JOINTS,
    LEG_JOINTS,
    LegPose,
    PassiveStiffness,
    run_dead_drop,
)

M_MODEL = 75.337
ROOT = Path(__file__).resolve().parents[1]
CACHED_AXIAL = ROOT / "results" / "opensim_fe" / "s1_h5_opensim.json"

#: 供慢测试复用 S3 脚本的坐标→状态路径解析
if str(ROOT / "scripts" / "opensim_fe") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts" / "opensim_fe"))

run_open = pytest.mark.skipif(
    os.environ.get("RUN_OPEN_FE") != "1",
    reason="慢（正动力学 ~10-20 s）；设 RUN_OPEN_FE=1 启用",
)


# --------------------------------------------------------------------------
# 接口：签名默认 = 旧行为（向后兼容）
# --------------------------------------------------------------------------
def test_run_dead_drop_signature_defaults() -> None:
    sig = inspect.signature(run_dead_drop)
    assert sig.parameters["posture"].default is None
    assert sig.parameters["passive"].default is None
    # 既有参数不动（向后兼容）
    assert sig.parameters["lock_leg_joints"].default is True
    assert sig.parameters["vx0_ms"].default == 0.0
    assert sig.parameters["vz0_ms"].default == 0.0
    assert sig.parameters["scale_mass_kg"].default is None


def test_leg_flex_joints_subset_of_leg_joints() -> None:
    assert set(LEG_FLEX_JOINTS).issubset(set(LEG_JOINTS))
    # 差集 = 非矢状面（髋内收/旋转），S3 不解锁
    assert set(LEG_JOINTS) - set(LEG_FLEX_JOINTS) == {
        "hip_adduction_r", "hip_adduction_l",
        "hip_rotation_r", "hip_rotation_l",
    }


# --------------------------------------------------------------------------
# LegPose —— 默认 0（= 旧行为）
# --------------------------------------------------------------------------
def test_legpose_default_is_zero() -> None:
    p = LegPose()
    assert p.knee_flexion_rad == 0.0
    assert p.hip_flexion_rad == 0.0
    assert p.ankle_rad == 0.0
    assert p.torso_lean_rad == 0.0  # 默认不动躯干（逐位一致）


def test_legpose_torso_lean() -> None:
    p = LegPose.flexed_landing(knee_deg=45.0, torso_lean_deg=-20.0)
    assert p.torso_lean_rad == pytest.approx(math.radians(-20.0))
    # as_default_value_dict 不含躯干（躯干由 _apply_leg_pose_locked 单独处理）
    assert "lumbar_extension" not in p.as_default_value_dict()


def test_legpose_flexed_landing_converts_degrees() -> None:
    p = LegPose.flexed_landing(knee_deg=30.0, hip_deg=15.0, ankle_deg=10.0)
    assert p.knee_flexion_rad == pytest.approx(math.radians(30.0))
    assert p.hip_flexion_rad == pytest.approx(math.radians(15.0))
    assert p.ankle_rad == pytest.approx(math.radians(10.0))


def test_legpose_default_value_dict_symmetric() -> None:
    d = LegPose.flexed_landing(knee_deg=45.0, hip_deg=20.0, ankle_deg=5.0).as_default_value_dict()
    assert d["knee_angle_r"] == d["knee_angle_l"] == pytest.approx(math.radians(45.0))
    assert d["hip_flexion_r"] == d["hip_flexion_l"] == pytest.approx(math.radians(20.0))
    assert d["ankle_angle_r"] == d["ankle_angle_l"] == pytest.approx(math.radians(5.0))
    assert set(d) == set(LEG_FLEX_JOINTS)


# --------------------------------------------------------------------------
# PassiveStiffness —— 默认关闭；预设开启
# --------------------------------------------------------------------------
def test_passive_default_inactive() -> None:
    assert PassiveStiffness().is_active() is False


def test_passive_nominal_and_soft_active() -> None:
    assert PassiveStiffness.nominal_landing().is_active() is True
    assert PassiveStiffness.soft_landing().is_active() is True


def test_passive_nominal_values_labeled_assumed() -> None:
    """nominal 为一组正刚度/阻尼；数值全部 > 0（ASSUMED，注释在 docstring）。"""
    p = PassiveStiffness.nominal_landing()
    for v in (
        p.knee_stiffness, p.hip_stiffness, p.ankle_stiffness,
        p.knee_damping, p.hip_damping, p.ankle_damping,
        p.knee_limit_stiffness, p.hip_limit_stiffness, p.ankle_limit_stiffness,
    ):
        assert v > 0.0
    assert p.transition_rad > 0.0


# --------------------------------------------------------------------------
# 错误契约（在 os.Model 载入**之前**抛出 → 无需积分）
# --------------------------------------------------------------------------
def test_posture_and_passive_mutually_exclusive() -> None:
    with pytest.raises(ValueError):
        run_dead_drop(
            None,  # 显式冲突在触碰 grf 之前抛出
            posture=LegPose(),
            passive=PassiveStiffness.nominal_landing(),
        )


def test_passive_inactive_raises() -> None:
    with pytest.raises(ValueError):
        run_dead_drop(None, passive=PassiveStiffness())


def test_legpose_is_frozen() -> None:
    p = LegPose()
    with pytest.raises(Exception):
        p.knee_flexion_rad = 1.0  # frozen dataclass


# --------------------------------------------------------------------------
# HARD GATE —— 默认路径数值复现既有轴向结果（慢，默认跳过）
# --------------------------------------------------------------------------
@run_open
def test_dead_drop_default_reproduces_cached_axial() -> None:
    from climbing.coupling.joint_loads import subtalar_reaction
    from climbing.coupling.opensim_grf import ground_reaction

    assert CACHED_AXIAL.is_file(), f"缺少轴向回归基准 {CACHED_AXIAL}"
    cached = json.loads(CACHED_AXIAL.read_text(encoding="utf-8"))

    grf = ground_reaction(height_m=5.0, mass_kg=M_MODEL)
    assert grf.impulse_ns == pytest.approx(cached["impulse_ns"], rel=1e-9)
    fall = run_dead_drop(grf)
    jr = subtalar_reaction(fall, grf, side="r")
    assert jr.peak_vertical_n == pytest.approx(
        cached["subtalar_peak_vertical_n"], rel=1e-6)
    assert jr.peak_force_n == pytest.approx(
        cached["subtalar_peak_force_n"], rel=1e-6)


@run_open
def test_posture_zero_bitwise_reproduces_default() -> None:
    """``posture=LegPose()``（全 0）必须与默认**逐位一致**（默认契约）。"""
    from climbing.coupling.joint_loads import subtalar_reaction
    from climbing.coupling.opensim_grf import ground_reaction

    grf = ground_reaction(height_m=5.0, mass_kg=M_MODEL)
    fall_d = run_dead_drop(grf)
    fall_p0 = run_dead_drop(grf, posture=LegPose())
    jr_d = subtalar_reaction(fall_d, grf, side="r")
    jr_p0 = subtalar_reaction(fall_p0, grf, side="r")
    assert jr_p0.peak_vertical_n == pytest.approx(jr_d.peak_vertical_n, rel=1e-12)
    assert jr_p0.peak_force_n == pytest.approx(jr_d.peak_force_n, rel=1e-12)


@run_open
def test_path_a_flexed_lock_changes_loads_and_joints_stay_locked() -> None:
    """Path (a)：锁在膝屈 45°；关节反力变化，且关节坐标严格 = 45°（锁定）。"""
    from climbing.coupling.joint_loads import subtalar_reaction
    from climbing.coupling.opensim_grf import ground_reaction

    grf = ground_reaction(height_m=5.0, mass_kg=M_MODEL)
    fall = run_dead_drop(grf, posture=LegPose.flexed_landing(knee_deg=45.0))
    # 关节轨迹里 knee_angle_r 应恒为 45°（锁定 = 不随时间变）
    tab = fall.states
    labels = [tab.getColumnLabel(j) for j in range(tab.getNumColumns())]
    path = [l for l in labels if l.endswith("/knee_angle_r/value") and "_beta" not in l][0]
    v = np.asarray(tab.getDependentColumnAtIndex(labels.index(path)).to_numpy())
    assert np.allclose(v, math.radians(45.0), atol=1e-9)
    # 反力为正、量级与施加 GRF 可比
    jr = subtalar_reaction(fall, grf, side="r")
    assert jr.peak_force_n > 0.3 * grf.f_right_n.max()


@run_open
def test_path_b_joints_flex_finitely_and_stay_in_range() -> None:
    """Path (b)：解锁 + nominal 被动刚度；关节**有限**弯曲且不越解剖范围。"""
    from climbing.coupling.opensim_grf import ground_reaction

    grf = ground_reaction(height_m=5.0, mass_kg=M_MODEL)
    fall = run_dead_drop(grf, passive=PassiveStiffness.nominal_landing())

    tab = fall.states
    labels = [tab.getColumnLabel(j) for j in range(tab.getNumColumns())]
    ranges = {
        "knee_angle_r": (-0.0, 2.0944),
        "hip_flexion_r": (-0.52359878, 2.0943951),
        "ankle_angle_r": (-0.6981317, 0.52359878),
    }
    for coord, (lo, hi) in ranges.items():
        path = [l for l in labels if l.endswith(f"/{coord}/value") and "_beta" not in l]
        assert path, f"缺少状态列 {coord}"
        v = np.asarray(tab.getDependentColumnAtIndex(labels.index(path[0])).to_numpy())
        # 有限弯曲：有位移（> 0）但远小于范围宽
        excursion = float(max(abs(v.min()), abs(v.max())))
        assert excursion > 0.0, f"{coord} 完全没动（被动刚度路径没生效？）"
        assert v.min() >= lo - 1e-6, f"{coord} 下越界 ({v.min()} < {lo})"
        assert v.max() <= hi + 1e-6, f"{coord} 上越界 ({v.max()} > {hi})"
        assert excursion <= 0.2, f"{coord} 位移 {excursion:.3f} rad 过大（疑似飞）"
