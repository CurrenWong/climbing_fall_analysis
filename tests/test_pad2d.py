"""P3 —— 抱石 2D 模型 pad2d 的回归测试。

四个不变量（来自 ``docs/方案更新_v2.md:896-902``）：

1. **角度 → 0 收敛回纯压缩模型** —— 默认 ``Posture2D()`` 与
   ``pad.simulate_boulder_fall`` 的 1D 结果 ``np.array_equal``
   （**不是** ``pytest.approx``）。这是 pad.py 的"硬门槛"约定。
2. **摩擦两端行为可解释** —— ``mu_foot=0`` ⇒ ``x_foot_m`` 保持 0、
   ``theta_rad`` 不被摩擦扰动；``mu_foot=1e9`` ⇒ ``|vx_foot| ≈ 0``
   （粘）。两条都需要注释解释**为什么**是这样。
3. **内翻 vs 外翻相对关系稳定** —— ``theta0 = ±X`` 的 peak rotational
   quantity 在 3+ 高度上**排序一致**（不假设哪一侧更重）。
4. **峰值时刻随高度单调提前** —— peak force 时刻距接触起点的
   延迟 dt_after_fall 随高度**单调非增**。

外加：
- 1D bit-identicality 回归（与不变量 1 同源，但用更多姿势覆盖）
- 下游消费者兼容性（``injury.py`` / ``bone.py`` / ``metrics.py``
  必须能消费 ``BoulderFallResult2D`` 而无感知差异）
- 接口契约：``z_torso[0] == -0.85``、shape 匹配等

约定（沿用 ``tests/test_nonvertical_s5.py`` 的"硬门槛"风格）
------------------------------------------------------------
* bit-identicality 必须用 ``np.array_equal`` —— ``pytest.approx``
  会把"对的不等"和"错的近似"混在一起，丢掉这条契约。
* 每个不变量都要有**正交**断言：只测一个量，不顺手测别的。
* Chinese 注释作为 house style —— 与 pad.py / tests/test_climbing.py
  保持一致。
"""

from __future__ import annotations

import inspect

import numpy as np
import pytest

from climbing import G
from climbing.coupling.ankle_ligament import (
    ATFL_FAILURE_STRAIN,
    ATFL_INVERSION_FAIL_MOMENT_NO_PRELOAD_NM,
)
from climbing.pad import (
    BoulderFallResult,
    CrashPad,
    POSTURES,
    simulate_boulder_fall,
)
from climbing.pad2d import (
    LEG_EXT0_M,
    BoulderFallResult2D,
    Posture2D,
    sign_reg,
    simulate_boulder_fall_2d,
)


# ==========================================================================
# 工具：取"冲击窗口内的峰值时刻"（避开反弹重峰污染）
# ==========================================================================
def impact_window_peak_time(result) -> float:
    """取 ``pad_force_n`` 在冲击窗口内的峰值时刻（s）。

    直接用 ``np.argmax`` 会拿到反弹重峰（h≥1 m 时常见，比初始峰高），
    反而**丢失**初始冲击的物理意义。pad.py 在 ``result.meta`` 里给了
    完整的"首次接触 → 首次分离"窗口，本函数用它来定位初始峰。
    """
    t = result.t_s
    f = result.pad_force_n
    win = result.meta.get("impact_window_s")
    assert win is not None, "result.meta['impact_window_s'] 缺失"
    t0, t1 = float(win[0]), float(win[1])
    sl = (t >= t0) & (t <= t1)
    if not sl.any():
        return float("nan")
    return float(t[sl][int(np.argmax(f[sl]))])


def contact_to_peak_time(result) -> float:
    """从首次接触到峰值力所经历的时间（s）。

    等价于 "冲击窗口内的 argmax 时刻 − 首次接触时刻"。这是物理上的
    "接触-峰延迟"，独立于自由落体时长 sqrt(2h/g)。

    [T22] Table 3 的口径：**这个量**随高度单调下降（接触更"硬"、
    更快达到峰值），而不是绝对峰值时刻。
    """
    t = result.t_s
    win = result.meta.get("impact_window_s")
    assert win is not None
    first_contact = float(win[0])
    return impact_window_peak_time(result) - first_contact


# ==========================================================================
# sign_reg：Coulomb 摩擦的正则化（避免 RK45 刚性问题）
# ==========================================================================
class TestSignReg:
    def test_positive_input_gives_positive(self):
        assert sign_reg(1.0) > 0.0
        assert sign_reg(0.1) > 0.0

    def test_negative_input_gives_negative(self):
        assert sign_reg(-1.0) < 0.0
        assert sign_reg(-0.1) < 0.0

    def test_at_zero_returns_zero(self):
        # sign_reg(0) = 0 / sqrt(eps²) = 0 —— 与 sign(0) 一致
        assert sign_reg(0.0) == 0.0

    def test_bounded_by_unit(self):
        """正则化 sign 的绝对值永远 ≤ 1。"""
        for v in [-1e6, -100.0, -10.0, -1e-3, 0.0, 1e-3, 10.0, 100.0, 1e6]:
            assert -1.0 <= sign_reg(v) <= 1.0

    def test_smooth_at_zero(self):
        """``|v| ≪ eps`` 范围内 ≈ v/eps，无不连续（用极小 v 保证线性区）。"""
        v = 1e-6                # ≪ eps = 1e-3 ⇒ v² 可忽略
        # 渐近展开 sign_reg(v) = (v/eps)·(1 − (v/eps)²/2 + …)
        # ⇒ 相对误差 ≈ (v/eps)²/2 = 5e-7。用 rel=1e-9 是**过紧**的（实测差 5e-7）。
        assert sign_reg(v) == pytest.approx(v / 1e-3, rel=1e-5)
        assert sign_reg(-v) == pytest.approx(-v / 1e-3, rel=1e-5)

    def test_custom_eps(self):
        """eps 可调。"""
        v = 0.5
        # eps 很小时 sign_reg 接近 sign，但 |v| = 0.5 >> 0.001 ⇒ 几乎为 1
        assert sign_reg(v, eps=1e-6) == pytest.approx(1.0, rel=1e-6)
        assert sign_reg(-v, eps=1e-6) == pytest.approx(-1.0, rel=1e-6)

    def test_array_input(self):
        """支持 numpy 数组（ODE RHS 内部使用）。"""
        v = np.array([-2.0, -0.5, 0.0, 0.5, 2.0])
        out = sign_reg(v)
        assert out.shape == v.shape
        # 严格保序
        assert np.all(np.diff(out) > 0)


# ==========================================================================
# Posture2D：默认值 + frozen + 可哈希
# ==========================================================================
class TestPosture2D:
    def test_all_defaults_are_zero(self):
        """所有字段默认 0.0 —— 这是默认与 1D bit-identical 的前提。"""
        p = Posture2D()
        for f in (
            "theta0_rad", "omega0_rad_s", "lateral_offset_m",
            "foot_x_lever_m", "mu_foot", "I_body_kgm2", "c_rot",
        ):
            assert getattr(p, f) == 0.0, f"{f} 默认值非零：{getattr(p, f)}"

    def test_is_frozen(self):
        """frozen dataclass —— 防止被意外修改。"""
        from dataclasses import FrozenInstanceError
        p = Posture2D()
        with pytest.raises(FrozenInstanceError):
            p.theta0_rad = 0.1

    def test_hashable(self):
        """frozen ⇒ 可哈希。可作缓存 key、集合元素、字典 key。"""
        p = Posture2D(mu_foot=0.5, theta0_rad=0.1)
        hash(p)  # 不抛错即通过

    def test_default_equals_all_explicit_zeros(self):
        """默认构造 ≡ 所有字段显式给 0.0。"""
        p1 = Posture2D()
        p2 = Posture2D(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        assert p1 == p2

    def test_default_does_not_equal_nonzero(self):
        p = Posture2D()
        assert p != Posture2D(mu_foot=0.1)
        assert p != Posture2D(theta0_rad=0.01)
        assert p != Posture2D(foot_x_lever_m=0.05)


# ==========================================================================
# 不变量 1 —— 角度→0 收敛回纯压缩模型（bit-identicality，硬门槛）
# ==========================================================================
class TestInvariant1BitIdentical:
    """默认 ``Posture2D()`` 必须给出与 1D 完全相同的轨迹。

    硬门槛：用 ``np.array_equal`` 而**不是** ``pytest.approx``。
    任何近似等价都说明 1D 子系统被改写了。
    """

    @pytest.mark.parametrize("posture", [
        "controlled-drop", "feet-first-stiff", "butt-impact",
        "flat-flop", "tuck-roll", "head-first",
    ])
    def test_z_foot_bit_identical(self, posture):
        r1d = simulate_boulder_fall(height_m=3.0, mass_kg=80.0, posture=posture)
        r2d = simulate_boulder_fall_2d(
            height_m=3.0, mass_kg=80.0, posture=posture,
            posture2d=Posture2D(),
        )
        assert np.array_equal(r2d.z_foot_m, r1d.z_foot_m)

    @pytest.mark.parametrize("posture", [
        "controlled-drop", "feet-first-stiff", "head-first",
    ])
    def test_z_torso_bit_identical(self, posture):
        r1d = simulate_boulder_fall(height_m=3.0, mass_kg=80.0, posture=posture)
        r2d = simulate_boulder_fall_2d(
            height_m=3.0, mass_kg=80.0, posture=posture,
            posture2d=Posture2D(),
        )
        assert np.array_equal(r2d.z_torso_m, r1d.z_torso_m)

    @pytest.mark.parametrize("posture", [
        "controlled-drop", "feet-first-stiff", "flat-flop", "head-first",
    ])
    def test_pad_force_bit_identical(self, posture):
        r1d = simulate_boulder_fall(height_m=3.0, mass_kg=80.0, posture=posture)
        r2d = simulate_boulder_fall_2d(
            height_m=3.0, mass_kg=80.0, posture=posture,
            posture2d=Posture2D(),
        )
        assert np.array_equal(r2d.pad_force_n, r1d.pad_force_n)

    @pytest.mark.parametrize("posture", [
        "controlled-drop", "feet-first-stiff", "butt-impact",
    ])
    def test_pad_compression_bit_identical(self, posture):
        r1d = simulate_boulder_fall(height_m=3.0, mass_kg=80.0, posture=posture)
        r2d = simulate_boulder_fall_2d(
            height_m=3.0, mass_kg=80.0, posture=posture,
            posture2d=Posture2D(),
        )
        assert np.array_equal(r2d.pad_compression_m, r1d.pad_compression_m)

    def test_scalar_peaks_match_exactly(self):
        """bit-identicality 不止时程，标量也必须完全相等（== 而非 approx）。"""
        for posture in ["controlled-drop", "feet-first-stiff", "head-first"]:
            r1d = simulate_boulder_fall(height_m=3.0, mass_kg=80.0, posture=posture)
            r2d = simulate_boulder_fall_2d(
                height_m=3.0, mass_kg=80.0, posture=posture,
                posture2d=Posture2D(),
            )
            assert r2d.peak_force_n == r1d.peak_force_n
            assert r2d.peak_accel_torso_g == r1d.peak_accel_torso_g
            assert r2d.peak_accel_leg_g == r1d.peak_accel_leg_g
            assert r2d.max_compression_m == r1d.max_compression_m
            assert r2d.energy_into_pad_j == r1d.energy_into_pad_j
            assert r2d.energy_into_flex_j == r1d.energy_into_flex_j
            assert r2d.energy_residual_j == r1d.energy_residual_j
            assert r2d.impact_speed_ms == r1d.impact_speed_ms
            assert r2d.impact_energy_j == r1d.impact_energy_j
            assert r2d.bottomed_out == r1d.bottomed_out

    def test_2d_arrays_are_exactly_zero(self):
        """默认 2D 字段全零 ⇒ ``x_foot`` / ``theta`` / ``omega`` 全 0。"""
        r2d = simulate_boulder_fall_2d(
            height_m=3.0, mass_kg=80.0, posture="controlled-drop",
            posture2d=Posture2D(),
        )
        assert np.array_equal(r2d.x_foot_m, np.zeros_like(r2d.t_s))
        assert np.array_equal(r2d.theta_rad, np.zeros_like(r2d.t_s))
        assert np.array_equal(r2d.omega_rad_s, np.zeros_like(r2d.t_s))
        assert r2d.peak_theta_rad == 0.0
        assert r2d.peak_omega_rad_s == 0.0

    def test_z_torso_starts_at_minus_leg_ext0(self):
        """回归：``z_torso[0] == -leg_ext0`` 恰为 -0.85（与 1D 共守的契约）。"""
        r2d = simulate_boulder_fall_2d(
            height_m=2.0, mass_kg=80.0, posture="controlled-drop",
            posture2d=Posture2D(),
        )
        assert r2d.z_torso_m[0] == -LEG_EXT0_M
        assert r2d.z_torso_m[0] == -0.85
        assert r2d.z_torso_m[0] < r2d.z_foot_m[0], "躯干初始应在脚上方"

    def test_meta_dict_preserved(self):
        """``BoulderFallResult2D.meta`` 必须包含 ``injury/bone`` 所需的所有键。"""
        r2d = simulate_boulder_fall_2d(
            height_m=3.0, mass_kg=80.0, posture="controlled-drop",
            posture2d=Posture2D(),
        )
        for key in ("posture_obj", "m_up", "m_low", "impact_window_s",
                    "head_is_contact", "f_flex_peak_n", "on_pad",
                    "pad_thickness_m", "k_flex", "x_flex_max"):
            assert key in r2d.meta, f"meta 缺少 injury/bone 必需的键 {key!r}"

    def test_meta_isolated_from_r1d(self):
        """2D 结果的 meta 不应与 r1d 共享同一个 dict 对象（避免下游误改）。"""
        r2d = simulate_boulder_fall_2d(
            height_m=3.0, mass_kg=80.0, posture="controlled-drop",
            posture2d=Posture2D(),
        )
        r1d = simulate_boulder_fall(height_m=3.0, mass_kg=80.0, posture="controlled-drop")
        # 修改 2D 的 meta 不应影响 r1d
        r2d.meta["_test_marker"] = 42
        assert "_test_marker" not in r1d.meta


# ==========================================================================
# 不变量 2 —— 摩擦两端行为可解释
# ==========================================================================
class TestInvariant2FrictionLimits:
    """``mu_foot`` 两端极限的行为契约。"""

    def test_zero_mu_no_lateral_motion(self):
        """``mu_foot=0`` ⇒ 横向力恒为 0 ⇒ ``x_f`` 永远保持 0。

        ``f_pad_x = -mu * f_pad * sign_reg(vx_f)`` 在 ``mu=0`` 下恒 0，
        无论 ``vx_f``、``f_pad`` 是多少。这是默认 1D 情形的退化：
        1D 模型里压根没有横向自由度，2D 模型里它必须保持静止。
        """
        r2d = simulate_boulder_fall_2d(
            height_m=3.0, mass_kg=80.0, posture="controlled-drop",
            posture2d=Posture2D(
                mu_foot=0.0,
                I_body_kgm2=50.0,        # 打开转动以确认与横向独立
                theta0_rad=0.1,
                foot_x_lever_m=0.02,    # 同上
            ),
        )
        assert np.array_equal(r2d.x_foot_m, np.zeros_like(r2d.t_s)), \
            "mu=0 时脚不应有横向位移"

    def test_zero_mu_does_not_perturb_theta(self):
        """``mu_foot=0`` ⇒ 横向力严格为 0 ⇒ ``theta`` 与摩擦反馈无关。

        本模型的 ``omega_dot`` 表达式里**没有** ``vx_f`` 项（见
        ``pad2d.py`` 中 ``rhs_2d``），所以横向力无法耦合到转动。
        这里跑两次（一次 mu=0、一次 mu=1 但 foot_x_lever=0）确认
        ``theta_rad`` 完全一致 —— 这是横向/转动**无耦合**的实证。
        """
        common_kwargs = dict(
            height_m=3.0, mass_kg=80.0, posture="controlled-drop",
        )
        r_mu0 = simulate_boulder_fall_2d(
            **common_kwargs,
            posture2d=Posture2D(
                mu_foot=0.0, foot_x_lever_m=0.0,
                I_body_kgm2=50.0, theta0_rad=0.05,
            ),
        )
        r_mu1 = simulate_boulder_fall_2d(
            **common_kwargs,
            posture2d=Posture2D(
                mu_foot=1.0, foot_x_lever_m=0.0,
                I_body_kgm2=50.0, theta0_rad=0.05,
            ),
        )
        assert np.array_equal(r_mu0.theta_rad, r_mu1.theta_rad)

    def test_infinite_mu_foot_sticks(self):
        """``mu_foot=1e9`` ⇒ 脚"粘"在原地（``|vx_foot| ≈ 0``）。

        数学上：在 ``vx_f=0`` 时 ``sign_reg(0)=0`` ⇒ ``f_pad_x=0`` ⇒
        ``vx_f_dot=0`` ⇒ ``vx_f`` 永远为 0 ⇒ 脚不会横向移动。

        这是个**平凡**的"stick"：本 2D 模型里横向 / 转动**无耦合**，
        也没有别的横向驱动（lateral_offset_m=0、theta→x 通道未建模），
        所以"无限摩擦"等同于"无摩擦+无驱动"。spec 显式要求这条断言
        —— 把这个平凡事实**钉在测试里**，防止将来添加横向-转动耦合
        时悄悄打破它。
        """
        r2d = simulate_boulder_fall_2d(
            height_m=3.0, mass_kg=80.0, posture="controlled-drop",
            posture2d=Posture2D(
                mu_foot=1e9,            # 物理上无穷大
                I_body_kgm2=50.0,      # 打开转动以便观察
                theta0_rad=0.05,
                foot_x_lever_m=0.02,    # 横向-转动通道仍关（无 vx_f→omega）
            ),
        )
        assert np.all(np.abs(r2d.x_foot_m) < 1e-12), \
            f"mu=1e9 时脚最大位移 {np.max(np.abs(r2d.x_foot_m)):.2e} m 过大"


# ==========================================================================
# 不变量 3 —— 内翻 vs 外翻相对关系稳定
# ==========================================================================
class TestInvariant3InversionEversion:
    """``theta0 = ±X`` (内翻 vs 外翻)，``foot_x_lever_m ≠ 0`` 时
    转动激励应有**稳定**的符号敏感不对称。

    关键：本测试**不**预设哪一侧更重，只预设"排序在所有高度上一致"。
    """

    @pytest.fixture(scope="class")
    def asymmetry_data(self):
        """跨 3 个高度的 +X / -X peak_omega 比较。

        ``scope="class"`` 一次性计算三个高度，避免每个测试重复跑 6 次
        仿真（~10 s/次）。pytest 8+ 推荐用 ``@classmethod`` 替代
        ``scope="class"`` 的实例方法，但本测试只读结果不修改状态，
        实例方法已足。
        """
        return _asymmetry_data()

    def test_sign_sensitive_asymmetry_exists(self, asymmetry_data):
        """每一高度上 ``peak_omega(+X)`` 和 ``peak_omega(-X)`` 必须不同。

        若两者相同 ⇒ ``foot_x_lever_m`` 没有引入符号敏感的不对称，
        模型退化成关于 theta 的对称方程。
        """
        for r in asymmetry_data:
            assert r["peak_omega_pos"] != r["peak_omega_neg"], \
                f"h={r['height_m']}: 内翻/外翻完全对称 —— foot_x_lever_m 无效"

    def test_asymmetry_direction_stable_across_heights(self, asymmetry_data):
        """3 个高度上 ``peak_omega(+X)`` vs ``peak_omega(-X)`` 的
        **排序方向**一致。

        spec 明确禁止断言物理方向（"内翻更轻"还是"外翻更轻"），
        只断言稳定。若某个高度上 +X 更大、另一高度上 -X 更大，
        说明模型有数值噪声 / 不稳定性，需要排查。
        """
        diffs = [r["peak_omega_pos"] - r["peak_omega_neg"] for r in asymmetry_data]
        all_pos = all(d > 0 for d in diffs)
        all_neg = all(d < 0 for d in diffs)
        assert all_pos or all_neg, (
            f"排序方向不一致: heights={[r['height_m'] for r in asymmetry_data]}, "
            f"diffs={diffs} —— 内翻/外翻不对称不稳定"
        )

    def test_zero_lever_arm_no_asymmetry(self):
        """``foot_x_lever_m=0`` 时内翻/外翻必须完全对称（这是健全性检查）。"""
        common_2d = dict(
            foot_x_lever_m=0.0,        # 关闭符号敏感项
            I_body_kgm2=50.0,
            mu_foot=0.5,
        )
        for h in [1.0, 2.0, 3.0]:
            r_pos = simulate_boulder_fall_2d(
                height_m=h, mass_kg=80.0, posture="controlled-drop",
                posture2d=Posture2D(theta0_rad=+0.1, **common_2d),
            )
            r_neg = simulate_boulder_fall_2d(
                height_m=h, mass_kg=80.0, posture="controlled-drop",
                posture2d=Posture2D(theta0_rad=-0.1, **common_2d),
            )
            assert r_pos.peak_omega_rad_s == r_neg.peak_omega_rad_s, \
                f"h={h}: foot_x_lever_m=0 时仍存在符号差异"

    def test_sign_sensitive_asymmetry_exists(self, asymmetry_data):
        """每一高度上 ``peak_omega(+X)`` 和 ``peak_omega(-X)`` 必须不同。

        若两者相同 ⇒ ``foot_x_lever_m`` 没有引入符号敏感的不对称，
        模型退化成关于 theta 的对称方程。
        """
        for r in asymmetry_data:
            assert r["peak_omega_pos"] != r["peak_omega_neg"], \
                f"h={r['height_m']}: 内翻/外翻完全对称 —— foot_x_lever_m 无效"

    def test_asymmetry_direction_stable_across_heights(self, asymmetry_data):
        """3 个高度上 ``peak_omega(+X)`` vs ``peak_omega(-X)`` 的
        **排序方向**一致。

        spec 明确禁止断言物理方向（"内翻更轻"还是"外翻更轻"），
        只断言稳定。若某个高度上 +X 更大、另一高度上 -X 更大，
        说明模型有数值噪声 / 不稳定性，需要排查。
        """
        diffs = [r["peak_omega_pos"] - r["peak_omega_neg"] for r in asymmetry_data]
        all_pos = all(d > 0 for d in diffs)
        all_neg = all(d < 0 for d in diffs)
        assert all_pos or all_neg, (
            f"排序方向不一致: heights={[r['height_m'] for r in asymmetry_data]}, "
            f"diffs={diffs} —— 内翻/外翻不对称不稳定"
        )

    def test_zero_lever_arm_no_asymmetry(self):
        """``foot_x_lever_m=0`` 时内翻/外翻必须完全对称（这是健全性检查）。"""
        common_2d = dict(
            foot_x_lever_m=0.0,        # 关闭符号敏感项
            I_body_kgm2=50.0,
            mu_foot=0.5,
        )
        for h in [1.0, 2.0, 3.0]:
            r_pos = simulate_boulder_fall_2d(
                height_m=h, mass_kg=80.0, posture="controlled-drop",
                posture2d=Posture2D(theta0_rad=+0.1, **common_2d),
            )
            r_neg = simulate_boulder_fall_2d(
                height_m=h, mass_kg=80.0, posture="controlled-drop",
                posture2d=Posture2D(theta0_rad=-0.1, **common_2d),
            )
            assert r_pos.peak_omega_rad_s == r_neg.peak_omega_rad_s, \
                f"h={h}: foot_x_lever_m=0 时仍存在符号差异"


# ==========================================================================
# 不变量 3 用的 fixture helper —— 顶层函数（pytest 8+ 推荐）
# ==========================================================================
def _asymmetry_data():
    """跨 3 个高度的 +X / -X peak_omega 比较（模块级函数，避开
    pytest 8 对 ``scope="class"`` 实例方法 fixture 的 deprecation 警告）。"""
    theta0 = 0.1                     # ~5.7°
    common_2d = dict(
        foot_x_lever_m=0.02,         # 触发符号敏感不对称
        I_body_kgm2=50.0,
        mu_foot=0.5,
    )
    heights = [1.0, 2.0, 3.0]
    rows = []
    for h in heights:
        r_pos = simulate_boulder_fall_2d(
            height_m=h, mass_kg=80.0, posture="controlled-drop",
            posture2d=Posture2D(theta0_rad=+theta0, **common_2d),
        )
        r_neg = simulate_boulder_fall_2d(
            height_m=h, mass_kg=80.0, posture="controlled-drop",
            posture2d=Posture2D(theta0_rad=-theta0, **common_2d),
        )
        rows.append(dict(
            height_m=h,
            peak_omega_pos=r_pos.peak_omega_rad_s,
            peak_omega_neg=r_neg.peak_omega_rad_s,
            peak_theta_pos=r_pos.peak_theta_rad,
            peak_theta_neg=r_neg.peak_theta_rad,
        ))
    return rows
class TestInvariant4PeakTimeMonotone:
    """[T22] Table 3：峰值力时刻**距接触起点**的延迟随高度单调非增。

    说明：spec 字面写的是 ``argmax(pad_force_n) time``，但绝对时刻
    包含 ``sqrt(2h/g)`` 自由落体段，它随高度**单调递增**（更高 ⇒ 更长
    自由落体），只有减去这段后的"接触-峰延迟"才单调非增。下文
    ``contact_to_peak_time`` 用 ``meta['impact_window_s']`` 把反弹重峰
    也排除掉（直接 ``argmax`` 在 h ≥ 1 m 时会被反弹重峰带跑）。
    """

    @pytest.mark.parametrize("heights", [
        [0.5, 1.0, 1.5, 2.0],
        [1.0, 2.0, 3.0, 5.0],
        [2.0, 3.0, 4.0, 5.0],
    ])
    def test_contact_to_peak_time_monotone_non_increasing(self, heights):
        """接触-峰延迟 dt_after_fall 在 heights 序列上单调非增。"""
        dts = []
        for h in heights:
            r = simulate_boulder_fall_2d(
                height_m=h, mass_kg=80.0, posture="controlled-drop",
                posture2d=Posture2D(),
            )
            dts.append(contact_to_peak_time(r))
        for i in range(1, len(heights)):
            assert dts[i] <= dts[i - 1] + 1e-9, (
                f"h={heights[i]} m 的接触-峰延迟 {dts[i]*1e3:.2f} ms "
                f"晚于 h={heights[i-1]} m 的 {dts[i-1]*1e3:.2f} ms"
            )

    def test_2d_peak_time_matches_1d_at_default(self):
        """``theta0=0``（默认） ⇒ 2D 峰值时刻 ≡ 1D 峰值时刻（bit-identicality 的直接推论）。"""
        for h in [0.5, 1.5, 3.0, 5.0, 8.0]:
            r1d = simulate_boulder_fall(height_m=h, mass_kg=80.0, posture="controlled-drop")
            r2d = simulate_boulder_fall_2d(
                height_m=h, mass_kg=80.0, posture="controlled-drop",
                posture2d=Posture2D(),
            )
            t1d = impact_window_peak_time(r1d)
            t2d = impact_window_peak_time(r2d)
            assert t1d == t2d, f"h={h}: 1D 峰值时刻 {t1d} ≠ 2D 峰值时刻 {t2d}"


# ==========================================================================
# 消费者兼容性：BoulderFallResult2D 必须能被 injury/bone/metrics 直接消费
# ==========================================================================
class TestConsumersCompatibility:
    """BoulderFallResult2D 是 BoulderFallResult 子类，下游消费者不应感知差异。"""

    def test_injury_assess_works(self):
        from climbing.injury import assess_boulder_fall
        r = simulate_boulder_fall_2d(
            height_m=3.0, mass_kg=80.0, posture="feet-first-stiff",
            posture2d=Posture2D(),
        )
        rep = assess_boulder_fall(r, "feet-first-stiff")
        assert rep.verdict in {"低", "中", "高", "极高"}
        assert rep.peak_pressure_kpa > 0.0
        assert rep.hic >= 0.0
        # injury 报的峰值力应等于 result 的峰值力（位单位换算后）
        assert rep.peak_force_kn == r.peak_force_n / 1e3

    def test_bone_assess_works(self):
        """``bone.assess_sites`` 在 ``meta['impact_window_s']`` 缺失时抛错
        —— BoulderFallResult2D 必须继承这个字段。
        """
        from climbing.bone import assess_sites
        r = simulate_boulder_fall_2d(
            height_m=3.0, mass_kg=77.0, on_pad=False, posture="feet-first-stiff",
            posture2d=Posture2D(),
        )
        v = assess_sites(r)
        assert v.worst is not None
        # 必须不抛错（冲击窗口被正确继承）
        assert v.peak_pad_n > 0.0

    def test_metrics_hic_works(self):
        from climbing.metrics import hic_from_result
        r = simulate_boulder_fall_2d(
            height_m=3.0, mass_kg=80.0, posture="head-first",
            posture2d=Posture2D(),
        )
        h = hic_from_result(r)
        assert h >= 0.0
        # head-first 时 primary_accel_g 应取 leg（meta['head_is_contact'] = True）
        assert r.meta["head_is_contact"] is True

    def test_inheritance_isinstance(self):
        """isinstance(r, BoulderFallResult) 必须为 True —— 消费者期望这点。"""
        r = simulate_boulder_fall_2d(
            height_m=3.0, mass_kg=80.0, posture="controlled-drop",
            posture2d=Posture2D(),
        )
        assert isinstance(r, BoulderFallResult2D)
        assert isinstance(r, BoulderFallResult)

    def test_1d_subclass_result_still_compatible(self):
        """1D BoulderFallResult 的现有消费者调用方式不被打破。"""
        r = simulate_boulder_fall_2d(
            height_m=2.0, mass_kg=80.0, posture="controlled-drop",
            posture2d=Posture2D(),
        )
        # 1D 属性都还在
        for attr in ("z_foot_m", "z_torso_m", "pad_force_n", "pad_compression_m",
                     "peak_force_n", "peak_accel_torso_g", "peak_accel_leg_g",
                     "max_compression_m", "bottomed_out", "meta"):
            assert hasattr(r, attr), f"1D 字段 {attr!r} 不见了"


# ==========================================================================
# 接口契约
# ==========================================================================
class TestSignature:
    def test_default_height(self):
        r = simulate_boulder_fall_2d()
        assert r.height_m == 3.0

    def test_default_mass(self):
        r = simulate_boulder_fall_2d()
        assert r.mass_kg == 80.0

    def test_default_posture2d_is_none(self):
        """``posture2d=None`` 是公开 API 的默认值，调用方不必显式构造 Posture2D()。"""
        sig = inspect.signature(simulate_boulder_fall_2d)
        assert sig.parameters["posture2d"].default is None

    def test_signature_accepts_posture2d_kwarg(self):
        sig = inspect.signature(simulate_boulder_fall_2d)
        assert "posture2d" in sig.parameters

    def test_signature_mirrors_pad(self):
        """1D 与 2D 入口的前 11 个参数名必须一致 —— 调用方可以在
        知道全部 1D 字段后**只加一个** ``posture2d=`` 即可升级到 2D。"""
        sig1d = inspect.signature(simulate_boulder_fall)
        sig2d = inspect.signature(simulate_boulder_fall_2d)
        common = ("height_m", "mass_kg", "posture", "pad", "on_pad", "g",
                  "dt_sample", "t_max", "max_ode_step", "rtol", "atol")
        for name in common:
            assert name in sig1d.parameters, f"1D 签名缺 {name!r}"
            assert name in sig2d.parameters, f"2D 签名缺 {name!r}"
            # 默认值也必须一致（"镜像" = 1D 默认值在 2D 仍然有效）
            assert sig1d.parameters[name].default == sig2d.parameters[name].default, \
                f"{name} 默认值在 1D/2D 间漂移"


class TestResultShape:
    def test_2d_arrays_match_t_s_shape(self):
        r = simulate_boulder_fall_2d(
            height_m=3.0, mass_kg=80.0, posture="controlled-drop",
            posture2d=Posture2D(),
        )
        assert r.x_foot_m.shape == r.t_s.shape
        assert r.theta_rad.shape == r.t_s.shape
        assert r.omega_rad_s.shape == r.t_s.shape

    def test_2d_with_nonzero_params_match_shape(self):
        """非默认 2D 参数下，shape 仍须与 t_s 对齐。"""
        r = simulate_boulder_fall_2d(
            height_m=3.0, mass_kg=80.0, posture="controlled-drop",
            posture2d=Posture2D(
                theta0_rad=0.1, foot_x_lever_m=0.02,
                I_body_kgm2=50.0, mu_foot=0.5,
            ),
        )
        assert r.x_foot_m.shape == r.t_s.shape
        assert r.theta_rad.shape == r.t_s.shape
        assert r.omega_rad_s.shape == r.t_s.shape

    def test_result_is_dataclass(self):
        from dataclasses import is_dataclass
        r = simulate_boulder_fall_2d(
            height_m=3.0, mass_kg=80.0, posture="controlled-drop",
            posture2d=Posture2D(),
        )
        assert is_dataclass(r)


class TestNonzero2DParametersActuallyAnimate:
    """健全性：非默认 2D 参数下，``x_foot`` / ``theta`` / ``omega``
    必须真的有非零演化 —— 不然"2D"是空壳。"""

    def test_theta_oscillates_under_gravity(self):
        """``I_body > 0, theta0 != 0`` ⇒ theta 应在重力矩驱动下**继续演化**。

        注意："演化"包括两部分 —— 接触段受 ``f_pad * d_lever`` 驱动，
        脱离接触后受 ``-m*g*L_arm*sin(theta)`` 摆动。``peak_theta_rad``
        应显著偏离初值。
        """
        r = simulate_boulder_fall_2d(
            height_m=3.0, mass_kg=80.0, posture="controlled-drop",
            posture2d=Posture2D(
                theta0_rad=0.3,        # 大幅初始倾角
                I_body_kgm2=50.0,
                mu_foot=0.5,
                foot_x_lever_m=0.02,
            ),
        )
        assert r.peak_theta_rad > 0.3 + 1e-3, \
            f"theta 没有演化（max={r.peak_theta_rad:.4f}, 初值 0.3）"

    def test_mu_zero_no_lateral_even_with_nonzero_other(self):
        """``mu_foot=0`` ⇒ ``x_f`` 永远保持初值 ``lateral_offset_m``（若有）。"""
        r = simulate_boulder_fall_2d(
            height_m=3.0, mass_kg=80.0, posture="controlled-drop",
            posture2d=Posture2D(
                mu_foot=0.0,
                lateral_offset_m=0.05,
                I_body_kgm2=50.0,      # 其它都打开
                theta0_rad=0.05, foot_x_lever_m=0.02,
            ),
        )
        # lateral_offset_m 给定 x_f(0)=0.05，mu=0 ⇒ 永远保持
        assert np.array_equal(r.x_foot_m, np.full_like(r.t_s, 0.05))


# ==========================================================================
# P2a/P2b —— 踝姿势角 + 滑动/卡住 + 韧带 / 骨失效
# ==========================================================================
# 复现 [T22]（Li Z et al., Forensic Sci Res 2022;7(3):518-527）。全部为
# **后处理派生量**，默认（β0=0 且 beta_from_theta=0）逐位归零、不改 1D 轨迹。
#
# 标定常量（见 ``src/climbing/pad2d.py``，另一 agent 的
# ``scripts/phase_p1_montecarlo.py`` 与之一致）：
#   mu_slide      = 0.4    —— |sinβ|>0.4（≈β>23.6°）⇒ 滑动
#   d_lateral     = 0.030 m （内翻，β≥0 的力矩臂）
#   d_medial      = 0.045 m （外翻，β<0）—— > lateral ⇒ V9「内翻轻于外翻」
#   ligament_r_mm = 22.0   （ATFL 力臂，[R6] 默认）
_ANKLE_D = dict(
    mu_slide=0.4,
    ankle_d_lateral_m=0.030,
    ankle_d_medial_m=0.045,
    ligament_r_mm=22.0,
)


def _ankle_run(height_m, beta0_deg, *, sign=1.0, beta_from_theta=0.0, **over):
    """跑一次带 P2 踝参数的 2D 仿真（default pad / controlled-drop）。"""
    params = dict(_ANKLE_D)
    params.update(over)
    return simulate_boulder_fall_2d(
        height_m=height_m, mass_kg=80.0, posture="controlled-drop",
        posture2d=Posture2D(
            ankle_beta0_rad=float(np.deg2rad(beta0_deg)) * float(sign),
            beta_from_theta=beta_from_theta,
            **params,
        ),
    )


_V6_CACHE: dict = {}


def _v6_grid():
    """[T22] Table 4 的 (h, β) 3×3 网格；模块内缓存避免重复 9 次求解。"""
    if not _V6_CACHE:
        for h in (5.0, 10.0, 20.0):
            for a in (10.0, 20.0, 30.0):
                _V6_CACHE[(h, a)] = _ankle_run(h, a)
    return _V6_CACHE


# ==========================================================================
# P2a —— 踝姿势角 β(t) = β0 + beta_from_theta·θ(t)
# ==========================================================================
class TestP2aAnklePose:
    """P2a：β 的精确映射 + β0 符号约定（>0 内翻 / <0 外翻）。"""

    def test_beta_follows_theta_exactly(self):
        """``β(t) = β0 + beta_from_theta·θ(t)`` 在采样栅格上逐位成立。"""
        p2d = Posture2D(
            theta0_rad=0.1, I_body_kgm2=50.0, mu_foot=0.5, foot_x_lever_m=0.02,
            ankle_beta0_rad=0.2, beta_from_theta=1.5,
            mu_slide=0.4, ankle_d_lateral_m=0.03, ankle_d_medial_m=0.045,
        )
        r = simulate_boulder_fall_2d(
            height_m=3.0, mass_kg=80.0, posture="controlled-drop", posture2d=p2d,
        )
        # θ(t) 必须真的非零，否则这条断言是空的
        assert not np.array_equal(r.theta_rad, np.zeros_like(r.t_s))
        expected = p2d.ankle_beta0_rad + p2d.beta_from_theta * r.theta_rad
        assert np.array_equal(r.ankle_beta_rad, expected), \
            "β(t) 必须逐位等于 β0 + beta_from_theta·θ(t)"

    def test_peak_beta_is_max_abs(self):
        r = _ankle_run(3.0, 20.0)
        assert r.peak_ankle_beta_rad == float(np.max(np.abs(r.ankle_beta_rad)))

    def test_beta_sign_convention(self):
        """β0>0 = 内翻、β0<0 = 外翻 ⇒ β(t) 时程符号跟随 β0。"""
        r_inv = _ankle_run(3.0, 20.0, sign=+1.0)
        r_ev = _ankle_run(3.0, 20.0, sign=-1.0)
        assert np.all(r_inv.ankle_beta_rad > 0.0)
        assert np.all(r_ev.ankle_beta_rad < 0.0)

    def test_sign_sensitive_asymmetry_stable(self):
        """同 |β|=20°，符号敏感：跨 3 高度 ``peak_inversion_moment_nmm``
        的排序一致（= 不变量 ③ 在踝层的体现）。"""
        heights = (2.0, 3.0, 5.0)
        diffs = []
        for h in heights:
            rp = _ankle_run(h, 20.0, sign=+1.0)
            rn = _ankle_run(h, 20.0, sign=-1.0)
            assert rp.peak_inversion_moment_nmm != rn.peak_inversion_moment_nmm, \
                f"h={h}: 内翻/外翻力矩完全对称 —— 符号不敏感"
            diffs.append(rn.peak_inversion_moment_nmm - rp.peak_inversion_moment_nmm)
        # d_medial > d_lateral ⇒ 外翻(−) 力矩更大；排序在 3 高度上一致。
        assert all(d > 0.0 for d in diffs), f"不对称方向不稳定: {diffs}"


# ==========================================================================
# P2b —— 滑动 vs 卡住
# ==========================================================================
class TestP2bSliding:
    """``sliding(t) = |sinβ(t)| > mu_slide``（Coulomb 锥）两端与边界。"""

    def test_mu_slide_zero_always_sliding(self):
        r = _ankle_run(3.0, 20.0, mu_slide=0.0)
        assert r.ankle_mode == "sliding"
        assert r.ankle_slide_frac == 1.0

    def test_mu_slide_huge_always_stuck(self):
        r = _ankle_run(3.0, 20.0, mu_slide=1e9)
        assert r.ankle_mode == "stuck"
        assert r.ankle_slide_frac == 0.0

    def test_boundary_is_sin_beta_equals_mu_slide(self):
        """边界恰为 ``|sinβ| = mu_slide``：β=30° ⇒ |sinβ|≈0.5。

        ``mu_slide=0.5`` ⇒ 不滑动（严格 ``>``）；略小 ⇒ 滑动。
        """
        r_stuck = _ankle_run(3.0, 30.0, mu_slide=0.5)
        r_slide = _ankle_run(3.0, 30.0, mu_slide=0.4999)
        assert r_stuck.ankle_mode == "stuck"
        assert r_slide.ankle_mode == "sliding"
        assert r_stuck.ankle_slide_frac == 0.0
        assert r_slide.ankle_slide_frac == 1.0

    def test_sliding_lowers_bone_utilization(self):
        """同 (h, β)：滑动卸载 ⇒ 骨利用率 < 卡住时的全额。"""
        h, beta = 5.0, 30.0
        r_stuck = _ankle_run(h, beta, mu_slide=1e9)
        r_slide = _ankle_run(h, beta, mu_slide=0.0)
        assert r_stuck.ankle_bone_utilization > r_slide.ankle_bone_utilization


# ==========================================================================
# V6 —— [T22] Table 4 的 (高度 × 角度) 定性结构
# ==========================================================================
class TestP2bV6Table4:
    """扫描 heights=(5,10,20) × angles=(10,20,30)°，复现 [T22] Table 4：

    ① 低高度 + 大角度 → 无骨折（滑动）
    ② 低高度 + 小角度 → 骨折
    ③ 高高度 → 无论角度都骨折
    ④ 存在"角度—高度"转变曲线
    """

    def test_table4_qualitative_structure(self):
        grid = _v6_grid()
        expect = {
            (5.0, 10.0): True,    # 5 m/10° 骨折（[T22] 腓骨+胫骨中下段）
            (5.0, 20.0): True,    # 5 m/20° 骨折
            (5.0, 30.0): False,   # 5 m/30° 滑动，力被卸掉 ⇒ 无骨折
            (10.0, 10.0): True,
            (10.0, 20.0): True,
            (10.0, 30.0): True,
            (20.0, 10.0): True,
            (20.0, 20.0): True,
            (20.0, 30.0): True,
        }
        for key, exp in expect.items():
            got = bool(grid[key].ankle_fracture)
            assert got == exp, (
                f"h={key[0]}m β={key[1]}°: 期望骨折={exp} 得到={got} "
                f"(util={grid[key].ankle_bone_utilization:.3f}, "
                f"mode={grid[key].ankle_mode})"
            )

    def test_large_angle_at_5m_is_sliding_no_fracture(self):
        r = _v6_grid()[(5.0, 30.0)]
        assert r.ankle_mode == "sliding"
        assert not r.ankle_fracture

    def test_transition_curve_exists(self):
        """至少一个相邻 (h,β) 骨折状态翻转 ⇒ 存在角度—高度转变曲线。"""
        grid = _v6_grid()
        flips = 0
        for h in (5.0, 10.0, 20.0):
            for a0, a1 in ((10.0, 20.0), (20.0, 30.0)):
                if grid[(h, a0)].ankle_fracture != grid[(h, a1)].ankle_fracture:
                    flips += 1
        for a in (10.0, 20.0, 30.0):
            for h0, h1 in ((5.0, 10.0), (10.0, 20.0)):
                if grid[(h0, a)].ankle_fracture != grid[(h1, a)].ankle_fracture:
                    flips += 1
        assert flips >= 1, "网格中不存在骨折翻转 ⇒ 无转变曲线"


# ==========================================================================
# V9 —— [T22] 内翻损伤轻于外翻（d_medial > d_lateral）
# ==========================================================================
class TestP2bV9:
    """V9：``ankle_d_medial_m > ankle_d_lateral_m`` 编码「内翻轻于外翻」。"""

    def test_eversion_worse_than_inversion(self):
        """同 |β|=20°、跨 3 高度：外翻(−) 的应变与骨利用率都 > 内翻(+)。"""
        for h in (2.0, 3.0, 5.0):
            r_inv = _ankle_run(h, 20.0, sign=+1.0)
            r_ev = _ankle_run(h, 20.0, sign=-1.0)
            assert r_ev.peak_ligament_strain > r_inv.peak_ligament_strain, \
                f"h={h}: 外翻应变应 > 内翻"
            assert r_ev.ankle_bone_utilization > r_inv.ankle_bone_utilization, \
                f"h={h}: 外翻骨利用率应 > 内翻"

    def test_ordering_stable_across_angles_and_heights(self):
        for a in (10.0, 30.0):
            for h in (2.0, 3.0, 5.0):
                r_inv = _ankle_run(h, a, sign=+1.0)
                r_ev = _ankle_run(h, a, sign=-1.0)
                assert r_ev.peak_ligament_strain > r_inv.peak_ligament_strain

    def test_equal_arms_no_asymmetry(self):
        """两侧力臂相等 ⇒ 无不对称（健全性检查，证明不对称源自力臂差）。"""
        over = dict(ankle_d_lateral_m=0.03, ankle_d_medial_m=0.03)
        rp = _ankle_run(3.0, 20.0, sign=+1.0, **over)
        rn = _ankle_run(3.0, 20.0, sign=-1.0, **over)
        assert rp.peak_ligament_strain == rn.peak_ligament_strain
        assert rp.ankle_bone_utilization == rn.ankle_bone_utilization


# ==========================================================================
# V6/[B25] —— 垫刚度 ⇒ 韧带应变非单调
# ==========================================================================
class TestP2bPadStiffnessNonMonotonic:
    """[B25]：垫刚度 ↑ ⇒ 韧带应变**非单调**（存在内部极大）。

    机制：``ε ∝ peak(f_pad(t)·|sinβ(t)|)``。刚度↑ 一方面抬高峰值力，另一方面
    使峰值力**提前**（不变量 ④）并改变它与踝旋转 β(t) 振荡峰的相位，于是
    **乘积**在中等刚度处出现内极大 —— 软垫（力小）与硬垫（相位错开）两侧更低。
    这正是 [B25]「提高刚度可能预防过度旋后，但高能场景可能无效甚至更糟」。
    """

    @staticmethod
    def _pads():
        # 平台应力 / 压实应力按同一因子 f 缩放 ⇒ 单调整体刚度族。
        fs = (10.0, 20.0, 50.0, 100.0)
        pads = [
            CrashPad(s_plateau_pa=35000.0 * f,
                     s_dense_pa=365000.0 * float(np.sqrt(f)))
            for f in fs
        ]
        return fs, pads

    def test_non_monotonic_local_max_exists(self):
        fs, pads = self._pads()
        p2d = Posture2D(
            theta0_rad=0.30, I_body_kgm2=50.0, foot_x_lever_m=0.02, mu_foot=0.5,
            ankle_beta0_rad=0.20, beta_from_theta=1.0, mu_slide=0.4,
            ankle_d_lateral_m=0.030, ankle_d_medial_m=0.045, ligament_r_mm=22.0,
        )
        eps = []
        for pad in pads:
            r = simulate_boulder_fall_2d(
                height_m=3.0, mass_kg=80.0, posture="controlled-drop",
                pad=pad, posture2d=p2d,
            )
            eps.append(r.peak_ligament_strain)
        # 存在 k1<k2<k3 使 ε(k2)>ε(k1) 且 ε(k2)>ε(k3)（局部极大 ⇒ 非单调）
        assert any(
            eps[j] > eps[j - 1] and eps[j] > eps[j + 1]
            for j in range(1, len(fs) - 1)
        ), f"ε 未出现非单调局部极大: f={fs}, ε={eps}"


# ==========================================================================
# 载荷分配系数（load-partition factor）—— ``ankle_load_share``
# ==========================================================================
# 契约：
#   1. ``ankle_load_share = 1.0``（默认）⇒ 与今日行为**逐位一致**（regression
#      锁），包括 bit-identicality 不变量（TestInvariant1BitIdentical）。
#   2. ``ankle_load_share = 0.0`` ⇒ 踝所有派生量精确为 0 / False / "none"。
#   3. ``peak_ligament_strain`` 与 ``ankle_bone_utilization`` 在
#      ``ankle_load_share`` 单调变化时也**单调变化**（线性缩放）。
#
# 【不要混淆】 ``ankle_load_share`` 不影响任何 1D 轨迹
# （不在 ``rhs_2d``），只缩放 ``_derive_ankle_state`` 中的 ``f_pad``。
class TestLoadPartition:
    """载荷分配系数 ``ankle_load_share`` 的契约测试。"""

    def test_default_share_equals_one(self):
        """默认 ``ankle_load_share == 1.0`` —— regression 锁。"""
        p = Posture2D()
        assert p.ankle_load_share == 1.0

    @pytest.mark.parametrize("posture", [
        "controlled-drop", "feet-first-stiff", "butt-impact",
        "flat-flop", "tuck-roll", "head-first",
    ])
    def test_share_one_reproduces_old_numbers_exactly(self, posture):
        """``ankle_load_share = 1.0`` ⇒ 与 ``Posture2D()``（不显式设
        ankle_load_share）等价：peak_ligament_strain / bone_util / 状态全部
        ``np.array_equal``（不是 ``pytest.approx``）。

        这是 regression 锁：lock 当前数字确保未来默认不变。
        """
        common = dict(
            ankle_beta0_rad=np.deg2rad(20.0),
            beta_from_theta=0.0,
            mu_slide=0.4, ankle_d_lateral_m=0.030, ankle_d_medial_m=0.045,
            ligament_r_mm=22.0,
        )
        r_default = simulate_boulder_fall_2d(
            height_m=3.0, mass_kg=80.0, posture=posture,
            posture2d=Posture2D(**common),
        )
        r_share1 = simulate_boulder_fall_2d(
            height_m=3.0, mass_kg=80.0, posture=posture,
            posture2d=Posture2D(ankle_load_share=1.0, **common),
        )
        assert np.array_equal(r_default.peak_ligament_strain,
                              r_share1.peak_ligament_strain)
        assert np.array_equal(r_default.ankle_bone_utilization,
                              r_share1.ankle_bone_utilization)
        assert r_default.peak_inversion_moment_nmm == r_share1.peak_inversion_moment_nmm
        assert r_default.ankle_sprain == r_share1.ankle_sprain
        assert r_default.ankle_fracture == r_share1.ankle_fracture

    def test_share_zero_zeroes_all_ankle_outputs(self):
        """``ankle_load_share = 0.0`` ⇒ 踝所有派生量 = 0 / False / 'none'。

        与 Posture2D() 默认（β=0）的「全 0」输出**结构上等价**——
        这是 ``f_pad_ankle = share * f_pad = 0`` ⇒ M = F_lat = F_lig = 0
        ⇒ 后续所有派生量 0。
        """
        r = simulate_boulder_fall_2d(
            height_m=3.0, mass_kg=80.0, posture="controlled-drop",
            posture2d=Posture2D(
                ankle_beta0_rad=np.deg2rad(20.0),    # 非零 β
                mu_slide=0.4, ankle_d_lateral_m=0.030,
                ankle_d_medial_m=0.045, ligament_r_mm=22.0,
                ankle_load_share=0.0,
            ),
        )
        assert np.array_equal(r.ankle_inversion_moment_nmm,
                              np.zeros_like(r.t_s))
        assert r.peak_inversion_moment_nmm == 0.0
        assert r.peak_ligament_strain == 0.0
        assert r.ankle_bone_utilization == 0.0
        assert r.ankle_fracture is False
        assert r.ankle_sprain is False

    @pytest.mark.parametrize("share", [1.0, 0.8, 0.6, 0.477, 0.4, 0.2])
    def test_peak_ligament_strain_monotonic_in_share(self, share):
        """``peak_ligament_strain`` 随 ``ankle_load_share`` 单调非增。

        ``F_lig = |M|/r ∝ share·f_pad`` ⇒ ``ε = F_lig/(k·L)`` 也线性。
        单调性来源：``f_pad`` 本身 ≥ 0（数值上）且 share ∈ [0, 1] ⇒ 乘积单调。
        """
        shares = [1.0, 0.8, 0.6, 0.477, 0.4, 0.2]
        strains = []
        for s in shares:
            r = simulate_boulder_fall_2d(
                height_m=3.0, mass_kg=80.0, posture="controlled-drop",
                posture2d=Posture2D(
                    ankle_beta0_rad=np.deg2rad(20.0),
                    beta_from_theta=0.0,
                    mu_slide=0.4, ankle_d_lateral_m=0.030,
                    ankle_d_medial_m=0.045, ligament_r_mm=22.0,
                    ankle_load_share=s,
                ),
            )
            strains.append(r.peak_ligament_strain)
        for i in range(1, len(strains)):
            assert strains[i] <= strains[i - 1] + 1e-12, (
                f"peak_ligament_strain 不单调: shares={shares}, strains={strains}"
            )
        # share=1.0 必须严格 > share=0.0（只要 β ≠ 0）
        r1 = simulate_boulder_fall_2d(
            height_m=3.0, mass_kg=80.0, posture="controlled-drop",
            posture2d=Posture2D(
                ankle_beta0_rad=np.deg2rad(20.0),
                mu_slide=0.4, ankle_d_lateral_m=0.030,
                ankle_d_medial_m=0.045, ligament_r_mm=22.0,
                ankle_load_share=1.0,
            ),
        )
        r0 = simulate_boulder_fall_2d(
            height_m=3.0, mass_kg=80.0, posture="controlled-drop",
            posture2d=Posture2D(
                ankle_beta0_rad=np.deg2rad(20.0),
                mu_slide=0.4, ankle_d_lateral_m=0.030,
                ankle_d_medial_m=0.045, ligament_r_mm=22.0,
                ankle_load_share=0.0,
            ),
        )
        assert r1.peak_ligament_strain > r0.peak_ligament_strain

    @pytest.mark.parametrize("share", [1.0, 0.8, 0.6, 0.477, 0.4, 0.2])
    def test_ankle_bone_utilization_monotonic_in_share(self, share):
        """``ankle_bone_utilization`` 随 ``ankle_load_share`` 单调非增。

        ``bone_util = max(σ_bend + τ_shear) / σ_c`` 中 ``σ_bend ∝ |M|·c/I ∝ share``，
        ``τ_shear ∝ |F_lat|/A ∝ share``。线性 ⇒ 单调。
        """
        shares = [1.0, 0.8, 0.6, 0.477, 0.4, 0.2]
        utils = []
        for s in shares:
            r = simulate_boulder_fall_2d(
                height_m=3.0, mass_kg=80.0, posture="controlled-drop",
                posture2d=Posture2D(
                    ankle_beta0_rad=np.deg2rad(20.0),
                    beta_from_theta=0.0,
                    mu_slide=0.4, ankle_d_lateral_m=0.030,
                    ankle_d_medial_m=0.045, ligament_r_mm=22.0,
                    ankle_load_share=s,
                ),
            )
            utils.append(r.ankle_bone_utilization)
        for i in range(1, len(utils)):
            assert utils[i] <= utils[i - 1] + 1e-12, (
                f"ankle_bone_utilization 不单调: shares={shares}, utils={utils}"
            )
        # 边界：share=1.0 ≫ share=0.0（只要 β ≠ 0）
        r1 = simulate_boulder_fall_2d(
            height_m=3.0, mass_kg=80.0, posture="controlled-drop",
            posture2d=Posture2D(
                ankle_beta0_rad=np.deg2rad(20.0),
                mu_slide=0.4, ankle_d_lateral_m=0.030,
                ankle_d_medial_m=0.045, ligament_r_mm=22.0,
                ankle_load_share=1.0,
            ),
        )
        r0 = simulate_boulder_fall_2d(
            height_m=3.0, mass_kg=80.0, posture="controlled-drop",
            posture2d=Posture2D(
                ankle_beta0_rad=np.deg2rad(20.0),
                mu_slide=0.4, ankle_d_lateral_m=0.030,
                ankle_d_medial_m=0.045, ligament_r_mm=22.0,
                ankle_load_share=0.0,
            ),
        )
        assert r1.ankle_bone_utilization > r0.ankle_bone_utilization

    def test_share_does_not_affect_1d_trajectory(self):
        """``ankle_load_share`` 不影响 1D 轨迹 —— bit-identicality 锁。

        在任何 share 下，``t_s`` / ``z_foot_m`` / ``z_torso_m`` /
        ``pad_force_n`` 等**逐位**一致（``np.array_equal``）。
        """
        common = dict(
            ankle_beta0_rad=np.deg2rad(20.0),
            mu_slide=0.4, ankle_d_lateral_m=0.030,
            ankle_d_medial_m=0.045, ligament_r_mm=22.0,
        )
        for s in (1.0, 0.477, 0.0):
            r = simulate_boulder_fall_2d(
                height_m=3.0, mass_kg=80.0, posture="controlled-drop",
                posture2d=Posture2D(ankle_load_share=s, **common),
            )
            assert r.peak_force_n > 0.0                          # sanity

        r1 = simulate_boulder_fall_2d(
            height_m=3.0, mass_kg=80.0, posture="controlled-drop",
            posture2d=Posture2D(ankle_load_share=1.0, **common),
        )
        r0477 = simulate_boulder_fall_2d(
            height_m=3.0, mass_kg=80.0, posture="controlled-drop",
            posture2d=Posture2D(ankle_load_share=0.477, **common),
        )
        r0 = simulate_boulder_fall_2d(
            height_m=3.0, mass_kg=80.0, posture="controlled-drop",
            posture2d=Posture2D(ankle_load_share=0.0, **common),
        )
        # 1D 轨迹逐位一致（不进入 rhs_2d）
        for attr in ("t_s", "z_foot_m", "z_torso_m", "pad_force_n",
                     "pad_compression_m", "peak_force_n"):
            assert np.array_equal(getattr(r1, attr), getattr(r0477, attr)), \
                f"share=0.477 与 share=1.0 在 {attr} 上不一致"
            assert np.array_equal(getattr(r1, attr), getattr(r0, attr)), \
                f"share=0.0 与 share=1.0 在 {attr} 上不一致"


# ==========================================================================
# P2 默认逐位一致 + 新输出归零
# ==========================================================================
class TestP2DefaultsBitIdentical:
    """全默认 ``Posture2D()`` ⇒ 1D 逐位一致 + 全部新踝输出 0/False/'none'。"""

    @pytest.mark.parametrize("height_m", [2.0, 5.0])
    def test_all_1d_fields_bit_identical(self, height_m):
        for posture in ("controlled-drop", "feet-first-stiff", "head-first"):
            r1d = simulate_boulder_fall(
                height_m=height_m, mass_kg=80.0, posture=posture)
            r2d = simulate_boulder_fall_2d(
                height_m=height_m, mass_kg=80.0, posture=posture,
                posture2d=Posture2D())
            for attr in ("t_s", "z_foot_m", "z_torso_m", "pad_compression_m",
                         "contact_area_m2", "pad_force_n", "accel_leg_g",
                         "accel_torso_g"):
                assert np.array_equal(getattr(r2d, attr), getattr(r1d, attr)), \
                    f"{attr} 非逐位一致 (h={height_m}, posture={posture})"

    def test_new_ankle_outputs_are_zero(self):
        r = simulate_boulder_fall_2d(
            height_m=3.0, mass_kg=80.0, posture="controlled-drop",
            posture2d=Posture2D())
        assert np.array_equal(r.ankle_beta_rad, np.zeros_like(r.t_s))
        assert np.array_equal(r.ankle_inversion_moment_nmm, np.zeros_like(r.t_s))
        assert r.peak_ankle_beta_rad == 0.0
        assert r.peak_inversion_moment_nmm == 0.0
        assert r.peak_ligament_strain == 0.0
        assert r.ankle_slide_frac == 0.0
        assert r.ankle_mode == "none"
        assert r.ankle_bone_utilization == 0.0
        assert r.ankle_fracture is False
        assert r.ankle_sprain is False

    def test_meta_mirrors_ankle_scalars(self):
        r = simulate_boulder_fall_2d(
            height_m=3.0, mass_kg=80.0, posture="controlled-drop",
            posture2d=Posture2D())
        for key in ("peak_ankle_beta_rad", "peak_inversion_moment_nmm",
                    "peak_ligament_strain", "ankle_slide_frac", "ankle_mode",
                    "ankle_bone_utilization", "ankle_fracture", "ankle_sprain",
                    "ankle_beta0_rad"):
            assert key in r.meta, f"meta 缺踝键 {key!r}"
        # 非默认下 meta 与属性一致（duck-typed 消费者口径）
        rn = _ankle_run(3.0, 20.0)
        assert rn.meta["ankle_mode"] == rn.ankle_mode
        assert rn.meta["peak_inversion_moment_nmm"] == rn.peak_inversion_moment_nmm
        assert rn.meta["peak_ligament_strain"] == rn.peak_ligament_strain
        assert rn.meta["ankle_bone_utilization"] == rn.ankle_bone_utilization
        assert rn.meta["ankle_sprain"] == rn.ankle_sprain
        assert rn.meta["ankle_fracture"] == rn.ankle_fracture
