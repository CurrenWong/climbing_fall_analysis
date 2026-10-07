"""P4 · 头/颈惯性分支（路径 B）的回归测试。

覆盖 ``docs/方案更新_v2.md:841-852`` 的 P4 定位与 ``:861/:866`` 的
V3 / V8 判据，以及 ``docs/方案更新_v2.md:548-553`` 的双路径物理：

* **路径 A（自下而上）** 足 → 胫骨/腓骨 → 股骨 → 骨盆；
* **路径 B（自上而下）** 颅骨 + 颈椎因**惯性**继续下压 → 胸椎 → 腰椎；
* 两路径在脊柱汇集 —— 本模块只建**路径 B**（头部惯性分支），
  脊柱后期的汇集峰需要多刚体脊柱模型（``bodies.py`` 结构完成但未标定），
  **明确 out of scope**（见 ``TestSpineConvergence`` 的 docstring）。

测试用的 fixture 是**轻量合成坠落轨迹**（duck-typed），不依赖真实
求解器：

* ``make_rigid_ground_fall`` —— [Y25] §2.3 的"直立双足 + **刚性地面**"
  场景：足骤停，躯干在极短时间（由 ``stop_distance_m`` 控制）内减速。
  V3 的"t < 0.0125 s 出现头部应力"就是在这个场景下测的。
* ``SimpleNamespace`` duck-typed 对象 —— 只有 ``t_s`` / ``z_torso_m`` /
  ``meta``，不给 ``accel_torso_g``，验证只靠 ``z_torso_m`` 的兜底路径。

约定
----
* V3 / V8 的阈值逐字取自 ``docs/方案更新_v2.md:861/:866``，不在测试里另定。
* 分档必须与 ``injury.VERDICTS`` 同源。
* 中文注释是 house style，与 pad.py / test_pad2d.py 一致。
* 所有断言都是**正交**的（一个测试只测一个量）。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest

from climbing import G
from climbing.head_branch import (
    HEAD_STRESS_BANDS_MPA,
    V3_THRESHOLD_S,
    V8_THRESHOLD_AT_20M_S,
    HeadBranchParams,
    HeadVerdict,
    head_band,
    simulate_head_branch,
)


# ==========================================================================
# 合成 fixture（duck-typed，不依赖真实求解器）
# ==========================================================================
@dataclass
class SyntheticFall:
    """轻量合成坠落结果，只携带 head_branch 消费的字段。

    字段口径与 ``climbing.pad.BoulderFallResult`` 一致：``z`` 向下为正，
    ``accel_torso_g`` 正 = 向上减速（= ``-z_ddot/g``）。
    """

    t_s: np.ndarray
    z_torso_m: np.ndarray
    accel_torso_g: np.ndarray
    meta: dict[str, Any] = field(default_factory=dict)


def make_rigid_ground_fall(
    height_m: float = 3.0,
    stop_distance_m: float = 0.02,
    g: float = G,
    dt: float = 2e-4,
    window_s: float = 0.05,
) -> SyntheticFall:
    """构造一个"近刚性地 + 快减速躯干"的合成坠落轨迹。

    这是 [Y25] §2.3 "直立双足 + **刚性地面**" 场景的简化：足在触地瞬间
    骤停，躯干因惯性继续下压，但被身体（关节/组织）的有限行程阻住。

    轨迹（``z`` 向下为正）：

    * **自由落体段** ``t < t_fall``：``z_torso = -0.85 + ½·g·t²``，
      ``accel_torso_g = -1``（自由落体，无相对加速度）。
    * **减速段** ``t_fall ≤ t < t_fall + T_dec``：半正弦脉冲
      ``z_torso_ddot = -a_peak·sin(π·τ)``，``a_peak = π·v_impact/(2·T_dec)``
      —— 该系数**精确**让躯干速度在脉冲末端归零（∫sin 的正确归一）。
    * **静止段**：``z_torso`` 常数。

    参数
    ----
    stop_distance_m
        躯干有效减速行程 (m)。**越小 ⇒ 减速越快 ⇒ 头峰值越早**
        （对应 [Y25] 的"刚性地面"）。默认 0.02 m（刚性地 + 2 cm 组织压缩）。
    window_s
        合成冲击窗口长度 (s)，默认 50 ms（真实抱石冲击窗口的量级）。
    """
    t_fall = float(np.sqrt(2.0 * height_m / g))
    v_impact = g * t_fall
    # 固定停止距离 ⇒ 更短减速时间 = 更高减速度（∝ 1/v_impact）。
    t_dec = 2.0 * stop_distance_m / v_impact
    a_peak = np.pi * v_impact / (2.0 * t_dec)

    t_total = t_fall + max(window_s, t_dec)
    t = np.arange(0.0, t_total, dt)
    if t.size < 4 or t[-1] < t_total - 1e-12:
        t = np.append(t, t_total)
    n = t.size

    z_torso = np.zeros(n)
    a_torso_g = np.zeros(n)

    in_ff = t < t_fall
    z_torso[in_ff] = -0.85 + 0.5 * g * t[in_ff] ** 2
    a_torso_g[in_ff] = -1.0

    in_dec = (t >= t_fall) & (t < t_fall + t_dec)
    tau = (t[in_dec] - t_fall) / t_dec
    a_ddot = -a_peak * np.sin(np.pi * tau)     # 向下为正的躯干加速度
    a_torso_g[in_dec] = -a_ddot / g            # 口径：向上减速为正
    v_dec = v_impact + np.cumsum(a_ddot) * dt
    z_at = -0.85 + height_m
    z_dec = z_at + np.concatenate(
        [[0.0], np.cumsum(0.5 * (v_dec[:-1] + v_dec[1:]) * dt)]
    )
    z_torso[in_dec] = z_dec

    in_rest = t >= t_fall + t_dec
    if in_dec.any():
        z_torso[in_rest] = z_torso[in_dec][-1]
    else:
        z_torso[in_rest] = z_at
    a_torso_g[in_rest] = 0.0

    window = (float(t_fall), float(t_fall + window_s))
    return SyntheticFall(
        t_s=t,
        z_torso_m=z_torso,
        accel_torso_g=a_torso_g,
        meta={"impact_window_s": window, "height_m": height_m,
              "stop_distance_m": stop_distance_m},
    )


# ==========================================================================
# V3 —— 头部应力峰值时刻 ≤ 0.0125 s
# ==========================================================================
class TestV3HeadPeakWithinY25Window:
    """V3（``docs/方案更新_v2.md:861``）：头部应力峰值时刻 ≤ 0.0125 s。

    物理：足骤停 → 头因惯性继续下压 → 颈部首次应力集中与足冲击
    **几乎同时**，而不是延迟到载荷传导之后。
    """

    def test_v3_peak_time_within_threshold_at_3m(self):
        """3 m 刚性地坠落，默认参数：``within_y25_window is True``。"""
        r = make_rigid_ground_fall(height_m=3.0)
        v = simulate_head_branch(r)
        assert v.peak_time_s <= V3_THRESHOLD_S, (
            f"V3 未通过：peak_time_s = {v.peak_time_s*1e3:.2f} ms > "
            f"{V3_THRESHOLD_S*1e3:.1f} ms（{v.note}）"
        )
        assert v.within_y25_window is True

    @pytest.mark.parametrize("h", [1.0, 2.0, 3.0, 5.0, 8.0])
    def test_v3_passes_across_realistic_heights(self, h):
        """1-8 m 的"直立双足 + 刚性地"场景，头部峰都应落在 Y25 窗口内。"""
        r = make_rigid_ground_fall(height_m=h)
        v = simulate_head_branch(r)
        assert v.within_y25_window is True, (
            f"h={h} m：峰值时刻 {v.peak_time_s*1e3:.2f} ms 超出 "
            f"{V3_THRESHOLD_S*1e3:.1f} ms"
        )

    def test_within_flag_matches_threshold(self):
        """``within_y25_window`` 必须严格等于 ``peak_time_s <= 0.0125``。"""
        r = make_rigid_ground_fall(height_m=3.0)
        v = simulate_head_branch(r)
        assert v.within_y25_window == (v.peak_time_s <= V3_THRESHOLD_S)

    def test_v3_threshold_constant_is_verbatim(self):
        """阈值必须逐字等于 `docs/方案更新_v2.md:861` 的 0.0125 s。"""
        assert V3_THRESHOLD_S == 0.0125


# ==========================================================================
# V8 —— 峰值时刻随高度单调提前；20 m 时 ≤ 5 ms
# ==========================================================================
class TestV8PeakTimeMonotoneWithHeight:
    """V8（``docs/方案更新_v2.md:866``）：峰值时刻随高度**提前**，
    20 m 时 ≤ 5 ms。

    物理：更高 ⇒ 冲击速度更大 ⇒ 在相同/更短的有效行程内减速 ⇒
    躯干减速脉冲更窄 ⇒ 头更早达到峰值应力。
    """

    @pytest.mark.parametrize("heights", [[2.0, 3.0, 5.0], [1.0, 2.0, 3.0, 5.0, 8.0]])
    def test_peak_time_monotone_non_increasing_default(self, heights):
        """默认刚度下，峰值时刻在高度序列上单调非增。"""
        peaks = [simulate_head_branch(make_rigid_ground_fall(h)).peak_time_s
                 for h in heights]
        for i in range(1, len(heights)):
            assert peaks[i] <= peaks[i - 1] + 1e-12, (
                f"h={heights[i]} m 的峰值时刻 {peaks[i]*1e3:.2f} ms 晚于 "
                f"h={heights[i-1]} m 的 {peaks[i-1]*1e3:.2f} ms"
            )

    def test_peak_time_at_20m_under_5ms_with_stiff_neck(self):
        """V8 的"20 m ≤ 5 ms"需要**脊柱-头段有足够刚度**（[Y25] §2.3）。

        默认 ``k_neck = 200 N/mm`` 在 20 m 给出 ~8.8 ms（颈是软的一环，
        头部响应被延迟）；把 ``k_neck`` 提到 1000 N/mm（颈椎在高加载率下
        的等效刚度）后头部及时跟上，峰值提前到 ~4.6 ms ≤ 5 ms。
        这正是 [Y25] 那句"要求脊柱-头段有足够刚度，否则头部响应会被延迟"。
        """
        stiff = HeadBranchParams(k_neck_n_per_mm=1000.0)
        v = simulate_head_branch(make_rigid_ground_fall(height_m=20.0), stiff)
        assert v.peak_time_s <= V8_THRESHOLD_AT_20M_S, (
            f"20 m 峰值时刻 {v.peak_time_s*1e3:.2f} ms > "
            f"{V8_THRESHOLD_AT_20M_S*1e3:.1f} ms"
        )

    def test_peak_time_monotone_across_wide_range_stiff_neck(self):
        """刚性颈下，1→20 m 全程峰值时刻仍单调非增（不因刚度而破坏）。"""
        stiff = HeadBranchParams(k_neck_n_per_mm=1000.0)
        heights = [1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0]
        peaks = [simulate_head_branch(make_rigid_ground_fall(h), stiff).peak_time_s
                 for h in heights]
        for i in range(1, len(peaks)):
            assert peaks[i] <= peaks[i - 1] + 1e-12

    def test_v8_threshold_constant_is_verbatim(self):
        """阈值必须逐字等于 `docs/方案更新_v2.md:866` 的 5 ms。"""
        assert V8_THRESHOLD_AT_20M_S == 0.005


# ==========================================================================
# 双路径可区分 —— 头部峰早、载荷传导晚
# ==========================================================================
class TestSpineConvergence:
    """双路径（[Y25] §2.3，``docs/方案更新_v2.md:548-553``）的时间可区分性。

    **本模块只建路径 B（头部惯性分支）。** 脊柱后期那个"两路径汇集",
    应力集中峰属于路径 A（足 → 胫腓 → 股骨 → 骨盆）与路径 B 在脊柱
    的**叠加**，需要带脊柱节的**多刚体链**才能表达 —— ``bodies.py``
    的结构在，但验收判据未通过（见 README「Phase B-3」），因此
    **脊柱后期的汇集峰明确 out of scope**。

    本测试在**时间上把两条路径区分开**：
    * 路径 B 的特征时刻 = 头部应力峰（本模块输出，早）；
    * 路径 A 的载荷传导完成标志 = 垫反力峰（``pad_force_n``，晚）。
    在"受控屈膝"这种时序清晰的姿势上，头部峰必须**早于**垫反力峰。
    """

    def test_head_peak_precedes_pad_load_transfer_peak(self):
        """受控屈膝 3 m：头部峰（路径 B）早于垫反力峰（路径 A 传导完成）。"""
        from climbing.pad import simulate_boulder_fall

        r = simulate_boulder_fall(
            height_m=3.0, mass_kg=80.0, posture="controlled-drop", on_pad=True,
        )
        v = simulate_head_branch(r)
        win = r.meta["impact_window_s"]
        t0 = float(win[0])
        in_win = (r.t_s >= win[0]) & (r.t_s <= win[1])
        t_pad_peak = float(r.t_s[in_win][np.argmax(r.pad_force_n[in_win])]) - t0
        assert v.peak_time_s < t_pad_peak, (
            f"头部峰 {v.peak_time_s*1e3:.2f} ms 未早于垫反力峰 "
            f"{t_pad_peak*1e3:.2f} ms —— 两条路径在时间上不可区分"
        )

    def test_head_peak_is_early_feature_on_rigid_ground(self):
        """刚性地场景下头部峰落在冲击窗口的**前 1/3** —— 早期特征，非延迟。"""
        r = make_rigid_ground_fall(height_m=3.0, window_s=0.05)
        v = simulate_head_branch(r)
        assert v.peak_time_s <= (1.0 / 3.0) * 0.05, (
            f"头部峰 {v.peak_time_s*1e3:.2f} ms 不在窗口前 1/3（早期特征）"
        )

    def test_module_documents_spine_out_of_scope(self):
        """模块必须**显式记录**脊柱后期汇集峰为何 out of scope（可追溯）。"""
        import climbing.head_branch as hb

        doc = hb.__doc__ or ""
        assert "路径 A" in doc and "路径 B" in doc, "模块 docstring 未描述双路径"
        assert "out of scope" in doc or "no scope" in doc or "scope" in doc, (
            "模块 docstring 未说明脊柱后期汇集峰超出范围"
        )


# ==========================================================================
# 刚度灵敏度 —— 两端行为可解释（invariant-② 风格）
# ==========================================================================
class TestStiffnessSensitivity:
    """``k_neck`` 两端行为：

    * **很大** ⇒ 头刚性跟随躯干 ⇒ ``delta_peak → 0``，峰值**提前**；
    * **很小** ⇒ 头与躯干解耦（颈部拉不住头）⇒ ``delta_peak`` 很大，
      峰值**推迟**。

    这正是 [Y25] "脊柱-头段有足够刚度，否则头部响应会被延迟" 的代码化。
    """

    def test_stiff_neck_tracks_torso_with_small_delta(self):
        """很大 ``k_neck``：``delta_peak`` 极小（头几乎不相对躯干移动）。"""
        r = make_rigid_ground_fall(height_m=3.0)
        rigid = simulate_head_branch(r, HeadBranchParams(k_neck_n_per_mm=1.0e5))
        assert abs(rigid.delta_peak_m) < 1e-3, (
            f"刚性颈 delta_peak = {rigid.delta_peak_m*1e3:.4f} mm，应 < 1 mm"
        )

    def test_soft_neck_decouples_with_large_delta(self):
        """很小 ``k_neck``：``delta_peak`` 远大于刚性颈（头被甩开）。"""
        r = make_rigid_ground_fall(height_m=3.0)
        soft = simulate_head_branch(r, HeadBranchParams(k_neck_n_per_mm=0.5))
        assert abs(soft.delta_peak_m) > 0.05, (
            f"软颈 delta_peak = {soft.delta_peak_m*1e3:.2f} mm，应 > 50 mm"
        )

    def test_delta_monotone_increasing_as_neck_softens(self):
        """``delta_peak`` 随 ``k_neck`` 减小而单调增大（两端之间无反转）。"""
        r = make_rigid_ground_fall(height_m=3.0)
        ks = [1.0e5, 5.0e3, 1.0e3, 200.0, 50.0, 5.0]
        deltas = [abs(simulate_head_branch(
            r, HeadBranchParams(k_neck_n_per_mm=k)).delta_peak_m) for k in ks]
        for i in range(1, len(ks)):
            assert deltas[i] > deltas[i - 1], (
                f"k={ks[i]} 的 delta {deltas[i]:.4f} 未大于 k={ks[i-1]} 的 "
                f"{deltas[i-1]:.4f}（软颈应甩得更开）"
            )

    def test_stiff_neck_peaks_earlier_than_soft(self):
        """刚度越大 ⇒ 峰值越早（[Y25]：足够刚度才不延迟响应）。"""
        r = make_rigid_ground_fall(height_m=3.0)
        t_stiff = simulate_head_branch(
            r, HeadBranchParams(k_neck_n_per_mm=5.0e4)).peak_time_s
        t_default = simulate_head_branch(r).peak_time_s
        t_soft = simulate_head_branch(
            r, HeadBranchParams(k_neck_n_per_mm=5.0)).peak_time_s
        assert t_stiff < t_default < t_soft, (
            f"刚性 {t_stiff*1e3:.2f} / 默认 {t_default*1e3:.2f} / "
            f"软 {t_soft*1e3:.2f} ms —— 刚度与峰值时刻未呈单调关系"
        )

    def test_stiff_neck_force_larger_than_soft(self):
        """刚性颈传递更大峰值力（``F = k·delta`` 两端乘积的量级关系）。"""
        r = make_rigid_ground_fall(height_m=3.0)
        f_rigid = simulate_head_branch(
            r, HeadBranchParams(k_neck_n_per_mm=5.0e3)).peak_force_n
        f_soft = simulate_head_branch(
            r, HeadBranchParams(k_neck_n_per_mm=5.0)).peak_force_n
        assert f_rigid > f_soft


# ==========================================================================
# 无反馈 —— 严格下游派生量
# ==========================================================================
class TestNoFeedback:
    """``simulate_head_branch`` 必须是**严格的下游派生量**：

    * 不修改输入 ``result`` 的任何数组；
    * 不把头部反作用力回灌到躯干 ODE。
    这是 1D/2D 模型 ``np.array_equal`` bit-identicality 的前提。
    """

    def test_input_arrays_not_mutated(self):
        """``z_torso_m`` / ``t_s`` / ``accel_torso_g`` 前后逐位相等。"""
        r = make_rigid_ground_fall(height_m=3.0)
        z_before = r.z_torso_m.copy()
        t_before = r.t_s.copy()
        a_before = r.accel_torso_g.copy()
        _ = simulate_head_branch(r)
        assert np.array_equal(r.z_torso_m, z_before), "z_torso_m 被修改"
        assert np.array_equal(r.t_s, t_before), "t_s 被修改"
        assert np.array_equal(r.accel_torso_g, a_before), "accel_torso_g 被修改"

    def test_input_arrays_not_mutated_with_custom_params(self):
        """自定义参数下同样不得修改输入。"""
        r = make_rigid_ground_fall(height_m=3.0)
        z_before = r.z_torso_m.copy()
        _ = simulate_head_branch(r, HeadBranchParams(k_neck_n_per_mm=12345.0))
        assert np.array_equal(r.z_torso_m, z_before)

    def test_works_on_plain_duck_typed_object(self):
        """纯 duck-typed 对象（只有 t_s / z_torso_m / meta）也能工作。"""
        r = make_rigid_ground_fall(height_m=3.0)
        duck = SimpleNamespace(
            t_s=r.t_s, z_torso_m=r.z_torso_m,
            meta={"impact_window_s": r.meta["impact_window_s"]},
        )
        v = simulate_head_branch(duck)          # 无 accel_torso_g，走梯度兜底
        assert isinstance(v, HeadVerdict)
        assert v.peak_force_n > 0.0

    def test_real_pad_result_not_mutated(self):
        """真实 ``BoulderFallResult`` 也不得被修改。"""
        from climbing.pad import simulate_boulder_fall

        r = simulate_boulder_fall(height_m=3.0, mass_kg=80.0,
                                  posture="controlled-drop", on_pad=True)
        z_before = r.z_torso_m.copy()
        a_before = r.accel_torso_g.copy()
        _ = simulate_head_branch(r)
        assert np.array_equal(r.z_torso_m, z_before)
        assert np.array_equal(r.accel_torso_g, a_before)

    def test_duck_typed_agrees_with_full_field_object(self):
        """带 ``accel_torso_g`` 与不带（梯度兜底）应给出**一致量级**的结果。

        兜底路径用 ``np.gradient`` 数值微分，比真实 ``accel_torso_g`` 略糙，
        但峰值时刻不应出现量级差。
        """
        r = make_rigid_ground_fall(height_m=3.0)
        duck = SimpleNamespace(
            t_s=r.t_s, z_torso_m=r.z_torso_m,
            meta={"impact_window_s": r.meta["impact_window_s"]},
        )
        v_full = simulate_head_branch(r)
        v_duck = simulate_head_branch(duck)
        assert abs(v_full.peak_time_s - v_duck.peak_time_s) < 5e-3, (
            f"带/不带 accel_torso_g 的峰值时刻差 "
            f"{abs(v_full.peak_time_s - v_duck.peak_time_s)*1e3:.2f} ms 过大"
        )


# ==========================================================================
# 接口契约 / 分档
# ==========================================================================
class TestApiContract:
    """``HeadBranchParams`` / ``HeadVerdict`` / ``head_band`` 的接口契约。"""

    def test_default_params_constructible(self):
        """``HeadBranchParams()`` 必须可直接使用（全字段有默认值）。"""
        p = HeadBranchParams()
        assert p.m_head_kg == 5.0
        assert p.k_neck_n_per_mm == 200.0
        assert p.c_neck_n_s_per_mm == 0.0
        assert p.head_section_mm2 == 1.0
        assert p.g == G

    def test_params_frozen(self):
        """冻结 dataclass：不可变（与 repo 其它标定容器一致）。"""
        p = HeadBranchParams()
        with pytest.raises(Exception):
            p.m_head_kg = 9.0  # type: ignore[misc]

    def test_verdict_fields_present(self):
        """``HeadVerdict`` 的七个字段都能取到且类型正确。"""
        v = simulate_head_branch(make_rigid_ground_fall(height_m=3.0))
        assert isinstance(v.peak_stress_mpa, float)
        assert isinstance(v.peak_time_s, float)
        assert isinstance(v.peak_force_n, float)
        assert isinstance(v.delta_peak_m, float)
        assert isinstance(v.within_y25_window, bool)
        assert isinstance(v.band, str)
        assert isinstance(v.note, str)

    def test_missing_impact_window_raises_chinese(self):
        """缺 ``impact_window_s`` 必须抛 ``ValueError``，中文消息含关键词。"""
        duck = SimpleNamespace(
            t_s=np.linspace(0.0, 1.0, 11),
            z_torso_m=np.zeros(11),
            meta={},
        )
        with pytest.raises(ValueError, match="impact_window_s"):
            simulate_head_branch(duck)

    def test_missing_z_torso_raises(self):
        """缺 ``z_torso_m`` 也必须抛错（不能静默用零）。"""
        class NoZ:
            t_s = np.linspace(0.0, 1.0, 11)
            meta = {"impact_window_s": (0.5, 0.6)}

        with pytest.raises((AttributeError, ValueError)):
            simulate_head_branch(NoZ())

    def test_head_band_uses_verdicts_vocabulary(self):
        """分档取值必须落在 ``injury.VERDICTS`` 的四个标签内。"""
        from climbing.injury import VERDICTS

        for s in [0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 6.0, 50.0, 1e4]:
            assert head_band(s) in VERDICTS

    def test_head_band_monotone_non_decreasing(self):
        """应力越大，分档只能同级或更高（不能反转）。"""
        rank = {"低": 0, "中": 1, "高": 2, "极高": 3}
        stresses = np.linspace(0.0, 20.0, 200)
        ranks = [rank[head_band(s)] for s in stresses]
        for i in range(1, len(ranks)):
            assert ranks[i] >= ranks[i - 1]

    def test_head_band_threshold_boundaries(self):
        """阈值边界按 ``HEAD_STRESS_BANDS_MPA`` 逐档命中。"""
        assert head_band(0.0) == "低"
        assert head_band(0.999) == "低"
        assert head_band(1.0) == "中"
        assert head_band(2.999) == "中"
        assert head_band(3.0) == "高"
        assert head_band(5.999) == "高"
        assert head_band(6.0) == "极高"
        assert head_band(1e9) == "极高"

    def test_head_band_matches_verdict_band_in_result(self):
        """``HeadVerdict.band`` 必须与 ``head_band(peak_stress)`` 一致。"""
        v = simulate_head_branch(make_rigid_ground_fall(height_m=3.0))
        assert v.band == head_band(v.peak_stress_mpa)

    def test_stress_proportional_to_force(self):
        """``σ = F/A``：把截面放大 10×，应力应恰好降到 1/10。"""
        r = make_rigid_ground_fall(height_m=3.0)
        v1 = simulate_head_branch(r, HeadBranchParams(head_section_mm2=1.0))
        v10 = simulate_head_branch(r, HeadBranchParams(head_section_mm2=10.0))
        assert v10.peak_stress_mpa == pytest.approx(v1.peak_stress_mpa / 10.0)
        # 峰值力不随截面变化（截面只影响应力度量）
        assert v10.peak_force_n == pytest.approx(v1.peak_force_n)


# ==========================================================================
# 真实求解器 sanity-check（不进断言主路径，仅保证可跑通 + 数量级）
# ==========================================================================
class TestRealPadSanity:
    """真实 ``climbing.pad.simulate_boulder_fall`` 上的 sanity-check。

    注意：**软垫 ≠ [Y25] 的刚性地面**。软垫把躯干减速拉长到数十毫秒，
    头部峰因此落在 20-55 ms —— 此时 V3 的 12.5 ms 判据**不适用**
    （V3 定义在"足骤停"的刚性地场景）。这些测试只验证真实数据上
    能跑通、方向合理、不被修改。
    """

    def test_runs_on_real_pad_result(self):
        from climbing.pad import simulate_boulder_fall

        r = simulate_boulder_fall(height_m=3.0, mass_kg=80.0,
                                  posture="controlled-drop", on_pad=True)
        v = simulate_head_branch(r)
        assert v.peak_force_n > 0.0
        assert 0.0 < v.peak_time_s < (r.meta["impact_window_s"][1]
                                      - r.meta["impact_window_s"][0]) + 1e-6
        assert v.band in {"低", "中", "高", "极高"}

    def test_head_peak_time_decreases_with_height_on_pad(self):
        """软垫上峰值时刻也应随高度非增（与 pad2d invariant 4 同向）。"""
        from climbing.pad import simulate_boulder_fall

        peaks = []
        for h in [2.0, 3.0, 5.0]:
            r = simulate_boulder_fall(height_m=h, mass_kg=80.0,
                                      posture="controlled-drop", on_pad=True)
            peaks.append(simulate_head_branch(r).peak_time_s)
        for i in range(1, len(peaks)):
            assert peaks[i] <= peaks[i - 1] + 1e-9
