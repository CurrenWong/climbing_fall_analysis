"""climbing 包的回归测试。

这些用例多数是为**锁死已修掉的 bug** 而写的 —— 每一个都对应一次
"输出看起来正常、但物理是错的"的静默失效。详见各用例注释。
"""

from __future__ import annotations

import numpy as np
import pytest

from climbing import G
from climbing.rope import RopeModel, simulate_rope_fall, default_rope
from climbing.belay import DEVICES, belayer_side_force, get_device
from climbing.pad import (
    POSTURES, CrashPad, Posture, hard_surface, simulate_boulder_fall,
)
from climbing.metrics import hic, hic_from_result, classify_g, classify_hic
from climbing.injury import assess_boulder_fall, PRESSURE_BANDS, BONE_BANDS


# ==========================================================================
# 缓存：抱石仿真是这套测试里的大头（单次 ~1-2 s），同一场景被多个用例引用。
# ==========================================================================
import functools


@functools.lru_cache(maxsize=None)
def boulder(posture: str, height: float = 3.0, mass: float = 80.0, on_pad: bool = True):
    b = simulate_boulder_fall(height_m=height, mass_kg=mass,
                              posture=posture, on_pad=on_pad)
    b.hic = hic_from_result(b)
    return b


# ==========================================================================
# metrics.hic
# ==========================================================================
class TestHIC:
    def test_constant_accel_known_value(self):
        """恒定加速度时 HIC 应等于 w * a^2.5（w 取窗口长度）。"""
        t = np.linspace(0, 0.5, 5000)
        a = np.full_like(t, 20.0)
        # 经典 HIC = (t2-t1) * (avg|a|)^2.5，15 ms 与 50 ms 窗口取最大。
        # 滑窗用步长扫描（step = n//8），靠近数组末端的窗口被截断，
        # 所以不要求严格相等。
        assert hic(t, a) == pytest.approx(0.05 * 20**2.5, rel=1e-3)

    def test_negative_accel_gives_nonzero(self):
        """回归：带符号的 a 取 2.5 次幂 -> NaN -> HIC 静默返回 0。

        自由落体段 a = -1g 曾经让整个 HIC 变成 0，看起来像"都不受伤"。
        """
        t = np.linspace(0, 0.5, 5000)
        assert hic(t, -np.ones_like(t)) > 0.0
        assert not np.isnan(hic(t, -np.ones_like(t)))

    def test_magnitude_equivalent(self):
        """HIC 用的是合矢量模长，±a 应当给出同样的结果。"""
        t = np.linspace(0, 0.5, 5000)
        assert hic(t, np.full_like(t, 30.0)) == pytest.approx(
            hic(t, -np.full_like(t, 30.0)), rel=1e-9)

    def test_rejects_nonfinite(self):
        t = np.linspace(0, 0.1, 100)
        with pytest.raises(ValueError):
            hic(t, np.array([np.nan] * 100))

    def test_rejects_bad_time(self):
        with pytest.raises(ValueError):
            hic(np.zeros(10), np.ones(10))


# ==========================================================================
# 姿势定义
# ==========================================================================
class TestPostures:
    def test_all_m_low_frac_strictly_inside_unit_interval(self):
        """回归：flat-flop / butt-impact 曾用 m_low_frac=1.00。

        那会让 m_up = 0，躯干方程退化成除以 1e-6，躯干一路自由落体
        （实测 2 s 掉 20 m），而峰值力看起来"正常"。
        """
        for name, p in POSTURES.items():
            assert 0.0 < p.m_low_frac < 1.0, f"{name} 的 m_low_frac={p.m_low_frac} 越界"

    def test_split_mass_raises_on_degenerate(self):
        p = Posture(name="bad", m_low_frac=1.0, k_flex=1e5, x_flex_max=0.1,
                    contact_schedule=((0, 0.1), (1, 0.2)), collapse_ref_m=0.2,
                    first_contact="x")
        with pytest.raises(ValueError, match="m_low_frac"):
            p.split_mass(80.0)

    def test_split_mass_positive(self):
        for name, p in POSTURES.items():
            m_low, m_up = p.split_mass(80.0)
            assert m_low > 0 and m_up > 0
            assert m_low + m_up == pytest.approx(80.0)

    def test_exactly_one_head_posture(self):
        heads = [n for n, p in POSTURES.items() if p.first_contact_is_head]
        assert heads == ["head-first"]

    def test_contact_area_monotone(self):
        for name, p in POSTURES.items():
            ks = [k for k, _ in p.contact_schedule]
            assert ks == sorted(ks), f"{name} 接触面积进度非单调"
            areas = [a for _, a in p.contact_schedule]
            assert areas == sorted(areas), f"{name} 接触面积非递增"
            assert min(areas) > 0


# ==========================================================================
# 绳索本构
# ==========================================================================
class TestRopeModel:
    def test_zero_strain_zero_tension(self):
        r = default_rope()
        assert r.tension(0.0) == 0.0
        assert r.tension(-0.1) == 0.0          # 绳松弛

    def test_tension_monotone_in_strain(self):
        r = default_rope()
        eps = np.linspace(0.0, 0.8, 200)
        t = r.tension(eps)
        assert np.all(np.diff(t) > 0)

    def test_eps_at_tension_inverts_tension(self):
        r = default_rope()
        for tn in [2000.0, 5000.0, 8000.0, 10000.0]:
            assert r.tension(r.eps_at_tension(tn)) == pytest.approx(tn, rel=1e-9)

    def test_stiffness_positive(self):
        assert default_rope().stiffness > 0


class TestRopeFall:
    @pytest.mark.parametrize("ff", [1.1, 1.25, 1.5, 1.77, 2.0])
    def test_peak_monotone_in_fall_factor(self, ff):
        """回归：峰值曾取全程 max，被回弹重绷污染，FF=1.5 反而 > FF=1.77。

        纯幂律 + 单向阻尼的绳几乎不耗能，3 s 内反复重绷，
        第二次重绷的瞬时冲击速度高于第一次。

        只对 FF > 1 要求严格递增：FF <= 1 时绳在释放瞬间就绷紧，
        自由落差为 0，峰值本来就与 FF 无关（见下一个用例）。
        """
        r = simulate_rope_fall(mass_kg=80, rope_length_m=2.0,
                              fall_factor=ff, t_hold_n=1e9)
        ref = simulate_rope_fall(mass_kg=80, rope_length_m=2.0,
                                 fall_factor=ff - 0.1, t_hold_n=1e9)
        assert r.peak_tension_n > ref.peak_tension_n, f"FF={ff} 峰值非单调"

    @pytest.mark.parametrize("ff", [0.3, 0.5, 0.75, 1.0])
    def test_peak_independent_of_ff_below_one(self, ff):
        """FF <= 1：没有自由落体段，绳一开始就吃力，峰值与 FF 无关。"""
        base = simulate_rope_fall(mass_kg=80, rope_length_m=2.0,
                                  fall_factor=ff, t_hold_n=1e9)
        for other in [0.3, 0.6, 1.0]:
            r = simulate_rope_fall(mass_kg=80, rope_length_m=2.0,
                                   fall_factor=other, t_hold_n=1e9)
            assert r.peak_tension_n == pytest.approx(base.peak_tension_n, rel=1e-9)

    def test_uiaa_criterion_at_standard_condition(self):
        """UIAA 101: 80 kg / FF=1.77 下峰值力应落在实测 8-9 kN 附近。"""
        r = simulate_rope_fall(mass_kg=80, rope_length_m=2.0,
                               fall_factor=1.77, t_hold_n=1e9)
        assert r.is_uiaa_pass
        assert 6.5 < r.peak_tension_kn < 11.0, r.peak_tension_kn
        # max_elongation_pct 的单位是 **百分数**，不是小数
        assert 15.0 < r.max_elongation_pct < 45.0, r.max_elongation_pct

    def test_taut_speed_matches_free_fall(self):
        """绳绷紧瞬间的速度必须等于自由落体速度。"""
        ff, l0 = 1.77, 2.0
        r = simulate_rope_fall(mass_kg=80, rope_length_m=l0,
                               fall_factor=ff, t_hold_n=1e9)
        assert r.free_fall_m == pytest.approx((ff - 1.0) * l0)
        assert np.abs(r.v_ms).max() == pytest.approx(
            np.sqrt(2 * G * r.free_fall_m), rel=0.05)

    def test_ff_le_1_has_no_free_fall(self):
        """FF<=1 时绳在释放瞬间就绷紧，没有自由落体段。"""
        for ff in [0.3, 0.5, 0.9]:
            r = simulate_rope_fall(mass_kg=80, rope_length_m=2.0,
                                   fall_factor=ff, t_hold_n=1e9)
            assert r.free_fall_m == 0.0

    def test_device_peak_is_capped_by_t_hold(self):
        """回归：旧默认峰值 5.58 kN 够不到任何 t_hold，保护器对比表四行完全相同。"""
        peaks = {}
        for name, d in DEVICES.items():
            r = simulate_rope_fall(mass_kg=80, rope_length_m=2.0, fall_factor=1.77,
                                   device=name, t_hold_n=d.t_hold_n)
            peaks[name] = r.peak_tension_n
            assert r.peak_tension_n <= d.t_hold_n + 1e-6 or d.t_hold_n > 1e8
        # 低保持力的设备必须真的被限幅
        assert peaks["Grigri"] == pytest.approx(4000.0, rel=1e-3)
        assert peaks["ATC-braked"] > peaks["Grigri"]

    def test_sliding_device_pays_out_rope(self):
        for name in ["Grigri", "ATC-nobrake", "Reverso"]:
            d = get_device(name)
            r = simulate_rope_fall(mass_kg=80, rope_length_m=2.0, fall_factor=1.77,
                                   device=name, t_hold_n=d.t_hold_n)
            assert r.max_payout_m > 0.0, f"{name} 应放绳"
            assert r.device_sat_time_s > 0.0

    def test_locked_device_does_not_pay_out(self):
        for name in ["ATC-braked", "Figure8-locked"]:
            d = get_device(name)
            r = simulate_rope_fall(mass_kg=80, rope_length_m=2.0, fall_factor=1.77,
                                   device=name, t_hold_n=d.t_hold_n)
            assert r.max_payout_m == pytest.approx(0.0, abs=1e-9)

    def test_peak_sits_inside_declared_window(self):
        """峰值必须落在 meta 声明的窗口内（窗口修复的回归测试）。"""
        r = simulate_rope_fall(mass_kg=80, rope_length_m=2.0, fall_factor=1.5,
                               t_hold_n=1e9)
        lo, hi = r.meta["peak_window_s"]
        assert lo <= r.t_s[int(np.argmax(r.tension_n))] <= hi

    def test_accel_is_in_g_units(self):
        """回归：绳侧 accel_g 一直是 m/s²/g，pad 侧曾经漏除，两边差 9.8 倍。"""
        r = simulate_rope_fall(mass_kg=80, rope_length_m=2.0, fall_factor=1.77,
                               t_hold_n=1e9)
        assert 1.0 < r.peak_accel_g < 60.0     # 若是 m/s² 会是几百


class TestBelay:
    def test_euler_amplification_increases_with_wrap_angle(self):
        prev = 0.0
        for th in [0.2, 0.6, 1.2, 2.0, 3.0]:
            f = belayer_side_force(8000.0, mu=0.20, theta_rad=th)
            assert f > prev
            prev = f

    def test_euler_matches_closed_form(self):
        import math
        f = belayer_side_force(8000.0, mu=0.20, theta_rad=1.2)
        assert f == pytest.approx(8000.0 * math.exp(0.20 * 1.2))

    def test_list_input_supported(self):
        out = belayer_side_force([1000.0, 2000.0], mu=0.2, theta_rad=1.0)
        assert out[1] == pytest.approx(2 * out[0])

    def test_unknown_device_raises(self):
        with pytest.raises(KeyError):
            get_device("nope")


# ==========================================================================
# 软垫本构
# ==========================================================================
class TestCrashPad:
    def test_stress_monotone(self):
        pad = CrashPad()
        xi = np.linspace(0, 0.9, 300)
        assert np.all(np.diff(pad.stress_pa(xi)) > 0)

    def test_zero_at_zero_strain(self):
        assert CrashPad().stress_pa(0.0) == pytest.approx(0.0)

    def test_matches_open_cell_pu_order_of_magnitude(self):
        """重标后的本构应落在 30 kg/m³ 开放孔 PU 的公开量级内。"""
        pad = CrashPad()
        for xi, lo, hi in [(0.10, 10e3, 20e3), (0.30, 30e3, 45e3),
                           (0.50, 60e3, 120e3), (0.70, 150e3, 300e3)]:
            s = float(pad.stress_pa(xi))
            assert lo <= s <= hi, f"xi={xi}: {s:.0f} Pa 不在 [{lo:.0f}, {hi:.0f}]"

    def test_hard_surface_is_linear(self):
        """回归：hard_surface 曾用饱和指数曲线，几厘米压缩下比"硬"值软 10 倍。"""
        h = hard_surface(1.0e7)
        for comp in [0.01, 0.05, 0.10]:
            assert float(h.stress_pa(comp)) == pytest.approx(1.0e7 * comp, rel=1e-6)

    def test_hard_surface_stiffer_than_pad(self):
        pad = CrashPad()
        hard = hard_surface()
        for xi in [0.05, 0.1, 0.2]:
            assert float(hard.stress_pa(xi)) > float(pad.stress_pa(xi))

    def test_capacity_positive(self):
        assert CrashPad().capacity_j(0.5) > 0.0


# ==========================================================================
# 抱石落地动力学
# ==========================================================================
class TestBoulderFall:
    def test_torso_does_not_free_fall(self):
        """回归：m_up=0 时躯干 2 s 掉 20 m（z_torso 0.85 -> 20.46）。"""
        for name in POSTURES:
            b = boulder(name)
            # 躯干下落量不可能超过落差 + 一点软垫压缩
            assert b.z_torso_m.max() < 3.0 + 0.6, f"{name}: {b.z_torso_m.max():.2f} m"

    def test_torso_starts_above_feet(self):
        """回归：z_torso 初值曾写成 +leg_ext0，把躯干放到了脚下方。

        z 是"自释放点起的向下位移"，向下为正。躯干在双足**之上**，
        所以它的 z 必须比脚**小**：z_torso = -leg_ext0。
        """
        b = simulate_boulder_fall(height_m=2.0, mass_kg=80, posture="controlled-drop")
        assert b.z_torso_m[0] == pytest.approx(-0.85)
        assert b.z_torso_m[0] < b.z_foot_m[0], "躯干初始应高于脚（z 向下为正）"
        # 腿长 = z_foot - z_torso，初始应为 0.85 m
        assert (b.z_foot_m[0] - b.z_torso_m[0]) == pytest.approx(0.85)

    def test_newton_third_law_com_acceleration(self):
        """回归：f_flex 是内力，曾被加进下段质量的"向上"方向。

        合外力变成 f_pad + 2*f_flex，系统凭空获得能量。
        正确的检验：质心加速度必须满足 m*a_com = f_pad - m*g。
        """
        b = boulder("controlled-drop")
        m = 80.0
        # 从结果数组重建两个质点的合外力，检查与质心加速度一致
        comp = b.pad_compression_m > 0
        idx = np.where(comp)[0]
        assert idx.size > 10
        mid = idx[len(idx) // 2]
        f_pad = b.pad_force_n[mid]
        a_leg_m = b.accel_leg_g[mid] * G
        a_torso_m = b.accel_torso_g[mid] * G
        # 质心减速度（向上为正）应等于 f_pad/m - g
        a_com = (a_leg_m + a_torso_m) / 2.0
        assert a_com == pytest.approx(f_pad / m - G, rel=0.25)

    def test_accel_in_g_units(self):
        """回归：accel_*_g 曾存 m/s²，与绳侧差 9.8 倍。"""
        b = boulder("controlled-drop")
        assert 1.0 < b.peak_accel_torso_g < 200.0

    def test_energy_absorbed_is_positive_and_bounded(self):
        """回归：np.trapezoid(f, x) 对弹性体净积分为 0 -> e_pad = -0.098 J。"""
        for name in POSTURES:
            b = boulder(name)
            assert b.energy_into_pad_j > 0.0, f"{name} 垫吸能为负"
            assert b.energy_into_flex_j > 0.0, f"{name} 屈曲吸能为负"
            tot = b.energy_into_pad_j + b.energy_into_flex_j
            # 允许略超 m*g*h（压缩过程中重力还在做功），但不能离谱
            assert tot < 1.5 * b.impact_energy_j, f"{name} 吸能超过输入"

    def test_flexion_within_posture_limit(self):
        for name in POSTURES:
            b = boulder(name)
            xf = np.maximum(0.0, 0.85 - (b.z_foot_m - b.z_torso_m))
            assert np.max(xf) < 0.85, f"{name} 腿被压穿（屈曲量 > 腿长）"

    def test_pad_reduces_pressure_vs_hard_ground(self):
        """软垫的核心作用：降低峰值接触压力。"""
        for name in POSTURES:
            on = boulder(name, on_pad=True)
            off = boulder(name, on_pad=False)
            pa = np.max(on.pad_force_n / np.maximum(on.contact_area_m2, 1e-9))
            pb = np.max(off.pad_force_n / np.maximum(off.contact_area_m2, 1e-9))
            assert pa < pb, f"{name}: 软垫压力 {pa:.0f} 未低于硬地 {pb:.0f}"

    def test_thicker_pad_is_safer(self):
        prev = None
        for th in [0.10, 0.20, 0.30]:
            b = simulate_boulder_fall(height_m=3.0, mass_kg=80,
                                      posture="controlled-drop",
                                      pad=CrashPad(thickness_m=th))
            p = np.max(b.pad_force_n / np.maximum(b.contact_area_m2, 1e-9))
            if prev is not None:
                assert p < prev, f"垫厚 {th}m 压力未下降"
            prev = p

    def test_impact_speed_matches_free_fall(self):
        for h in [0.5, 1.5, 3.0]:
            b = simulate_boulder_fall(height_m=h, mass_kg=80, posture="tuck-roll")
            assert b.impact_speed_ms == pytest.approx(np.sqrt(2 * G * h), rel=1e-6)

    def test_head_first_uses_head_for_hic(self):
        """回归：HIC 曾一律取躯干，头朝下算出 82 显示"全场最安全"，排名反了。"""
        b = boulder("head-first")
        assert b.meta["head_is_contact"] is True
        assert np.shares_memory(np.asarray(b.primary_accel_g),
                                np.asarray(b.accel_leg_g))
        assert b.peak_primary_g == pytest.approx(b.peak_accel_leg_g)


# ==========================================================================
# 损伤评估
# ==========================================================================
class TestInjury:
    def test_pressure_is_force_over_area(self):
        b = boulder("toe-point")
        r = assess_boulder_fall(b, "toe-point")
        p = b.pad_force_n / np.maximum(b.contact_area_m2, 1e-9)
        assert r.peak_pressure_kpa == pytest.approx(np.max(p) / 1e3, rel=1e-9)

    def test_small_contact_area_gives_high_pressure(self):
        """前脚掌点地的接触面积远小于平拍 -> 压力反而更高。"""
        toe = assess_boulder_fall(boulder("toe-point"), "toe-point")
        flop = assess_boulder_fall(boulder("flat-flop"), "flat-flop")
        assert toe.peak_pressure_kpa > flop.peak_pressure_kpa

    def test_good_technique_beats_bad_technique(self):
        """技术动作的核心收益：团身/屈膝应显著优于平拍。"""
        good = assess_boulder_fall(boulder("tuck-roll"), "tuck-roll")
        bad = assess_boulder_fall(boulder("flat-flop"), "flat-flop")
        assert good.peak_accel_g < 0.6 * bad.peak_accel_g
        assert good.hic < bad.hic

    def test_spine_uses_inertial_load_not_damping(self):
        """回归：脊柱载荷曾用 f_flex（含虚高阻尼项），虚高到 44-48 kN。"""
        b = boulder("feet-first-stiff")
        r = assess_boulder_fall(b, "feet-first-stiff")
        spine = [p for p in r.parts if p.name == "脊柱"][0]
        m_up = b.meta["m_up"]
        assert spine.metric == pytest.approx(
            m_up * G * b.peak_accel_torso_g / 1e3, rel=1e-6)

    def test_all_parts_have_mechanism(self):
        r = assess_boulder_fall(boulder("butt-impact"), "butt-impact")
        assert len(r.parts) >= 4
        for p in r.parts:
            assert p.band in {"低", "中", "高", "极高"}
            assert p.mechanism, f"{p.name} 缺少机制说明"

    def test_verdict_takes_worst_part(self):
        r = assess_boulder_fall(boulder("toe-point"), "toe-point")
        assert r.verdict in {"低", "中", "高", "极高"}

    def test_bands_are_ordered(self):
        for table in (PRESSURE_BANDS, BONE_BANDS):
            for a, b in zip(table, table[1:]):
                assert a.hi == b.lo, "分档之间必须连续不留空隙"


class TestClassify:
    def test_classify_g_bands(self):
        assert classify_g(5) == "低"
        assert classify_g(15) == "中"
        assert classify_g(30) == "高"
        assert classify_g(60) == "极高"

    def test_classify_hic_bands(self):
        assert classify_hic(100) == "低"
        assert classify_hic(500) == "中"
        assert classify_hic(800) == "高"
        assert classify_hic(1500) == "极高"

# ==========================================================================
# 静默失效：积分时长不足导致"撞击前就结束"
# ==========================================================================
# 早期版本 ``t_max`` 固定 2.0 s，而自由落体到接触需 t = √(2h/g)：
#   h = 20 m  ->  t = 2.019 s  >  2.0 s
# 于是 h ≥ 19.6 m 时积分在**撞击发生之前**就结束，模型静默返回
#   peak_force_n = 0 kN、max_compression_m = 0
# 不抛异常、字段齐全、数值看着像"这么高反而很安全"。
# 属于 阶段总结.md §3.1 那类静默失效：越危险的情形输出越无害。
#
# 现修复：t_max 自适应为 max(2.0, √(2h/g) + 1.0)；
# 并加一道守卫 —— 全程未接触且 h>0 时直接抛 RuntimeError。
class TestIntegrationHorizon:
    def test_high_fall_is_not_silently_zero(self):
        """h=25 m（自由落体 2.26 s）必须有非零峰值。"""
        r = simulate_boulder_fall(height_m=25.0, mass_kg=77.0,
                                  on_pad=False, posture="feet-first-stiff")
        assert r.peak_force_n > 0.0, "25 m 坠落报 0 峰值 —— 积分时长不足"
        assert r.max_compression_m > 0.0

    def test_peak_force_monotone_across_horizon(self):
        """跨过旧的 19.6 m 断点，峰值必须继续单调上升。"""
        hs = [16.0, 18.0, 20.0, 22.0]
        peaks = [
            simulate_boulder_fall(height_m=h, mass_kg=77.0, on_pad=False,
                                  posture="feet-first-stiff").peak_force_n
            for h in hs
        ]
        for h, pk in zip(hs, peaks):
            assert pk > 0.0, f"h={h} m 报 0 峰值"
        diffs = np.diff(peaks)
        assert np.all(diffs > 0.0), f"峰值未单调上升: {peaks}"

    def test_insufficient_t_max_raises(self):
        """显式给一个不够的 t_max，必须抛错而不是静默返回 0。"""
        with pytest.raises(RuntimeError, match="仍未发生接触"):
            simulate_boulder_fall(height_m=30.0, mass_kg=77.0, on_pad=False,
                                  posture="feet-first-stiff", t_max=1.0)


# ==========================================================================
# 分部位骨骼阈值（方案 v2 · P0）
# ==========================================================================
class TestBoneSites:
    def test_two_ends_weaker_than_shaft(self):
        """[Y25] Table 1 的核心关系：骨端（松质骨）弱于骨干（皮质骨）。"""
        from climbing.bone import MATERIAL_STRENGTH_MPA as M
        for ends, shaft in (("tibia_ends", "tibia_shaft"),
                            ("fibula_ends", "fibula_shaft")):
            se, ss = M[ends], M[shaft]
            assert se[0] < ss[0] and se[1] < ss[1], f"{ends} 应弱于 {shaft}"

    def test_site_threshold_unit_conversion(self):
        """1 MPa = 1 N/mm²，故阈值(N) = σ(MPa) × A(mm²)。"""
        from climbing.bone import BONE_SITES
        for s in BONE_SITES.values():
            assert s.force_threshold_n == pytest.approx(s.sigma_mpa * s.area_mm2)
            assert s.threshold_kn == pytest.approx(s.force_threshold_n / 1e3)

    def test_mode_selects_right_strength(self):
        """压缩模式必须选压缩强度（更大的那个）。"""
        from climbing.bone import BONE_SITES, MATERIAL_STRENGTH_MPA as M
        for s in BONE_SITES.values():
            st, sc = M[s.key]
            want = sc if s.mode == "compression" else st
            assert s.sigma_mpa == want

    def test_assess_sites_requires_impact_window(self):
        """缺失冲击窗口信息时必须报错，不能对全程取峰值。"""
        from climbing.bone import assess_sites

        class Fake:
            t_s = np.linspace(0, 1, 11)
            pad_force_n = np.zeros(11)
            accel_torso_g = np.zeros(11)
            height_m = 3.0
            mass_kg = 80.0
            posture = "x"
            meta: dict = {}

        with pytest.raises(ValueError, match="impact_window_s"):
            assess_sites(Fake())

    def test_high_fall_fractures_more_than_low(self):
        """利用率必须随高度上升（分部位评估的方向性检查）。"""
        from climbing.bone import assess_sites
        lo = assess_sites(simulate_boulder_fall(
            height_m=3.0, mass_kg=77.0, on_pad=False, posture="feet-first-stiff"))
        hi = assess_sites(simulate_boulder_fall(
            height_m=25.0, mass_kg=77.0, on_pad=False, posture="feet-first-stiff"))
        assert hi.worst.utilization > lo.worst.utilization

    def test_no_fracture_at_low_height_on_pad(self):
        """1 m 落在软垫上，不应有任何部位超过名义阈值。"""
        from climbing.bone import assess_sites
        v = assess_sites(simulate_boulder_fall(
            height_m=1.0, mass_kg=77.0, on_pad=True, posture="controlled-drop"))
        assert not v.fractured_sites, v.summary()


# ==========================================================================
# 场景库（方案 v2 · P1）
# ==========================================================================
from climbing.scenarios import (                                     # noqa: E402
    DIMENSIONS, HEIGHT_BANDS, LANDING_TO_POSTURE, MASS_CLIP_KG, MASS_MODEL,
    NAMED_SCENARIOS, OBSERVED_INJURY_LOCATION, OBSERVED_INJURY_TYPE,
    SURFACE_TO_ON_PAD, WALL_HEIGHT_M, ScenarioSampler, coverage_report,
)


class TestScenarioTables:
    def test_every_dimension_sums_to_one(self):
        """抄录校验：每个维度含 unknown 的占比之和应为 1（±1%）。

        回归目标：初稿把 [B25] **Table 4 的下肢列**当成全体分布抄了进来
        （脚先 87% 而非 73%、无旋转 39% 而非 30%）。表若混了，和会明显偏离 1。
        """
        for key, dim in DIMENSIONS.items():
            tot = dim.raw_total()
            assert abs(tot - 1.0) <= 0.01, f"{key} 占比和 = {tot:.3f}"

    def test_labels_unique_and_nonempty(self):
        for key, dim in DIMENSIONS.items():
            labels = [c.label for c in dim.cats]
            assert len(labels) == len(set(labels)), f"{key} 有重复标签"
            assert all(labels)
        # 除 landing_surface 外都应有 unknown 项：Table 3 的落点表面四项
        # 已经和为 1.00，没有 unknown 列；其余维度都有。
        without_unknown = [k for k, d in DIMENSIONS.items()
                           if not any(c.is_unknown for c in d.cats)]
        assert without_unknown == ["landing_surface"], without_unknown

    def test_ci_brackets_point_estimate(self):
        for key, dim in DIMENSIONS.items():
            for c in dim.cats:
                if c.ci is None:
                    continue
                lo, hi = c.ci
                assert lo < hi, f"{key}/{c.label} CI 反了"
                assert lo - 0.02 <= c.prob <= hi + 0.02, \
                    f"{key}/{c.label}: 点估计 {c.prob} 不在 CI {c.ci} 内"

    def test_landing_probability_drops_unknown_and_renormalizes(self):
        """剔除 unknown 后必须重新归一（否则采样会静默少抽）。"""
        labels, p = DIMENSIONS["rotation"].probs(drop_unknown=True)
        assert "unknown" not in labels
        assert p.sum() == pytest.approx(1.0)
        # [B25] Table 3：无旋转 30/92
        assert dict(zip(labels, p))["without"] == pytest.approx(0.30 / 0.92, rel=1e-9)

    def test_b25_table3_landing_is_73_percent_feet_first(self):
        """口径回归：脚先落地是 73%（Table 3），不是 87%（Table 4 下肢列）。"""
        lp = {c.label: c.prob for c in DIMENSIONS["landing_position"].cats}
        assert lp["standing_on_feet"] + lp["on_feet_leaning"] == pytest.approx(0.73)

    def test_no_rotation_is_thirty_percent(self):
        rot = {c.label: c.prob for c in DIMENSIONS["rotation"].cats}
        assert rot["without"] == pytest.approx(0.30)
        assert 1 - rot["without"] - rot["unknown"] == pytest.approx(0.62, abs=0.01)


class TestPostureMapping:
    def test_every_landing_category_is_mapped(self):
        """映射必须完备：每个落地类别都要能落到一个 POSTURES 姿势上。"""
        for cat in DIMENSIONS["landing_position"].cats:
            if cat.is_unknown:
                continue
            assert cat.label in LANDING_TO_POSTURE, f"{cat.label} 没有映射"

    def test_mapped_postures_exist_and_weights_sum_to_one(self):
        for lab, m in LANDING_TO_POSTURE.items():
            tot = 0.0
            for name, w in m.postures:
                assert name in POSTURES, f"{lab} 指向不存在的姿势 {name}"
                assert w > 0
                tot += w
            assert tot == pytest.approx(1.0), f"{lab} 权重和 {tot}"

    def test_every_surface_has_on_pad_decision(self):
        for cat in DIMENSIONS["landing_surface"].cats:
            assert cat.label in SURFACE_TO_ON_PAD, f"{cat.label} 没有落点判定"

    def test_gaps_between_pads_are_not_soft(self):
        """垫与墙之间 / 两垫之间都是硬缝隙 —— 不能当成软垫。"""
        assert SURFACE_TO_ON_PAD["between_pads_and_wall"] is False
        assert SURFACE_TO_ON_PAD["between_2_pads"] is False
        assert SURFACE_TO_ON_PAD["pads"] is True

    def test_ligament_risk_flagged_only_for_tilted_landing(self):
        """踝旋后机制只标在「在脚上且倾斜」上（[B25] 明确写的机制）。"""
        flagged = [k for k, m in LANDING_TO_POSTURE.items() if m.ligament_risk]
        assert flagged == ["on_feet_leaning"]

    def test_height_bands_inside_wall_and_ordered(self):
        """高度带必须互不重叠、递增，且上界不超过抱石墙高。"""
        bounds = sorted(HEIGHT_BANDS.values(), key=lambda x: x[0])
        for (lo, hi) in bounds:
            assert 0 < lo < hi <= WALL_HEIGHT_M + 1e-9
        for (_, hi_a), (lo_b, _) in zip(bounds, bounds[1:]):
            assert hi_a == pytest.approx(lo_b), "高度带之间必须首尾相接"

    def test_mass_model_probabilities_sum_to_one(self):
        assert sum(v[0] for v in MASS_MODEL.values()) == pytest.approx(1.0)
        for _, mu, sd in MASS_MODEL.values():
            assert MASS_CLIP_KG[0] < mu < MASS_CLIP_KG[1]
            assert sd > 0


class TestScenarioSampler:
    def test_seed_is_reproducible(self):
        a = ScenarioSampler(seed=42).sample_many(50)
        b = ScenarioSampler(seed=42).sample_many(50)
        assert a == b
        c = ScenarioSampler(seed=43).sample_many(50)
        assert a != c

    def test_samples_match_measured_distribution(self):
        """大样本下 MC 频率应贴合 [B25] Table 3（4σ 二项带）。"""
        n = 8000
        rows = ScenarioSampler(seed=1234).sample_many(n)
        for key, dim in DIMENSIONS.items():
            labels, p = dim.probs(drop_unknown=True)
            for lab, pr in zip(labels, p):
                if key == "landing_surface" and lab == "bump_into_someone":
                    continue                     # p 太小，4σ 带退化，跳过
                hits = sum(1 for r in rows if getattr(r, key) == lab)
                f = hits / n
                se = float(np.sqrt(pr * (1 - pr) / n))
                assert abs(f - pr) <= 4 * se + 0.01, \
                    f"{key}/{lab}: MC {f:.4f} vs 实测 {pr:.4f}"

    def test_unknown_never_sampled(self):
        rows = ScenarioSampler(seed=5).sample_many(2000)
        for r in rows:
            for k, v in r.to_dict().items():
                assert not v.startswith("unknown"), f"{k} 抽到了 unknown"

    def test_simulation_spec_is_feedable(self):
        s = ScenarioSampler(seed=9)
        for _ in range(300):
            spec = s.sample_simulation()
            assert spec.posture in POSTURES
            band = HEIGHT_BANDS[spec.scenario.height_band]
            assert band[0] <= spec.height_m <= band[1]
            assert MASS_CLIP_KG[0] <= spec.mass_kg <= MASS_CLIP_KG[1]
            assert isinstance(spec.on_pad, bool)

    def test_flags_are_consistent(self):
        rows = ScenarioSampler(seed=77).sample_many(500)
        for r in rows:
            assert r.rotation_free == (r.rotation == "without")
            assert r.foot_first == (
                r.landing_position in ("standing_on_feet", "on_feet_leaning"))
            assert r.upper_limb_mechanism == (
                r.start_position == "leaning_backward")
            assert r.on_pad == SURFACE_TO_ON_PAD[r.landing_surface]

    def test_weight_exact_is_conservative(self):
        """「结论可信」必须同时排除旋转、韧带机制、上肢机制。"""
        s = ScenarioSampler(seed=11)
        n = 4000
        for _ in range(n):
            sp = s.sample_simulation()
            assert sp.weight_exact == (
                sp.rotation_representable and not sp.ligament_dominant
                and not sp.upper_limb_mechanism)
            # 可信子集必然是无旋转的
            if sp.weight_exact:
                assert sp.scenario.rotation_free


class TestCoverage:
    def test_report_contains_key_numbers(self):
        txt = coverage_report()
        assert "P(无旋转) = 0.326" in txt
        assert "与 [B25] 0.73" not in txt          # 防呆：不该出现 87%
        assert "0.163" in txt                       # 骨折×下肢
        assert "踝（40%" in txt or "踝（40%，" in txt

    def test_observed_tables_sum_to_one(self):
        assert sum(v[0] for v in OBSERVED_INJURY_LOCATION.values()) \
            == pytest.approx(1.0, abs=0.01)
        assert sum(v[0] for v in OBSERVED_INJURY_TYPE.values()) \
            == pytest.approx(1.0, abs=0.01)

    def test_ankle_is_the_top_body_part_but_unmodelled(self):
        """[B25] 第一大部位是踝（40%），而模型部位表里没有踝。"""
        from climbing.bone import BONE_SITES
        assert OBSERVED_INJURY_TYPE["sprain"][0] == pytest.approx(0.36)
        assert max(OBSERVED_INJURY_TYPE["fracture"][1]) == 0.276
        assert not any("踝" == s.name_cn for s in BONE_SITES.values())

    def test_only_one_named_scenario_is_rotation_free(self):
        """[B25] Figure 4 的三个实测联合场景里，只有一个无旋转。"""
        free = [ns for ns in NAMED_SCENARIOS if ns.rotation_free]
        assert len(free) == 1
        assert free[0].share == pytest.approx(0.17)
        assert free[0].injury_location == "lower_limb"

    def test_fracture_share_used_for_coverage_matches_table(self):
        """覆盖度用的骨折占比必须来自 [B25] Table 1，不能手写。"""
        frac, ci, nn = OBSERVED_INJURY_TYPE["fracture"]
        assert frac == pytest.approx(68 / 301, abs=0.01)
        assert nn == 68


# ==========================================================================
# 快速求解器预设的等价性（P1 蒙特卡洛用）
# ==========================================================================
# 蒙特卡洛要跑 1000 次。pad.py 的默认求解设置单次要 ~5-6 s（t_max 固定
# max(2.0, t_fall+1.0)、max_ode_step=1e-4、rtol=1e-8），1000 次要 90 分钟。
# P1 因此用了收紧的 t_max + 放宽的步长/容差，单次降到 ~0.26 s。
#
# **放宽求解精度必须有测试守住等价性**，否则就是拿"跑得快"换"算得错"：
# 这是 阶段总结.md §3.7 那次教训（B-3 的 f_iso×12 实验）的同类风险。
class TestFastSolverPreset:
    @staticmethod
    def _fast_kw(h: float) -> dict:
        return dict(t_max=float(np.sqrt(2 * h / G)) + 0.25,
                    max_ode_step=1e-3, rtol=1e-6, atol=1e-9)

    @pytest.mark.parametrize("posture,h,on_pad", [
        ("feet-first-stiff", 2.20, True),
        ("one-leg-awkward", 3.35, False),
    ])
    def test_peak_force_and_utilisation_match_strict(self, posture, h, on_pad):
        from climbing.bone import assess_sites
        common = dict(height_m=h, mass_kg=75.0, posture=posture, on_pad=on_pad)
        fast = simulate_boulder_fall(**common, **self._fast_kw(h))
        strict = simulate_boulder_fall(**common)

        assert fast.peak_force_n == pytest.approx(
            strict.peak_force_n, rel=1e-3), "峰值力在快速预设下变了"

        vf = {s.key: s.utilization for s in assess_sites(fast).sites}
        vs = {s.key: s.utilization for s in assess_sites(strict).sites}
        for k in vf:
            assert vf[k] == pytest.approx(vs[k], rel=2e-3, abs=2e-4), \
                f"{k} 的利用率在快速预设下变了"

    def test_truncated_horizon_still_captures_contact(self):
        """t_max = t_fall + 0.25 必须够到撞击，否则守卫会抛错。"""
        r = simulate_boulder_fall(height_m=4.5, mass_kg=90.0,
                                  posture="feet-first-stiff",
                                  **self._fast_kw(4.5))
        assert r.peak_force_n > 0
        assert r.max_compression_m > 0

    def test_horizon_that_misses_contact_raises(self):
        """余量不够时守卫必须抛错 —— 不许静默给 0。"""
        h = 4.5
        with pytest.raises(RuntimeError, match="仍未发生接触"):
            simulate_boulder_fall(height_m=h, mass_kg=90.0,
                                  posture="feet-first-stiff",
                                  t_max=float(np.sqrt(2 * h / G)) - 0.05,
                                  max_ode_step=1e-3, rtol=1e-6, atol=1e-9)

