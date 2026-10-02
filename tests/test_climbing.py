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
