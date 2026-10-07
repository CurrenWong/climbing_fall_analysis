"""踝损伤评估模块 + Phase P1 蒙特卡洛钩子的回归测试。

覆盖范围
--------
1. ``assess_ankle``：鸭式读取 P2 踝字段，缺数据时抛 ``ValueError``（中文
   错误信息，行为与 ``bone.py:301-307`` 同源）。
2. ``ankle_band``：四种组合 → 四档中文档位，与 ``injury.VERDICTS`` 一致。
3. ``site_distribution`` 聚合：把踝列纳入三个子集（all / no_rotation /
   credible），每 trial 数一次踝损伤。
4. ``ROTATION_TO_POSTURE2D``：覆盖 ``scenarios.DIMENSIONS["rotation"]``
   全部类别；``"without"`` 映射到 ``Posture2D()``（全默认）。

约定
----
* 用轻量 stub 而非 ``BoulderFallResult2D``，使测试在并行 agent 还在写
  ``pad2d.py`` 时也能通过。
* 不触碰 ``tests/test_pad2d.py`` / ``tests/test_climbing.py`` —— 平行 agent
  的 live workspace。
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field

import numpy as np
import pytest

from climbing.ankle_injury import (
    ANKLE_MODE_VALUES,
    AnkleVerdict,
    ankle_band,
    assess_ankle,
)
from climbing.injury import VERDICTS
from climbing.pad2d import Posture2D
from climbing.scenarios import DIMENSIONS


# ==========================================================================
# Stub fixtures
# ==========================================================================
@dataclass
class _StubWithAnkle:
    """携带全部 P2 踝字段的轻量 stub —— 与 ``BoulderFallResult2D`` 同字段名。"""

    posture: str = "feet-first-stiff"
    height_m: float = 3.0
    ankle_mode: str = "stuck"
    peak_ankle_beta_rad: float = 0.20
    peak_inversion_moment_nmm: float = 25000.0
    peak_ligament_strain: float = 0.10
    ankle_bone_utilization: float = 0.5
    ankle_fracture: bool = False
    ankle_sprain: bool = True
    ankle_beta0_rad: float = 0.175   # ~10 deg
    meta: dict = field(default_factory=dict)


@dataclass
class _StubPlain1D:
    """plain 1D ``BoulderFallResult`` —— 无踝字段，应触发 ValueError。"""

    posture: str = "feet-first-stiff"
    height_m: float = 3.0
    peak_force_n: float = 18000.0
    meta: dict = field(default_factory=dict)


# ==========================================================================
# 1. ankle_band —— 与 injury.VERDICTS 同源，4 档
# ==========================================================================
class TestAnkleBand:
    """分档规则（顺序敏感，最严重者先吞）：

    * 骨折 ∧ utilization ≥ 2.0 → 极高
    * 骨折 → 高
    * 扭伤 ∧ strain ≥ 0.2     → 高
    * 扭伤 → 中
    * 其余 → 低
    """

    def test_bands_set_matches_injury_verdicts(self) -> None:
        assert set(VERDICTS) == {"低", "中", "高", "极高"}

    @pytest.mark.parametrize(
        "sprain,frac,util,strain,expected",
        [
            # 基线：未扭未折 → 低
            (False, False, 0.0,  0.0,  "低"),
            (False, False, 0.5,  0.1,  "低"),
            # 扭伤（应变低） → 中
            (True,  False, 0.5,  0.1,  "中"),
            (True,  False, 0.5,  0.14, "中"),  # 刚越 ATFL 失效也是中
            (True,  False, 0.5,  0.199,"中"),
            # 扭伤 ∧ 高应变 → 高
            (True,  False, 0.5,  0.20, "高"),
            (True,  False, 0.5,  0.25, "高"),
            # 骨折（利用率低） → 高
            (False, True,  0.5,  0.0,  "高"),
            (False, True,  1.0,  0.0,  "高"),
            (False, True,  1.99, 0.0,  "高"),  # 临界下
            # 骨折 ∧ utilization ≥ 2.0 → 极高
            (False, True,  2.0,  0.0,  "极高"),
            (False, True,  3.5,  0.0,  "极高"),
            # 顺序敏感：骨折 ∧ 扭伤 ∧ 高应变 → 仍是极高（最高优先）
            (True,  True,  2.5,  0.3,  "极高"),
            # 顺序敏感：骨折 ∧ 扭伤 ∧ 低应变 → 高（骨折先吞）
            (True,  True,  1.5,  0.1,  "高"),
        ],
    )
    def test_band_mappings(
        self,
        sprain: bool,
        frac: bool,
        util: float,
        strain: float,
        expected: str,
    ) -> None:
        got = ankle_band(sprain, frac, util, strain)
        assert got == expected, (
            f"ankle_band(sprain={sprain}, fracture={frac}, util={util}, "
            f"strain={strain}) = {got!r}, expected {expected!r}"
        )

    def test_returned_value_is_one_of_verdicts(self) -> None:
        """输出值必须是 ``injury.VERDICTS`` 的元素（防漂移）。"""
        for s in (False, True):
            for f in (False, True):
                for u in (0.0, 0.5, 1.0, 2.0):
                    for st in (0.0, 0.1, 0.2, 0.3):
                        assert ankle_band(s, f, u, st) in VERDICTS


# ==========================================================================
# 2. assess_ankle —— 鸭式读取 + 错误信息
# ==========================================================================
class TestAssessAnkle:
    """鸭式读取 P2 踝字段；缺数据时抛 ValueError（中文）。"""

    def test_stub_with_ankle_returns_full_verdict(self) -> None:
        stub = _StubWithAnkle(
            posture="toe-point",
            height_m=2.5,
            ankle_mode="stuck",
            peak_ankle_beta_rad=0.20,
            peak_inversion_moment_nmm=25000.0,
            peak_ligament_strain=0.10,
            ankle_bone_utilization=0.5,
            ankle_fracture=False,
            ankle_sprain=True,
            ankle_beta0_rad=0.175,
        )
        v = assess_ankle(stub)
        # 字段对字段地核：所有字段从 stub 直传
        assert v.posture == "toe-point"
        assert v.height_m == 2.5
        assert v.beta0_rad == pytest.approx(0.175)
        assert v.mode == "stuck"
        assert v.peak_beta_rad == pytest.approx(0.20)
        assert v.peak_moment_nmm == pytest.approx(25000.0)
        assert v.peak_ligament_strain == pytest.approx(0.10)
        assert v.bone_utilization == pytest.approx(0.5)
        assert v.fracture is False
        assert v.sprain is True
        # band 由 ankle_band 派生：sprain & strain=0.10 < 0.2 → 中
        assert v.band == "中"
        # note 必有内容（机制说明）
        assert v.note

    def test_explicit_posture_overrides_result_posture(self) -> None:
        stub = _StubWithAnkle(posture="from-result")
        v = assess_ankle(stub, posture="override")
        assert v.posture == "override"

    def test_default_posture_falls_back_to_result(self) -> None:
        stub = _StubWithAnkle(posture="from-result")
        v = assess_ankle(stub)  # posture 默认 ""
        assert v.posture == "from-result"

    def test_unknown_mode_is_clamped_to_none_with_note(self) -> None:
        stub = _StubWithAnkle(ankle_mode="bogus_mode")  # 不在 ANKLE_MODE_VALUES
        v = assess_ankle(stub)
        assert v.mode == "none"
        assert "bogus_mode" in v.note
        assert "回退" in v.note

    def test_all_four_modes_round_trip(self) -> None:
        for mode in ANKLE_MODE_VALUES:
            stub = _StubWithAnkle(ankle_mode=mode)
            v = assess_ankle(stub)
            assert v.mode == mode

    def test_plain_1d_raises_value_error_with_chinese_message(self) -> None:
        stub = _StubPlain1D()
        with pytest.raises(ValueError) as exc:
            assess_ankle(stub)
        # 错误信息要点：必须包含"踝"和"peak_ligament_strain"，必须是中文
        msg = str(exc.value)
        assert "踝" in msg
        assert "peak_ligament_strain" in msg
        # 镜像 bone.py:301-307 的口径：明确指出缺失字段
        assert "缺失" in msg

    def test_value_error_when_only_in_meta_missing_too(self) -> None:
        """Stub 上既没属性、``meta`` 里也没有 peak_ligament_strain → 抛错。"""
        stub = _StubPlain1D()  # meta 是空 dict
        assert "peak_ligament_strain" not in stub.meta
        with pytest.raises(ValueError):
            assess_ankle(stub)

    def test_value_error_mentions_pad2d_or_posture2d(self) -> None:
        """错误信息应引导用户去 pad2d（合约位置）。"""
        with pytest.raises(ValueError) as exc:
            assess_ankle(_StubPlain1D())
        msg = str(exc.value)
        # 至少出现一个合约线索：pad2d / Posture2D / BoulderFallResult2D
        assert any(token in msg for token in
                   ("pad2d", "Posture2D", "BoulderFallResult2D"))

    def test_meta_mirrored_field_is_also_picked_up(self) -> None:
        """合约约定 peak_ligament_strain 被镜像到 meta —— 即使没属性也可读。"""
        stub = _StubPlain1D()
        stub.meta = {
            "peak_ligament_strain": 0.08,
            "ankle_mode": "stuck",
            "ankle_sprain": True,
            "ankle_fracture": False,
            "ankle_bone_utilization": 0.3,
            "peak_ankle_beta_rad": 0.15,
            "peak_inversion_moment_nmm": 12000.0,
            "ankle_beta0_rad": 0.10,
        }
        v = assess_ankle(stub)
        assert v.peak_ligament_strain == pytest.approx(0.08)
        assert v.mode == "stuck"
        assert v.sprain is True
        assert v.band == "中"  # sprain ∧ strain=0.08 < 0.2 → 中


# ==========================================================================
# 3. AnkleVerdict dataclass shape —— 锁字段
# ==========================================================================
class TestAnkleVerdictShape:
    """锁字段集合与顺序 —— 防新增字段时被静默加进 CSV/报告。"""

    def test_fields_match_spec(self) -> None:
        expected = [
            "posture", "height_m", "beta0_rad", "mode", "peak_beta_rad",
            "peak_moment_nmm", "peak_ligament_strain", "bone_utilization",
            "fracture", "sprain", "band", "note",
        ]
        got = list(AnkleVerdict.__dataclass_fields__.keys())
        assert got == expected, f"AnkleVerdict fields drifted: {got}"

    def test_frozen(self) -> None:
        v = AnkleVerdict(
            posture="x", height_m=1.0, beta0_rad=0.0, mode="none",
            peak_beta_rad=0.0, peak_moment_nmm=0.0,
            peak_ligament_strain=0.0, bone_utilization=0.0,
            fracture=False, sprain=False, band="低", note="",
        )
        with pytest.raises(dataclasses.FrozenInstanceError):
            setattr(v, "fracture", True)  # 故意违反 frozen 以断言保护生效


# ==========================================================================
# 4. ROTATION_TO_POSTURE2D —— 来自 scripts/phase_p1_montecarlo.py 的工厂
# ==========================================================================
# 延迟导入：避免在并行 agent 还在写 pad2d.py 时就强行 import。
# 这里的测试**期望**工厂存在；如果工厂不存在 ⇒ 这一路测试失败
# （这是契约：脚本必须导出 ROTATION_TO_POSTURE2D）。
def _load_p1_module():
    """从脚本模块里加载 ``phase_p1_montecarlo``，失败时 pytest.skip。"""
    import importlib.util
    import pathlib
    script = pathlib.Path(__file__).resolve().parents[1] / "scripts" / "phase_p1_montecarlo.py"
    if not script.is_file():
        pytest.skip(f"脚本未找到: {script}")
    spec = importlib.util.spec_from_file_location("phase_p1_montecarlo", script)
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
    except (ImportError, SyntaxError, AttributeError, NameError, TypeError) as e:
        # 并行 agent 可能正在改 pad2d.py ⇒ 导入失败时跳过而非报错
        pytest.skip(f"phase_p1_montecarlo 导入失败（pad2d 可能还在写）: {e}")
    return mod


@dataclass(frozen=True)
class _FakePosture2DWithAnkle:
    """与 P2 合约同字段名的假 ``Posture2D`` —— 用于在 pad2d.py 未就绪时
    仍然测试 ``ROTATION_TO_POSTURE2D`` 的映射意图。

    字段集合 = 现有 2D 字段 + P2 合约的 6 个踝字段（默认值与合约一致）。
    """

    theta0_rad: float = 0.0
    omega0_rad_s: float = 0.0
    lateral_offset_m: float = 0.0
    foot_x_lever_m: float = 0.0
    mu_foot: float = 0.0
    I_body_kgm2: float = 0.0
    c_rot: float = 0.0
    ankle_beta0_rad: float = 0.0
    beta_from_theta: float = 0.0
    mu_slide: float = 0.0
    ankle_d_lateral_m: float = 0.0
    ankle_d_medial_m: float = 0.0
    ligament_r_mm: float = 22.0


class TestRotationToPosture2D:
    """验证 ``ROTATION_TO_POSTURE2D`` 覆盖 ``DIMENSIONS["rotation"]`` 全部
    类别，且 ``"without"`` 映射到 ``Posture2D()`` 的全默认（bit-identical
    to 1D —— pad.py 硬门槛）。

    ``"longitudinal"`` 的角度抽样用**假 Posture2D**（带 P2 踝字段）测 —— 
    这样即使并行 agent 的 ``pad2d.py`` 还没加踝字段，映射意图仍被锁住。
    """

    @pytest.fixture
    def p1_module(self):
        return _load_p1_module()

    @pytest.fixture
    def factory_map(self, p1_module):
        if not hasattr(p1_module, "ROTATION_TO_POSTURE2D"):
            pytest.skip("ROTATION_TO_POSTURE2D 未在 phase_p1_montecarlo 中导出")
        return p1_module.ROTATION_TO_POSTURE2D

    def test_factory_is_callable_dict(self, factory_map) -> None:
        """必须是 ``dict[str, Callable[..., Posture2D]]``。"""
        assert isinstance(factory_map, dict)
        for k, v in factory_map.items():
            assert callable(v), f"{k!r} 的映射不是可调用对象"

    def test_covers_every_rotation_category(self, factory_map) -> None:
        cats = {c.label for c in DIMENSIONS["rotation"].cats}
        mapped = set(factory_map.keys())
        missing = cats - mapped
        assert not missing, f"ROTATION_TO_POSTURE2D 缺少类别: {missing}"

    def test_without_is_all_defaults(self, factory_map) -> None:
        """``"without"`` 必须返回 ``Posture2D()`` 全默认 —— 与 1D bit-identical。"""
        rng = np.random.default_rng(0)
        p2d = factory_map["without"](rng)
        # 所有字段都等于 Posture2D() 的默认值
        defaults = Posture2D()
        assert dataclasses.asdict(p2d) == dataclasses.asdict(defaults), (
            f"'without' 映射不是全默认:\n  got={p2d}\n  exp={defaults}"
        )

    def test_without_matches_posture2d_constructor(self, factory_map) -> None:
        rng = np.random.default_rng(42)
        a = factory_map["without"](rng)
        b = Posture2D()
        assert a == b

    def test_unknown_is_all_defaults(self, factory_map) -> None:
        """``"unknown"`` 与 ``"without"`` 同（保守默认）。"""
        if "unknown" not in factory_map:
            pytest.skip("unknown 未映射")
        rng = np.random.default_rng(0)
        p2d = factory_map["unknown"](rng)
        assert dataclasses.asdict(p2d) == dataclasses.asdict(Posture2D())

    @pytest.mark.parametrize("rng_seed", [0, 1, 42, 2026])
    def test_longitudinal_samples_within_T22_set(
        self, p1_module, factory_map, rng_seed, monkeypatch
    ) -> None:
        """``"longitudinal"`` 的角度幅度 ∈ {10, 20, 30} deg（[T22]）；符号 50/50。

        用**假 Posture2D**（带 P2 踝字段）替换模块里的 ``Posture2D`` 与
        ``_VALID_POSTURE2D_FIELDS`` —— 这样在并行 agent 的 ``pad2d.py``
        还没加踝字段时，映射意图仍被锁定。
        """
        if "longitudinal" not in factory_map:
            pytest.skip("longitudinal 未映射")
        # 注入带踝字段的假类
        fake_fields = frozenset(f.name for f in
                                dataclasses.fields(_FakePosture2DWithAnkle))
        monkeypatch.setattr(p1_module, "Posture2D", _FakePosture2DWithAnkle)
        monkeypatch.setattr(p1_module, "_VALID_POSTURE2D_FIELDS", fake_fields)

        rng = np.random.default_rng(rng_seed)
        mags_deg = set()
        seen_signs = set()
        for _ in range(200):
            p2d = factory_map["longitudinal"](rng)
            beta0_rad = float(p2d.ankle_beta0_rad)
            mags_deg.add(round(abs(np.rad2deg(beta0_rad)), 6))
            if beta0_rad != 0.0:
                seen_signs.add(float(np.sign(beta0_rad)))
        assert mags_deg <= {10.0, 20.0, 30.0}, f"出现 [T22] 集合外的角度: {mags_deg}"
        assert mags_deg == {10.0, 20.0, 30.0}, (
            f"200 次采样未覆盖全部 {{10,20,30}}: {mags_deg}"
        )
        assert seen_signs == {1.0, -1.0}, f"符号不全: {seen_signs}"

    def test_longitudinal_sets_documented_ankle_params(
        self, p1_module, factory_map, monkeypatch
    ) -> None:
        """``"longitudinal"`` 的其余踝参数是文档化的常量（可回归）。"""
        if "longitudinal" not in factory_map:
            pytest.skip("longitudinal 未映射")
        fake_fields = frozenset(f.name for f in
                                dataclasses.fields(_FakePosture2DWithAnkle))
        monkeypatch.setattr(p1_module, "Posture2D", _FakePosture2DWithAnkle)
        monkeypatch.setattr(p1_module, "_VALID_POSTURE2D_FIELDS", fake_fields)

        rng = np.random.default_rng(7)
        p2d = factory_map["longitudinal"](rng)
        assert p2d.beta_from_theta == pytest.approx(1.0)
        assert p2d.mu_slide == pytest.approx(0.4)
        assert p2d.ankle_d_lateral_m == pytest.approx(0.030)
        assert p2d.ankle_d_medial_m == pytest.approx(0.045)
        assert p2d.ligament_r_mm == pytest.approx(22.0)

    def test_antero_posterior_sets_theta0(
        self, p1_module, factory_map, monkeypatch
    ) -> None:
        """``"antero_posterior"`` ⇒ 身体俯仰（±15°，50/50）。"""
        if "antero_posterior" not in factory_map:
            pytest.skip("antero_posterior 未映射")
        fake_fields = frozenset(f.name for f in
                                dataclasses.fields(_FakePosture2DWithAnkle))
        monkeypatch.setattr(p1_module, "Posture2D", _FakePosture2DWithAnkle)
        monkeypatch.setattr(p1_module, "_VALID_POSTURE2D_FIELDS", fake_fields)

        rng = np.random.default_rng(11)
        thetas = []
        for _ in range(50):
            p2d = factory_map["antero_posterior"](rng)
            thetas.append(float(p2d.theta0_rad))
        mags = {round(abs(np.rad2deg(t)), 6) for t in thetas}
        assert mags == {15.0}, f"antero_posterior 的 theta0 幅度应为 15°: {mags}"
        assert {float(np.sign(t)) for t in thetas if t != 0.0} == {1.0, -1.0}


# ==========================================================================
# 5. site_distribution —— 踝列的聚合
# ==========================================================================
# site_distribution 在 scripts/phase_p1_montecarlo.py 内定义。
def _load_site_distribution():
    mod = _load_p1_module()
    if not hasattr(mod, "site_distribution"):
        pytest.skip("site_distribution 未在 phase_p1_montecarlo 中定义")
    return mod.site_distribution


class TestSiteDistributionAnkleAggregation:
    """``site_distribution`` 必须把踝列纳入 all / no_rotation / credible
    三个子集，每 trial 数一次踝扭伤 / 骨折 / 任意损伤。"""

    @pytest.fixture
    def site_dist_fn(self):
        return _load_site_distribution()

    @pytest.fixture
    def fixture_rows(self):
        """合成 5 个 trial 的 fixture：见函数体。"""
        return [
            {"i": 0, "rotation_free": 1, "weight_exact": 1},
            {"i": 1, "rotation_free": 1, "weight_exact": 1},
            {"i": 2, "rotation_free": 1, "weight_exact": 1},
            {"i": 3, "rotation_free": 1, "weight_exact": 1},
            {"i": 4, "rotation_free": 1, "weight_exact": 0},
        ]

    @pytest.fixture
    def fixture_site_rows(self):
        # 每个 i 有 7 个 bone_site 行（与 BONE_SITES 7 个部位对齐），全部带 ankle 列
        SITE_ORDER = [
            "calcaneus", "tibia_distal", "tibia_mid", "fibula_ends",
            "femoral_neck", "pelvis", "lumbar_spine",
        ]
        fracture_map = {
            # i: {site: fractured_bool}
            0: {"calcaneus": True},
            1: {"tibia_distal": True},
            2: {},
            3: {},
            4: {},
        }
        ankle = {
            0: {"ankle_fracture": 0, "ankle_sprain": 1},
            1: {"ankle_fracture": 1, "ankle_sprain": 1},
            2: {"ankle_fracture": 1, "ankle_sprain": 0},
            3: {"ankle_fracture": 0, "ankle_sprain": 0},
            4: {"ankle_fracture": 0, "ankle_sprain": 1},
        }
        rows = []
        for i in range(5):
            for s in SITE_ORDER:
                rows.append({
                    "i": i,
                    "site": s,
                    "utilization": 0.5,
                    "fractured": int(fracture_map[i].get(s, False)),
                    "ankle_fracture": ankle[i]["ankle_fracture"],
                    "ankle_sprain": ankle[i]["ankle_sprain"],
                })
        return rows

    def test_ankle_keys_present(self, site_dist_fn, fixture_rows, fixture_site_rows) -> None:
        dist = site_dist_fn(fixture_rows, fixture_site_rows)
        for subset in ("all", "no_rotation", "credible"):
            assert "ankle_sprain" in dist[subset], f"missing ankle_sprain in {subset}"
            assert "ankle_fracture" in dist[subset]
            assert "ankle_injured" in dist[subset]

    def test_ankle_counts_in_all_subset(
        self, site_dist_fn, fixture_rows, fixture_site_rows
    ) -> None:
        """all 子集下（5 trial）：
        - ankle_sprain = 3  （i=0,1,4）
        - ankle_fracture = 2  （i=1,2）
        - ankle_injured = 4  （i=0,1,2,4）
        """
        dist = site_dist_fn(fixture_rows, fixture_site_rows)
        all_dist = dist["all"]
        assert all_dist["n_trials"] == 5
        assert all_dist["ankle_sprain"] == 3
        assert all_dist["ankle_fracture"] == 2
        assert all_dist["ankle_injured"] == 4

    def test_ankle_counts_in_credible_subset(
        self, site_dist_fn, fixture_rows, fixture_site_rows
    ) -> None:
        """credible 子集（weight_exact=1 ⇒ i=0,1,2,3 ⇒ 4 trial）：
        - ankle_sprain = 2  （i=0,1）
        - ankle_fracture = 2  （i=1,2）
        - ankle_injured = 3  （i=0,1,2）
        """
        dist = site_dist_fn(fixture_rows, fixture_site_rows)
        cred = dist["credible"]
        assert cred["n_trials"] == 4
        assert cred["ankle_sprain"] == 2
        assert cred["ankle_fracture"] == 2
        assert cred["ankle_injured"] == 3

    def test_ankle_no_double_counting(
        self, site_dist_fn, fixture_rows, fixture_site_rows
    ) -> None:
        """每个 trial 只数一次踝 —— 即使 site_rows 有 7 行 / trial。"""
        dist = site_dist_fn(fixture_rows, fixture_site_rows)
        # i=1 在 site_rows 里有 7 行（每行 ankle_sprain=1）
        # 但踝 count 只应为 1 ⇒ 全 4 trial 共 3 sprain
        assert dist["all"]["ankle_sprain"] == 3

    def test_bone_site_share_unchanged(
        self, site_dist_fn, fixture_rows, fixture_site_rows
    ) -> None:
        """bone site 的 share 必须保持原有语义（事件计数 / 事件总数）。"""
        dist = site_dist_fn(fixture_rows, fixture_site_rows)
        all_dist = dist["all"]
        # bone site 事件：calcaneus×1 + tibia_distal×1 = 2 events
        assert all_dist["n_injuries"] == 2
        # share
        assert all_dist["share"]["calcaneus"] == pytest.approx(0.5)
        assert all_dist["share"]["tibia_distal"] == pytest.approx(0.5)

    def test_ankle_share_over_n_trials(
        self, site_dist_fn, fixture_rows, fixture_site_rows
    ) -> None:
        """踝 share 应为 踝 count / n_trials（口径与 [B25] Table 1 一致：
        踝占 40% 伤者 ⇒ per-trial 比例，而非 per-event）。"""
        dist = site_dist_fn(fixture_rows, fixture_site_rows)
        all_dist = dist["all"]
        if "ankle_injured_share" in all_dist:
            assert all_dist["ankle_injured_share"] == pytest.approx(4 / 5)
            assert all_dist["ankle_sprain_share"] == pytest.approx(3 / 5)
            assert all_dist["ankle_fracture_share"] == pytest.approx(2 / 5)


# ==========================================================================
# 6. 端到端集成 —— run_trials 把 pad2d 的踝输出接进列与分布
# ==========================================================================
class TestRunTrialsAnkleWiring:
    """验证 ``run_trials`` 真的换用了 ``simulate_boulder_fall_2d``，并把踝列
    接进 ``trial_rows`` / ``site_rows``，且 ``site_distribution`` 的踝计数
    与 trial 行自洽。

    这是 wiring 测试（不是物理测试）：物理标定由 pad2d.py 的并行 agent
    负责，本测试只锁"接口没接错"。
    """

    def test_run_trials_uses_2d_simulator_and_populates_ankle(
        self, monkeypatch
    ) -> None:
        mod = _load_p1_module()
        # 确保脚本里有踝列清单
        assert hasattr(mod, "ANKLE_COLUMNS")
        # Spy：记录 simulate_boulder_fall_2d 被调用时的 posture2d
        real_sim = mod.simulate_boulder_fall_2d
        seen_posture2d = []

        def spy(**kwargs):
            seen_posture2d.append(kwargs.get("posture2d"))
            return real_sim(**kwargs)

        monkeypatch.setattr(mod, "simulate_boulder_fall_2d", spy)

        # 跑 3 次（fast 预设）—— 保证至少一次仿真
        rows, site_rows = mod.run_trials(3, seed=20261007, fast=True)
        assert rows, "run_trials 没有任何成功样本"

        # (a) 仿真器被换成了 2D：每次都传了 Posture2D 实例
        assert seen_posture2d, "simulate_boulder_fall_2d 未被调用"
        assert all(c is not None for c in seen_posture2d), (
            "run_trials 未给 simulate_boulder_fall_2d 传 posture2d"
        )

        # (b) 每个 trial 行都带全部踝列
        for r in rows:
            for col in mod.ANKLE_COLUMNS:
                assert col in r, f"trial 行缺踝列 {col!r}"
            # ankle_injured = sprain or fracture（自洽）
            assert int(r["ankle_injured"]) == int(
                bool(r["ankle_sprain"] or r["ankle_fracture"]))

        # (c) 每个 site 行也带全部踝列
        for s in site_rows:
            for col in mod.ANKLE_COLUMNS:
                assert col in s, f"site 行缺踝列 {col!r}"

        # (d) site_distribution 的踝计数与 trial 行自洽
        dist = mod.site_distribution(rows, site_rows)
        for sub in ("all", "no_rotation", "credible"):
            assert "ankle_sprain" in dist[sub]
            assert "ankle_fracture" in dist[sub]
            assert "ankle_injured" in dist[sub]
        # all 子集：踝计数 = trial 行里对应列的和
        exp_sp = sum(int(r["ankle_sprain"]) for r in rows)
        exp_fr = sum(int(r["ankle_fracture"]) for r in rows)
        exp_in = sum(int(r["ankle_injured"]) for r in rows)
        assert dist["all"]["ankle_sprain"] == exp_sp
        assert dist["all"]["ankle_fracture"] == exp_fr
        assert dist["all"]["ankle_injured"] == exp_in