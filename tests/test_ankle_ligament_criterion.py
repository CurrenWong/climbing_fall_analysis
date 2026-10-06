"""R6 —— 踝韧带扭伤（sprain）判据回归测试。

覆盖新模块 ``climbing.coupling.ankle_ligament`` 与脚本产物
``results/opensim_fe/ankle_ligament_criterion.json``。

**硬门槛（默认行为不变）**：
* 新模块在 ``M_inv=0``（β=0）时全部 risk=0、``sprain=False``；
* 既有 ``ankle_supination`` 行为不变（``supination_load`` 默认仍为零）；
* 脚本的 ``default_regression`` 必须 PASS 且相对误差 ≤ 1e-12。

快测试只读既有 JSON 产物（毫秒级），不重跑 OpenSim / FEBio。
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from climbing.coupling.ankle_ligament import (
    ATFL_FAILURE_STRAIN,
    ATFL_INVERSION_FAIL_ANGLE_2KN_DEG,
    ATFL_INVERSION_FAIL_ANGLE_DEG,
    ATFL_INVERSION_FAIL_MOMENT_2KN_PRELOAD_NM,
    ATFL_INVERSION_FAIL_MOMENT_NO_PRELOAD_HIGH_NM,
    ATFL_INVERSION_FAIL_MOMENT_NO_PRELOAD_NM,
    ATFL_PEAK_FORCE_STANCE_N,
    ATFL_RESTING_LENGTH_MM,
    ATFL_STIFFNESS_N_PER_MM,
    ATFL_ULTIMATE_LOAD_N,
    CFL_FAILURE_STRAIN,
    CFL_ULTIMATE_LOAD_N,
    DEFAULT_LIGAMENT_LEVER_ARM_MM,
    R_LIG_SENSITIVITY_MM,
    TOCHIGI_ATFL_NORMAL_STRAIN_MAX,
    TOCHIGI_CFL_NORMAL_STRAIN_MAX,
    TOCHIGI_TAUT_STRAIN,
    ankle_ligament_criterion,
    atfl_sprain_risk,
    cfl_sprain_risk,
    critical_inversion_angle_deg as ligament_critical_inversion_angle_deg,
    critical_inversion_moment_nmm,
    ligament_force_n,
    ligament_strain,
    preload_capacity_interp,
    sprain_risk,
)
from climbing.coupling.ankle_supination import supination_load

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results" / "opensim_fe"
ARTIFACT = RES / "ankle_ligament_criterion.json"
S5_JSON = RES / "nonvertical_s5.json"


# --------------------------------------------------------------------------
# 1. 参数常量（measured，文献）
# --------------------------------------------------------------------------
def test_ligament_ultimate_loads() -> None:
    assert ATFL_ULTIMATE_LOAD_N == pytest.approx(200.0)
    assert CFL_ULTIMATE_LOAD_N == pytest.approx(340.0)


def test_ligament_failure_strains() -> None:
    assert ATFL_FAILURE_STRAIN == pytest.approx(0.14)
    assert CFL_FAILURE_STRAIN == pytest.approx(0.13)


def test_joint_moment_thresholds_ordering() -> None:
    # 无预载下限 < 无预载上限（Funk < Parenteau）< 2 kN 预载
    assert (
        ATFL_INVERSION_FAIL_MOMENT_NO_PRELOAD_NM
        < ATFL_INVERSION_FAIL_MOMENT_NO_PRELOAD_HIGH_NM
        < ATFL_INVERSION_FAIL_MOMENT_2KN_PRELOAD_NM
    )


def test_joint_angle_thresholds_ordering() -> None:
    assert 33.0 <= ATFL_INVERSION_FAIL_ANGLE_DEG <= 40.0
    assert 33.0 <= ATFL_INVERSION_FAIL_ANGLE_2KN_DEG <= 42.0


def test_tochigi_reference_band() -> None:
    assert TOCHIGI_ATFL_NORMAL_STRAIN_MAX == pytest.approx(0.062)
    assert TOCHIGI_CFL_NORMAL_STRAIN_MAX == pytest.approx(0.021)
    assert TOCHIGI_TAUT_STRAIN == pytest.approx(0.010)
    # 正常带上限必须低于失效阈值（参考带 ≠ 失效阈值）
    assert TOCHIGI_ATFL_NORMAL_STRAIN_MAX < ATFL_FAILURE_STRAIN
    assert TOCHIGI_CFL_NORMAL_STRAIN_MAX < CFL_FAILURE_STRAIN


def test_bahr_stance_peak_is_half_ultimate() -> None:
    # 生理 sanity：步态峰值 ≈ 极限载荷一半（安全系数 ≈ 2）
    assert ATFL_PEAK_FORCE_STANCE_N == pytest.approx(0.5 * ATFL_ULTIMATE_LOAD_N, rel=0.3)


# --------------------------------------------------------------------------
# 2. 关节级判据（主口径）
# --------------------------------------------------------------------------
def test_atfl_default_zero_moment() -> None:
    c = atfl_sprain_risk(0.0)
    assert c.ligament == "ATFL"
    assert c.risk == 0.0
    assert c.sprain is False
    assert c.preload_model == "conservative"
    assert c.M_fail_nm == pytest.approx(ATFL_INVERSION_FAIL_MOMENT_NO_PRELOAD_NM)


def test_cfl_default_zero_moment() -> None:
    c = cfl_sprain_risk(0.0)
    assert c.ligament == "CFL"
    assert c.risk == 0.0
    assert c.sprain is False


def test_atfl_risk_definition() -> None:
    # risk = M_inv / M_fail；M_fail = 21 N·m = 21000 N·mm
    m_inv = 21000.0
    c = atfl_sprain_risk(m_inv)
    assert c.risk == pytest.approx(1.0)
    assert c.sprain is True
    assert c.M_inv_nm == pytest.approx(21.0)
    assert c.M_fail_nmm == pytest.approx(21000.0)


def test_atfl_risk_is_abs_and_linear() -> None:
    assert atfl_sprain_risk(-21000.0).risk == pytest.approx(1.0)
    assert atfl_sprain_risk(42000.0).risk == pytest.approx(2.0, rel=1e-12)


def test_atfl_conservative_vs_elevated() -> None:
    m_inv = 60000.0  # 60 N·m
    cons = atfl_sprain_risk(m_inv, preload_model="conservative")
    elev = atfl_sprain_risk(m_inv, preload_model="elevated")
    assert cons.risk > elev.risk
    assert cons.sprain is True
    assert elev.sprain is False  # 60 < 77
    assert elev.M_fail_nm == pytest.approx(ATFL_INVERSION_FAIL_MOMENT_2KN_PRELOAD_NM)


def test_atfl_explicit_m_fail_overrides_model() -> None:
    c = atfl_sprain_risk(30000.0, m_fail_nm=30.0)
    assert c.risk == pytest.approx(1.0)
    assert c.M_fail_nm == pytest.approx(30.0)


def test_sprain_risk_bad_preload_model_raises() -> None:
    with pytest.raises(ValueError):
        atfl_sprain_risk(1000.0, preload_model="bogus")  # type: ignore[arg-type]


# --------------------------------------------------------------------------
# 3. 预载承载力插值
# --------------------------------------------------------------------------
def test_preload_interp_anchors() -> None:
    assert preload_capacity_interp(0.0) == pytest.approx(
        ATFL_INVERSION_FAIL_MOMENT_NO_PRELOAD_NM)
    assert preload_capacity_interp(2000.0) == pytest.approx(
        ATFL_INVERSION_FAIL_MOMENT_2KN_PRELOAD_NM)


def test_preload_interp_midpoint_is_linear() -> None:
    lo = ATFL_INVERSION_FAIL_MOMENT_NO_PRELOAD_NM
    hi = ATFL_INVERSION_FAIL_MOMENT_2KN_PRELOAD_NM
    assert preload_capacity_interp(1000.0) == pytest.approx(0.5 * (lo + hi))


def test_preload_interp_increases_with_preload() -> None:
    ps = [0.0, 500.0, 1000.0, 1500.0, 2000.0]
    ms = [preload_capacity_interp(p) for p in ps]
    assert all(b >= a for a, b in zip(ms, ms[1:]))


def test_preload_interp_clamps_above_2kn() -> None:
    assert preload_capacity_interp(21000.0, extrapolation="clamp") == pytest.approx(
        ATFL_INVERSION_FAIL_MOMENT_2KN_PRELOAD_NM)
    # linear 外推：继续升高
    assert preload_capacity_interp(21000.0, extrapolation="linear") > (
        ATFL_INVERSION_FAIL_MOMENT_2KN_PRELOAD_NM)


def test_preload_interp_negative_clamps_to_no_preload() -> None:
    assert preload_capacity_interp(-100.0) == pytest.approx(
        ATFL_INVERSION_FAIL_MOMENT_NO_PRELOAD_NM)


def test_preload_interp_bad_extrapolation_raises() -> None:
    with pytest.raises(ValueError):
        preload_capacity_interp(5000.0, extrapolation="bogus")  # type: ignore[arg-type]


# --------------------------------------------------------------------------
# 4. 力 / 应变（r_lig assumed）
# --------------------------------------------------------------------------
def test_ligament_force_definition() -> None:
    assert ligament_force_n(4400.0, 22.0) == pytest.approx(200.0)
    assert ligament_force_n(-4400.0, 22.0) == pytest.approx(200.0)
    assert ligament_force_n(4400.0, 20.0) == pytest.approx(220.0)


def test_ligament_force_bad_r_raises() -> None:
    with pytest.raises(ValueError):
        ligament_force_n(1000.0, 0.0)
    with pytest.raises(ValueError):
        ligament_force_n(1000.0, -5.0)


def test_ligament_strain_definition() -> None:
    # ε = F/(k·L)；k=14, L=22 → 308 N ⇒ ε=1.0
    assert ligament_strain(308.0, stiffness_n_per_mm=14.0,
                           resting_length_mm=22.0) == pytest.approx(1.0)
    assert ligament_strain(154.0, stiffness_n_per_mm=14.0,
                           resting_length_mm=22.0) == pytest.approx(0.5)


def test_ligament_strain_zero_guard() -> None:
    assert ligament_strain(100.0, stiffness_n_per_mm=0.0, resting_length_mm=22.0) == 0.0
    assert ligament_strain(100.0, stiffness_n_per_mm=14.0, resting_length_mm=0.0) == 0.0


# --------------------------------------------------------------------------
# 5. 力/应变双判据
# --------------------------------------------------------------------------
def test_sprain_risk_default_zero() -> None:
    out = sprain_risk(0.0)
    assert out["F_lig_n"] == 0.0
    assert out["strain"] == 0.0
    assert out["sprain"] is False
    assert out["risk_force"] == 0.0
    assert out["risk_strain"] == 0.0


def test_sprain_risk_force_threshold_met() -> None:
    # M_inv = 200 N · 22 mm = 4400 N·mm → F_lig = 200 N
    out = sprain_risk(4400.0)
    assert out["F_lig_n"] == pytest.approx(200.0)
    assert out["risk_force"] == pytest.approx(1.0)
    assert out["sprain_force"] is True


def test_sprain_risk_strain_threshold() -> None:
    # ε_fail·k·L = 触发力；× r_lig = 触发 M_inv
    m_thr = (ATFL_FAILURE_STRAIN * ATFL_STIFFNESS_N_PER_MM
             * ATFL_RESTING_LENGTH_MM * DEFAULT_LIGAMENT_LEVER_ARM_MM)
    out = sprain_risk(m_thr)
    assert out["risk_strain"] == pytest.approx(1.0, rel=1e-9)
    assert out["sprain_strain"] is True
    # 略高于阈值仍扭伤（无浮点边界抖动）
    assert sprain_risk(m_thr * 1.01)["sprain"] is True


def test_sprain_risk_reports_self_consistency_diagnostic() -> None:
    out = sprain_risk(4400.0)
    # 力判据与应变判据的触发力不同（文献参数在 toe/失效区不自洽）
    assert out["force_trigger_n"] == pytest.approx(
        ATFL_FAILURE_STRAIN * ATFL_STIFFNESS_N_PER_MM * ATFL_RESTING_LENGTH_MM)
    assert out["force_threshold_n"] == pytest.approx(ATFL_ULTIMATE_LOAD_N)
    assert out["params_consistent"] is False


def test_sprain_risk_cfl_uses_cfl_params() -> None:
    out = sprain_risk(10000.0, ligament="CFL")
    assert out["ligament"] == "CFL"
    assert out["ultimate_load_n"] == pytest.approx(CFL_ULTIMATE_LOAD_N)
    assert out["failure_strain"] == pytest.approx(CFL_FAILURE_STRAIN)


def test_sprain_risk_bad_ligament_raises() -> None:
    with pytest.raises(ValueError):
        sprain_risk(1000.0, ligament="PTFL")  # type: ignore[arg-type]


def test_sprain_risk_monotone_in_moment() -> None:
    ms = [0.0, 500.0, 1000.0, 2000.0, 4000.0]
    fs = [sprain_risk(m)["F_lig_n"] for m in ms]
    assert all(b >= a for a, b in zip(fs, fs[1:]))


# --------------------------------------------------------------------------
# 6. 临界力矩 / 角
# --------------------------------------------------------------------------
def test_critical_inversion_moment_thresholds() -> None:
    assert critical_inversion_moment_nmm(preload_model="conservative") == pytest.approx(21000.0)
    assert critical_inversion_moment_nmm(preload_model="elevated") == pytest.approx(77000.0)


def test_critical_inversion_angle_bracket() -> None:
    fvec = (0.0, 21000.0, 0.0)
    b_cons = ligament_critical_inversion_angle_deg(fvec, preload_model="conservative")
    b_elev = ligament_critical_inversion_angle_deg(fvec, preload_model="elevated")
    assert b_cons is not None and 0.0 < b_cons < 60.0
    assert b_elev is not None and 0.0 < b_elev < 60.0
    # 高承载口径需要更大的内翻角才越阈
    assert b_elev > b_cons


def test_critical_inversion_angle_consistency() -> None:
    """在临界角处，M_inv 应恰好等于 M_fail（risk=1）。"""
    from climbing.coupling.ankle_supination import (
        frontal_lateral_force_n,
        inversion_moment_nmm,
    )

    fvec = (0.0, 21000.0, 0.0)
    b = ligament_critical_inversion_angle_deg(
        fvec, preload_model="conservative", lever_arm_mm=30.0)
    assert b is not None
    fl = abs(frontal_lateral_force_n(fvec, b))
    m = inversion_moment_nmm(fl, 30.0)
    assert m == pytest.approx(21000.0, rel=1e-6)


def test_critical_inversion_angle_none_when_below() -> None:
    # 极小力 → 60° 内不越阈
    assert ligament_critical_inversion_angle_deg(
        (0.0, 1.0, 0.0), preload_model="conservative") is None


def test_critical_inversion_angle_bad_vector_raises() -> None:
    with pytest.raises(ValueError):
        ligament_critical_inversion_angle_deg((0.0, 1.0))


# --------------------------------------------------------------------------
# 7. 单格联合判据
# --------------------------------------------------------------------------
def test_ankle_ligament_criterion_keys() -> None:
    d = ankle_ligament_criterion(50000.0)
    for key in ("m_inv_nmm", "m_inv_nm", "r_lig_mm", "preload_model",
                "atfl_joint", "cfl_joint", "atfl_force_strain",
                "cfl_force_strain", "r_lig_sensitivity"):
        assert key in d
    assert d["m_inv_nm"] == pytest.approx(d["m_inv_nmm"] / 1000.0)
    # r_lig 灵敏度含 20/25
    assert set(d["r_lig_sensitivity"].keys()) <= {20.0, 25.0}


def test_ankle_ligament_criterion_default_zero() -> None:
    d = ankle_ligament_criterion(0.0)
    assert d["atfl_joint"].sprain is False
    assert d["cfl_joint"].sprain is False
    assert d["atfl_force_strain"]["sprain"] is False
    assert d["cfl_force_strain"]["sprain"] is False


# --------------------------------------------------------------------------
# 8. 默认行为不变（硬门槛）
# --------------------------------------------------------------------------
def test_existing_supination_default_is_unchanged() -> None:
    """既有 ankle_supination 行为不变：β=0 ⇒ F_lat=M_inv=0。"""
    sl = supination_load((0.0, 21000.0, 0.0))
    assert sl.is_zero
    assert sl.f_lat_n == 0.0
    assert sl.m_inv_nmm == 0.0


def test_new_module_default_is_inactive() -> None:
    """新模块默认口径（conservative, M_inv=0）⇒ 所有判据零风险。"""
    assert atfl_sprain_risk(0.0).risk == 0.0
    assert cfl_sprain_risk(0.0).risk == 0.0
    assert sprain_risk(0.0)["sprain"] is False
    assert ankle_ligament_criterion(0.0)["atfl_joint"].sprain is False


def test_default_r_lig_constant() -> None:
    assert DEFAULT_LIGAMENT_LEVER_ARM_MM == pytest.approx(22.0)
    assert set(R_LIG_SENSITIVITY_MM) == {20.0, 25.0}


# --------------------------------------------------------------------------
# 9. 产物一致性（若已生成脚本产物）
# --------------------------------------------------------------------------
@pytest.mark.skipif(not ARTIFACT.is_file(), reason="ankle_ligament_criterion.json 未生成")
def test_artifact_default_regression_passes() -> None:
    data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    reg = data["default_regression"]
    assert reg["verdict"] == "PASS"
    assert reg["scenario_max_rel_err"] <= 1e-12
    assert reg["grf_max_rel_err"] <= 1e-9
    for name, ok in reg["checks"].items():
        assert ok is True, name


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="ankle_ligament_criterion.json 未生成")
def test_artifact_grid_and_preload_columns() -> None:
    data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    rows = data["scenario"]["rows"]
    heights = sorted({r["height_m"] for r in rows})
    angles = sorted({r["supination_deg"] for r in rows})
    preloads = sorted({r["preload_model"] for r in rows})
    assert heights == [2.0, 3.0, 3.5, 4.5]
    assert angles == [0.0, 5.0, 10.0, 15.0, 30.0]
    assert preloads == ["conservative", "elevated"]
    # 4 h × 5 β × 2 口径 = 40 格
    assert len(rows) == 40


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="ankle_ligament_criterion.json 未生成")
def test_artifact_zero_beta_is_safe_everywhere() -> None:
    data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    for r in data["scenario"]["rows"]:
        if r["supination_deg"] == 0.0:
            assert r["M_inv_nmm"] == 0.0
            assert r["atfl_joint_sprain"] is False
            assert r["cfl_joint_sprain"] is False
            assert r["atfl_fe_sprain"] is False
            assert r["atfl_force_n"] == 0.0
            assert r["atfl_strain"] == 0.0


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="ankle_ligament_criterion.json 未生成")
def test_artifact_conservative_sprains_at_beta5() -> None:
    data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    cons = [r for r in data["scenario"]["rows"]
            if r["preload_model"] == "conservative" and r["supination_deg"] >= 5.0]
    assert cons, "缺少 β≥5° 的保守口径格点"
    for r in cons:
        assert r["atfl_joint_sprain"] is True, r
        assert r["atfl_fe_sprain"] is True, r


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="ankle_ligament_criterion.json 未生成")
def test_artifact_elevated_matches_threshold() -> None:
    data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    for r in data["scenario"]["rows"]:
        if r["preload_model"] != "elevated":
            continue
        assert r["atfl_joint_sprain"] == (r["M_inv_nm"] >= 77.0), r


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="ankle_ligament_criterion.json 未生成")
def test_artifact_r_lig_sensitivity_all_sprain() -> None:
    data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    sens = data["r_lig_sensitivity"]
    assert set(sens.keys()) == {"20", "22", "25"}
    # β≥5° 在所有 r_lig 档位下都力/应变扭伤
    for r in sens["22"]:
        if r["supination_deg"] >= 5.0:
            for key in ("20", "22", "25"):
                match = [x for x in sens[key]
                         if x["height_m"] == r["height_m"]
                         and x["supination_deg"] == r["supination_deg"]]
                assert match and match[0]["atfl_fe_sprain"] is True


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="ankle_ligament_criterion.json 未生成")
def test_artifact_preload_curve_monotone() -> None:
    data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    curve = data["preload_capacity"]["curve"]
    ms = [c["M_fail_nm"] for c in curve]
    assert all(b >= a for a, b in zip(ms, ms[1:]))


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="ankle_ligament_criterion.json 未生成")
def test_artifact_clinical_comparison() -> None:
    data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    clin = data["clinical_comparison"]
    assert clin["heck2024_ankle_sprain_share_pct"] == pytest.approx(71.3)
    assert clin["heck2024_ankle_fracture_share_pct"] == pytest.approx(27.4)


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="ankle_ligament_criterion.json 未生成")
def test_artifact_matches_r1_f_vert() -> None:
    """脚本用的 F_vert 必须来自 R1 缓存（只读复用）。"""
    if not S5_JSON.is_file():
        pytest.skip("nonvertical_s5.json 未生成")
    s5 = json.loads(S5_JSON.read_text(encoding="utf-8"))
    art = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    cached = {float(m["height_m"]): float(m["grf_peak_right_n"])
              for m in s5["measured_heights"] if not m.get("on_pad", True)}
    for r in art["scenario"]["rows"]:
        assert r["F_vert_n"] == pytest.approx(cached[float(r["height_m"])], rel=1e-12)
