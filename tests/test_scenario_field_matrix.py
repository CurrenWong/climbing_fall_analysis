"""回归：``scripts/opensim_fe/scenario_field_matrix.py``（快测试，不跑 OpenSim）。

覆盖：
* 场景网格几何（3×3×4=36）与“只用 opt-in 旋钮”约束；
* 现场 7 运动学（Beurienne 2025）的 top-3 与损伤倾向；
* 弱子区域强度与任务给定值一致；
* 载荷分配倍率（胫/腓并联、中轴 mass-above）与 ``risk_1d_loadshare`` 一致；
* 纯轴向 σ / 弱子区域旗标的合成数值；
* 临床分桶覆盖 9 部位；空 analyses 下临床/R5 聚合不崩。

OpenSim 正动力学很慢，本测试**不**跑载荷链（那属于脚本主流程，缓存于
``temp/opensim_fe/scenario_field_matrix/``）。
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parents[1]
_MOD_PATH = _ROOT / "scripts" / "opensim_fe" / "scenario_field_matrix.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("scenario_field_matrix", _MOD_PATH)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


m = _load_module()

ALLOWED_SCENARIO_KEYS = {
    "tilt_deg", "split", "foot_forces", "vx0_ms", "vz0_ms",
    "run", "representable", "height_m",
}


# --------------------------------------------------------------------------
# 网格 / 现场运动学
# --------------------------------------------------------------------------
def test_grid_axes_are_3x3x4():
    assert m.GRID_AXES["start"] == ("vertical", "leaning-back", "leaning-forward")
    assert m.GRID_AXES["rotation"] == ("none", "longitudinal", "other")
    assert m.GRID_AXES["landing"] == (
        "upright-on-feet", "feet-leaning", "head-first", "buttocks")


def test_full_grid_is_36_and_covers_all_combos():
    g = m.enumerate_full_grid()
    assert len(g) == 36
    combos = {(r["start"], r["rotation"], r["landing"]) for r in g}
    assert len(combos) == 36
    assert {r["start"] for r in g} == set(m.GRID_AXES["start"])
    # 只有足先落地（双足/侧倾）可由足端载荷链表示
    for r in g:
        assert r["representable"] == (r["landing"] in ("upright-on-feet", "feet-leaning"))


def test_full_grid_only_uses_allowed_knobs_and_varies_them():
    g = m.enumerate_full_grid()
    for r in g:
        assert r["tilt_deg"] in (0.0, 25.0, -25.0)
        assert r["vx0_ms"] in (0.0, 1.2)
        assert r["vz0_ms"] == 0.0
        if r["representable"]:
            if r["rotation"] == "longitudinal":
                assert tuple(r["split"]) == (0.3, 0.7)
            elif r["rotation"] == "none":
                assert tuple(r["split"]) == (0.5, 0.5)
            else:
                assert tuple(r["split"]) == (0.2, 0.8)
            assert r["foot_forces"] is None
        else:
            # 接触不在足端 → 无 split/foot_forces
            assert r["split"] is None and r["foot_forces"] is None
    # 观测到的 7 个组合都被标记
    assert sum(1 for r in g if r["observed"]) == 7


def test_field_kinematics_top3_and_tendency():
    fk = m.FIELD_KINEMATICS
    assert len(fk) == 7
    top = [k for k in fk if k["freq_pct"] is not None]
    assert [k["kind"] for k in top] == ["K1", "K2", "K3"]
    assert [k["freq_pct"] for k in top] == [17.0, 14.0, 11.0]
    # K1/K2 → 下肢；K3（后仰+纵轴旋转+侧倾）→ 上肢（论文 §4）
    assert top[0]["injury_tendency"] == "下肢"
    assert top[1]["injury_tendency"] == "下肢"
    assert top[2]["injury_tendency"] == "上肢"


def test_scenarios_top3_plus_upperlimb_are_run():
    run = [s for s in m.SCENARIOS if s["run"]]
    kinds = {s["field_kind"] for s in run}
    assert {"K1", "K2", "K3"}.issubset(kinds), kinds
    assert len(run) >= 4, "任务要求 top-3 + 至少一个 leaning/upper-limb 场景"


def test_scenarios_use_only_optin_knobs():
    for s in m.SCENARIOS:
        assert "posture" not in s, "不得依赖 posture 旋钮"
        extra = set(s) - ALLOWED_SCENARIO_KEYS - {
            "id", "start", "rotation", "landing", "field_kind", "note"}
        assert not extra, f"{s['id']} 含不允许的字段 {extra}"
        if s["representable"] in (True, "partial"):
            assert (s["split"] is None) != (s["foot_forces"] is None), \
                "split 与 foot_forces 恰有其一（ground_reaction 语义）"
        else:
            assert s["split"] is None and s["foot_forces"] is None


def test_headfirst_and_buttocks_are_proxy_not_representable():
    by_kind = {s["field_kind"]: s for s in m.SCENARIOS}
    assert by_kind["K5"]["representable"] is False   # head-first
    assert by_kind["K6"]["representable"] is False   # buttocks
    assert by_kind["K5"]["run"] is False
    assert by_kind["K6"]["run"] is False


# --------------------------------------------------------------------------
# 强度 / 载荷分配
# --------------------------------------------------------------------------
def test_weak_subregion_strengths_match_task():
    assert m._sc("calcaneus") == 150.0
    assert m._sc("tibia_ends") == 70.0
    assert m._sc("fibula_ends") == 70.0
    assert m._sc("femoral_neck") == 80.0
    assert m._sc("pelvis") == 180.0
    assert m._sc("spine") == 150.0
    assert m._sc("skull") == 160.0


def test_region_table_matches_expected_9_parts():
    assert m.REGION_NAMES == (
        "足部", "胫骨", "腓骨", "股骨", "骨盆", "腰椎", "胸椎", "颈椎", "颅骨")
    bones = {r[1] for r in m.REGIONS}
    assert bones == {"calcaneus_r", "tibia_r", "fibula_r", "femur_r",
                     "R_HIPBONE", "L3", "T6", "C5", "parietal_r"}


def test_load_share_multipliers_from_repo_artifacts():
    if not (m.ROUTE1_JSON.is_file() and m.LOADSHARE_JSON.is_file()):
        pytest.skip("既有产物缺失（fracture_matrix 前置未生成）")
    st = m.load_static_inputs()
    mult = st["multipliers"]
    # 胫/腓并联分流之和 = 1
    assert abs(mult["tibia_r"] + mult["fibula_r"] - 1.0) < 1e-9
    # 单骨倍率
    for b in ("calcaneus_r", "femur_r", "R_HIPBONE"):
        assert mult[b] == 1.0
    # 中轴骨 mass-above 比严格小于 1
    for b in ("L3", "T6", "C5", "parietal_r"):
        assert 0.0 < mult[b] < 1.0


# --------------------------------------------------------------------------
# 纯轴向 σ / 旗标（合成数值）
# --------------------------------------------------------------------------
def _synthetic_static(a=100.0):
    mult = {b: 1.0 for _, b, _, _, _ in m.REGIONS}
    mult.update({"tibia_r": 0.9, "fibula_r": 0.1,
                 "L3": 0.4, "T6": 0.3, "C5": 0.1, "parietal_r": 0.095})
    return {
        "A_section_mm2": {b: a for _, b, _, _, _ in m.REGIONS},
        "sigma_c_stored_mpa": {b: 100.0 for _, b, _, _, _ in m.REGIONS},
        "multipliers": mult,
    }


def test_region_sigma_formula_and_weak_flag():
    static = _synthetic_static(a=100.0)
    loads = {j: [0.0] * 4 for j in m.JOINTS}
    loads["ankle_r"] = [100_000.0] * 4      # 胫 0.9→900 MPa；腓 0.1→100 MPa
    loads["subtalar_r"] = [15_000.0] * 4    # 跟骨 1.0→150 MPa
    sig = m.region_sigma_by_height(loads, static)
    assert sig["fibula_r"][0] == pytest.approx(0.1 * 100_000.0 / 100.0)   # 100 MPa
    assert sig["tibia_r"][0] == pytest.approx(0.9 * 100_000.0 / 100.0)    # 900 MPa
    assert sig["calcaneus_r"][0] == pytest.approx(150.0)
    flags = m.region_flagged_by_height(sig)
    assert flags["fibula_r"][0] == 1        # 100 ≥ 70
    assert flags["calcaneus_r"][0] == 1     # 150 ≥ 150（边界含等号）
    assert flags["femur_r"][0] == 0         # 0 MPa < 80
    assert m.first_flagged_height([2.0, 3.0, 3.5, 4.5], [0, 1, 1, 1]) == 3.0
    assert m.first_flagged_height([2.0, 3.0], [0, 0]) is None


def test_analyze_scenario_loads_returns_rows_and_flags():
    static = _synthetic_static(a=50.0)
    loads = {j: [0.0] * 4 for j in m.JOINTS}
    loads["ankle_r"] = [10_000.0, 20_000.0, 30_000.0, 40_000.0]
    a = m.analyze_scenario_loads(loads, static, [2.0, 3.0, 3.5, 4.5])
    assert len(a["rows"]) == 9
    row = next(r for r in a["rows"] if r["region_cn"] == "腓骨")
    # σ = 0.1*F/50 → 20,40,60,80 MPa；阈值 70 → 旗标高度 4.5
    assert row["sigma_at_heights_mpa"] == pytest.approx([20.0, 40.0, 60.0, 80.0])
    assert row["flagged_heights"] == [4.5]
    assert "腓骨" in a["flagged_regions"]
    assert "足部" not in a["flagged_regions"]  # subtalar 负荷为 0


# --------------------------------------------------------------------------
# 临床分桶 / 聚合
# --------------------------------------------------------------------------
def test_clinical_bucket_covers_all_regions():
    for cn in m.REGION_NAMES:
        assert m.bucket_of(cn) != "其他", cn
    assert m.bucket_of("足部") == "下肢(合计)"
    assert m.bucket_of("腰椎") == "脊柱"
    assert m.bucket_of("颅骨") == "颅骨"


def test_clinical_and_r5_aggregation_empty_analyses():
    cc = m.clinical_comparison({}, {})
    assert cc["n_run_scenarios"] == 0
    assert cc["lower_limb_share_modeled"] == 0.0
    r5 = m.r5_lower_falls({}, {})
    assert r5["trend"] is None
    assert len(r5["by_height"]) == len(m.HEIGHTS)


def test_clinical_targets_have_paper_numbers():
    t = m.CLINICAL_TARGETS
    assert t["lower_limb_pct_heck"] == 61.1
    assert t["ankle_pct_heck"] == 36.7
    assert t["knee_pct_heck"] == 16.8
    assert t["spine_pct_heck"] == 7.2
    assert t["ankle_sprain_rank"] == 1
    assert t["r5_low_fall_lower_limb_pct"] == 79.1
    assert t["r5_high_fall_lower_limb_pct"] == 58.3
