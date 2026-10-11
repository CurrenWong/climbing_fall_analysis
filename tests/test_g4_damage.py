"""tests/test_g4_damage.py -- G4 损伤演化最小路径的接口、不变量与诚实判定锁。

对应设计 ``docs/methods/G4_损伤演化设计.md``（预注册冻结）与实现
``src/climbing/coupling/bone_damage.py`` / ``scripts/opensim_fe/g4_damage_evolution.py`` /
``scripts/opensim_fe/bone_feb.py``（加法式）。**纯增量**：不改既有模块数值。

锁定的不变量
------------
1. 映射不变量：σ_y→∞ 或 E_t=E ⇒ 双线性律**退化为恒等**（冻结线弹性逐格复现）；
   双线性膝点 σ(ε) 与后处理映射 σ(σ_lin) 逐值正确。
2. ``read_materials`` **逐位不变**（加法字段默认 ``None`` 不改变既有构造）；
   ``read_material_damage`` 读卡 ``nums[4]/nums[5]``（非调参）。
3. rank-2 提升回归：损伤律把 PC2 解释方差比从冻结的 ~2.4e-4 抬到 ~8.8e-3
   （数值脱离 rank-1 天花板），但数据驱动 null 门如实为 **FAIL**（PC2 < null）。
4. 规则复用锁：转录的 ``a4_rank2_from_sigma`` 在**线性**轨迹上与冻结的
   ``s4_alignment_bc2.a4_from_matrix`` **逐值相等**（LP/PCA 规则未变）。
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "src", ROOT / "scripts" / "opensim_fe"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from climbing.coupling import bone_damage as bd  # noqa: E402
from scripts.opensim_fe import bone_feb  # noqa: E402
from scripts.opensim_fe import g4_damage_evolution as G  # noqa: E402

MASTER = ROOT / "model" / "AM50_V71_Occupant" / "main_THUMS_AM50_V71.k"
_BC2 = ROOT / "results" / "opensim_fe" / "fracture_matrix_bc2.json"
_H5 = ROOT / "temp" / "opensim_fe" / "fracture_matrix_bc2" / "tibia_r" / "tibia_r.hdf5"
need_data = pytest.mark.skipif(
    not (MASTER.is_file() and _BC2.is_file() and _H5.is_file()),
    reason="THUMS 卡 / 冻结 bc2 JSON / hdf5 未随仓库分发或未生成")


# --------------------------------------------------------------------------
# 1. 映射不变量
# --------------------------------------------------------------------------
def test_damage_map_identity_when_no_yield() -> None:
    """σ_y=+∞ ⇒ 永不屈服 ⇒ 恒等（冻结严格线弹性逐格复现）。"""
    x = np.array([1.0, 5.0, 40.0, 300.0])
    out = bd.damage_map(x, np.inf, 1234.0, 100.0)
    assert np.array_equal(out, x)


def test_damage_map_identity_when_tangent_equals_elastic() -> None:
    """E_t=E ⇒ 屈服后斜率=1 ⇒ 整体恒等（一个 E_t 参数的连续性端点）。"""
    x = np.array([1.0, 49.0, 50.0, 51.0, 300.0])
    out = bd.damage_map(x, 50.0, 100.0, 100.0)
    assert np.allclose(out, x)


def test_damage_map_bilinear_knee_values() -> None:
    """σ_y 以下恒等；σ_y 以上斜率 E_t/E。"""
    x = np.array([40.0, 50.0, 100.0])          # σ_lin
    out = bd.damage_map(x, 50.0, 20.0, 100.0)  # σ_y=50, E_t/E=0.2
    assert out[0] == pytest.approx(40.0)
    assert out[1] == pytest.approx(50.0)
    assert out[2] == pytest.approx(50.0 + 0.2 * (100.0 - 50.0))  # = 60


def test_bone_damage_law_sigma_matches_manual() -> None:
    law = bd.BoneDamageLaw(E_mpa=100.0, yield_stress_mpa=50.0, tangent_modulus_mpa=20.0)
    assert law.yield_strain == pytest.approx(0.5)
    assert law.tangent_ratio == pytest.approx(0.2)
    assert law.sigma(0.4) == pytest.approx(40.0)          # 弹性
    assert law.sigma(0.5) == pytest.approx(50.0)          # 膝点
    assert law.sigma(1.0) == pytest.approx(50.0 + 20.0 * 0.5)  # 塑性
    assert law.sigma_from_linear_stress(40.0) == pytest.approx(40.0)
    assert law.sigma_from_linear_stress(100.0) == pytest.approx(60.0)


def test_sigma_h_from_law_matches_manual() -> None:
    law = bd.BoneDamageLaw(E_mpa=18000.0, yield_stress_mpa=36.1, tangent_modulus_mpa=14800.0)
    load = np.array([1000.0, 5000.0, 25000.0])
    sig = bd.sigma_h_from_law(load, law, sigma_ref_mpa=50.0, load_ref_n=5000.0)
    lin = 50.0 * load / 5000.0
    manual = bd.damage_map(lin, 36.1, 14800.0, 18000.0)
    assert np.allclose(sig, manual)


# --------------------------------------------------------------------------
# 2. 分类与 null 工具
# --------------------------------------------------------------------------
def test_classify_categories() -> None:
    lin = np.linspace(1.0, 100.0, 50)
    heights = np.arange(1, 51, dtype=float)
    sat = bd.classify_step_progressive(
        lin, bd.BoneDamageLaw(100.0, 50.0, 0.0), heights=heights)
    assert sat["class"] == "saturation"
    grad = bd.classify_step_progressive(
        lin, bd.BoneDamageLaw(100.0, 50.0, 90.0), heights=heights)
    assert grad["class"] == "gradual"
    step = bd.classify_step_progressive(
        lin, bd.BoneDamageLaw(100.0, 50.0, 30.0), heights=heights)
    assert step["class"] == "step"
    elas = bd.classify_step_progressive(
        np.linspace(1.0, 10.0, 50), bd.BoneDamageLaw(100.0, 50.0, 0.0), heights=heights)
    assert elas["class"] == "elastic_no_knee"
    assert elas["yielded_in_range"] is False


def test_broken_stick_properties() -> None:
    bs = bd.broken_stick(9)
    assert bs.shape == (9,)
    assert bs.sum() == pytest.approx(1.0, abs=1e-12)
    assert np.all(np.diff(bs) < 0)      # 单调递减
    assert bs[1] == pytest.approx(0.20322, abs=1e-4)


def test_parallel_analysis_null_deterministic() -> None:
    scales = np.ones(6)
    p1, _ = bd.parallel_analysis_null(scales, n_samples=50, n_trials=200, seed=7)
    p2, _ = bd.parallel_analysis_null(scales, n_samples=50, n_trials=200, seed=7)
    assert np.array_equal(p1, p2)
    p3, _ = bd.parallel_analysis_null(scales, n_samples=50, n_trials=200, seed=8)
    assert not np.array_equal(p1, p3)


# --------------------------------------------------------------------------
# 3. read_materials 逐位不变 + read_material_damage
# --------------------------------------------------------------------------
@need_data
def test_read_materials_bit_identical() -> None:
    """既有构造未被新增字段改变（钉住逐值 + 新字段默认 None）。"""
    pids = [81000700, 89001501, 81001300]
    mats = bone_feb.read_materials(MASTER, pids)
    assert set(mats) == set(pids)
    m = mats[81000700]
    assert (m.pid, m.E_mpa, m.nu) == (81000700, 18000.0, 0.3)
    assert m.src_keyword == "*MAT_PIECEWISE_LINEAR_PLASTICITY"
    assert m.note == ("elastic approximation of MAT_PIECEWISE_LINEAR_PLASTICITY "
                      "(plastic/damage ignored)")
    assert m.src_card == [81000700.0, 2e-09, 18000.0, 0.3, 36.1, 14800.0, 0.0, 0.0]
    assert m.yield_stress_mpa is None and m.tangent_modulus_mpa is None
    assert mats[89001501].src_keyword == "*MAT_PLASTICITY_WITH_DAMAGE"
    assert mats[89001501].E_mpa == 13020.0


@need_data
def test_read_material_damage_fields() -> None:
    """read_material_damage 读卡 nums[4]/nums[5]，且既有字段与 read_materials 相同。"""
    pids = [81000700, 89001501, 81001300]
    base = bone_feb.read_materials(MASTER, pids)
    dmg = bone_feb.read_material_damage(MASTER, pids)
    assert set(base) == set(dmg)
    for pid in pids:
        b, d = base[pid], dmg[pid]
        assert (b.pid, b.E_mpa, b.nu, b.src_keyword, b.src_card, b.note) == \
               (d.pid, d.E_mpa, d.nu, d.src_keyword, d.src_card, d.note)
    assert dmg[81000700].yield_stress_mpa == pytest.approx(36.1)
    assert dmg[81000700].tangent_modulus_mpa == pytest.approx(14800.0)
    assert dmg[89001501].yield_stress_mpa == pytest.approx(80.0)
    assert dmg[89001501].tangent_modulus_mpa == pytest.approx(0.0)
    assert dmg[81001300].yield_stress_mpa == pytest.approx(224.0)


# --------------------------------------------------------------------------
# 4. 驱动：规则复用锁 + rank-2 提升 + 诚实判定
# --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def g4res() -> dict:
    return G.build()


@need_data
def test_transcribed_rule_equals_frozen_a4_from_matrix() -> None:
    """转录规则在线性轨迹上与冻结 ``a4_from_matrix`` 逐值相等（规则未变）。"""
    import json
    fm = json.loads(_BC2.read_text(encoding="utf-8"))
    heights = np.asarray(fm["heights_m"], dtype=float)
    yi = G.s4.route_yi("median")
    sig_lin = G.build_linear_trajectories(yi, fm)
    ref = G.s4.a4_from_matrix(yi, heights, fm)
    tr = G.a4_rank2_from_sigma(sig_lin)
    assert tr["separable"] == ref["separable"]
    assert tr["verdict"] == ref["verdict"]
    assert tr["explained_var_ratio"][0] == pytest.approx(
        ref["explained_var_ratio"][0], rel=1e-12)
    assert tr["explained_var_ratio"][1] == pytest.approx(
        ref["explained_var_ratio"][1], rel=1e-9)


@need_data
def test_rank2_improvement_and_honest_verdict(g4res: dict) -> None:
    r2 = g4res["rank2_damaged"]
    # 数值上脱离 rank-1 天花板（冻结基线 PC2 ≈ 2.42e-4）
    assert g4res["rank2_baseline"]["pc2_explained_ratio"] == pytest.approx(
        2.421e-4, rel=1e-2)
    assert r2["baseline_pc2_explained_ratio"] == pytest.approx(2.421e-4, rel=1e-2)
    assert r2["pc2_explained_ratio"] == pytest.approx(8.791e-3, rel=1e-2)
    assert r2["pc2_explained_ratio"] > 10.0 * r2["baseline_pc2_explained_ratio"]
    assert g4res["verdicts"]["break_rank1_ceiling_numeric"] is True
    # 但数据驱动 null 门如实为 FAIL（PC2 < 平行分析 95% null）
    assert r2["pc2_exceeds_pa_null_p95"] is False
    assert r2["rank2_holds"] is False
    assert g4res["verdicts"]["a4_rank2"] == "FAIL"
    assert g4res["verdicts"]["a4_lp_separation"] == "FAIL"


@need_data
def test_constant_saturated_parts_detected(g4res: dict) -> None:
    """脊柱 σ_y=80 < σ_c 且 h=1 即屈服 ⇒ 损伤后 σ 恒定 ⇒ 从 rank-2 检验剔除。"""
    r2 = g4res["rank2_damaged"]
    assert r2["n_active_parts"] == 6
    assert set(r2["constant_saturated_parts"]) == {"腰椎", "胸椎", "颈椎"}


@need_data
def test_a2_femur_thoracic_no_discrimination(g4res: dict) -> None:
    a2 = g4res["a2"]
    assert a2["femur"]["damage_status"] == "all_zero"
    assert a2["thoracic"]["damage_status"] == "all_zero"
    # 股骨基线亦为常值0（未获区分度）；胸椎由常值1退化为常值0
    assert a2["femur"]["line_status"] == "all_zero"
    assert a2["thoracic"]["line_status"] == "all_one"


@need_data
def test_driver_output_schema(g4res: dict) -> None:
    assert set(g4res["verdicts"]) == {
        "a4_rank2", "a4_lp_separation", "break_rank1_ceiling_numeric"}
    assert len(g4res["parts"]) == 9
    p = g4res["parts"][1]
    assert p["bone"] == "tibia_r" and p["pid"] == 81000700
    assert p["yield_stress_mpa"] == pytest.approx(36.1)
    assert p["tangent_ratio_Et_over_E"] == pytest.approx(14800.0 / 18000.0)
    assert len(p["sigma_h_linear_mpa"]) == 50
    assert len(p["sigma_h_damaged_mpa"]) == 50
    # 损伤后 σ ≤ 线性 σ（E_t/E ≤ 1；饱和端严格小于）
    lin = np.array([q["sigma_h_linear_mpa"] for q in g4res["parts"]])
    dmg = np.array([q["sigma_h_damaged_mpa"] for q in g4res["parts"]])
    assert np.all(dmg <= lin + 1e-9)
