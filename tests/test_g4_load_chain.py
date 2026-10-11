"""tests/test_g4_load_chain.py -- G4 载荷链顺应性杠杆的接口、不变量与诚实判定锁。

对应设计 ``docs/methods/G4_损伤演化设计.md``（§3 行4 / §6 D6）与实现
``scripts/opensim_fe/g4_load_chain_compliant.py``。**纯增量**、**不含慢 OpenSim 调用**
（OpenSim 载荷链结果从既有 ``load_chain_compliant.json`` 读取）。

锁定的不变量
------------
1. **载荷注入恒等**：``build_linear_trajectories_override(yi, fm, fm["loads_n"])`` 与
   冻结的 ``g4_damage_evolution.build_linear_trajectories(yi, fm)`` **逐值相等**
   （9 部位）——即"仅换载荷来源"这一参数化未改轨迹规则。
2. **损伤轨迹恒等**：同一载荷向量下，经同一 ``BoneDamageLaw`` 映射得到的损伤 σ(h)
   与最小路径 ``g4_damage_matrix.json`` 的 ``sigma_h_damaged_mpa`` **逐值相等**。
3. **新件可解析**：``load_chain_compliant.json`` / ``g4_loadchain_matrix.json`` /
   ``g4_loadchain_run.json`` schema 齐备（关节、50 高度、meta 记录预设/质量/命令）。
4. **冻结点未变**：``fracture_matrix_bc2.json`` 的 ``loads_n`` 逐值钉住（未被本杠杆改动）；
   ``g4_damage_*`` 最小路径产物仍在且数值未变。
5. **诚实判定锁**：顺应链的三路结果（rank-2 / LP / A2 股骨·胸椎 / §7 触发）**如实**落盘。
"""
from __future__ import annotations

import json
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
from scripts.opensim_fe import g4_load_chain_compliant as LC  # noqa: E402
from scripts.opensim_fe import paper_alignment as pa  # noqa: E402

MASTER = ROOT / "model" / "AM50_V71_Occupant" / "main_THUMS_AM50_V71.k"
_BC2 = ROOT / "results" / "opensim_fe" / "fracture_matrix_bc2.json"
_H5 = ROOT / "temp" / "opensim_fe" / "fracture_matrix_bc2" / "tibia_r" / "tibia_r.hdf5"
_DMG = ROOT / "results" / "opensim_fe" / "g4_damage_matrix.json"
_LC_JSON = ROOT / "results" / "opensim_fe" / "load_chain_compliant.json"
_LC_MATRIX = ROOT / "results" / "opensim_fe" / "g4_loadchain_matrix.json"
_LC_RUN = ROOT / "results" / "opensim_fe" / "g4_loadchain_run.json"
_LC_REPORT = ROOT / "results" / "opensim_fe" / "G4_LOADCHAIN_REPORT.md"

need_data = pytest.mark.skipif(
    not (MASTER.is_file() and _BC2.is_file() and _H5.is_file()),
    reason="THUMS 卡 / 冻结 bc2 JSON / hdf5 未随仓库分发或未生成")
need_loads = pytest.mark.skipif(
    not (_LC_JSON.is_file() and _LC_MATRIX.is_file() and _LC_RUN.is_file()),
    reason="顺应链产物未生成（先跑 scripts/opensim_fe/g4_load_chain_compliant.py）")


def _fm() -> dict:
    return json.loads(_BC2.read_text(encoding="utf-8"))


def _laws() -> dict:
    pid_by_bone = G.main_domain_pid(_fm())
    mats = bone_feb.read_material_damage(
        MASTER, sorted(set(pid_by_bone.values())))
    out = {}
    for cn, bone, _joint, _subs in pa.REGIONS:
        m = mats[pid_by_bone[bone]]
        out[cn] = bd.BoneDamageLaw(E_mpa=m.E_mpa, yield_stress_mpa=m.yield_stress_mpa,
                                   tangent_modulus_mpa=m.tangent_modulus_mpa)
    return out


# --------------------------------------------------------------------------
# 1/2. 载荷注入恒等 + 损伤轨迹恒等
# --------------------------------------------------------------------------
@need_data
def test_override_identity_with_frozen_loads() -> None:
    """传入冻结 ``loads_n`` 时，注入版轨迹与冻结 ``build_linear_trajectories`` 逐值相等。"""
    fm = _fm()
    yi = G.s4.route_yi("median")
    ref = G.build_linear_trajectories(yi, fm)
    got = LC.build_linear_trajectories_override(yi, fm, fm["loads_n"])
    assert set(got) == set(ref) == set(pa.REGION_NAMES)
    for cn in pa.REGION_NAMES:
        assert np.array_equal(np.asarray(got[cn]), np.asarray(ref[cn])), cn


@need_data
def test_identical_load_reproduces_damage_trajectory() -> None:
    """同一载荷向量 ⇒ 同一双线性映射 ⇒ 与最小路径 ``g4_damage_matrix.json`` 逐值相等。"""
    if not _DMG.is_file():
        pytest.skip("g4_damage_matrix.json 缺失（最小路径产物）")
    fm = _fm()
    yi = G.s4.route_yi("median")
    lin = LC.build_linear_trajectories_override(yi, fm, fm["loads_n"])
    laws = _laws()
    dmg_json = json.loads(_DMG.read_text(encoding="utf-8"))
    ref = {p["part_cn"]: np.asarray(p["sigma_h_damaged_mpa"], dtype=float)
           for p in dmg_json["parts"]}
    for cn in pa.REGION_NAMES:
        got = laws[cn].sigma_from_linear_stress(lin[cn])
        assert np.allclose(got, ref[cn], rtol=1e-12, atol=0.0), cn


# --------------------------------------------------------------------------
# 3. 新件可解析
# --------------------------------------------------------------------------
@need_loads
def test_load_chain_compliant_parses() -> None:
    d = json.loads(_LC_JSON.read_text(encoding="utf-8"))
    m = d["meta"]
    assert m["preset"] == LC.PRESET_NAME
    assert m["mass_kg"] == LC.MODEL_MASS_KG
    assert m["heights_m"] == list(range(1, 51))
    assert set(d["loads_n"]) == {"ankle_r", "hip_r"}
    for j, arr in d["loads_n"].items():
        assert len(arr) == 50 and np.all(np.isfinite(np.asarray(arr, float))), j
    assert "command" in m and "g4_load_chain_compliant.py" in m["command"]


@need_loads
def test_loadchain_matrix_and_run_parse() -> None:
    mat = json.loads(_LC_MATRIX.read_text(encoding="utf-8"))
    for key in ("meta", "heights_m", "parts", "verdicts",
                "a4_baseline_rank1", "a4_compliant_linear", "a4_compliant_damaged",
                "rank2_baseline_rigid", "rank2_damage_only", "rank2_compliant",
                "a2_compliant", "loads_n_compliant"):
        assert key in mat, key
    assert len(mat["parts"]) == 9
    p = mat["parts"][1]
    assert p["bone"] == "tibia_r" and p["pid"] == 81000700
    assert len(p["sigma_h_linear_compliant_mpa"]) == 50
    assert len(p["sigma_h_damaged_compliant_mpa"]) == 50

    run = json.loads(_LC_RUN.read_text(encoding="utf-8"))
    assert set(run["pc2_explained_ratio"]) == {
        "rank1_baseline", "damage_only", "damage_compliant"}
    assert set(run["rank2_holds"]) == {"damage_only", "damage_compliant"}
    assert set(run["a4_lp_separation"]) == {
        "rank1_baseline", "damage_only", "damage_compliant"}
    assert set(run["a2_femur_thoracic"]["damage_compliant"]) == {"股骨", "胸椎"}
    assert _LC_REPORT.is_file() and _LC_REPORT.stat().st_size > 0


# --------------------------------------------------------------------------
# 4. 冻结点未变（只读输入 + 最小路径产物）
# --------------------------------------------------------------------------
@need_data
def test_frozen_inputs_pinned() -> None:
    fm = _fm()
    # 冻结（刚性腿）载荷链逐值钉住 —— 本杠杆**不改**冻结输入。
    assert fm["loads_n"]["ankle_r"][0] == pytest.approx(15454.5, rel=1e-4)
    assert fm["loads_n"]["ankle_r"][4] == pytest.approx(25266.2, rel=1e-4)
    assert fm["loads_n"]["hip_r"][0] == pytest.approx(9783.8, rel=1e-4)
    assert fm["loads_n"]["hip_r"][49] == pytest.approx(57941.3, rel=1e-4)
    # 最小路径产物仍在且数值未变（未被本杠杆覆盖）。
    if _DMG.is_file():
        d = json.loads(_DMG.read_text(encoding="utf-8"))
        assert d["verdicts"]["a4_rank2"] == "FAIL"
        assert d["rank2_damaged"]["pc2_explained_ratio"] == pytest.approx(
            8.791e-3, rel=1e-3)


@need_loads
def test_compliant_loads_differ_from_frozen() -> None:
    """顺应链**确实**改变了载荷（非恒等）——否则杠杆未生效。"""
    fm = _fm()
    d = json.loads(_LC_JSON.read_text(encoding="utf-8"))
    for j in ("ankle_r", "hip_r"):
        a = np.asarray(d["loads_n"][j], float)
        b = np.asarray(fm["loads_n"][j], float)
        assert not np.allclose(a, b)


# --------------------------------------------------------------------------
# 5. 诚实判定锁（三路结果）
# --------------------------------------------------------------------------
@need_loads
def test_honest_verdicts_locked() -> None:
    run = json.loads(_LC_RUN.read_text(encoding="utf-8"))
    v = run["verdicts"]
    # A4 rank-2 与 LP 均如实 FAIL；§7 触发条件未成立。
    assert v["a4_rank2"] == "FAIL"
    assert v["a4_lp_separation"] == "FAIL"
    assert v["section7_trigger_met"] is False
    assert run["rank2_holds"]["damage_only"] is False
    assert run["rank2_holds"]["damage_compliant"] is False
    assert run["a4_lp_separation"]["damage_compliant"] == "FAIL"
    # 数值上仍脱离 rank-1 天花板（PC2 抬升），但未超数据驱动 null。
    assert (run["pc2_explained_ratio"]["damage_compliant"]
            > 10.0 * run["pc2_explained_ratio"]["rank1_baseline"])
    # A2：股骨在顺应链下获阶跃区分度；胸椎仍无（冻结 median + 封顶）。
    a2 = run["a2_femur_thoracic"]["damage_compliant"]
    assert a2["股骨"] == "complete_separation"
    assert a2["胸椎"] == "all_zero"
