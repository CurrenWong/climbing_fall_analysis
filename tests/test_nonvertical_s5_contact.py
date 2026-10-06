"""Phase-S5-contact 波4（跟骨跖面真实接触 + 摩擦）结构/契约回归。

**不调用 FEBio、不跑 OpenSim**（毫秒级），只锁死：

* CLI / argparse 表面（`--out-json` / `--out-md` / `--fast`）；
* 覆盖保护：产物已存在 → `main()` 返回 **2**（不覆盖既有 `results/opensim_fe/*`）；
* JSON schema 形状：**11 个顶层键**齐全、**无 `sigma_over_sigma_law` 键**（禁用键）；
* scope-limit 绑定字段（§4 条 7）；
* G1 = `N/A_minimum_version`、G7 = `RETIRED_by_contract_条6-8`（二者均不省略）；
* 报告 13 节 `## 0.`–`## 12.`；
* 默认 `plantar_bc="fixed"` + `contact_mu` opt-in 默认不变（防回归）。

慢 FE 数值矩阵在 `scripts/opensim_fe/nonvertical_s5_contact.py`，一次跑完并写报告。
"""

from __future__ import annotations

import inspect
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts" / "opensim_fe") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts" / "opensim_fe"))

import nonvertical_s5_contact as S5C  # noqa: E402

JSON_PATH = ROOT / "results" / "opensim_fe" / "nonvertical_s5_contact.json"
MD_PATH = ROOT / "results" / "opensim_fe" / "NONVERTICAL_S5_CONTACT_REPORT.md"

#: §1 定义的 11 个顶层键（顺序 + 集合）。
EXPECTED_TOP_KEYS = (
    "meta", "scope_limit", "interface_regression", "axial_regression",
    "scenario_matrix", "contact_recipe", "runs", "regression_gates",
    "sensitivity", "assumptions", "repro",
)


# --------------------------------------------------------------------------
# 1. CLI / argparse + 覆盖保护
# --------------------------------------------------------------------------
def test_output_basenames_fixed_by_design() -> None:
    assert S5C.OUT_JSON.name == "nonvertical_s5_contact.json"
    assert S5C.OUT_MD.name == "NONVERTICAL_S5_CONTACT_REPORT.md"
    # 不能占用 supination 的 basename
    assert S5C.OUT_JSON.name != "nonvertical_s5.json"


def test_cli_surface_accepts_out_json_out_md_fast() -> None:
    """`main` 的 argparse 接受 --out-json / --out-md / --fast（--help 退出码 0）。"""
    with pytest.raises(SystemExit) as exc:
        S5C.main(["--help"])
    assert exc.value.code == 0


def test_overwrite_guard_returns_2(tmp_path: Path) -> None:
    """产物已存在 → 返回 2，且不触碰（内容不变）。"""
    j = tmp_path / "nonvertical_s5_contact.json"
    m = tmp_path / "NONVERTICAL_S5_CONTACT_REPORT.md"
    j.write_text("{}", encoding="utf-8")
    m.write_text("# x\n", encoding="utf-8")
    rc = S5C.main(["--out-json", str(j), "--out-md", str(m)])
    assert rc == 2
    assert j.read_text(encoding="utf-8") == "{}"
    assert m.read_text(encoding="utf-8") == "# x\n"


def test_overwrite_guard_returns_2_if_only_md_exists(tmp_path: Path) -> None:
    m = tmp_path / "NONVERTICAL_S5_CONTACT_REPORT.md"
    m.write_text("# x\n", encoding="utf-8")
    rc = S5C.main(["--out-json", str(tmp_path / "new.json"), "--out-md", str(m)])
    assert rc == 2


# --------------------------------------------------------------------------
# 2. JSON schema 形状（读已生成产物；缺失则跳过）
# --------------------------------------------------------------------------
pytestmark_json = pytest.mark.skipif(not JSON_PATH.is_file(),
                                     reason="nonvertical_s5_contact.json 未生成")


@pytestmark_json
def test_json_top_level_has_exactly_11_keys() -> None:
    data = json.loads(JSON_PATH.read_text(encoding="utf-8"))
    assert set(data.keys()) == set(EXPECTED_TOP_KEYS)
    assert len(data) == len(EXPECTED_TOP_KEYS)


@pytestmark_json
def test_json_has_no_sigma_over_sigma_law_key_anywhere() -> None:
    """§1.9 / §5.1：任何 σ/σ_law 键都不得存在（含嵌套）。"""
    text = JSON_PATH.read_text(encoding="utf-8")
    assert "sigma_over_sigma_law" not in text
    assert "sigma_contact_over_sigma_law" not in text
    assert "sigma_law" not in text

    def _walk(o):
        if isinstance(o, dict):
            for k, v in o.items():
                assert "sigma_over_sigma_law" not in str(k)
                _walk(v)
        elif isinstance(o, list):
            for v in o:
                _walk(v)

    _walk(json.loads(text))


@pytestmark_json
def test_scope_limit_binds_contract_clause_7() -> None:
    sl = json.loads(JSON_PATH.read_text(encoding="utf-8"))["scope_limit"]
    assert "条 7" in sl["contract_clause"]
    # 四个 machine-checkable 的三元组必须都在，且 value 未翻转
    assert sl["no_quotable_sigma_ratio"]["value"] is True
    assert sl["citation_policy"]["value"] == "QUOTABLE-only"
    assert sl["recipe_constraint"]["value"] == "minimum_version_only"
    assert sl["fe_output_constraint"]["value"] == "absolute_quantities_only"


@pytestmark_json
def test_gates_g1_and_g7_codes_present_and_distinct() -> None:
    gates = json.loads(JSON_PATH.read_text(encoding="utf-8"))["regression_gates"]
    assert gates["g1_pad_stress_reproduced"]["verdict"] == "N/A_minimum_version"
    assert gates["g7_penetration_window"]["verdict"] == "RETIRED_by_contract_条6-8"
    # 二者语义不同，不得共用同一个码
    assert (gates["g1_pad_stress_reproduced"]["verdict"]
            != gates["g7_penetration_window"]["verdict"])


@pytestmark_json
def test_runs_and_scenario_matrix_present() -> None:
    data = json.loads(JSON_PATH.read_text(encoding="utf-8"))
    assert isinstance(data["runs"], list) and data["runs"]
    axes = data["scenario_matrix"]["axes"]
    assert axes["mu"] == [0.0, 0.3, 0.6]
    assert "hard_surface" in axes["pad_scenario"]
    # 每条 run 必须带 run_verdict
    for r in data["runs"]:
        assert r["run_verdict"] in ("PASS", "FAIL", "SKIPPED")
        assert r["plantar_bc"] in ("fixed", "spring", "contact")


@pytestmark_json
def test_report_has_thirteen_sections() -> None:
    text = MD_PATH.read_text(encoding="utf-8")
    for i in range(13):
        assert f"## {i}." in text, f"缺 ## {i}. 节"


# --------------------------------------------------------------------------
# 3. 默认逐位不变（防回归；不跑 FEBio）
# --------------------------------------------------------------------------
def test_thums_plantar_bc_default_fixed() -> None:
    import thums_feb

    sig = inspect.signature(thums_feb.build_thums_feb)
    assert sig.parameters["plantar_bc"].default == "fixed"


def test_thums_contact_mu_optin_default_is_recipe_mu() -> None:
    """`contact_mu` 是 opt-in；默认 = plantar_bc.DEFAULT_CONTACT_MU（不改现有默认）。"""
    import thums_feb

    from climbing.coupling.plantar_bc import DEFAULT_CONTACT_MU

    sig = inspect.signature(thums_feb.build_thums_feb)
    assert sig.parameters["contact_mu"].default == DEFAULT_CONTACT_MU


def test_recipe_is_new_penalty_not_auglag() -> None:
    """单配方：NEW PENALTY（禁止 OLD AUGLAG；§5.1 #1）。"""
    assert S5C.RECIPE["laugon"] == "PENALTY"
    assert S5C.RECIPE["penalty"] == 0.1
    assert S5C.RECIPE["feb_section_order"] == ["Loads", "Boundary", "Contact"]
