"""tests/test_g8_fibula_or.py -- G8 腓骨 1D 载荷分享接入 G3 的回归锁。

锁住（见 `docs/methods/G8_腓骨OR_1D载荷分享设计.md`）：
* 数据源换算：`σ_fib(h) = risk_by_height(h) × σc_row(160)`，**非 ×70**（勘误）。
* 混合矩阵：只换腓骨行，其余 8 行逐值不变；原 inputs 不被就地改动。
* 其余 8 部位 OR 与 G3 主设定逐值相同（DoD-2）。
* 腓骨从 `all_one` 变**渐进**（status=ok、OR 有限）——机制解锁。
* A3 ≥2/3、`logistic_fit` 接口稳定、极限退化收敛。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

_ROOT = Path(__file__).resolve().parents[1]
for _p in (_ROOT, _ROOT / "src", _ROOT / "scripts" / "opensim_fe"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from scripts.opensim_fe import g3_subject_variability as g3  # noqa: E402
from scripts.opensim_fe import g8_fibula_or_share as g8      # noqa: E402
from scripts.opensim_fe import paper_alignment as pa         # noqa: E402


@pytest.fixture(scope="module")
def inputs():
    return g3.load_inputs()


@pytest.fixture(scope="module")
def fib05():
    return g8.load_fibula_1d("rows_with_share_lit05")


@pytest.fixture(scope="module")
def primary_result(inputs):
    return g8.run_share(inputs, g8.load_fibula_1d("rows_with_share_lit05"))


# ---------------------------------------------------------------------------
def test_height_grid_matches(inputs, fib05):
    chk = g8.height_grid_check(inputs, fib05)
    assert chk["match"] is True
    assert chk["g3_heights_n"] == 50
    assert chk["fibula_1d_heights_n"] == 50


def test_fibula_sigma_is_risk_times_row_sigma_c(fib05):
    # 关键：σ_fib = risk × σc_row(160)，非 ×70（设计 §3 勘误锁）
    assert fib05["sigma_c_row_mpa"] == 160.0
    assert np.allclose(fib05["sigma_fib_mpa"], fib05["risk_by_height"] * 160.0)
    assert not np.allclose(fib05["sigma_fib_mpa"], fib05["risk_by_height"] * 70.0)
    assert fib05["risk_at_5m"] == pytest.approx(0.0775, abs=1e-3)


def test_mixed_inputs_only_fibula_row_changes(inputs, fib05):
    base = np.asarray(inputs["sigma_nom"], dtype=float).copy()
    mixed, idx = g8.build_mixed_inputs(inputs, fib05["sigma_fib_mpa"])
    assert inputs["bones"][idx] == "fibula_r"
    out = np.asarray(mixed["sigma_nom"], dtype=float)
    for k in range(base.shape[0]):
        if k == idx:
            continue
        assert np.array_equal(base[k], out[k]), f"row {k} changed"
    assert np.allclose(out[idx], fib05["sigma_fib_mpa"])
    # 原 inputs 未被就地修改
    assert np.array_equal(np.asarray(inputs["sigma_nom"], dtype=float), base)


def test_fibula_becomes_progressive_not_all_one(primary_result):
    rec = primary_result["fibula_rec"]
    assert rec["status"] == "ok"
    ov = rec["OR"]
    assert ov is not None and np.isfinite(ov)
    # 去奇异：本轮实测 OR≈1.152（5% 分担）
    assert 1.10 <= ov <= 1.20


def test_other8_or_unchanged(inputs, primary_result):
    idx = list(inputs["bones"]).index("fibula_r")
    unc = g8.unchanged_check(inputs, idx, primary_result["primary"])
    assert unc["available"] is True
    assert unc["all_equal"] is True
    assert unc["max_abs_or_diff"] <= 1e-9


def test_a3_still_at_least_two_of_three(primary_result):
    a3 = primary_result["a3"]
    assert a3["n_within"] >= 2
    assert a3["verdict"] == "PASS"


def test_logistic_fit_interface_stable():
    x = np.arange(1, 51, dtype=float)
    y = (x >= 25).astype(int)
    rec = pa.logistic_fit(x, y)
    assert rec["status"] in {"ok", "all_zero", "all_one", "complete_separation"}


def test_degenerate_converges_to_base(inputs, fib05):
    mixed, _ = g8.build_mixed_inputs(inputs, fib05["sigma_fib_mpa"])
    _, p = g3.simulate(g3.OFFICIAL_SEED, 120, g3.COV_DEGENERATE,
                       g3.COV_DEGENERATE, g3.COV_DEGENERATE, mixed,
                       mass_in_load=False)
    base = g3.deterministic_base(mixed["sigma_nom"], mixed["sigma_c"])
    assert np.array_equal((p >= 0.5).astype(int), base)


def test_verdict_records_honest_fail(inputs, primary_result):
    v = g8.verdict_text(primary_result)
    # 本轮实测：OR≈1.152 < 靶窗下沿 1.177 ⇒ DoD-1 如实 FAIL（未回解）
    assert v["fibula_status"] == "ok"
    assert v["fibula_within"] is False
    assert v["fibula_window"][0] == pytest.approx(1.1774, abs=1e-3)
    assert v["dod3_pass"] is True
