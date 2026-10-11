"""tests/test_g3_subject_variability.py -- G3 被试级变异统计层的接口与不变量锁。

对应设计 `docs/methods/G3_被试级变异_统计层设计.md`（预注册）与脚本
`scripts/opensim_fe/g3_subject_variability.py`。**纯增量**：不改任何既有模块。

锁定的不变量
------------
1. lognormal 约定（mean-1）：E[X]=1；COV→0 ⇒ (0,0)。
2. 输入口径：部位顺序 = ``paper_alignment.REGION_NAMES``；σc = route-乙 弱区。
3. 极限退化（DoD-2）：COV=0（+ 质量不入载荷）逐格复现确定性基准；
   COV=0.001 经多数判定复现。
4. 单调性（DoD-3）：每部位 p[h] 随 h 非降（对任意 n 成立）。
5. 种子复现：同 seed 逐位可复现；异 seed 不同。
6. 复用锁：估计器就是 ``paper_alignment.logistic_fit`` 本体（不换、不包）；
   A3 靶与判定规则 = ``s4_alignment_bc2`` 原件。
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

from scripts.opensim_fe import g3_subject_variability as G  # noqa: E402
from scripts.opensim_fe import paper_alignment as pa  # noqa: E402
from scripts.opensim_fe import s4_alignment_bc2 as s4  # noqa: E402

SEED = G.OFFICIAL_SEED
N_SMALL = 50


@pytest.fixture(scope="module")
def inputs() -> G.Inputs:
    return G.load_inputs()


# --------------------------------------------------------------------------
# 1. lognormal 约定
# --------------------------------------------------------------------------
@pytest.mark.parametrize("cov", [0.0, 0.05, 0.10, 0.15, 0.30])
def test_lognormal_params_documented_convention(cov: float):
    mu, sig = G.log_params(cov)
    # sigma_ln = sqrt(ln(1+c^2))
    assert sig == pytest.approx(np.sqrt(np.log1p(cov * cov)), abs=1e-15)
    # mu = -sigma^2/2  ⇒ E[X] = exp(mu + sigma^2/2) = 1
    assert mu == pytest.approx(-0.5 * sig * sig, abs=1e-15)
    assert np.exp(mu + 0.5 * sig * sig) == pytest.approx(1.0, abs=1e-12)


def test_lognormal_cov_zero_is_identity():
    assert G.log_params(0.0) == (0.0, 0.0)
    # np.random.lognormal(0, 0) == 1
    rng = np.random.default_rng([1, 2])
    assert np.all(rng.lognormal(0.0, 0.0, size=5) == 1.0)


# --------------------------------------------------------------------------
# 2. 输入口径
# --------------------------------------------------------------------------
def test_inputs_site_order_and_weak_sigma_c(inputs: G.Inputs):
    assert inputs["sites"] == list(pa.REGION_NAMES)
    # route-乙 弱区 σc（design §2 第 6 行）
    expect = {"足部": 150.0, "胫骨": 70.0, "腓骨": 70.0, "股骨": 80.0,
              "骨盆": 180.0, "腰椎": 150.0, "胸椎": 150.0, "颈椎": 150.0,
              "颅骨": 160.0}
    assert inputs["weak_sigma_c_mpa"] == expect
    assert inputs["sigma_nom"].shape == (9, 50)
    assert inputs["heights"].shape == (50,)


def test_deterministic_base_is_threshold_indicator(inputs: G.Inputs):
    base = G.deterministic_base(inputs["sigma_nom"], inputs["sigma_c"])
    assert base.shape == (9, 50)
    assert np.array_equal(base.astype(bool),
                          inputs["sigma_nom"] >= inputs["sigma_c"][:, None])


# --------------------------------------------------------------------------
# 3. 极限退化（DoD-2）
# --------------------------------------------------------------------------
def test_degenerate_limit_cov_zero_bit_exact(inputs: G.Inputs):
    base = G.deterministic_base(inputs["sigma_nom"], inputs["sigma_c"])
    _y, p = G.simulate(SEED, N_SMALL, 0.0, 0.0, 0.0, inputs, mass_in_load=False)
    assert np.array_equal(p.astype(int), base)


def test_degenerate_limit_cov_near_zero_majority_matches(inputs: G.Inputs):
    base = G.deterministic_base(inputs["sigma_nom"], inputs["sigma_c"])
    _y, p = G.simulate(SEED, 200, G.COV_DEGENERATE, G.COV_DEGENERATE,
                       G.COV_DEGENERATE, inputs, mass_in_load=False)
    assert np.array_equal((p >= 0.5).astype(int), base)


# --------------------------------------------------------------------------
# 4. 单调性（DoD-3）
# --------------------------------------------------------------------------
def test_probability_monotone_non_decreasing(inputs: G.Inputs):
    _y, p = G.simulate(SEED, 120, G.COV_B_PRIMARY, G.COV_E_PRIMARY,
                       G.COV_A_PRIMARY, inputs, mass_in_load=True)
    assert np.all(np.diff(p, axis=1) >= -1e-12)


# --------------------------------------------------------------------------
# 5. 种子复现 / RNG 纪律
# --------------------------------------------------------------------------
def test_seed_reproducible_bit_identical(inputs: G.Inputs):
    y1, p1 = G.simulate(SEED, 30, G.COV_B_PRIMARY, G.COV_E_PRIMARY,
                        G.COV_A_PRIMARY, inputs, mass_in_load=True)
    y2, p2 = G.simulate(SEED, 30, G.COV_B_PRIMARY, G.COV_E_PRIMARY,
                        G.COV_A_PRIMARY, inputs, mass_in_load=True)
    assert np.array_equal(y1, y2) and np.array_equal(p1, p2)


def test_different_seed_changes_outcome(inputs: G.Inputs):
    y1, _ = G.simulate(SEED, 30, G.COV_B_PRIMARY, G.COV_E_PRIMARY,
                       G.COV_A_PRIMARY, inputs, mass_in_load=True)
    y2, _ = G.simulate(SEED + 1, 30, G.COV_B_PRIMARY, G.COV_E_PRIMARY,
                       G.COV_A_PRIMARY, inputs, mass_in_load=True)
    assert not np.array_equal(y1, y2)


def test_rng_uses_seed_subject_pair_convention(inputs: G.Inputs):
    # phase_p1_montecarlo.py:353 约定：default_rng([seed, subject])
    _ = inputs  # fixture 仅用于模块级 load_inputs 预热
    fac = G.sample_factors(SEED, 3, G.COV_B_PRIMARY, G.COV_E_PRIMARY,
                           G.COV_A_PRIMARY, 9, mass_in_load=True)
    mb, sb = G.log_params(G.COV_B_PRIMARY)
    rng = np.random.default_rng([SEED, 0])
    assert fac["b"][0] == pytest.approx(rng.lognormal(mb, sb), rel=1e-12)


# --------------------------------------------------------------------------
# 6. 复用锁：logistic_fit 本体 + A3 规则原件
# --------------------------------------------------------------------------
def test_logistic_fit_is_reused_unchanged():
    assert G.pa.logistic_fit is pa.logistic_fit
    assert G.pa is pa


def test_logistic_fit_interface_and_or_definition():
    x = np.array([1.0, 2, 3, 4, 5, 6])
    y = np.array([0, 1, 0, 1, 1, 1])          # 非完全分离
    r = pa.logistic_fit(x, y)
    assert r["status"] == "ok"
    assert set(["OR", "beta", "se", "ci95", "n", "n_fracture"]).issubset(r)
    assert r["OR"] == pytest.approx(np.exp(r["beta"]), rel=1e-12)   # OR=exp(beta)
    assert len(r["ci95"]) == 2 and r["ci95"][0] <= r["OR"] <= r["ci95"][1]


def test_logistic_fit_degenerate_statuses():
    x = np.arange(1, 11, dtype=float)
    assert pa.logistic_fit(x, np.zeros(10, int))["status"] == "all_zero"
    assert pa.logistic_fit(x, np.ones(10, int))["status"] == "all_one"
    assert pa.logistic_fit(x, (x >= 5).astype(int))["status"] == "complete_separation"


def test_a3_targets_unchanged():
    assert G.A3_TARGETS == {"腓骨": 1.682, "骨盆": 1.236, "颅骨": 1.576}
    assert G.A3_TARGETS == dict(s4.A3_TARGETS)


def test_a3_rule_unchanged_30pct_2of3():
    # 2/3 落在 ±30% → PASS
    lg: dict[str, G.OrRecord] = {
        "腓骨": {"status": "ok", "OR": 1.682 * 1.10},
        "骨盆": {"status": "ok", "OR": 1.236 * 0.90},
        "颅骨": {"status": "all_one"},
    }
    v = G.a3_verdict(lg)
    assert v["n_within"] == 2 and v["verdict"] == "PASS"
    # 0/3 → FAIL
    lg2: dict[str, G.OrRecord] = {cn: {"status": "ok", "OR": 5.0}
                                  for cn in ("腓骨", "骨盆", "颅骨")}
    assert G.a3_verdict(lg2)["verdict"] == "FAIL"


def test_a3_verdict_reused_from_s4():
    # G3 的 a3_verdict 直接委托 s4 原件：用同一输入交叉核对结果一致
    lg: dict[str, G.OrRecord] = {"腓骨": {"status": "ok", "OR": 1.682},
                                 "骨盆": {"status": "ok", "OR": 1.236},
                                 "颅骨": {"status": "ok", "OR": 1.576}}
    assert G.a3_verdict(lg) == s4.a3_verdict(lg)
    assert G.s4.A3_TARGETS == s4.A3_TARGETS


# --------------------------------------------------------------------------
# 7. Wilson CI 健全性
# --------------------------------------------------------------------------
def test_wilson_ci_bounds():
    # k = 计数（n=400）
    lo, hi = G.wilson_ci(np.array([0.0, 200.0, 400.0]), 400.0)
    assert np.all(lo >= 0.0) and np.all(hi <= 1.0)
    assert lo[1] < 0.5 < hi[1]
    # k=0 / k=n 的区间仍合法（不在端点塌缩）
    assert 0.0 <= lo[0] <= hi[0] < 0.5
    assert 0.5 < lo[2] <= hi[2] <= 1.0
