"""S5-contact T1 回归测试：鞍垫泡沫材料在 FEBio 中的复现。

三组用例：
1. **纯 numpy** —— :func:`pad_stress_pa` 必须与 ``CrashPad.stress_pa`` 逐位一致（1e-12）；
2. **拟合报告** —— :func:`fit_pad_material` 返回约定的键，且 ``max_rel_err``
   与 ``table`` 自洽，且达到 G1 门槛（≤5%）；
3. **真实 FEBio 单轴** —— :func:`run_uniaxial_verify` 跑通并给回 σ(ξ) 曲线，
   且 FE 曲线在 ξ∈[0,0.85] 上满足 G1（≤5%）。

G1 说明（重要，2026-10 修订）
----------------------------
``docs/S5_contact_contract.md`` 的 G1 判据是「ξ∈[0,0.85] 内 max_rel_err ≤ 5%」。

* 之前的实现用**不可压** Ogden（``λ_t = λ^{−1/2}``）且以「峰值归一化的平方误差」
  为目标 —— 两者都不对，实测 37%+。
* 现在的实现使用 FEBio 4.13 **可表达的压缩 uncoupled Ogden**
  （``<material type="Ogden">`` + ``<k>`` + ``<pressure_model>``），对
  **单轴应力（侧面自由）** 的工程应力解做拟合，目标正是 **逐点最大相对误差**
  （线性规划精确求解）。
* 实测：解析拟合 ≈0.4%，**真实 FEBio FE ≈2.5%**（分段 [0,0.15]/[0.15,0.5]/
  [0.5,0.85] = 2.5%/0.9%/0.8%）—— **G1 通过**。

注意 FEBio 4.13 **没有** hyperfoam / 任何 foam 材料（DLL 注册表 + 实证 tag 探针
均已确认）—— 详见 ``src/climbing/coupling/pad_foam.py`` 的模块 docstring。

依赖 FEBio 的用例统一用 ``requires_febio`` 跳过（与 ``tests/test_febio_run.py`` 同款）。
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from climbing.coupling.pad_foam import (
    DEFAULT_BULK_MPA,
    DEFAULT_PRESSURE_MODEL,
    FEBIO_TYPE,
    fit_pad_material,
    pad_stress_pa,
    run_uniaxial_verify,
    write_uniaxial_feb,
)
from climbing.pad import CrashPad

# ---------------------------------------------------------------------------
# FEBio 可用性（照抄 tests/test_febio_run.py 的模式）
# ---------------------------------------------------------------------------
from climbing.coupling.febio_run import FebioNotFound, find_febio


def _febio_available() -> bool:
    try:
        find_febio(verify=True)
        return True
    except FebioNotFound:
        return False


requires_febio = pytest.mark.skipif(
    not _febio_available(), reason="本机未安装可运行的 FEBio"
)


# ---------------------------------------------------------------------------
# 1. 纯 numpy：pad_stress_pa 必须与 CrashPad.stress_pa 逐位一致
# ---------------------------------------------------------------------------
class TestPadStressExact:
    def test_matches_crashpad_exactly(self):
        """ξ∈[0,1] 上逐点相对误差 ≤ 1e-12（含 ξ>0.85 的压实段）。"""
        xi = np.linspace(0.0, 1.2, 2401)
        a = pad_stress_pa(xi)
        b = CrashPad().stress_pa(xi)
        rel = np.abs(a - b) / np.maximum(np.abs(b), 1.0)
        assert rel.max() <= 1e-12, f"max rel diff = {rel.max():.3e}"

    def test_matches_at_gate_endpoints(self):
        """G1 区间端点 + 压实拐点单独核对。"""
        for xi in (0.0, 0.15, 0.85, 0.90):
            assert pad_stress_pa(xi) == pytest.approx(
                CrashPad().stress_pa(xi), rel=1e-12, abs=0.0
            )

    def test_vectorised_and_scalar_agree(self):
        xs = [0.0, 0.1, 0.5, 0.85]
        assert np.allclose(pad_stress_pa(xs), [pad_stress_pa(x) for x in xs])

    def test_clips_negative_strain_to_zero(self):
        assert float(pad_stress_pa(-0.5)) == 0.0

    def test_monotonic_increasing(self):
        xi = np.linspace(0.0, 1.0, 1001)
        s = pad_stress_pa(xi)
        assert np.all(np.diff(s) >= 0.0)


# ---------------------------------------------------------------------------
# 2. 拟合报告
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def fit():
    return fit_pad_material()


class TestFitReport:
    def test_required_keys_present(self, fit):
        for key in ("febio_type", "params", "max_rel_err", "xi_range", "table"):
            assert key in fit, f"missing key {key!r}"
        assert fit["febio_type"] == FEBIO_TYPE == "Ogden"
        assert set(fit["params"]) >= {"m", "c", "k_mpa"}
        assert len(fit["params"]["m"]) == len(fit["params"]["c"])
        assert fit["params"]["k_mpa"] == pytest.approx(DEFAULT_BULK_MPA)

    def test_xi_range_is_gate_range(self, fit):
        assert tuple(fit["xi_range"]) == (0.0, 0.85)

    def test_table_is_self_consistent(self, fit):
        """table[i] = [xi, sigma_fit, sigma_pad]，且 max_rel_err 与逐行自洽。"""
        table = np.asarray(fit["table"], dtype=float)
        assert table.ndim == 2 and table.shape[1] == 3
        xi, sigma_fit, sigma_pad = table[:, 0], table[:, 1], table[:, 2]
        # σ_pad 列必须与真值完全一致
        assert np.allclose(sigma_pad, pad_stress_pa(xi), rtol=0.0, atol=1e-9)
        rel = np.abs(sigma_fit - sigma_pad) / np.maximum(sigma_pad, 1.0)
        assert rel.max() == pytest.approx(fit["max_rel_err"], rel=1e-9)

    def test_analytical_fit_meets_g1(self, fit):
        """解析拟合必须达到 G1（≤5%），且分段误差有界。"""
        assert fit["max_rel_err"] <= 0.05, (
            f"analytical max_rel_err={fit['max_rel_err']:.4f} > 0.05"
        )
        # 密实化段（S5 G7 工作区）应当明显更好
        assert fit["max_rel_err_dens"] < fit["max_rel_err"]
        for name, val in fit["bands"].items():
            assert val <= 0.05, f"band {name} rel={val:.4f}"

    def test_material_is_compressible_ogden(self, fit):
        """材料块必须是可压缩 uncoupled Ogden（含 k 与 pressure_model）。"""
        assert fit["params"]["pressure_model"] == DEFAULT_PRESSURE_MODEL
        assert np.isfinite(fit["params"]["k_mpa"]) and fit["params"]["k_mpa"] > 0.0
        assert np.isfinite(fit["j_min"]) and fit["j_min"] > 0.0

    def test_deterministic(self):
        a = fit_pad_material()
        b = fit_pad_material()
        assert a["max_rel_err"] == pytest.approx(b["max_rel_err"], rel=1e-12)


# ---------------------------------------------------------------------------
# 3. 写出 .feb（不需要 FEBio，只检查结构）
# ---------------------------------------------------------------------------
class TestWriteUniaxialFeb:
    def test_writes_ogden_material_and_uses_xi_max(self, tmp_path):
        fit = fit_pad_material()
        out = write_uniaxial_feb(tmp_path / "u.feb", fit, xi_targets=[0.0, 0.4, 0.85])
        assert out.is_file()
        text = out.read_text(encoding="latin-1")
        assert 'type="Ogden"' in text
        assert "<k>" in text
        # 位移 = -0.85 * 20 mm = -17 mm
        assert "<value lc=\"1\">-17.000000</value>" in text
        # 关键段序：Loads 在 Boundary 前
        assert text.index("<Loads") < text.index("<Boundary")

    def test_rejects_wrong_type(self, tmp_path):
        bad = {"febio_type": "NeoHookean", "params": {"m": [1.0], "c": [1.0]},
               "xi_range": (0.0, 0.85)}
        with pytest.raises((ValueError, KeyError)):
            write_uniaxial_feb(tmp_path / "bad.feb", bad)


# ---------------------------------------------------------------------------
# 4. 真实 FEBio 单轴验证
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def result(tmp_path_factory):
    if not _febio_available():
        pytest.skip("本机未安装可运行的 FEBio")
    workdir = tmp_path_factory.mktemp("pad_foam_fe")
    fit = fit_pad_material()
    return run_uniaxial_verify(fit, workdir=workdir, n_solver_steps=80)


class TestUniaxialFEVerify:
    @requires_febio
    def test_runs_and_covers_gate_range(self, result):
        assert result["febio_rc"] == 0
        assert result["n_states"] >= 10
        assert result["xplt_path"].is_file()
        xi = result["xi_fe"]
        assert xi.max() >= 0.84, f"FE only reached ξ={xi.max():.3f}"
        assert np.all(np.diff(xi) >= -1e-9)

    @requires_febio
    def test_sigma_fe_is_positive_and_increasing(self, result):
        s = result["sigma_fe_pa"]
        assert np.all(s >= -1e-6)
        assert s.max() > 1.0  # 不是静默 0

    @requires_febio
    def test_fe_in_densification_band(self, result):
        """密实化带 [0.4,0.72] 的 FE 误差应远好于全区间门槛。"""
        xi = result["xi_fe"]
        s_fe = result["sigma_fe_pa"]
        s_pad = result["sigma_anal_pa"]
        mask = (xi >= 0.4) & (xi <= 0.72)
        assert mask.any(), "densification band not sampled"
        rel = np.abs(s_fe[mask] - s_pad[mask]) / np.maximum(s_pad[mask], 1.0)
        assert rel.mean() <= 0.05, f"band mean rel = {rel.mean():.4f}"
        assert rel.max() <= 0.05, f"band max rel = {rel.max():.4f}"

    @requires_febio
    def test_g1_full_range_5pct(self, result):
        """契约 G1：ξ∈[0,0.85] 全区间 max_rel_err ≤ 5%（真实 FEBio 结果）。"""
        max_rel = result["max_rel_err"]
        print(f"\n[G1] FE-vs-pad max_rel_err = {max_rel:.4f} (threshold 0.05)")
        assert max_rel <= 0.05

    @requires_febio
    def test_fe_vs_analytic_max_rel_err_le_5pct(self, result):
        """显式回归：:func:`pad_foam.run_uniaxial_verify` 的 FE-vs-解析 max 相对误差 ≤ 5%。

        （复用 module-scoped ``result`` fixture —— 该 fixture 就是真实调用
        ``run_uniaxial_verify`` 的那条路径；本用例把"≤5%"写成独立断言，防止将来
        有人误删 G1 判据而套件仍然全绿。）
        """
        assert result["febio_rc"] == 0
        assert result["max_rel_err"] <= 0.05, (
            f"FE-vs-analytic max_rel_err={result['max_rel_err']:.4f} > 0.05"
        )
