"""S5-contact 回归测试：``pad_contact``（垫-压头 滑动接触 + 库仑摩擦）。

两组用例：

1. **纯 deck 结构**（不需要 FEBio）—— :func:`build_pad_contact_feb` 产出的
   ``.feb`` 必须满足契约 §5 的硬规则：两个固体域、真实 tri3 facet 的
   ``<Surface>``、``<SurfacePair>`` 在 ``<Surface>`` 之后、段序
   ``Loads -> Boundary -> Contact``、含 ``<fric_coeff>``、**不含**
   ``<var type="contact traction"/>``；以及 ``run_pad_contact`` 必需的 API。
2. **真实 FEBio**（缺则 skip）—— 跑一个小接触算例，断言：

   * G3 守恒：接触面反力合力 == 驱动面反力（≤2%）；
   * G4 摩擦锥：``|F_t| <= mu*F_n``（所有 output state）；
   * 接触确实发生（indentation > 0）。

依赖 FEBio 的用例统一用 ``requires_febio`` 跳过（与 ``tests/test_febio_run.py`` 同款）。
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from climbing.coupling.febio_run import FebioNotFound, find_febio
from climbing.coupling.pad_contact import (
    DEFAULT_MU,
    _combine_domains,
    _fast_meshes,
    build_pad_contact_feb,
    indenter_block_tet,
    run_pad_contact,
    tri_area,
)
from climbing.coupling.pad_foam import fit_pad_material
from climbing.coupling.pad_mesh import pad_block_tet


def _febio_available() -> bool:
    try:
        find_febio(verify=True)
        return True
    except FebioNotFound:
        return False


requires_febio = pytest.mark.skipif(
    not _febio_available(), reason="本机未安装可运行的 FEBio"
)

# 相对误差阈值（G3）
REL_TOL = 0.02


# ---------------------------------------------------------------------------
# 1. deck 结构（不依赖 FEBio）
# ---------------------------------------------------------------------------
class TestDeckStructure:
    @pytest.fixture
    def deck_text(self, tmp_path):
        pad = pad_block_tet(3, 3, 2, (300.0, 300.0, 200.0))
        ind = indenter_block_tet(2, 2, 1, (200.0, 300.0, 5.0),
                                 offset_mm=(50.0, 0.0, 200.5))
        out = build_pad_contact_feb(tmp_path / "c.feb", pad_mesh=pad, indenter_mesh=ind,
                                    mu=0.6, penalty=1000.0,
                                    load={"mode": "displacement", "sinkage_mm": 10.0})
        assert isinstance(out, Path) and out.is_file()
        return out.read_text(encoding="latin-1")

    def test_two_solid_domains(self, deck_text):
        assert deck_text.count("<SolidDomain") == 2
        assert 'name="pad_tets"' in deck_text
        assert 'name="indenter_tets"' in deck_text

    def test_section_order_loads_boundary_contact(self, deck_text):
        """契约 §5：Loads -> Boundary -> Contact（顺序错 = unrecognized tag）。"""
        i_loads = deck_text.index("<Loads")
        i_bnd = deck_text.index("<Boundary")
        i_ct = deck_text.index("<Contact>")
        assert i_loads < i_bnd < i_ct

    def test_surfacepair_after_surface(self, deck_text):
        """契约：SurfacePair 必须定义在 Surface 之后（User Manual §3.6.8）。"""
        last_surf = deck_text.rindex("<Surface ")
        pair = deck_text.index("<SurfacePair")
        assert pair > last_surf
        # 关标签也要在 pair 之前
        assert deck_text.rindex("</Surface>") < pair

    def test_real_tri_facets(self, deck_text):
        """接触面必须是真实 tri3 facet（不是 quad / 合成四边形）。"""
        assert "<tri3" in deck_text
        # 接触面片数 = pad top 18 + indenter top 18 + indenter bottom 8
        assert deck_text.count("<tri3") >= 30

    def test_friction_coefficient_present(self, deck_text):
        assert "<fric_coeff>0.6</fric_coeff>" in deck_text
        assert 'type="sliding-elastic"' in deck_text

    def test_no_contact_traction_var(self, deck_text):
        """★ 不能出现 contact traction —— 它会让 pyfebio 的 xplt->hdf5 整体报废。"""
        assert "contact traction" not in deck_text.lower()

    def test_plotfile_has_reaction_forces(self, deck_text):
        assert "reaction forces" in deck_text
        assert 'var type="displacement"' in deck_text

    def test_material_types(self, deck_text):
        assert 'type="Ogden"' in deck_text          # 垫
        assert 'type="isotropic elastic"' in deck_text  # 压头


class TestApiSurface:
    def test_exports_present(self):
        import climbing.coupling.pad_contact as pc

        for name in ("build_pad_contact_feb", "run_pad_contact", "indenter_block_tet"):
            assert name in pc.__all__

    def test_indenter_block_tet_offset(self):
        ind = indenter_block_tet(1, 1, 1, (10.0, 10.0, 5.0),
                                 offset_mm=(5.0, 0.0, 100.0))
        assert ind["nodes"][:, 2].min() == pytest.approx(100.0)
        assert ind["nodes"][:, 2].max() == pytest.approx(105.0)
        # surfaces schema
        assert set(ind["surfaces"]) == {"top", "bottom", "sides"}

    def test_load_validation(self, tmp_path):
        pad = pad_block_tet(1, 1, 1, (10.0, 10.0, 10.0))
        ind = indenter_block_tet(1, 1, 1, (10.0, 10.0, 5.0), offset_mm=(0.0, 0.0, 10.5))
        with pytest.raises(ValueError):
            build_pad_contact_feb(tmp_path / "bad.feb", pad_mesh=pad, indenter_mesh=ind,
                                  load={"mode": "displacement"})  # 缺 sinkage_mm
        with pytest.raises(ValueError):
            build_pad_contact_feb(tmp_path / "bad2.feb", pad_mesh=pad, indenter_mesh=ind,
                                  mu=-1.0)

    def test_full_domain_index_accounting(self):
        pad, ind = _fast_meshes(pad_shape=(3, 3, 2), ind_shape=(2, 2, 1))
        comb = _combine_domains(pad, ind)
        assert comb["n_pad"] == len(pad["nodes"])
        assert comb["tets_b"].min() >= comb["n_pad"]
        # 接触面积 = 压头底面 200x300
        assert tri_area(comb["nodes"], comb["ind_surfaces"]["top"]) == pytest.approx(60000.0)


# ---------------------------------------------------------------------------
# 2. 真实 FEBio：G3 + G4
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def result(tmp_path_factory):
    if not _febio_available():
        pytest.skip("本机未安装可运行的 FEBio")
    workdir = tmp_path_factory.mktemp("pad_contact_fe")
    pad, ind = _fast_meshes(pad_shape=(3, 3, 2), ind_shape=(2, 2, 1))
    # Convergence regime for the fitted-Ogden + tet4 contact: penalty=1,
    # auto_penalty=0 (enforced in ``pad_contact.py::_emit_contact``), and
    # ``n_solver_steps ≥ 1000``.  See
    # ``temp/opensim_fe/s5_t7/probe28_pensweep.py`` for the convergence map.
    return run_pad_contact(
        workdir / "c.feb", pad_mesh=pad, indenter_mesh=ind,
        mu=DEFAULT_MU, pad_mat=fit_pad_material(), penalty=1.0,
        load={"mode": "displacement", "sinkage_mm": 30.0, "lateral_mm": (0.3, 0.0)},
        n_solver_steps=1000,
    )


@requires_febio
class TestContactFE:
    def test_runs_rc_zero(self, result):
        assert result["febio_rc"] == 0
        assert result["xplt_path"].is_file()
        assert result["n_states"] >= 5

    def test_contact_engaged(self, result):
        """接触必须真的发生（压入量 > 0）。"""
        assert result["indentation_mm"] > 5.0, result["indentation_mm"]

    def test_g3_conservation(self, result):
        """G3：接触面反力合力 == 驱动面反力（≤2%）。"""
        f_contact = abs(float(result["pad_top_reaction"][-1, 2]))
        f_drive = abs(float(result["indenter_reaction"][-1, 2]))
        assert f_drive > 0
        rel = abs(f_contact - f_drive) / f_drive
        assert rel <= REL_TOL, (
            f"G3 守恒失败：contact={f_contact:.4e} drive={f_drive:.4e} rel={rel:.4f}"
        )

    def test_g4_friction_cone(self, result):
        """G4：任意 output state 上 |F_t| <= mu*F_n。"""
        assert result["friction_ok"], result["f_n_t_n"]
        f_n, f_t = result["f_n_t_n"].T
        assert np.all(f_t <= DEFAULT_MU * f_n + 1e-6 * np.maximum(f_n, 1.0))
        # 摩擦必须真的被激励（侧向位移 0.3mm）
        assert result["max_f_t_over_mu_f_n"] > 0.01

    def test_boundary_gap_positive_required_note(self):
        """文档性：0 gap 不收敛 → 模块默认 gap>0（这里只断言常量）。"""
        from climbing.coupling.pad_contact import DEFAULT_GAP_MM

        assert DEFAULT_GAP_MM > 0.0


# ---------------------------------------------------------------------------
# 3. AUGLAG 增广控制（收敛配方，见 docs/S5_contact_convergence.md）
# ---------------------------------------------------------------------------
class TestAuglagControls:
    @staticmethod
    def _deck(tmp_path, **kw):
        pad = pad_block_tet(3, 3, 2, (300.0, 300.0, 200.0))
        ind = indenter_block_tet(2, 2, 1, (200.0, 300.0, 5.0),
                                 offset_mm=(50.0, 0.0, 200.5))
        out = build_pad_contact_feb(
            tmp_path / "a.feb", pad_mesh=pad, indenter_mesh=ind, mu=0.6,
            penalty=1.0, load={"mode": "displacement", "sinkage_mm": 10.0}, **kw)
        return out.read_text(encoding="latin-1")

    def test_default_emits_no_aug_tags(self, tmp_path):
        """不传 aug_controls ⇒ 输出与旧版逐字一致（无增广子标签）。"""
        txt = self._deck(tmp_path, laugon="AUGLAG")
        block = txt[txt.index("<contact"):txt.index("</contact>")]
        for tag in ("gaptol", "minaug", "maxaug", "smooth_aug"):
            assert f"<{tag}>" not in block

    def test_convergent_controls_emitted(self, tmp_path):
        from climbing.coupling.pad_contact import AUGLAG_CONVERGENT_CONTROLS

        txt = self._deck(tmp_path, laugon="AUGLAG", tolerance=0.005,
                         aug_controls=AUGLAG_CONVERGENT_CONTROLS)
        block = txt[txt.index("<contact"):txt.index("</contact>")]
        assert block.count("<tolerance>") == 1
        assert "<tolerance>0.005</tolerance>" in block
        assert "<gaptol>0.001</gaptol>" in block
        assert "<minaug>3</minaug>" in block
        assert "<maxaug>200</maxaug>" in block
        assert "<smooth_aug>1</smooth_aug>" in block

    def test_aug_controls_rejected_for_penalty_laugon(self, tmp_path):
        from climbing.coupling.pad_contact import AUGLAG_CONVERGENT_CONTROLS

        with pytest.raises(ValueError, match="AUGLAG"):
            self._deck(tmp_path, laugon="PENALTY",
                       aug_controls=AUGLAG_CONVERGENT_CONTROLS)
