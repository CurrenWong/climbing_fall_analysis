"""波3 · S5-contact 最小版跖面 BC（``plantar_bc="contact"``）结构回归。

不调用 FEBio（毫秒级），只锁死：

* ``"contact"`` 被 :data:`climbing.coupling.plantar_bc.BC_KINDS` 接受；
* 默认路径（``kind="fixed"``）**逐位不变**——不追加任何 surface / node_set / contact；
* ``kind="contact"`` 产出一个合法接触 deck（``Surface`` + ``SurfacePair`` +
  ``sliding-elastic`` ``<Contact>`` + ``<fric_coeff>`` + 刚性对偶面）；
* ``build_calcaneus_feb`` 的 ``plantar_bc`` 默认仍是 ``"fixed"``（不改签名）；
* ``build_thums_feb`` 的 ``plantar_bc`` 默认仍是 ``"fixed"``。

数值/收敛验收（G2/G3/G4/G5/G6）属 wave4（泡沫实体）与 ``pad_contact`` 集成，
不在本文件范围——契约 §4 条 8：最小版不在 G7 阻塞内。
"""
from __future__ import annotations

import inspect
import sys
from pathlib import Path

import numpy as np
import pytest

from climbing.coupling.plantar_bc import (
    BC_KINDS,
    DEFAULT_CONTACT_MU,
    DEFAULT_CONTACT_PENALTY,
    apply_plantar_bc,
)
from climbing.coupling.febio_model import build_calcaneus_feb

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts" / "opensim_fe") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts" / "opensim_fe"))


# ---------------------------------------------------------------------------
# 合成网格（与 test_nonvertical_s4 同款立方体：有效 tet，跖面在底面）
# ---------------------------------------------------------------------------
def _cube_mesh() -> dict:
    nodes = np.array(
        [
            [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
            [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],
        ],
        dtype=float,
    ) * 10.0
    tets = np.array([[0, 1, 2, 6], [0, 2, 3, 6], [0, 3, 7, 6], [0, 7, 4, 6], [0, 4, 5, 6]])
    return {
        "nodes": nodes,
        "tets": tets,
        "surfaces": {
            "subtalar_joint": np.array([[4, 5, 6], [4, 6, 7]]),
            "plantar": np.array([[0, 1, 2], [0, 2, 3]]),
        },
        "node_sets": {
            "subtalar_joint": np.array([4, 5, 6, 7]),
            "plantar": np.array([0, 1, 2, 3]),
        },
    }


def _minimal_model(mesh: dict):
    """一个只够 ``apply_plantar_bc`` 用的 pyfebio Model。"""
    from pyfebio import mesh as fmesh
    from pyfebio.model import Model

    model = Model()
    nodes = np.asarray(mesh["nodes"], dtype=np.float64)
    tets = np.asarray(mesh["tets"], dtype=np.int64)
    model.mesh_.add_node_domain(fmesh.numpy_to_nodes(nodes, name="calcaneus"))
    model.mesh_.add_element_domain(
        fmesh.numpy_to_elements(tets + 1, "tet4", name="calcaneus")
    )
    plantar = np.asarray(mesh["node_sets"]["plantar"], dtype=np.int64)
    model.mesh_.add_node_set(
        fmesh.NodeSet(name="plantar", text=",".join(map(str, plantar + 1)))
    )
    return model, nodes, plantar


# ==========================================================================
# 1. BC_KINDS / apply_plantar_bc —— opt-in 接受 + 默认不变
# ==========================================================================
def test_bc_kinds_includes_contact() -> None:
    assert "contact" in BC_KINDS
    # 追加，不是替换：四种历史 kind 仍在，且顺序不变
    assert BC_KINDS[:4] == ("fixed", "roller", "roller_free", "spring")
    assert BC_KINDS[-1] == "contact"


def test_apply_plantar_bc_contact_accepted() -> None:
    mesh = _cube_mesh()
    model, nodes, plantar = _minimal_model(mesh)
    bcs: list = []
    meta = apply_plantar_bc(
        model, nodes=nodes, plantar_local=plantar,
        plantar_tris=np.asarray(mesh["surfaces"]["plantar"]), bcs=bcs,
        kind="contact",
    )
    assert meta["kind"] == "contact"
    assert meta["n_plantar_nodes"] == int(len(plantar))
    assert meta["contact"] is not None
    assert meta["mu"] == DEFAULT_CONTACT_MU
    assert meta["penalty"] == DEFAULT_CONTACT_PENALTY
    # 1 个 pad 固定 BC（追加进 bcs）
    assert len(bcs) == 1
    # mesh 里多了 plantar Surface + pad_top Surface + SurfacePair
    surf_names = [s.name for s in model.mesh_.surfaces]
    assert "plantar" in surf_names
    assert "pad_top" in surf_names
    assert [p.name for p in model.mesh_.surface_pairs] == ["plantar_pad_pair"]
    # pad 固定节点集存在
    assert "pad_floor_fixed" in [ns.name for ns in model.mesh_.node_sets]


def test_apply_plantar_bc_default_fixed_unchanged() -> None:
    """默认 ``fixed``：只 append 1 个 BC，不加 surface / node_set / contact。"""
    mesh = _cube_mesh()
    model, nodes, plantar = _minimal_model(mesh)
    bcs: list = []
    meta = apply_plantar_bc(
        model, nodes=nodes, plantar_local=plantar,
        plantar_tris=np.asarray(mesh["surfaces"]["plantar"]), bcs=bcs,
    )
    assert meta == {"kind": "fixed", "n_plantar_nodes": int(len(plantar))}
    assert "contact" not in meta
    assert len(bcs) == 1
    assert model.mesh_.surfaces == []
    assert model.mesh_.surface_pairs == []


def test_apply_plantar_bc_contact_rejects_empty_tris() -> None:
    mesh = _cube_mesh()
    model, nodes, plantar = _minimal_model(mesh)
    with pytest.raises(ValueError):
        apply_plantar_bc(
            model, nodes=nodes, plantar_local=plantar,
            plantar_tris=np.empty((0, 3), dtype=np.int64), bcs=[],
            kind="contact",
        )


def test_apply_plantar_bc_contact_recipe_kwargs_flow_through() -> None:
    """opt-in 配方 kwargs 真的改写 contact 块。"""
    mesh = _cube_mesh()
    model, nodes, plantar = _minimal_model(mesh)
    meta = apply_plantar_bc(
        model, nodes=nodes, plantar_local=plantar,
        plantar_tris=np.asarray(mesh["surfaces"]["plantar"]), bcs=[],
        kind="contact", contact_mu=0.3, contact_penalty=0.5,
        contact_search_radius=12.0, contact_tolerance=0.02,
    )
    assert meta["mu"] == pytest.approx(0.3)
    assert meta["penalty"] == pytest.approx(0.5)
    assert meta["search_radius"] == pytest.approx(12.0)
    assert meta["tolerance"] == pytest.approx(0.02)
    cc = meta["contact"].all_contact_interfaces[0]
    assert cc.fric_coeff == pytest.approx(0.3)
    assert cc.penalty == pytest.approx(0.5)


# ==========================================================================
# 2. build_calcaneus_feb —— 发射的 deck 含 contact 块；默认不含
# ==========================================================================
class TestCalcaneusDeck:
    def test_default_has_no_contact_block(self, tmp_path) -> None:
        out = build_calcaneus_feb(_cube_mesh(), tmp_path / "fixed.feb")
        text = out.read_text(encoding="ISO-8859-1")
        assert "<Contact>" not in text
        assert "sliding-elastic" not in text

    def test_contact_emits_contact_block(self, tmp_path) -> None:
        # pressure path (use_rigid=False) emits an explicit <Loads> section, so the
        # required Loads → Boundary → Contact order is exercised end-to-end.
        out = build_calcaneus_feb(
            _cube_mesh(), tmp_path / "contact.feb",
            plantar_bc="contact", use_rigid=False,
        )
        text = out.read_text(encoding="ISO-8859-1")
        assert "<Contact>" in text
        assert "<SurfacePair" in text
        assert 'type="sliding-elastic"' in text
        assert f"<fric_coeff>{DEFAULT_CONTACT_MU}</fric_coeff>" in text
        assert f"<penalty>{DEFAULT_CONTACT_PENALTY}</penalty>" in text
        # 段序（契约 §5）：Contact 必须晚于 Loads / Boundary。
        # （step 内 pyfebio 序列化为 Boundary → Loads → Contact；FEBio 4.13 接受，
        # 既有 calcaneus deck 亦如此，见 test_opensim_fe/test_febio_run。）
        assert text.index("<Loads>") < text.index("<Contact>")
        assert text.index("<Boundary>") < text.index("<Contact>")
        # 对偶面是真实单元面（pad_top），不是自由四边形
        assert 'name="pad_top"' in text

    def test_contact_emits_after_boundary_rigid_path(self, tmp_path) -> None:
        """幽灵刚体路径：contact 仍在 Boundary 之后（Loads 由 Rigid 承载）。"""
        out = build_calcaneus_feb(
            _cube_mesh(), tmp_path / "contact_rigid.feb", plantar_bc="contact"
        )
        text = out.read_text(encoding="ISO-8859-1")
        assert text.index("<Boundary>") < text.index("<Contact>")

    def test_contact_deck_parses_section_order(self, tmp_path) -> None:
        """SurfacePair 必须晚于最后一个 Surface（FEBio User Manual §3.6.8）。"""
        out = build_calcaneus_feb(
            _cube_mesh(), tmp_path / "contact_order.feb", plantar_bc="contact"
        )
        text = out.read_text(encoding="ISO-8859-1")
        assert text.rindex("</Surface>") < text.index("<SurfacePair")


# ==========================================================================
# 3. 签名默认值 —— 默认逐位复现旧行为
# ==========================================================================
def test_calcaneus_plantar_bc_default_fixed() -> None:
    sig = inspect.signature(build_calcaneus_feb)
    assert sig.parameters["plantar_bc"].default == "fixed"


def test_thums_plantar_bc_default_fixed() -> None:
    import thums_feb

    sig = inspect.signature(thums_feb.build_thums_feb)
    assert sig.parameters["plantar_bc"].default == "fixed"
