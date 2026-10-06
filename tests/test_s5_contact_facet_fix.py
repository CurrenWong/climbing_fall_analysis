"""S5-contact facet-type fix (Defect A) + G7 solver recipe (Defect B) — 结构回归。

不调用 FEBio（毫秒级），只锁死：

**Defect A — 跖面 contact surface 必须是 quad4**
- 真实 THUMS 网格的跖面属于 CORT hex8 域（191/191 节点在 CORT），
  而 ``mesh['surfaces']['plantar']`` 把 quad4 切成 292 tri3 ⇒ 挂在 hex8 域上
  FEBio 报 ``292 invalid facets`` ⇒ 接触不承载。详见
  ``docs/S5_contact_facet_fix.md`` §1。
- 新增 ``plantar_quads`` opt-in kwarg：传入 ``(K, 4)`` ⇒ 跖面 Surface 发射
  ``<quad4>``；未传或 ``None`` ⇒ 走 tri3 兼容路径（合成 cube + 旧版测试）。
- **防静默退化**：传入 ``plantar_quads`` 但行数为 0 ⇒ 显式 ``ValueError``，
  绝不回退到 tri3。

**Defect B — 求解器 / 加载配置必须用 G7 配方**
- G7 配方（``docs/S5_contact_udg.md`` §7 实测验证）：
  ``time_steps=2400`` / ``step_size=1/2400`` / ``cutback=0.125`` /
  ``max_retries=20`` / ``opt_iter=15`` / BFGS ``qn_method max_ups=10`` /
  ``reform_each_time_step=1`` / load-curve ramp。
- 新增 opt-in kwargs ``g7_solver_recipe: bool=False`` 与
  ``load_ramp_time_s: float|None=None``；默认 ⇒ 字节级保留旧版控制块。
"""
from __future__ import annotations

import inspect
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts" / "opensim_fe") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts" / "opensim_fe"))


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _cube_mesh() -> dict:
    """tet4 cube (与 test_plantar_contact / test_nonvertical_s4 同款)。"""
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


def _fake_thums_mesh() -> dict:
    """最小 THUMS 双域网格 (CORT hex8 + SPON tet4)；用 ``_build_thums_plantar_quads``
    重建跖面 quad4 facet。"""
    nodes = np.array(
        [
            [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],   # 0..3  bottom (y=0)
            [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],   # 4..7  top (y=1)
        ],
        dtype=float,
    )
    cort = np.array([[0, 1, 2, 3, 4, 5, 6, 7]], dtype=np.int64)        # 1 hex8
    spon = np.array([[0, 1, 2, 5]], dtype=np.int64)                   # 1 tet4 (interior)
    # boundary polys: bottom quad4 (0,1,2,3) on cort [hex bottom face],
    # plus side tris (on cort [hex sides]) + top quad (on cort [hex top face]).
    # The bottom quad4 has normal.y == -1.0 (plantar).
    polys = [
        (0, 1, 2, 3),    # bottom (cort, owner idx=0)   plantar
        (4, 5, 6, 7),    # top    (cort, owner idx=0)   subtalar joint
        (0, 1, 5, 4),    # front  (cort, owner idx=0)
        (1, 2, 6, 5),    # right  (cort, owner idx=0)
        (2, 3, 7, 6),    # back   (cort, owner idx=0)
        (3, 0, 4, 7),    # left   (cort, owner idx=0)
    ]
    owners = np.array([0, 0, 0, 0, 0, 0], dtype=np.int64)
    normals = np.array([
        [0, -1, 0],   # bottom
        [0,  1, 0],   # top
        [0,  0, -1],  # front
        [1,  0, 0],   # right
        [0,  0,  1],  # back
        [-1, 0, 0],   # left
    ], dtype=np.float64)
    return {
        "nodes": nodes,
        "nodes_global": np.arange(1, 9, dtype=np.int64),
        "elements_cort": cort,
        "elements_spon": spon,
        "elements": [cort[0], spon[0]],
        "surfaces": {
            # top quad4 split into 2 tri3 (the pressure path needs tri3 for joint_area_mm2)
            "subtalar_joint": np.array([[4, 5, 6], [4, 6, 7]], dtype=np.int64),
            "plantar": np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int64),
            "achilles": np.empty((0, 3), dtype=np.int64),
        },
        "node_sets": {
            "subtalar_joint": np.array([4, 5, 6, 7], dtype=np.int64),
            "plantar": np.array([0, 1, 2, 3], dtype=np.int64),
            "achilles": np.empty(0, dtype=np.int64),
        },
        "boundary_polys": polys,
        "boundary_owners": owners,
        "boundary_normals": normals,
        "boundary_centroids": np.array(
            [[0.5, 0.0, 0.5], [0.5, 1.0, 0.5], [0.5, 0.5, 0.0],
             [1.0, 0.5, 0.5], [0.5, 0.5, 1.0], [0.0, 0.5, 0.5]],
            dtype=np.float64,
        ),
        "boundary_faces": np.array(
            [[0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7],
             [0, 1, 5], [0, 5, 4], [1, 2, 6], [1, 6, 5],
             [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7]],
            dtype=np.int64,
        ),
        "element_centroids": np.array([[0.5, 0.5, 0.5], [0.5, 0.25, 0.75]]),
        "joint_center_mm": (0.0, 0.0, 0.0),
        "n_nodes": 8, "n_elements_cort": 1, "n_elements_spon": 1,
        "n_nodes_cort": 8, "n_nodes_spon": 4,
        "volume_cort_mm3": 1.0, "volume_spon_mm3": 0.125,
        "n_boundary_faces": 6, "n_boundary_tris": 12,
        "n_joint_faces": 1, "n_plantar_faces": 1, "n_achilles_faces": 0,
        "n_joint_tris": 2, "n_plantar_tris": 2, "n_achilles_tris": 0,
        "n_joint_nodes": 4, "n_plantar_nodes": 4, "n_achilles_nodes": 0,
    }


# ===========================================================================
# 1. Defect A — 跖面 Surface 必须是 quad4（真实 THUMS 路径）
# ===========================================================================
def test_thums_plantar_quads_rebuild_from_boundary_polys() -> None:
    """``_build_thums_plantar_quads`` 从 boundary_polys 重建跖面 quad4。

    191/191 跖面节点在 CORT hex8 域（实测），所以跖面命中 hex 面的 4-tuple。
    """
    import thums_feb as tf
    quads = tf._build_thums_plantar_quads(_fake_thums_mesh())
    assert quads.ndim == 2 and quads.shape[1] == 4
    # fake mesh has exactly one hex8 bottom face → exactly 1 quad4
    assert len(quads) == 1
    np.testing.assert_array_equal(quads[0], [0, 1, 2, 3])


def test_thums_plantar_quads_rejects_non_thums_mesh() -> None:
    """非 THUMS 双域网格（缺 elements_cort 或 boundary_polys）显式抛错。"""
    import thums_feb as tf
    cube = _cube_mesh()
    with pytest.raises(ValueError, match="不是 THUMS 两域"):
        tf._build_thums_plantar_quads(cube)


def test_apply_plantar_bc_quad4_path_emits_quad4_not_tri3() -> None:
    """``plantar_quads`` opt-in ⇒ 跖面 Surface 含 ``<quad4>``，不含 ``<tri3>``。"""
    from climbing.coupling.plantar_bc import apply_plantar_bc

    mesh = _cube_mesh()
    model, nodes, plantar = _minimal_model(mesh)
    quads = np.array([[0, 1, 2, 3]], dtype=np.int64)
    meta = apply_plantar_bc(
        model, nodes=nodes, plantar_local=plantar,
        plantar_tris=mesh["surfaces"]["plantar"],
        bcs=[],
        kind="contact",
        contact_plantar_quads=quads,
    )
    assert meta["plantar_facet_kind"] == "quad4"
    assert meta["n_surface_quads"] == 1
    assert meta["n_surface_tris"] == 2  # tri3 still recorded for back-compat
    # 只看 plantar Surface（不是 pad_top）
    surf_names = [s.name for s in model.mesh_.surfaces]
    assert "plantar" in surf_names
    plantar_surf = next(s for s in model.mesh_.surfaces if s.name == "plantar")
    # pyfebio Surface 暴露 all_quad4 / all_tri3 / all_tri6 列表
    assert len(plantar_surf.all_quad4) == 1, \
        f"plantar Surface 应有 1 个 quad4，实际 {len(plantar_surf.all_quad4)}"
    assert len(plantar_surf.all_tri3) == 0, \
        f"plantar Surface 应有 0 个 tri3，实际 {len(plantar_surf.all_tri3)}"


def test_apply_plantar_bc_quad4_rejects_empty_quads() -> None:
    """防静默退化：``plantar_quads`` 非 None 但行数 0 ⇒ 显式 ``ValueError``。"""
    from climbing.coupling.plantar_bc import apply_plantar_bc

    mesh = _cube_mesh()
    model, nodes, plantar = _minimal_model(mesh)
    with pytest.raises(ValueError, match="plantar_quads 为空"):
        apply_plantar_bc(
            model, nodes=nodes, plantar_local=plantar,
            plantar_tris=mesh["surfaces"]["plantar"],
            bcs=[],
            kind="contact",
            contact_plantar_quads=np.empty((0, 4), dtype=np.int64),
        )


def test_apply_plantar_bc_quad4_rejects_wrong_shape() -> None:
    """``plantar_quads`` 形状不对 → ``ValueError``。"""
    from climbing.coupling.plantar_bc import apply_plantar_bc

    mesh = _cube_mesh()
    model, nodes, plantar = _minimal_model(mesh)
    with pytest.raises(ValueError, match="(K,4)"):
        apply_plantar_bc(
            model, nodes=nodes, plantar_local=plantar,
            plantar_tris=mesh["surfaces"]["plantar"],
            bcs=[],
            kind="contact",
            contact_plantar_quads=np.array([[0, 1, 2]], dtype=np.int64),  # 3 列，不是 4
        )


def test_apply_plantar_bc_no_quads_falls_back_to_tri3() -> None:
    """``plantar_quads=None``（默认）⇒ 走 tri3 路径（合成 / 旧版兼容）。"""
    from climbing.coupling.plantar_bc import apply_plantar_bc

    mesh = _cube_mesh()
    model, nodes, plantar = _minimal_model(mesh)
    meta = apply_plantar_bc(
        model, nodes=nodes, plantar_local=plantar,
        plantar_tris=mesh["surfaces"]["plantar"],
        bcs=[],
        kind="contact",
    )
    assert meta["plantar_facet_kind"] == "tri3"
    assert meta["n_surface_quads"] == 0
    assert meta["n_surface_tris"] == 2


# ===========================================================================
# 2. Defect B — G7 求解器 / load-curve 配方
# ===========================================================================
def test_thums_feb_signature_g7_default_false() -> None:
    """``g7_solver_recipe`` 默认 ``False`` ⇒ G0 字节级隔离。"""
    import thums_feb as tf
    sig = inspect.signature(tf.build_thums_feb)
    assert sig.parameters["g7_solver_recipe"].default is False
    assert sig.parameters["load_ramp_time_s"].default is None


def test_thums_feb_default_fixed_unchanged() -> None:
    """``plantar_bc`` 默认仍是 ``"fixed"``。"""
    import thums_feb as tf
    sig = inspect.signature(tf.build_thums_feb)
    assert sig.parameters["plantar_bc"].default == "fixed"


def test_thums_feb_g7_recipe_emits_control_block() -> None:
    """``g7_solver_recipe=True`` ⇒ time_steps=2400, max_retries=20,
    cutback=0.125, opt_iter=15, dtmax=1/2400, BFGS max_ups=10,
    reform_each_time_step=1。"""
    import thums_feb as tf
    mesh = _fake_thums_mesh()
    out = Path("temp/s5_contact_facet_fix_test_g7.feb").resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    tf.build_thums_feb(
        mesh, out,
        use_rigid=False,
        plantar_bc="contact",
        contact_mu=0.6,
        g7_solver_recipe=True,
    )
    text = out.read_text(encoding="ISO-8859-1")
    assert "time_steps>2400</time_steps>" in text, "G7 time_steps=2400 not found"
    assert "max_retries>20</max_retries>" in text, "G7 max_retries=20 not found"
    assert "cutback>0.125</cutback>" in text, "G7 cutback=0.125 not found"
    assert "opt_iter>15</opt_iter>" in text, "G7 opt_iter=15 not found"
    # dtmax=1/2400 ≈ 4.166e-4
    assert "dtmax>0.0004166666666666667</dtmax>" in text, "G7 dtmax=1/2400 not found"
    # BFGS max_ups=10 + reform_each_time_step=1 是 SolidSolver 的 pyfebio 默认
    # （也是 FEBio 4.13 的内置默认）；``_fix_febio413`` 把 <solver> 块换成
    # ``<solver/>`` ⇒ FEBio 仍读默认 = G7 配方要求值。
    assert "<solver/>" in text, "G7 应保留 <solver/> 占位（默认值即 G7）"
    # 进一步：seed 一下 SolidSolver 默认值断言
    from pyfebio import control as ctl
    ss = ctl.SolidSolver()
    assert ss.qn_method.max_ups == 10, "pyfebio SolidSolver 默认 BFGS max_ups 应为 10"
    assert ss.reform_each_time_step == 1, "pyfebio SolidSolver 默认 reform_each_time_step 应为 1"


def test_thums_feb_default_no_g7_in_control_block() -> None:
    """默认（``g7_solver_recipe=False``）⇒ step 级 ``<Control>`` 是 time_steps=1,
    step_size=1.0；**不含任何 G7 值**（max_retries=20 / cutback=0.125 / dtmax=1/2400）。

    默认 path 与 wave-3 baseline 一致：pyfebio 仍发射 top-level + step 两份
    ``<Control>``（top-level 是 FEBio 默认值，step 级是调用方输入）。G7 只改 step 级。
    """
    import thums_feb as tf
    mesh = _fake_thums_mesh()
    out = Path("temp/s5_contact_facet_fix_test_default.feb").resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    tf.build_thums_feb(mesh, out, use_rigid=False, plantar_bc="fixed")
    text = out.read_text(encoding="ISO-8859-1")
    assert "time_steps>1</time_steps>" in text
    assert "step_size>1.0</step_size>" in text
    # 任何 G7 配方值都不得出现
    assert "max_retries>20</max_retries>" not in text, "默认 path 不应有 G7 max_retries=20"
    assert "cutback>0.125</cutback>" not in text, "默认 path 不应有 G7 cutback=0.125"
    assert "dtmax>0.0004166666666666667</dtmax>" not in text, "默认 path 不应有 G7 dtmax"


def test_thums_feb_load_ramp_emits_load_data() -> None:
    """``load_ramp_time_s=0.5`` ⇒ ``<LoadData><load_controller id=1>`` 出现，
    ramp points (0,0)→(0.5,1)。FEBio 4.13 只接受 id=1（单曲线 deck）。"""
    import thums_feb as tf
    mesh = _fake_thums_mesh()
    out = Path("temp/s5_contact_facet_fix_test_ramp.feb").resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    tf.build_thums_feb(
        mesh, out,
        use_rigid=False,
        plantar_bc="contact",
        contact_mu=0.6,
        load_ramp_time_s=0.5,
    )
    text = out.read_text(encoding="ISO-8859-1")
    assert "<LoadData>" in text, "G7 ramp 未发射 <LoadData>"
    assert 'load_controller id="1"' in text, "G7 ramp lc id=1 未出现"
    assert "0.0,0.0" in text and "0.5,1.0" in text, "G7 ramp points (0,0)→(0.5,1) 未出现"
    # PressureLoad 必须 reference lc=1
    assert 'lc="1"' in text, "ramp 曲线未被任何 load 引用"


def test_thums_feb_rigid_load_references_ramp_lc() -> None:
    """``use_rigid=True`` + ramp ⇒ RigidForceLoad 的 value 带 ``lc="1"`` 引用。"""
    import thums_feb as tf
    mesh = _fake_thums_mesh()
    out = Path("temp/s5_contact_facet_fix_test_rigid_ramp.feb").resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    tf.build_thums_feb(
        mesh, out,
        use_rigid=True,
        plantar_bc="contact",
        contact_mu=0.6,
        load_ramp_time_s=1.0,
    )
    text = out.read_text(encoding="ISO-8859-1")
    # rigid_load 的 value 应带 lc="1"（否则 load controller 是 unreferenced）
    import re
    m = re.search(r'<rigid_load[^>]*>.*?</rigid_load>', text, re.S)
    assert m is not None, "rigid_load 未发射"
    assert 'lc="1"' in m.group(0), "RigidForceLoad.value 必须引用 lc=1"


def test_thums_feb_default_no_load_data_block() -> None:
    """默认（``load_ramp_time_s=None``）⇒ 不写 ``<LoadData>`` ⇒ 字节级保留旧版。"""
    import thums_feb as tf
    mesh = _fake_thums_mesh()
    out = Path("temp/s5_contact_facet_fix_test_default2.feb").resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    tf.build_thums_feb(mesh, out, use_rigid=False, plantar_bc="fixed")
    text = out.read_text(encoding="ISO-8859-1")
    assert "<LoadData>" not in text, "默认不应发射 <LoadData>"


def test_thums_feb_real_plantar_is_quad4_in_deck() -> None:
    """真实 THUMS 网格（``elements_cort`` + ``boundary_polys``）→ 自动走 quad4 路径，
    跖面 Surface 不含 ``<tri3>``。"""
    import thums_feb as tf
    real_npz = ROOT / "temp" / "opensim_fe" / "thums_calcaneus" / "calcaneus_r_calcnframe.npz"
    if not real_npz.is_file():
        pytest.skip("real THUMS mesh not present")
    mesh = tf.load_thums_mesh(real_npz)
    out = Path("temp/s5_contact_facet_fix_test_real.feb").resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    tf.build_thums_feb(
        mesh, out,
        use_rigid=False,
        plantar_bc="contact",
        contact_mu=0.6,
        g7_solver_recipe=True,
        load_ramp_time_s=1.0,
    )
    text = out.read_text(encoding="ISO-8859-1")
    # 找 plantar Surface 块
    import re
    sm = re.search(r'<Surface name="plantar">(.*?)</Surface>', text, re.S)
    assert sm is not None, "no plantar Surface found"
    inner = sm.group(1)
    assert "<quad4" in inner, "plantar Surface 应有 quad4"
    assert "<tri3" not in inner, "plantar Surface 不应有 tri3（hex8 域不接受）"
    # 146 个 quad4（与 boundary_polys plantar_mask 命中数一致）
    n_quad = len(re.findall(r"<quad4\s", inner))
    assert n_quad == 146, f"plantar quad4 count {n_quad} ≠ expected 146"


def test_thums_feb_default_plantar_is_tri3_in_deck() -> None:
    """非 contact 路径（``plantar_bc='fixed'``）⇒ 不应有 quad4 / tri3 facet emission
    from plantar（因为 fixed 路径只发 1 个 BCZeroDisplacement，不建 surface）。"""
    import thums_feb as tf
    mesh = _fake_thums_mesh()
    out = Path("temp/s5_contact_facet_fix_test_fixed.feb").resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    tf.build_thums_feb(mesh, out, use_rigid=False, plantar_bc="fixed")
    text = out.read_text(encoding="ISO-8859-1")
    # fixed 路径下不应发射 <Surface name="plantar">
    assert '<Surface name="plantar"' not in text
