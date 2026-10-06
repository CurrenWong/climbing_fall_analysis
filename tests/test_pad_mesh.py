"""S5-contact 回归测试：矩形 crash-pad 体网格生成器（``pad_block_tet``）。

锁死的契约（见 ``docs/S5_contact_contract.md`` §3.2 / §4）：

* 返回 dict 与 ``climbing.coupling.meshing.vtp_to_tet`` **同 schema**；
* 单位 mm；``surfaces`` 至少含 ``top``（+Z 接触面）、``bottom``（−Z 固定面）、
  ``sides``（四个侧面合并）；
* 裸方块体积 == ``sx*sy*sz``（误差 < 0.1%）；全部 tet 有**正** Jacobian；
  ``jacobian_frac_lt_0_3 == 0.0``；
* 粗网格（3×3×2）与细网格（8×8×6）**都**通过。

只依赖 numpy + pytest（``out_msh`` 用例可选 meshio）；不启动 gmsh、不跑 FEBio。
"""
from __future__ import annotations

import numpy as np
import pytest

from climbing.coupling.pad_mesh import pad_block_tet

# ---------------------------------------------------------------------------
# 常量 / 参数化
# ---------------------------------------------------------------------------
SIZE_MM: tuple[float, float, float] = (300.0, 300.0, 200.0)
SX, SY, SZ = SIZE_MM
EXPECTED_VOLUME_MM3 = SX * SY * SZ  # 18_000_000 mm³

#: 契约要求的三个具名面（top = 脚接触面 +Z，bottom = 固定面 −Z）。
REQUIRED_SURFACES = ("top", "bottom", "sides")

#: 契约 §3.2 要求返回的 schema 键（与 meshing.vtp_to_tet 的公共子集）。
REQUIRED_KEYS = (
    "nodes",
    "tets",
    "surfaces",
    "node_sets",
    "volume_mm3",
    "surface_area_mm2",
    "min_jacobian",
    "jacobian_frac_lt_0_3",
    "method",
    "n_nodes",
    "n_tets",
)

#: (nx, ny, nz) —— 粗 / 细两档，都必须通过。
COARSE = (3, 3, 2)
FINE = (8, 8, 6)
SHAPES = [
    pytest.param(COARSE, id="coarse-3x3x2"),
    pytest.param(FINE, id="fine-8x8x6"),
]


# ---------------------------------------------------------------------------
# 几何工具（纯 numpy，不依赖任何 mesh 库）
# ---------------------------------------------------------------------------
def _signed_tet_volumes(nodes: np.ndarray, tets: np.ndarray) -> np.ndarray:
    """每个 tet 的**有向**体积（det(b−a, c−a, d−a) / 6）。正 = 取向正确。"""
    q = nodes[tets]
    return (
        np.einsum(
            "ij,ij->i",
            q[:, 1] - q[:, 0],
            np.cross(q[:, 2] - q[:, 0], q[:, 3] - q[:, 0]),
        )
        / 6.0
    )


def _face_areas(nodes: np.ndarray, faces: np.ndarray) -> np.ndarray:
    p = nodes[faces]
    return 0.5 * np.linalg.norm(
        np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=1
    )


def _face_normals(nodes: np.ndarray, faces: np.ndarray) -> np.ndarray:
    p = nodes[faces]
    n = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    return n / np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-30)


def _expected_side_nodes(nx: int, ny: int, nz: int) -> int:
    """四个侧面的唯一节点数 = 全部节点 − (x,y) 内部节点（保留全部 z 层）。"""
    return (nx + 1) * (ny + 1) * (nz + 1) - (nx - 1) * (ny - 1) * (nz + 1)


# ---------------------------------------------------------------------------
# 1. schema / dtype
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("shape", SHAPES)
def test_schema_matches_meshing(shape: tuple[int, int, int]) -> None:
    """契约 §3.2：返回 dict 必须含 meshing.vtp_to_tet 同 schema 的键。"""
    nx, ny, nz = shape
    r = pad_block_tet(*shape, SIZE_MM)

    missing = set(REQUIRED_KEYS) - set(r)
    assert not missing, f"缺少 meshing 同 schema 键：{missing}"

    nodes, tets = r["nodes"], r["tets"]
    assert nodes.dtype == np.float64 and nodes.ndim == 2 and nodes.shape[1] == 3
    assert tets.dtype == np.int64 and tets.ndim == 2 and tets.shape[1] == 4
    assert r["n_nodes"] == len(nodes)
    assert r["n_tets"] == len(tets)
    assert isinstance(r["method"], str) and r["method"]

    # 连接关系是合法的 0-based 节点索引；每个节点都被引用（无孤立节点）
    assert tets.min() >= 0 and tets.max() < len(nodes)
    assert len(np.unique(tets)) == len(nodes)
    # 结构化计数
    assert len(nodes) == (nx + 1) * (ny + 1) * (nz + 1)
    assert len(tets) == 6 * nx * ny * nz


@pytest.mark.parametrize("shape", SHAPES)
def test_required_surfaces_exist_and_are_consistent(
    shape: tuple[int, int, int],
) -> None:
    """契约 §3.2：``surfaces``／``node_sets`` 至少含 top / bottom / sides。"""
    r = pad_block_tet(*shape, SIZE_MM)
    for name in REQUIRED_SURFACES:
        assert name in r["surfaces"], f"surfaces 缺 {name!r}"
        assert name in r["node_sets"], f"node_sets 缺 {name!r}"
        faces = r["surfaces"][name]
        nodes_ = r["node_sets"][name]
        assert faces.ndim == 2 and faces.shape[1] == 3
        assert np.issubdtype(faces.dtype, np.integer)
        assert nodes_.ndim == 1 and nodes_.dtype == np.int64
        assert faces.shape[0] > 0 and nodes_.shape[0] > 0
        assert faces.min() >= 0 and faces.max() < r["n_nodes"]
        assert nodes_.min() >= 0 and nodes_.max() < r["n_nodes"]
        # node_set 必须恰好是 faces 用到的节点的并集
        assert np.array_equal(nodes_, np.unique(faces))


# ---------------------------------------------------------------------------
# 2. 体积
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("shape", SHAPES)
def test_volume_matches_box_within_0_1_percent(shape: tuple[int, int, int]) -> None:
    """契约 §4：``volume_mm3`` == 300·300·200 mm³，误差 < 0.1%。"""
    r = pad_block_tet(*shape, SIZE_MM)
    err_pct = abs(r["volume_mm3"] - EXPECTED_VOLUME_MM3) / EXPECTED_VOLUME_MM3 * 100.0
    assert err_pct < 0.1, f"体积误差 {err_pct:.4g}% ≥ 0.1%"
    assert r["reference_volume_mm3"] == pytest.approx(EXPECTED_VOLUME_MM3, rel=1e-12)
    assert r["volume_error_pct"] < 0.1
    # 有向体积之和也必须等于盒子体积
    signed = _signed_tet_volumes(r["nodes"], r["tets"])
    assert abs(float(signed.sum()) - EXPECTED_VOLUME_MM3) / EXPECTED_VOLUME_MM3 < 1e-3


# ---------------------------------------------------------------------------
# 3. Jacobian / 单元质量
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("shape", SHAPES)
def test_all_tets_have_positive_jacobian(shape: tuple[int, int, int]) -> None:
    """契约红线：所有 tet 必须有**正**有向体积（一致定向）。"""
    r = pad_block_tet(*shape, SIZE_MM)
    sv = _signed_tet_volumes(r["nodes"], r["tets"])
    assert np.all(sv > 0), (
        f"{shape}: {(sv <= 0).sum()} 个 tet 非正 Jacobian (min={sv.min():.3e})"
    )
    assert r["min_jacobian"] > 0.0
    assert r["min_jacobian"] <= 1.0 + 1e-12  # 上界（单位立方体理想 tet = 1）


@pytest.mark.parametrize("shape", SHAPES)
def test_jacobian_frac_lt_0_3_is_zero(shape: tuple[int, int, int]) -> None:
    """契约：裸方块 ``jacobian_frac_lt_0_3`` 恒为 0（结构化拆分保证）。"""
    r = pad_block_tet(*shape, SIZE_MM)
    assert r["jacobian_frac_lt_0_3"] == 0.0, (
        f"{shape}: 有 {r['jacobian_frac_lt_0_3']*100:.2f}% 的 tet scaled-Jacobian < 0.3"
    )


# ---------------------------------------------------------------------------
# 4. 具名面：节点数 / 面积 / 三角数与盒子尺寸一致
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("shape", SHAPES)
def test_surface_node_counts(shape: tuple[int, int, int]) -> None:
    """三个面的**去重节点数**与 (nx, ny, nz) 闭式一致。"""
    nx, ny, nz = shape
    r = pad_block_tet(*shape, SIZE_MM)
    assert len(np.unique(r["surfaces"]["top"])) == (nx + 1) * (ny + 1)
    assert len(np.unique(r["surfaces"]["bottom"])) == (nx + 1) * (ny + 1)
    assert len(np.unique(r["surfaces"]["sides"])) == _expected_side_nodes(nx, ny, nz)
    # node_sets 与 surfaces 的 unique 一致
    for name in REQUIRED_SURFACES:
        assert np.array_equal(r["node_sets"][name], np.unique(r["surfaces"][name]))


@pytest.mark.parametrize("shape", SHAPES)
def test_surface_triangle_counts(shape: tuple[int, int, int]) -> None:
    """三角数与结构化面片数一致（每 quad 切 2 三角）。"""
    nx, ny, nz = shape
    r = pad_block_tet(*shape, SIZE_MM)
    assert r["surfaces"]["top"].shape[0] == 2 * nx * ny
    assert r["surfaces"]["bottom"].shape[0] == 2 * nx * ny
    assert r["surfaces"]["sides"].shape[0] == 4 * nz * (nx + ny)


@pytest.mark.parametrize("shape", SHAPES)
def test_surface_areas_match_box(shape: tuple[int, int, int]) -> None:
    """面积与盒子解析面积一致：
    top = bottom = sx·sy；sides = 2·sx·sz + 2·sy·sz。
    """
    r = pad_block_tet(*shape, SIZE_MM)
    top_area = float(_face_areas(r["nodes"], r["surfaces"]["top"]).sum())
    bot_area = float(_face_areas(r["nodes"], r["surfaces"]["bottom"]).sum())
    side_area = float(_face_areas(r["nodes"], r["surfaces"]["sides"]).sum())

    assert top_area == pytest.approx(SX * SY, rel=1e-9)
    assert bot_area == pytest.approx(SX * SY, rel=1e-9)
    assert side_area == pytest.approx(2 * SX * SZ + 2 * SY * SZ, rel=1e-9)
    # 总面积 = 2(sx·sy) + 2(sx·sz) + 2(sy·sz) = 420_000 mm²
    assert r["surface_area_mm2"] == pytest.approx(
        2 * SX * SY + 2 * SX * SZ + 2 * SY * SZ, rel=1e-9
    )


@pytest.mark.parametrize("shape", SHAPES)
def test_top_is_plus_z_bottom_is_minus_z(shape: tuple[int, int, int]) -> None:
    """top = +Z（脚接触面）、bottom = −Z（固定面），法向朝外。"""
    r = pad_block_tet(*shape, SIZE_MM)
    nodes, s = r["nodes"], r["surfaces"]

    # 几何位置
    assert np.allclose(nodes[s["top"]][:, :, 2], SZ)
    assert np.allclose(nodes[s["bottom"]][:, :, 2], 0.0)
    # x/y 落在盒子截面内
    for name in ("top", "bottom"):
        c = nodes[s[name]].mean(axis=1)
        assert c[:, 0].min() >= 0.0 and c[:, 0].max() <= SX
        assert c[:, 1].min() >= 0.0 and c[:, 1].max() <= SY

    # 外法向符号
    n_top = _face_normals(nodes, s["top"])
    n_bot = _face_normals(nodes, s["bottom"])
    n_side = _face_normals(nodes, s["sides"])
    assert np.all(n_top[:, 2] > 0.99)
    assert np.all(n_bot[:, 2] < -0.99)
    assert np.all(np.abs(n_side[:, 2]) < 1e-9)


@pytest.mark.parametrize("shape", SHAPES)
def test_boundary_normals_are_outward(shape: tuple[int, int, int]) -> None:
    """模块自带的 ``boundary_normals`` 与几何外法向一致
    （顺序：top → sides → bottom）。
    """
    r = pad_block_tet(*shape, SIZE_MM)
    faces = r["boundary_faces"]
    normals = r["boundary_normals"]
    assert faces.shape == (
        r["surfaces"]["top"].shape[0]
        + r["surfaces"]["sides"].shape[0]
        + r["surfaces"]["bottom"].shape[0],
        3,
    )
    assert normals.shape == (len(faces), 3)
    geom = _face_normals(r["nodes"], faces)
    assert np.allclose(normals, geom, atol=1e-12)


@pytest.mark.parametrize("shape", SHAPES)
def test_boundary_faces_cover_closed_shell(shape: tuple[int, int, int]) -> None:
    """闭合外壳：每条边恰好被 2 个三角面共享（无裂缝 / 无内翻）。"""
    r = pad_block_tet(*shape, SIZE_MM)
    faces = np.concatenate(
        [r["surfaces"]["top"], r["surfaces"]["sides"], r["surfaces"]["bottom"]],
        axis=0,
    )
    edges: dict[tuple[int, int], int] = {}
    for a, b, c in faces:
        for u, v in ((a, b), (b, c), (c, a)):
            key = (int(min(u, v)), int(max(u, v)))
            edges[key] = edges.get(key, 0) + 1
    counts = np.array(list(edges.values()))
    assert counts.size > 0
    assert np.all(counts == 2), f"非闭合外壳：边共享数={sorted(set(counts.tolist()))}"


# ---------------------------------------------------------------------------
# 5. 粗 / 细两档都必须通过（显式重复，便于失败定位）
# ---------------------------------------------------------------------------
def test_coarse_mesh_passes() -> None:
    r = pad_block_tet(*COARSE, SIZE_MM)
    assert r["n_tets"] == 6 * 3 * 3 * 2
    assert r["min_jacobian"] > 0
    assert r["jacobian_frac_lt_0_3"] == 0.0
    assert abs(r["volume_mm3"] - EXPECTED_VOLUME_MM3) / EXPECTED_VOLUME_MM3 < 1e-3


def test_fine_mesh_passes() -> None:
    r = pad_block_tet(*FINE, SIZE_MM)
    assert r["n_tets"] == 6 * 8 * 8 * 6
    assert r["min_jacobian"] > 0
    assert r["jacobian_frac_lt_0_3"] == 0.0
    assert abs(r["volume_mm3"] - EXPECTED_VOLUME_MM3) / EXPECTED_VOLUME_MM3 < 1e-3


# ---------------------------------------------------------------------------
# 6. 非默认尺寸 / 节点范围 / 确定性
# ---------------------------------------------------------------------------
def test_non_default_size_scales() -> None:
    """自定义尺寸（200×400×100）：体积 / 面积随尺寸解析缩放。"""
    size = (200.0, 400.0, 100.0)
    nx, ny, nz = 4, 8, 2
    r = pad_block_tet(nx, ny, nz, size)
    sx, sy, sz = size
    assert r["volume_mm3"] == pytest.approx(sx * sy * sz, rel=1e-9)
    assert r["surface_area_mm2"] == pytest.approx(
        2 * (sx * sy + sx * sz + sy * sz), rel=1e-9
    )
    assert _face_areas(r["nodes"], r["surfaces"]["top"]).sum() == pytest.approx(
        sx * sy, rel=1e-9
    )
    assert _face_areas(r["nodes"], r["surfaces"]["sides"]).sum() == pytest.approx(
        2 * sx * sz + 2 * sy * sz, rel=1e-9
    )
    assert np.all(_signed_tet_volumes(r["nodes"], r["tets"]) > 0)
    assert r["jacobian_frac_lt_0_3"] == 0.0


def test_default_size_is_300_300_200() -> None:
    r = pad_block_tet(2, 2, 2)
    assert r["size_mm"] == (300.0, 300.0, 200.0)
    assert r["volume_mm3"] == pytest.approx(EXPECTED_VOLUME_MM3, rel=1e-12)


def test_nodes_span_the_box_exactly() -> None:
    """节点坐标范围必须是 [0, sx] × [0, sy] × [0, sz]（无外溢 / 无内缩）。"""
    r = pad_block_tet(5, 4, 3, SIZE_MM)
    n = r["nodes"]
    assert n[:, 0].min() == pytest.approx(0.0)
    assert n[:, 0].max() == pytest.approx(SX)
    assert n[:, 1].min() == pytest.approx(0.0)
    assert n[:, 1].max() == pytest.approx(SY)
    assert n[:, 2].min() == pytest.approx(0.0)
    assert n[:, 2].max() == pytest.approx(SZ)


def test_deterministic_across_calls() -> None:
    """同一输入两次调用必须逐位一致（结构化、无随机性）。"""
    a = pad_block_tet(*COARSE, SIZE_MM)
    b = pad_block_tet(*COARSE, SIZE_MM)
    assert np.array_equal(a["nodes"], b["nodes"])
    assert np.array_equal(a["tets"], b["tets"])
    assert np.array_equal(a["boundary_faces"], b["boundary_faces"])
    for name in REQUIRED_SURFACES:
        assert np.array_equal(a["surfaces"][name], b["surfaces"][name])


# ---------------------------------------------------------------------------
# 7. out_msh 落盘 + Physical Surface
# ---------------------------------------------------------------------------
def test_out_msh_written_and_readable(tmp_path) -> None:
    meshio = pytest.importorskip("meshio")
    out = tmp_path / "pad_coarse.msh"
    r = pad_block_tet(*COARSE, SIZE_MM, out_msh=str(out))
    assert out.is_file()
    assert r.get("msh_path") == str(out)

    m = meshio.read(str(out))
    assert m.points.shape[0] == r["n_nodes"]
    assert "tetra" in m.cells_dict
    assert len(m.cells_dict["tetra"]) == r["n_tets"]
    assert "triangle" in m.cells_dict
    # 具名 Physical Surface 三个
    assert {"top", "sides", "bottom"} <= set(m.field_data)


def test_out_msh_creates_parent_dir(tmp_path) -> None:
    out = tmp_path / "nested" / "deep" / "pad.msh"
    r = pad_block_tet(1, 1, 1, SIZE_MM, out_msh=str(out))
    assert out.is_file()
    assert r["msh_path"] == str(out)


# ---------------------------------------------------------------------------
# 8. 输入校验（失败必须显式抛错，不静默）
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("bad", [(0, 1, 1), (1, 0, 1), (1, 1, 0), (-1, 2, 2)])
def test_non_positive_counts_raise(bad: tuple) -> None:
    with pytest.raises(ValueError):
        pad_block_tet(*bad)


@pytest.mark.parametrize(
    "bad_size", [(0.0, 1.0, 1.0), (1.0, -1.0, 1.0), (1.0, 1.0, 0.0)]
)
def test_non_positive_size_raises(bad_size: tuple) -> None:
    with pytest.raises(ValueError):
        pad_block_tet(1, 1, 1, bad_size)


def test_wrong_size_arity_raises() -> None:
    with pytest.raises(ValueError):
        pad_block_tet(1, 1, 1, (1.0, 1.0))  # type: ignore[arg-type]


@pytest.mark.parametrize("bad", [(1.5, 2, 2), (2, 2.0, 2), (2, 2, "x")])
def test_non_integer_counts_raise(bad: tuple) -> None:
    with pytest.raises(ValueError):
        pad_block_tet(*bad)


# ---------------------------------------------------------------------------
# pad_block_hex（六面体路径）—— 见 docs/S5_contact_convergence.md §7
# 说明：hex8 的**单元**在 ξ≳0.5 反演（比 tet4 差），故未用于 S5 主线；
# 这里只锁住网格生成器本身的几何正确性（体积/面积/面朝向/单元数）。
# ---------------------------------------------------------------------------
class TestPadBlockHex:
    def test_geometry_exact(self) -> None:
        from climbing.coupling.pad_mesh import pad_block_hex

        m = pad_block_hex(6, 6, 4, (300.0, 300.0, 200.0))
        assert m["n_nodes"] == 7 * 7 * 5
        assert m["n_hexes"] == 6 * 6 * 4
        assert m["n_tets"] == 0
        assert abs(m["volume_error_pct"]) < 1e-6
        assert m["method"] == "structured-hex8"
        assert m["min_jacobian"] > 0.0

    def test_surfaces_are_outward_quad4(self) -> None:
        from climbing.coupling.pad_mesh import pad_block_hex

        m = pad_block_hex(3, 3, 2, (300.0, 300.0, 200.0))
        for name in ("top", "bottom", "sides"):
            assert m["surfaces"][name].shape[1] == 4, name
        pt = m["nodes"][m["surfaces"]["top"][0]]
        assert np.cross(pt[1] - pt[0], pt[2] - pt[0])[2] > 0.0      # +Z
        pb = m["nodes"][m["surfaces"]["bottom"][0]]
        assert np.cross(pb[1] - pb[0], pb[2] - pb[0])[2] < 0.0      # -Z

    def test_surface_area_matches_box(self) -> None:
        from climbing.coupling.pad_mesh import pad_block_hex

        m = pad_block_hex(2, 2, 2, (100.0, 100.0, 100.0))
        assert abs(m["surface_area_mm2"] - 60000.0) < 1e-6

    def test_node_set_sizes(self) -> None:
        from climbing.coupling.pad_mesh import pad_block_hex

        m = pad_block_hex(4, 4, 2, (200.0, 200.0, 100.0))
        assert len(m["node_sets"]["top"]) == 5 * 5
        assert len(m["node_sets"]["bottom"]) == 5 * 5