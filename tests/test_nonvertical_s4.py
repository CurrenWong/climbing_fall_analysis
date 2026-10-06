"""Phase-S4 非垂直落地：3D wrench + 跖面弹性支撑 的快速结构/数学回归。

不调用 FEBio（毫秒级）：锁死

* 默认参数（``subtalar_force=None, subtalar_moment=None``）→ ``.feb`` 逐位不变；
* 关节面几何 / 力矩→压力梯度 / 分带压力 的数学自洽（平面解析解）；
* opt-in 打开后确实新增了 wrench 载荷项。

慢 FE 数值回归（默认逐位复现改前 gauge 186.756）在
``scripts/opensim_fe/nonvertical_s4.py``，一次跑完并写报告。
"""

from __future__ import annotations

import inspect
import sys
from pathlib import Path

import numpy as np
import pytest

from climbing.coupling.febio_model import (
    build_calcaneus_feb,
    joint_face_frame,
    linear_pressure_bands,
    moment_to_pressure_gradient,
)

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts" / "opensim_fe") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts" / "opensim_fe"))


# --------------------------------------------------------------------------
# 合成网格
# --------------------------------------------------------------------------
def _cube_mesh() -> dict:
    """10 mm 立方体：顶面=关节面（2 三角），底面=跖面；5 个 tet。"""
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


def _plate_mesh(n: int = 8) -> dict:
    """10×10 平板：顶面 = n×n 网格三角（关节面），底面 4 角（跖面）。

    只用于纯几何/压力梯度数学测试（不需要 tet）。细网格让线性压力场能精确
    表达面内力矩。
    """

    def idx(i: int, j: int) -> int:
        return i * (n + 1) + j

    xs = np.linspace(0.0, 10.0, n + 1)
    top = np.array([[x, y, 10.0] for x in xs for y in xs], dtype=float)
    tris: list[list[int]] = []
    for i in range(n):
        for j in range(n):
            a, b, c, d = idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)
            tris.append([a, b, c])
            tris.append([a, c, d])
    tris_arr = np.asarray(tris, dtype=np.int64)
    bottom = np.array([[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]], dtype=float)
    nodes = np.vstack([top, bottom])
    nb = len(top)
    return {
        "nodes": nodes,
        "surfaces": {
            "subtalar_joint": tris_arr,
            "plantar": np.array([[nb, nb + 1, nb + 2], [nb, nb + 2, nb + 3]]),
        },
        "node_sets": {
            "subtalar_joint": np.arange(nb),
            "plantar": np.arange(nb, nb + 4),
        },
    }


# ==========================================================================
# signature defaults —— 默认逐位复现旧行为
# ==========================================================================
def test_build_calcaneus_signature_defaults() -> None:
    sig = inspect.signature(build_calcaneus_feb)
    assert sig.parameters["subtalar_force"].default is None
    assert sig.parameters["subtalar_moment"].default is None
    assert sig.parameters["moment_bands"].default == 6
    assert sig.parameters["plantar_bc"].default == "fixed"   # A7: 默认仍三向全固定
    assert sig.parameters["spring_k"].default == 1000.0
    assert sig.parameters["use_rigid"].default is True


def test_thums_build_signature_defaults() -> None:
    import thums_feb

    sig = inspect.signature(thums_feb.build_thums_feb)
    assert sig.parameters["subtalar_force"].default is None
    assert sig.parameters["subtalar_moment"].default is None
    assert sig.parameters["moment_bands"].default == 6
    assert sig.parameters["plantar_bc"].default == "fixed"


# ==========================================================================
# joint face frame —— 平面解析值
# ==========================================================================
def test_joint_face_frame_flat_square() -> None:
    frame = joint_face_frame(_plate_mesh(8), "subtalar_joint")
    assert frame["area_mm2"] == pytest.approx(100.0)
    assert frame["centroid_mm"] == pytest.approx([5.0, 5.0, 10.0])
    assert frame["normal_unit"] == pytest.approx([0.0, 0.0, 1.0])
    I = np.asarray(frame["second_moment_mm4"])
    assert I[0, 0] == pytest.approx(10.0 * 10.0**3 / 12.0, rel=1e-6)
    assert I[1, 1] == pytest.approx(10.0 * 10.0**3 / 12.0, rel=1e-6)
    assert abs(I[0, 1]) < 1e-6


def test_moment_gradient_zero_is_uniform() -> None:
    frame = joint_face_frame(_plate_mesh(8), "subtalar_joint")
    slope, direction, meta = moment_to_pressure_gradient((0.0, 0.0, 0.0), frame)
    assert slope == pytest.approx(0.0)
    assert np.allclose(direction, 0.0)
    assert meta["dropped_normal_axis_nm"] == pytest.approx(0.0)


def test_moment_gradient_flat_reproduces_inplane() -> None:
    frame = joint_face_frame(_plate_mesh(8), "subtalar_joint")
    # M=(1,0,0) N·m，对平面面精确可表达
    slope, direction, meta = moment_to_pressure_gradient((1.0, 0.0, 0.0), frame)
    assert meta["achievable_moment_nm"] == pytest.approx([1.0, 0.0, 0.0], abs=1e-6)
    assert meta["moment_residual_nm"] == pytest.approx(0.0, abs=1e-6)
    assert slope > 0.0
    assert abs(float(np.dot(direction, [0.0, 0.0, 1.0]))) < 1e-9


def test_linear_pressure_bands_force_conserved_and_moment() -> None:
    from pyfebio import mesh as fmesh
    from pyfebio.model import Model

    mesh = _plate_mesh(8)
    frame = joint_face_frame(mesh, "subtalar_joint")
    slope, direction, _ = moment_to_pressure_gradient((1.0, 0.0, 0.0), frame)
    bands, meta = linear_pressure_bands(
        Model(), fmesh, mesh, "subtalar_joint", frame,
        p0_mpa=7.0, slope=slope, direction=direction, n_bands=8,
    )
    # 合力精确守恒（梯度项零均值）：Σ P_b A_b == p0·A
    assert meta["net_normal_force_n"] == pytest.approx(7.0 * frame["area_mm2"], rel=1e-9)
    # 细网格 + 分带 → 实际力矩逼近目标
    assert meta["realized_moment_nm"] == pytest.approx([1.0, 0.0, 0.0], abs=0.05)
    assert 1 <= len(bands) <= 8


def test_wrench_bad_inputs_raise() -> None:
    mesh = _cube_mesh()
    with pytest.raises(ValueError):
        build_calcaneus_feb(mesh, "unused.feb", use_rigid=False, subtalar_force=(1.0, 2.0))
    with pytest.raises(ValueError):
        build_calcaneus_feb(
            mesh, "unused.feb", use_rigid=False, subtalar_moment=(float("nan"), 0.0, 0.0)
        )


# ==========================================================================
# .feb 输出：默认逐位不变 / opt-in 才新增 wrench 项
# ==========================================================================
pytest.importorskip("pyfebio")


def _build(tmp_path, name: str, **kw) -> Path:
    out = tmp_path / name / "m.feb"
    build_calcaneus_feb(_cube_mesh(), out, **kw)
    return out


def test_default_feb_bitwise_when_wrench_none(tmp_path) -> None:
    a = _build(tmp_path, "d1", use_rigid=False)
    b = _build(tmp_path, "d2", use_rigid=False, subtalar_force=None, subtalar_moment=None)
    assert a.read_bytes() == b.read_bytes()


def test_default_pressure_feb_has_single_uniform_load(tmp_path) -> None:
    a = _build(tmp_path, "d3", use_rigid=False)
    text = a.read_text(encoding="ISO-8859-1")
    assert text.count('type="pressure"') == 1
    assert 'type="traction"' not in text
    assert "rigid_moment" not in text
    assert "subtalar_joint_g" not in text


def test_moment_pressure_path_adds_gradient_and_traction(tmp_path) -> None:
    a = _build(
        tmp_path, "d4", use_rigid=False,
        subtalar_force=(0.0, -100.0, 0.0), subtalar_moment=(1.0, 0.0, 0.0),
    )
    text = a.read_text(encoding="ISO-8859-1")
    assert "subtalar_joint_g" in text          # 分带梯度面
    assert text.count('type="pressure"') > 1
    assert 'type="traction"' in text           # 切向牵引


def test_rigid_moment_path_adds_rigid_moment(tmp_path) -> None:
    a = _build(
        tmp_path, "d5", use_rigid=True,
        subtalar_force=(0.0, -100.0, 0.0), subtalar_moment=(1.0, 0.0, 0.0),
    )
    text = a.read_text(encoding="ISO-8859-1")
    assert 'type="rigid_moment"' in text
    assert 'type="rigid_force"' in text
