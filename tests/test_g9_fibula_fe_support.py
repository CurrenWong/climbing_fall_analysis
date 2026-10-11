"""tests/test_g9_fibula_fe_support.py -- G9 加法式 `lateral_support` 的回归锁。

锁住：
* 默认 `lateral_support=None` ⇒ `.feb` **无** `lat_support` 节点集/BC（逐位回归）。
* 给带 ⇒ 新增 `lat_support_k` 节点集 + **只约束横断面 2 自由度**（轴向 dof=0）。
* 空带 ⇒ 明确报错（不静默）。
* 返回元数据含 `lateral_support` / `lateral_support_bands`。

只用**合成小网格**（12 节点 / 2 hex8），**不调用 FEBio、不读 THUMS deck**（快速）。
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pytest

_ROOT = Path(__file__).resolve().parents[1]
for _p in (_ROOT, _ROOT / "src", _ROOT / "scripts" / "opensim_fe"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from bone_feb import (  # noqa: E402
    BoneFebError,
    BoneSpec,
    Material,
    SolidDomainSpec,
    build_bone_feb,
)


def _synth_spec() -> BoneSpec:
    mat = Material(pid=1, E_mpa=1000.0, nu=0.3, src_keyword="*MAT_ELASTIC",
                   src_card=[1, 0, 1000.0, 0.3])
    sd = SolidDomainSpec(pids=[1], mat_name="m", mesh_name="sol", material=mat)
    return BoneSpec(name="t", solid_domains=[sd])


def _synth_mesh() -> dict:
    # elongated box in Z with a mid layer (z=0,5,10) -> 12 nodes, 2 hex8
    nodes = np.array([
        [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
        [0, 0, 5], [1, 0, 5], [1, 1, 5], [0, 1, 5],
        [0, 0, 10], [1, 0, 10], [1, 1, 10], [0, 1, 10],
    ], dtype=float)
    return {
        "nodes": nodes,
        "by_solid": {"sol": {
            "tets": np.zeros((0, 4), np.int64),
            "hexes": np.array([[1, 2, 3, 4, 5, 6, 7, 8],
                               [5, 6, 7, 8, 9, 10, 11, 12]], np.int64)}},
        "by_shell": {},
    }


def _build(tmp_path: Path, name: str, **kw) -> tuple[str, dict]:
    feb = tmp_path / f"{name}.feb"
    meta = build_bone_feb(_synth_spec(), feb, load_n=10.0, mesh=_synth_mesh(),
                          align_pca=False, **kw)
    return feb.read_text(encoding="ISO-8859-1"), meta


def _norm(text: str) -> str:
    # the only run-to-run difference of identical builds is the plotfile name
    return re.sub(r'file="[^"]*\.xplt"', 'file="X.xplt"', text)


def test_default_has_no_lateral_support(tmp_path):
    txt, meta = _build(tmp_path, "a")
    assert "lat_support" not in txt
    assert meta["lateral_support"] is None
    assert meta["lateral_support_bands"] == []
    assert set(re.findall(r'set="([^"]+)"', txt)) == {"bottom_fix", "top_load"}


def test_default_path_identical_modulo_plotfile(tmp_path):
    t1, _ = _build(tmp_path, "a")
    t2, _ = _build(tmp_path, "c")
    assert _norm(t1) == _norm(t2)


def test_lateral_support_adds_transverse_only_bc(tmp_path):
    txt, meta = _build(tmp_path, "b", lateral_support=[(0.45, 0.55)])
    assert meta["lateral_support"] == [[0.45, 0.55]]
    assert meta["lateral_support_bands"] == ["lat_support_0"]
    assert "lat_support_0" in re.findall(r'set="([^"]+)"', txt)
    # mid-layer nodes 5..8 are in the support set
    m = re.search(r'<NodeSet name="lat_support_0">([^<]*)</NodeSet>', txt)
    assert m and set(m.group(1).split(",")) == {"5", "6", "7", "8"}
    # the lateral BC fixes only the two transverse dofs (load axis = Z here)
    bc = re.search(r'<bc type="zero displacement" node_set="lat_support_0">(.*?)</bc>',
                   txt, re.S)
    assert bc is not None
    assert "<x_dof>1</x_dof>" in bc.group(1)
    assert "<y_dof>1</y_dof>" in bc.group(1)
    assert "<z_dof>0</z_dof>" in bc.group(1)


def test_lateral_support_three_bands(tmp_path):
    # synthetic box has nodes only at z=0/5/10 -> pick 3 non-empty bands
    txt, meta = _build(tmp_path, "d",
                       lateral_support=[(0.0, 0.1), (0.45, 0.55), (0.9, 1.0)])
    assert meta["lateral_support_bands"] == ["lat_support_0", "lat_support_1",
                                             "lat_support_2"]
    assert txt.count('type="zero displacement"') == 4  # bottom_fix + 3 bands


def test_empty_band_raises(tmp_path):
    # z in [1.1, 3.9] contains no box node (nodes only at z=0/5/10)
    with pytest.raises(BoneFebError):
        _build(tmp_path, "e", lateral_support=[(0.11, 0.39)])


def test_load_axis_band_matches_axis_idx(tmp_path):
    # the synthetic box's longest bbox span is Z -> axis_idx == 2
    _, meta = _build(tmp_path, "f", lateral_support=[(0.45, 0.55)])
    assert meta["axis_idx"] == 2
