# -*- coding: utf-8 -*-
"""``climbing.coupling.loadchain`` 的单测 —— 闸门 10 的工具。

**关键：给检查器配「已知正/负对照」**（纪律 4）：
* **正对照**：同一张 → 同指纹；常数缩放 → 判为「常数缩放」而**不是**「同链」
* **负对照**：形状不同（比值随高度变）→ 判为「形状不同」而**不是**「常数缩放」
* **边界**：高度覆盖不同时，必须在**共有高度**上比（否则把同链误判为异链）
"""
from __future__ import annotations

import io
import json

import pytest

from climbing.coupling.loadchain import compare, fingerprint, table_from_json

CUR = {1.0: 18029.4, 2.0: 22723.586482, 3.0: 24286.3324, 5.0: 28520.6078}


def test_same_table_same_fingerprint():
    """正对照：同一张表 ⇒ 指纹相同。"""
    assert fingerprint(CUR) == fingerprint(dict(CUR))


def test_fingerprint_is_order_insensitive():
    items = list(CUR.items())
    rev = dict(reversed(items))
    assert fingerprint(CUR) == fingerprint(rev)


def test_fingerprint_absorbs_float_noise():
    """末位噪声（<6 位有效数字）不应改指纹。"""
    noisy = {h: v * (1 + 1e-12) for h, v in CUR.items()}
    assert fingerprint(CUR) == fingerprint(noisy)


def test_constant_scale_is_not_same_chain():
    """负对照 A：常数缩放 ⇒ `is_constant_scale` 且**非** same。"""
    c = compare(CUR, {h: v * 1.081 for h, v in CUR.items()})
    assert not c["same_fingerprint"]
    assert c["is_constant_scale"]
    assert c["ratio_min"] == pytest.approx(1.081, abs=1e-6)
    assert c["ratio_span_pct"] == pytest.approx(0.0, abs=1e-9)


def test_shape_change_detected():
    """负对照 B：形状变了（比值随高度变）⇒ `not is_constant_scale`。"""
    other = dict(CUR)
    other[5.0] = CUR[5.0] * 1.2          # 只动一个高度
    c = compare(CUR, other)
    assert not c["same_fingerprint"]
    assert not c["is_constant_scale"]
    assert c["ratio_span_pct"] > 1.0


def test_partial_height_coverage_is_same_chain():
    """边界：高度覆盖不同 ⇒ 必须在**共有高度**上比，否则误判。"""
    part = {1.0: CUR[1.0], 2.0: CUR[2.0]}          # 只覆盖 2 个高度
    c = compare(part, CUR)
    assert c["same_fingerprint"], "共有高度上同链却报异链 ⇒ 回退到整表取指纹了"
    assert c["n_common"] == 2


def test_table_from_json_rows_and_dict(tmp_path):
    """两种 JSON 形态（行列表 / {"rows": ...}）都要能读；_kN 自动归一到 N。"""
    rows = [{"h": 1, "act": 0, "subtalar_kN": 18.0294},
            {"h": 2, "act": 0, "subtalar_kN": 22.7236},
            {"h": 2, "act": 1, "subtalar_kN": 99.0}]
    for payload in (rows, {"rows": rows}):
        p = tmp_path / "a.json"
        p.write_text(json.dumps(payload), encoding="utf-8")
        t = table_from_json(p, act=0.0)
        assert set(t) == {1.0, 2.0}
        assert t[1.0] == pytest.approx(18029.4, rel=1e-12)
        assert t[2.0] == pytest.approx(22723.6, rel=1e-12)
    p1 = tmp_path / "b.json"
    p1.write_text(json.dumps([{"h": 1, "act": 0, "foo": 1}]), encoding="utf-8")
    with pytest.raises(ValueError):
        table_from_json(p1, act=0.0)


def test_archived_vs_current_is_shape_change():
    """**回归锚**（真实数据）：归档 sweep 与当前链是**形状不同**，不是常数缩放。

    若将来有人误把归档 json 当基线，本测会先炸。
    """
    import os
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]     # tests/ 就在仓库根下
    archived = root / "results" / "opensim_fe" / "s2_thums_height_sweep.json"
    current = root / "results" / "opensim_fe" / "LOADCHAIN_CURRENT.json"
    if not (archived.is_file() and current.is_file()):
        pytest.skip("归档 sweep 或当前链缓存缺失")
    a = table_from_json(archived, act=0.0)
    cur = {float(k): float(v)
           for k, v in json.load(io.open(current, encoding="utf-8")).items()}
    c = compare(a, cur)
    assert not c["same_fingerprint"]
    assert not c["is_constant_scale"], "归档是常数缩放？⇒ 说明链又变了，需复核"
    assert c["ratio_span_pct"] > 5.0
