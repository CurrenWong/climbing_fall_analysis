"""守卫测试：`m_theta_table` 的**契约**（见 `CANONICAL.md` §0）。

目的：让「格式错」**当场报错**，而不是静默算错。
背景（2026-10-17 事故）：传 `(x_rad_array, y_array)` 元组时，`np.asarray` 会把它
变成 2×N 数组 ⇒ `_mt[:, 0]` 取到**整行 x** ⇒ 再 `np.radians` 一次 ⇒ 角度被压到
0~0.04 rad ⇒ `np.interp` 对真实 θ(~0.5–1.8 rad) 全部撞到末端 ⇒ **力矩被钉成常数**
⇒ 整条链的力静默算错，并导致同一处返工三次。
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import climbing.pad as P  # noqa: E402

_KEYS = ("k_theta_nm_per_rad", "segment_l_m", "theta_max_rad",
         "m_ecc_max_nm", "ecc_rate_ref_ms", "m_theta_table")
_DEG = np.arange(0.0, 135.0 + 1e-9, 2.5)
_M = np.full_like(_DEG, 200.0)
#: 契约格式：`[[θ_deg, M_Nm], ...]`（2 维、两列、度数、升序）
GOOD = [[float(a), float(b)] for a, b in zip(_DEG, _M)]


def _construct(table):
    """把 table 塞进 D1 门并触发构造。格式非法应在**求解前**抛错。"""
    prev = {k: P.TWO_SEG_LEG.get(k) for k in _KEYS}
    try:
        P.TWO_SEG_LEG.update({"k_theta_nm_per_rad": 350.0, "segment_l_m": 0.425,
                              "theta_max_rad": 1.94, "m_ecc_max_nm": 540.0,
                              "ecc_rate_ref_ms": 1.0, "m_theta_table": table})
        P.simulate_boulder_fall(height_m=0.5, mass_kg=80.0, on_pad=False,
                                posture="anchored-self")
    finally:
        for k in _KEYS:
            P.TWO_SEG_LEG[k] = prev[k]


def test_reject_tuple_of_two_arrays():
    """★ 事故原型：`(x_rad, y)` 元组必须被拒（否则力矩被钉成常数）。"""
    with pytest.raises(ValueError):
        _construct((np.radians(_DEG), _M))


def test_reject_theta_not_ascending():
    with pytest.raises(ValueError):
        _construct([[10.0, 200.0], [5.0, 210.0], [20.0, 220.0]])


def test_reject_radians_passed_as_degrees():
    """契约要求**度数**；传弧度会静默算出偏小的力矩 ⇒ 必须拒。"""
    with pytest.raises(ValueError):
        _construct([[0.0, 200.0], [1000.0, 210.0], [2000.0, 220.0]])


def test_contract_shape_and_conversion_are_self_consistent():
    """正面对照：契约格式是 2 维两列，且「度数 → 弧度」后的插值符合预期。"""
    arr = np.asarray(GOOD, dtype=float)
    assert arr.ndim == 2 and arr.shape[1] == 2 and arr.shape[0] >= 2
    assert float(np.interp(1.0, np.radians(arr[:, 0]), arr[:, 1])) == pytest.approx(200.0, abs=3.0)
