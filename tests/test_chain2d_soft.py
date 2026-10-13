"""G6 · 载荷段软组织串联柔顺 —— 单元 + 旁路隔离 + 回归保护。

对应 ``docs/methods/G6_载荷段软组织串联柔顺_预注册.md`` §3/§4/§10。
**纯新增**：本文件只 import ``climbing.chain2d_soft`` 与只读的 ``chain2d``。

重点保护三件事：
1. **旁路隔离**：``SEGMENTS_6`` **未被修改**（默认路径逐位不变）；
2. **公式推导**（非手抄）：``k_eff`` 与串联式自洽（预注册表曾手抄错 10 倍）；
3. **判据观测量定义**不退化（首版用「足节点速度过零」无区分力）。
"""

from __future__ import annotations

import numpy as np
import pytest

from climbing.chain2d import SEGMENTS_6
from climbing.chain2d_soft import (
    ALPHA_ST,
    E_ST_PA,
    SEGMENTS_6_SOFT,
    arrest_stiffness_n_per_m,
    chain_effective_stiffness_n_per_m,
    segments_soft,
)

#: 冻结锚点（预注册 §2/§3）—— 测试只断言**不被改**。
A1_LO_KN, A1_HI_KN = 30.0, 170.0
A2_LO_UE, A2_HI_UE = 237.0, 1250.0


# --------------------------------------------------------------------------
# 旁路隔离：SEGMENTS_6 必须原封不动
# --------------------------------------------------------------------------
def test_base_segments_untouched():
    """默认段定义**不得**被本件修改（跑前冻结的解剖参数）。"""
    k = {s.name: s.k_axial_n_per_m for s in SEGMENTS_6}
    assert k["小腿-胫骨"] == pytest.approx(1.3e7)
    assert k["大腿"] == pytest.approx(8.0e6)
    assert k["足"] == pytest.approx(2.0e6)
    # 观测量的分母（骨应变）也须原封不动
    s = SEGMENTS_6[1]
    assert s.tissue_area_m2 == pytest.approx(3.20e-4)
    assert s.tissue_modulus_pa == pytest.approx(1.7e10)


def test_soft_segments_only_change_axial_stiffness():
    """旁路段只改 ``k_axial``，其余字段（含观测量分母）逐值不变。"""
    for a, b in zip(SEGMENTS_6, SEGMENTS_6_SOFT):
        assert b.name == a.name
        assert b.length_m == a.length_m
        assert b.mass_frac == a.mass_frac
        assert b.tissue_area_m2 == a.tissue_area_m2
        assert b.tissue_modulus_pa == a.tissue_modulus_pa
        assert b.contact_area_m2 == a.contact_area_m2
        assert b.k_axial_n_per_m <= a.k_axial_n_per_m      # 只会更软


def test_soft_is_softer_but_torso_less_affected():
    """软组织元件保留「远端骨主导」分层：躯干段受影响远小于胫骨段。"""
    ratio = {a.name: a.k_axial_n_per_m / b.k_axial_n_per_m
             for a, b in zip(SEGMENTS_6, SEGMENTS_6_SOFT)}
    assert ratio["小腿-胫骨"] > 100.0          # 远端骨主导 ⇒ 大幅变软
    assert ratio["骨盆-躯干"] < 30.0           # 近端 ⇒ 影响小
    assert ratio["足"] < 10.0


# --------------------------------------------------------------------------
# 公式推导（非手抄）
# --------------------------------------------------------------------------
def test_k_eff_derived_from_formula():
    """``k_eff = 1/(1/k_bone + 1/k_soft)``，``k_soft = E_st·α·A_cort/L``。"""
    for a, b in zip(SEGMENTS_6, SEGMENTS_6_SOFT):
        k_soft = E_ST_PA * ALPHA_ST * a.tissue_area_m2 / a.length_m
        expect = 1.0 / (1.0 / a.k_axial_n_per_m + 1.0 / k_soft)
        assert b.k_axial_n_per_m == pytest.approx(expect, rel=1e-12)


def test_foot_row_arithmetic_regression():
    """**回归保护**：预注册表 `足` 行曾算术错 10 倍（漏除 L）。

    足段：``k_soft = 10e6×3×1.50e-3/0.067 = 6.716e5``（**非** 6.72e6）。
    """
    foot = SEGMENTS_6_SOFT[0]
    k_soft = E_ST_PA * ALPHA_ST * SEGMENTS_6[0].tissue_area_m2 / SEGMENTS_6[0].length_m
    assert k_soft == pytest.approx(6.716e5, rel=1e-3)
    assert k_soft != pytest.approx(6.716e6, rel=1e-2)          # 旧的错值
    assert foot.k_axial_n_per_m == pytest.approx(5.028e5, rel=1e-3)


def test_anchors_frozen():
    assert E_ST_PA == pytest.approx(10.0e6)
    assert ALPHA_ST == pytest.approx(3.0)


def test_soft_variant_is_softer_than_base():
    kb = chain_effective_stiffness_n_per_m(SEGMENTS_6)
    ks = chain_effective_stiffness_n_per_m(SEGMENTS_6_SOFT)
    assert ks < kb / 10.0


def test_segments_soft_is_pure_function():
    """``segments_soft`` 不改入参（``SEGMENTS_6`` 不被就地修改）、可重复给同值。"""
    a = segments_soft()
    b = segments_soft()
    assert [s.k_axial_n_per_m for s in a] == [s.k_axial_n_per_m for s in b]
    src = tuple(s.k_axial_n_per_m for s in SEGMENTS_6)
    segments_soft()
    assert tuple(s.k_axial_n_per_m for s in SEGMENTS_6) == src


# --------------------------------------------------------------------------
# 判据观测量：arrest_stiffness 不得退化
# --------------------------------------------------------------------------
def test_arrest_stiffness_monotone_in_t_arrest():
    """``k = m·(π/(2T))²``：T 越小 k 越大；且不得为常数（退化检测）。"""
    ks = [arrest_stiffness_n_per_m(80.0, T) for T in (0.005, 0.01, 0.02, 0.05)]
    assert all(a > b for a, b in zip(ks, ks[1:]))     # 严格递减
    assert len(set(round(k, 3) for k in ks)) == len(ks)   # 有区分力（非退化）
    assert arrest_stiffness_n_per_m(80.0, 0.05) == pytest.approx(78.96e3, rel=2e-2)


def test_arrest_stiffness_invalid_input():
    assert np.isnan(arrest_stiffness_n_per_m(80.0, 0.0))
    assert np.isnan(arrest_stiffness_n_per_m(80.0, -1.0))
