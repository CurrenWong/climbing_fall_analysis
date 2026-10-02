"""Phase A-4 / B-4 —— 伤害指标与阈值判定。

**数据性质警告**：本模块的阈值来自公开的生物力学/临床文献**量级**，
是共识区间而非精确判据。任何单个数字都应视为"这个量级的力/加速度
值得关注"，而不是"超过就必然骨折"。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "hic",
    "hic_from_result",
    "HIC_BANDS",
    "RISK_G",
    "classify_hic",
    "classify_g",
    "rope_risk",
    "pad_risk",
    "summary_table",
]


# --------------------------------------------------------------------------
# 头部伤害判据 HIC (Head Injury Criterion)
# --------------------------------------------------------------------------
#: 滑动窗口时长 (s)。HIC 的标准定义用 15 ms 和 50 ms 两个窗口，取最大。
HIC_WINDOWS_S = (0.015, 0.050)


def hic(t_s, accel_g, windows_s=HIC_WINDOWS_S) -> float:
    """计算 Head Injury Criterion。

    公式（Maag et al. 的标准形式）::

        HIC = max over (t1,t2) of  (t2-t1) * ( (1/(t2-t1)) ∫ |a| dt )^2.5

    Parameters
    ----------
    t_s
        时间 (s)，等间隔。
    accel_g
        头部线性加速度 (g)。**归一化到 1 g** 后代入（经典 HIC 即如此定义）。

    Returns
    -------
    float
        HIC 值。典型：<240 低风险；>700 脑震荡关注；>1000 高风险。

    Notes
    -----
    **为什么取 |a|**：HIC 定义里的加速度是**合矢量**（resultant），
    恒为非负。本模型是 1D 竖直的，冲击段给出的 ``a`` 为正（向上减速），
    但自由落体段是 ``-1 g``。早期版本直接对带符号的 ``a`` 取 2.5 次幂，
    负数 → NaN → ``if val > best`` 对 NaN 恒为 False → **HIC 静默返回 0**，
    看起来"所有场景都不受伤"。取绝对值既符合 HIC 的原始定义，也消除了
    这个静默失效。

    1D 竖直模型只能近似"头部合矢量"：本模块的 ``accel_*_g`` 传的是
    躯干/腿部代表点，不是真正的头颈合加速度。因此这里的 HIC 用于
    **场景间相对比较**，不应作为绝对临床判据。
    """
    t = np.asarray(t_s, dtype=float)
    a = np.asarray(accel_g, dtype=float)
    if t.ndim != 1 or a.shape != t.shape:
        raise ValueError("t_s 与 accel_g 必须同形一维数组")
    if t.size < 2:
        return 0.0
    dt = float(np.median(np.diff(t)))
    if dt <= 0:
        raise ValueError("时间必须严格递增")

    # HIC 的被积量是合矢量模长 -> 取 |a|（见 docstring 的 Notes）
    a = np.abs(a)
    if not np.all(np.isfinite(a)):
        raise ValueError("accel_g 含 NaN/Inf，无法计算 HIC")

    # 累积积分，滑动窗口用向量化切片
    cum = np.concatenate([[0.0], np.cumsum((a[:-1] + a[1:]) * 0.5 * dt)])
    best = 0.0
    for w in windows_s:
        n = max(2, int(round(w / dt)))
        if n >= t.size:
            continue
        # 用步长扫描（步长 = n//4，兼顾精度与速度）
        step = max(1, n // 8)
        for i in range(0, t.size - n, step):
            dv = cum[i + n] - cum[i]
            val = (n * dt) * (dv / (n * dt)) ** 2.5
            if val > best:
                best = val
    return float(best)


def hic_from_result(result) -> float:
    """从 ``RopeFallResult`` / ``BoulderFallResult`` 提取并计算 HIC。

    优先级：``primary_accel_g``（结果对象自己知道"哪块身体先着地"）
    -> ``accel_torso_g`` -> ``accel_g``。

    绳索场景取躯干（悬挂者的代表点）；抱石场景取**第一接触部位**——
    头朝下落地时那是 6.4 kg 的头而不是 73.6 kg 的躯干，HIC 算错对象
    会让最危险的姿势显示成最安全。
    """
    t = result.t_s
    for attr in ("primary_accel_g", "accel_torso_g", "accel_g"):
        if hasattr(result, attr):
            a = getattr(result, attr)
            break
    else:
        raise AttributeError(
            f"{type(result).__name__} 没有可用的加速度时程"
            "（需要 primary_accel_g / accel_torso_g / accel_g 之一）"
        )
    return hic(t, np.asarray(a, dtype=float))


# --------------------------------------------------------------------------
# 分级
# --------------------------------------------------------------------------
HIC_BANDS = [
    (0.0, 240.0, "低"),
    (240.0, 700.0, "中"),
    (700.0, 1000.0, "高"),
    (1000.0, float("inf"), "极高"),
]

#: 全身峰值加速度 (g) 的量级分档（参考跌倒 / 汽车碰撞生物力学）
RISK_G = [
    (0.0, 10.0, "低"),
    (10.0, 20.0, "中"),
    (20.0, 40.0, "高"),
    (40.0, float("inf"), "极高"),
]


def _band(value: float, table) -> str:
    for lo, hi, label in table:
        if lo <= value < hi:
            return label
    return table[-1][2]


def classify_hic(v: float) -> str:
    return _band(v, HIC_BANDS)


def classify_g(v: float) -> str:
    return _band(v, RISK_G)


# --------------------------------------------------------------------------
# 绳索场景汇总
# --------------------------------------------------------------------------
def rope_risk(res, uiaa_limit_n: float = 12_000.0) -> dict:
    """把 ``RopeFallResult`` 汇总成一行风险判读。"""
    return {
        "device": res.device,
        "FF": res.fall_factor,
        "L0_m": res.rope_length_m,
        "peak_kN": res.peak_tension_kn,
        "peak_g": res.peak_accel_g,
        "drop_m": res.max_drop_m,
        "slide_m": res.total_rope_slide_m,
        "elong_pct": res.max_elongation_pct,
        "payout_m": res.max_payout_m,
        "uiaa_ok": res.peak_tension_n < uiaa_limit_n,
        "g_risk": classify_g(res.peak_accel_g),
    }


# --------------------------------------------------------------------------
# 抱石场景汇总
# --------------------------------------------------------------------------
def pad_risk(res) -> dict:
    """把 ``BoulderFallResult`` 汇总成一行风险判读。"""
    return {
        "posture": res.posture,
        "surface": res.pad_label,
        "h_m": res.height_m,
        "v_ms": res.impact_speed_ms,
        "E_kJ": res.impact_energy_j / 1e3,
        "peak_kN": res.peak_force_kn,
        "torso_g": res.peak_accel_torso_g,
        "leg_g": res.peak_accel_leg_g,
        "comp_cm": res.max_compression_m * 100,
        "bottomed": res.bottomed_out,
        "absorbed": res.absorbed_fraction,
        "residual_kJ": res.energy_residual_j / 1e3,
        "g_risk": classify_g(res.peak_accel_torso_g),
    }


def summary_table(rows: list[dict]) -> str:
    """把 list[dict] 渲染成等宽表格。"""
    if not rows:
        return "(空)"
    cols = list(rows[0].keys())
    widths = {
        c: max(len(str(c)), *(len(_fmt(r.get(c))) for r in rows)) for c in cols
    }
    head = "  ".join(str(c).ljust(widths[c]) for c in cols)
    sep = "  ".join("-" * widths[c] for c in cols)
    body = "\n".join(
        "  ".join(_fmt(r.get(c)).ljust(widths[c]) for c in cols) for r in rows
    )
    return f"{head}\n{sep}\n{body}"


def _fmt(v) -> str:
    if isinstance(v, bool):
        return "是" if v else "否"
    if isinstance(v, float):
        if v != v:              # NaN
            return "-"
        if abs(v) >= 1000:
            return f"{v:,.0f}"
        if abs(v) >= 10:
            return f"{v:.1f}"
        if abs(v) >= 1:
            return f"{v:.2f}"
        return f"{v:.3f}"
    return str(v)
