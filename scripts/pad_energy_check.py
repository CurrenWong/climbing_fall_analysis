"""`pad.py` 结果的有效性判据 —— 沿轨迹的机械能守恒审计。

为什么需要它
------------
`pad.py` 自带的 `absorbed_fraction`（以及 `energy_residual_j`）**不能**用来判「数值是否物理」：
其分子 `cumtrapz(f_flex·d(xf)/dt, t)` 里 `f_flex` **含阻尼项** `c_flex·rate`（不可恢复、
随身体继续下沉累积、**不受 `mgh` 约束**），而分母 `impact_energy_j = m·g·h` 只对应
**首次接触那一刻** ⇒ 不同物理量 + 不同时间基准 ⇒ `>1` 是预期行为。
详见 `docs/notes/PAD_GOTCHAS.md` §6。

本模块给一个**同时刻、同时间基准**的判据：

    E_mech(t) = KE(t) + PE_grav(t) + E_elastic(t)     必须 **非增**（阻尼只耗散）

- `KE       = ½·m_low·vzf² + ½·m_up·vzt²`，速度**从 `result.impact_speed_ms` 偏移起积**
  （⚠️ `result.t_s` **只覆盖冲击窗**，窗起点速度是 `v_impact = √(2gh)`，**不是 0** —— 漏了这一步
  检验就无效，2026-10-16 已踩过一次）。
- 加速度符号：`pad.py` 报的 `accel_*_g` 以「向上减速为正」⇒ `az_down = −g·accel_*_g`。
- `PE_grav  = −[m_low·g·(z_foot−z_foot0) + m_up·g·(z_torso−z_torso0)]`，参考点取窗首。
- `E_elastic = ∫F_pad·d(comp) + 分段弹簧能（含 `k_stop` 止段）`。

用法::

    from scripts.pad_energy_check import audit_energy, ENERGY_TOL_J
    a = audit_energy(result, posture, mass_kg)
    assert a["violation_j"] <= ENERGY_TOL_J
"""
from __future__ import annotations

import numpy as np

G = 9.80665
LEG_EXT0 = 0.85          # 与 pad.py 的 leg_ext0 一致（躯干 CoM 在双足之上的初始高差）
K_STOP = 3.0e5           # 与 pad.py 的 k_stop 一致
#: 容许的数值残差（相对 `mgh` 的 1%）。超此值 ⇒ 判「能量被创造」。
ENERGY_TOL_FRAC = 0.01


def elastic_energy_flex(xf: np.ndarray, k_flex: float, x_flex_max: float) -> np.ndarray:
    """分段弹簧能：`½k·x²`；超过 `x_flex_max` 后接 `k_stop` 的二次段。"""
    x = np.maximum(np.asarray(xf, dtype=float), 0.0)
    lo = np.minimum(x, x_flex_max)
    over = np.maximum(x - x_flex_max, 0.0)
    return (0.5 * k_flex * lo ** 2
            + k_flex * x_flex_max * over
            + 0.5 * K_STOP * over ** 2)


def _cumtrapz(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    out = np.zeros_like(y)
    if y.size > 1:
        out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))
    return out


def audit_energy(result, posture, mass_kg: float = 80.0) -> dict:
    """沿轨迹审计机械能。返回 `violation_j`（>0 ⇒ 能量被创造）等诊断量。

    ⚠️ **数组口径（2026-10-16 实测钉死，`temp/v1_arrdump.py`）**：
    - `result.t_s` 从 **t = 0（释放瞬间）** 起（**不是**接触时刻）⇒ 窗首速度是 **0**，
      **不要**用 `impact_speed_ms` 作偏移（那是接触瞬间的速度，用错会系统性偏差 ≈2·mgh）；
    - `z_foot_m` = 首次接触点**相对释放点**的下落 ⇒ 接触压缩量 = `max(0, z_foot − height_m)`
      （**不是** `z_foot − z_foot[0]`）；
    - `z_torso` 初值 = `−leg_ext0` ⇒ 重力势参考点取 t=0 时 PE(0)=0。
    """
    t = np.asarray(result.t_s, dtype=float)
    zf = np.asarray(result.z_foot_m, dtype=float)
    zt = np.asarray(result.z_torso_m, dtype=float)
    al = np.asarray(result.accel_leg_g, dtype=float)
    at = np.asarray(result.accel_torso_g, dtype=float)
    f_pad = np.asarray(result.pad_force_n, dtype=float)
    h = float(getattr(result, "height_m", 0.0))

    m_low, m_up = posture.split_mass(mass_kg)

    # 向下为正的加速度（pad.py 的 accel_*_g 以「向上减速为正」）；v(t=0)=0
    vzf = _cumtrapz(-G * al, t)
    vzt = _cumtrapz(-G * at, t)

    ke = 0.5 * m_low * vzf ** 2 + 0.5 * m_up * vzt ** 2
    pe = -(m_low * G * zf + m_up * G * (zt + LEG_EXT0))     # t=0 时为 0
    xf = np.maximum(0.0, LEG_EXT0 - (zf - zt))
    e_flex = elastic_energy_flex(xf, posture.k_flex, posture.x_flex_max)
    comp = np.maximum(0.0, zf - h)
    e_pad = _cumtrapz(f_pad, comp)

    e = ke + pe + e_flex + e_pad          # 必须 ≤ 0（t=0 时为 0，之后只被阻尼拉低）
    viol = float(np.max(e))
    mgh = mass_kg * G * max(1e-9, h)
    return {
        "violation_j": viol,
        "tolerance_j": ENERGY_TOL_FRAC * mgh,
        "pass": bool(viol <= ENERGY_TOL_FRAC * mgh),
        "e_mech_rel_max": viol,
        "e0_j": float(e[0]),
        "e_end_rel_j": float(e[-1]),
        "monotone": bool(np.all(np.diff(e) <= 1e-6 * max(1.0, mgh))),
        "ke0_j": float(ke[0]),
        "max_xf_m": float(xf.max()),
        "max_rate_xf": float(np.max(np.gradient(xf, t))) if t.size > 1 else 0.0,
        "max_comp_m": float(comp.max()),
        "n_samples": int(t.size),
    }
